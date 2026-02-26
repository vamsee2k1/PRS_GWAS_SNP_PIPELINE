from __future__ import annotations

from datetime import datetime
from typing import Any, Literal
from uuid import UUID

from pydantic import BaseModel, Field


JobMode = Literal["full", "variant_only", "vcf_interpretation", "gwas_summary"]
VariantDataMode = Literal["auto", "vcf_interpretation", "gwas_summary"]
JobStatus = Literal["queued", "dispatching", "running", "succeeded", "failed", "canceled"]
EventLevel = Literal["info", "warning", "error"]


class JobCreateRequest(BaseModel):
    name: str | None = None
    mode: JobMode = "variant_only"
    variant_data_mode: VariantDataMode = "vcf_interpretation"
    configfile: str = Field(..., description="Config file path on the pipeline VM/repo")
    output_dir: str = Field(..., description="Output directory for this run (relative to pipeline repo root)")
    cores: int = Field(default=4, ge=1, le=128)
    run_until: str | None = Field(default=None, description="Optional Snakemake --until target")
    snakemake_args: list[str] = Field(default_factory=list, description="Extra Snakemake CLI args")
    repo_ref: str | None = Field(default=None, description="Optional git ref to checkout on VM")
    input_payload: dict[str, Any] = Field(default_factory=dict, description="Opaque user inputs/metadata")
    labels: dict[str, str] = Field(default_factory=dict, description="UI labels/tags")


class JobEventRecord(BaseModel):
    id: int | None = None
    job_id: UUID
    level: EventLevel = "info"
    event_type: str
    message: str
    payload: dict[str, Any] = Field(default_factory=dict)
    created_at: datetime


class JobRecord(BaseModel):
    id: UUID
    name: str | None = None
    mode: JobMode
    variant_data_mode: VariantDataMode
    status: JobStatus
    configfile: str
    output_dir: str
    cores: int
    run_until: str | None = None
    snakemake_args: list[str] = Field(default_factory=list)
    repo_ref: str | None = None
    input_payload: dict[str, Any] = Field(default_factory=dict)
    labels: dict[str, str] = Field(default_factory=dict)
    dispatch_payload: dict[str, Any] = Field(default_factory=dict)
    result_payload: dict[str, Any] = Field(default_factory=dict)
    error_message: str | None = None
    callback_token: str
    created_at: datetime
    updated_at: datetime
    started_at: datetime | None = None
    finished_at: datetime | None = None


class JobResponse(BaseModel):
    id: UUID
    name: str | None = None
    mode: JobMode
    variant_data_mode: VariantDataMode
    status: JobStatus
    configfile: str
    output_dir: str
    cores: int
    run_until: str | None = None
    snakemake_args: list[str]
    repo_ref: str | None = None
    input_payload: dict[str, Any]
    labels: dict[str, str]
    dispatch_payload: dict[str, Any]
    result_payload: dict[str, Any]
    error_message: str | None = None
    created_at: datetime
    updated_at: datetime
    started_at: datetime | None = None
    finished_at: datetime | None = None


class JobListResponse(BaseModel):
    items: list[JobResponse]


class DispatchResponse(BaseModel):
    job: JobResponse
    dispatch: dict[str, Any]


class WorkerJobSpec(BaseModel):
    job_id: UUID
    configfile: str
    output_dir: str
    cores: int
    run_until: str | None = None
    snakemake_args: list[str] = Field(default_factory=list)
    repo_ref: str
    api_base_url: str


class WorkerEventRequest(BaseModel):
    token: str
    level: EventLevel = "info"
    event_type: str
    message: str
    payload: dict[str, Any] = Field(default_factory=dict)


class WorkerCallbackRequest(BaseModel):
    token: str
    status: Literal["running", "succeeded", "failed", "canceled"]
    message: str
    payload: dict[str, Any] = Field(default_factory=dict)
    results_zip_path: str | None = None
    results_zip_uri: str | None = None
    pipeline_log_path: str | None = None
    run_manifest_path: str | None = None


def to_job_response(job: JobRecord) -> JobResponse:
    data = job.model_dump(exclude={"callback_token"})
    return JobResponse(**data)
