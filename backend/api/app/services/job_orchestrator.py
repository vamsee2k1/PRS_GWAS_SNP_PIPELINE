from __future__ import annotations

import secrets
from datetime import timezone
from typing import Any
from uuid import UUID

from fastapi import HTTPException, status

from app.config import Settings
from app.repositories.base import JobRepository, utcnow
from app.schemas import (
    DispatchResponse,
    JobCreateRequest,
    JobRecord,
    JobResponse,
    WorkerCallbackRequest,
    WorkerEventRequest,
    WorkerJobSpec,
    to_job_response,
)
from app.services.aws_ssm_runner import JobRunner


class JobOrchestrator:
    def __init__(self, settings: Settings, repo: JobRepository, runner: JobRunner) -> None:
        self.settings = settings
        self.repo = repo
        self.runner = runner

    def create_job(self, req: JobCreateRequest) -> JobResponse:
        now = utcnow()
        callback_token = secrets.token_urlsafe(24)
        job = self.repo.create_job(
            {
                "name": req.name,
                "mode": req.mode,
                "variant_data_mode": req.variant_data_mode,
                "status": "queued",
                "configfile": req.configfile,
                "output_dir": req.output_dir,
                "cores": req.cores,
                "run_until": req.run_until,
                "snakemake_args": req.snakemake_args,
                "repo_ref": req.repo_ref or self.settings.pipeline_repo_ref,
                "input_payload": req.input_payload,
                "labels": req.labels,
                "dispatch_payload": {},
                "result_payload": {},
                "error_message": None,
                "callback_token": callback_token,
                "created_at": now,
                "updated_at": now,
            }
        )
        self.repo.append_event(
            {
                "job_id": job.id,
                "level": "info",
                "event_type": "job_created",
                "message": "Job created and queued.",
                "payload": {},
            }
        )
        return to_job_response(job)

    def list_jobs(self, limit: int = 50) -> list[JobResponse]:
        return [to_job_response(j) for j in self.repo.list_jobs(limit=limit)]

    def get_job(self, job_id: UUID) -> JobResponse:
        return to_job_response(self._require_job(job_id))

    def list_events(self, job_id: UUID):
        self._require_job(job_id)
        return self.repo.list_events(job_id)

    def dispatch_job(self, job_id: UUID) -> DispatchResponse:
        job = self._require_job(job_id)
        if job.status in {"running", "succeeded"}:
            raise HTTPException(status_code=status.HTTP_409_CONFLICT, detail=f"Job already {job.status}")

        job = self.repo.update_job(job_id, {"status": "dispatching", "error_message": None})
        self.repo.append_event(
            {
                "job_id": job_id,
                "level": "info",
                "event_type": "dispatch_requested",
                "message": "Dispatching job to runner backend.",
                "payload": {"runner_backend": self.settings.runner_backend},
            }
        )

        try:
            dispatch_payload = self.runner.dispatch(job)
        except Exception as exc:
            failed = self.repo.update_job(job_id, {"status": "failed", "error_message": str(exc), "finished_at": utcnow()})
            self.repo.append_event(
                {
                    "job_id": job_id,
                    "level": "error",
                    "event_type": "dispatch_failed",
                    "message": "Runner dispatch failed.",
                    "payload": {"error": str(exc)},
                }
            )
            raise HTTPException(status_code=500, detail=f"Dispatch failed: {exc}") from exc

        job = self.repo.update_job(job_id, {"dispatch_payload": dispatch_payload})
        self.repo.append_event(
            {
                "job_id": job_id,
                "level": "info",
                "event_type": "dispatched",
                "message": "Job dispatched to runner.",
                "payload": dispatch_payload,
            }
        )
        return DispatchResponse(job=to_job_response(job), dispatch=dispatch_payload)

    def get_worker_spec(self, job_id: UUID, token: str) -> WorkerJobSpec:
        job = self._require_job(job_id)
        self._validate_token(job, token)
        return WorkerJobSpec(
            job_id=job.id,
            configfile=job.configfile,
            output_dir=job.output_dir,
            cores=job.cores,
            run_until=job.run_until,
            snakemake_args=job.snakemake_args,
            repo_ref=job.repo_ref or self.settings.pipeline_repo_ref,
            api_base_url=self.settings.api_public_base_url,
        )

    def record_worker_event(self, job_id: UUID, req: WorkerEventRequest) -> None:
        job = self._require_job(job_id)
        self._validate_token(job, req.token)
        self.repo.append_event(
            {
                "job_id": job_id,
                "level": req.level,
                "event_type": req.event_type,
                "message": req.message,
                "payload": req.payload,
            }
        )

    def handle_worker_callback(self, job_id: UUID, req: WorkerCallbackRequest) -> JobResponse:
        job = self._require_job(job_id)
        self._validate_token(job, req.token)

        now = utcnow()
        changes: dict[str, Any] = {
            "status": req.status,
            "error_message": None if req.status != "failed" else req.message,
        }

        if req.status == "running":
            if job.started_at is None:
                changes["started_at"] = now
        else:
            if job.started_at is None:
                changes["started_at"] = now
            changes["finished_at"] = now

        result_payload = dict(job.result_payload or {})
        if req.results_zip_path:
            result_payload["results_zip_path"] = req.results_zip_path
        if req.results_zip_uri:
            result_payload["results_zip_uri"] = req.results_zip_uri
        if req.pipeline_log_path:
            result_payload["pipeline_log_path"] = req.pipeline_log_path
        if req.run_manifest_path:
            result_payload["run_manifest_path"] = req.run_manifest_path
        if req.payload:
            result_payload["worker_payload"] = req.payload
        changes["result_payload"] = result_payload

        job = self.repo.update_job(job_id, changes)
        self.repo.append_event(
            {
                "job_id": job_id,
                "level": "error" if req.status == "failed" else "info",
                "event_type": f"worker_callback_{req.status}",
                "message": req.message,
                "payload": req.payload,
            }
        )
        return to_job_response(job)

    def _require_job(self, job_id: UUID) -> JobRecord:
        job = self.repo.get_job(job_id)
        if not job:
            raise HTTPException(status_code=404, detail="Job not found")
        return job

    @staticmethod
    def _validate_token(job: JobRecord, token: str) -> None:
        if token != job.callback_token:
            raise HTTPException(status_code=401, detail="Invalid worker token")
