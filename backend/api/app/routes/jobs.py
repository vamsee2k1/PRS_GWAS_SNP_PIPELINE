from __future__ import annotations

from uuid import UUID

from fastapi import APIRouter, Depends, Query

from app.deps import get_orchestrator
from app.schemas import (
    DispatchResponse,
    JobCreateRequest,
    JobListResponse,
    JobResponse,
    WorkerCallbackRequest,
    WorkerEventRequest,
    WorkerJobSpec,
)
from app.services.job_orchestrator import JobOrchestrator

router = APIRouter(prefix="/api/v1", tags=["jobs"])


@router.post("/jobs", response_model=JobResponse)
def create_job(
    req: JobCreateRequest,
    orchestrator: JobOrchestrator = Depends(get_orchestrator),
) -> JobResponse:
    return orchestrator.create_job(req)


@router.get("/jobs", response_model=JobListResponse)
def list_jobs(
    limit: int = Query(default=50, ge=1, le=200),
    orchestrator: JobOrchestrator = Depends(get_orchestrator),
) -> JobListResponse:
    return JobListResponse(items=orchestrator.list_jobs(limit=limit))


@router.get("/jobs/{job_id}", response_model=JobResponse)
def get_job(job_id: UUID, orchestrator: JobOrchestrator = Depends(get_orchestrator)) -> JobResponse:
    return orchestrator.get_job(job_id)


@router.get("/jobs/{job_id}/events")
def get_job_events(job_id: UUID, orchestrator: JobOrchestrator = Depends(get_orchestrator)):
    return {"items": [e.model_dump() for e in orchestrator.list_events(job_id)]}


@router.post("/jobs/{job_id}/dispatch", response_model=DispatchResponse)
def dispatch_job(job_id: UUID, orchestrator: JobOrchestrator = Depends(get_orchestrator)) -> DispatchResponse:
    return orchestrator.dispatch_job(job_id)


@router.get("/worker/jobs/{job_id}/spec", response_model=WorkerJobSpec, tags=["worker"])
def worker_get_spec(
    job_id: UUID,
    token: str = Query(...),
    orchestrator: JobOrchestrator = Depends(get_orchestrator),
) -> WorkerJobSpec:
    return orchestrator.get_worker_spec(job_id, token)


@router.post("/worker/jobs/{job_id}/events", tags=["worker"])
def worker_event(
    job_id: UUID,
    req: WorkerEventRequest,
    orchestrator: JobOrchestrator = Depends(get_orchestrator),
):
    orchestrator.record_worker_event(job_id, req)
    return {"ok": True}


@router.post("/worker/jobs/{job_id}/callback", response_model=JobResponse, tags=["worker"])
def worker_callback(
    job_id: UUID,
    req: WorkerCallbackRequest,
    orchestrator: JobOrchestrator = Depends(get_orchestrator),
) -> JobResponse:
    return orchestrator.handle_worker_callback(job_id, req)
