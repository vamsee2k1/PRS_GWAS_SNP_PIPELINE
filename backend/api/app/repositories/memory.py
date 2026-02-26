from __future__ import annotations

from typing import Any
from uuid import UUID, uuid4

from app.repositories.base import JobRepository, utcnow
from app.schemas import JobEventRecord, JobRecord


class MemoryJobRepository(JobRepository):
    def __init__(self) -> None:
        self._jobs: dict[UUID, JobRecord] = {}
        self._events: list[JobEventRecord] = []
        self._next_event_id = 1

    def create_job(self, payload: dict[str, Any]) -> JobRecord:
        now = utcnow()
        row = {
            "id": payload.get("id") or uuid4(),
            "created_at": payload.get("created_at") or now,
            "updated_at": payload.get("updated_at") or now,
            "started_at": payload.get("started_at"),
            "finished_at": payload.get("finished_at"),
            "dispatch_payload": payload.get("dispatch_payload", {}),
            "result_payload": payload.get("result_payload", {}),
            "error_message": payload.get("error_message"),
            **payload,
        }
        job = JobRecord(**row)
        self._jobs[job.id] = job
        return job

    def get_job(self, job_id: UUID) -> JobRecord | None:
        return self._jobs.get(job_id)

    def list_jobs(self, limit: int = 50) -> list[JobRecord]:
        items = sorted(self._jobs.values(), key=lambda j: j.created_at, reverse=True)
        return items[:limit]

    def update_job(self, job_id: UUID, changes: dict[str, Any]) -> JobRecord:
        current = self._jobs[job_id]
        merged = current.model_copy(update={**changes, "updated_at": changes.get("updated_at") or utcnow()})
        self._jobs[job_id] = merged
        return merged

    def append_event(self, payload: dict[str, Any]) -> JobEventRecord:
        event = JobEventRecord(
            id=self._next_event_id,
            created_at=payload.get("created_at") or utcnow(),
            **payload,
        )
        self._events.append(event)
        self._next_event_id += 1
        return event

    def list_events(self, job_id: UUID) -> list[JobEventRecord]:
        items = [e for e in self._events if e.job_id == job_id]
        return sorted(items, key=lambda e: e.created_at)
