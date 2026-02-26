from __future__ import annotations

from abc import ABC, abstractmethod
from datetime import datetime
from typing import Any
from uuid import UUID

from app.schemas import JobEventRecord, JobRecord


class JobRepository(ABC):
    @abstractmethod
    def create_job(self, payload: dict[str, Any]) -> JobRecord:
        raise NotImplementedError

    @abstractmethod
    def get_job(self, job_id: UUID) -> JobRecord | None:
        raise NotImplementedError

    @abstractmethod
    def list_jobs(self, limit: int = 50) -> list[JobRecord]:
        raise NotImplementedError

    @abstractmethod
    def update_job(self, job_id: UUID, changes: dict[str, Any]) -> JobRecord:
        raise NotImplementedError

    @abstractmethod
    def append_event(self, payload: dict[str, Any]) -> JobEventRecord:
        raise NotImplementedError

    @abstractmethod
    def list_events(self, job_id: UUID) -> list[JobEventRecord]:
        raise NotImplementedError


def utcnow() -> datetime:
    from datetime import timezone

    return datetime.now(timezone.utc)
