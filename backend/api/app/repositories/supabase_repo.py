from __future__ import annotations

from typing import Any
from uuid import UUID

from app.config import Settings
from app.repositories.base import JobRepository, utcnow
from app.schemas import JobEventRecord, JobRecord


class SupabaseJobRepository(JobRepository):
    def __init__(self, settings: Settings) -> None:
        try:
            from supabase import create_client
        except Exception as exc:  # pragma: no cover - import path depends on installed extras
            raise RuntimeError("supabase package is required for SupabaseJobRepository") from exc

        self.settings = settings
        self._client = create_client(settings.supabase_url, settings.supabase_service_role_key)

    def _jobs_table(self):
        return self._client.table(self.settings.supabase_jobs_table)

    def _events_table(self):
        return self._client.table(self.settings.supabase_job_events_table)

    def create_job(self, payload: dict[str, Any]) -> JobRecord:
        now = utcnow()
        row = {
            **payload,
            "created_at": payload.get("created_at") or now.isoformat(),
            "updated_at": payload.get("updated_at") or now.isoformat(),
        }
        res = self._jobs_table().insert(row).execute()
        data = (res.data or [None])[0]
        if not data:
            raise RuntimeError("Supabase insert returned no job row")
        return JobRecord(**data)

    def get_job(self, job_id: UUID) -> JobRecord | None:
        res = self._jobs_table().select("*").eq("id", str(job_id)).limit(1).execute()
        data = (res.data or [None])[0]
        return JobRecord(**data) if data else None

    def list_jobs(self, limit: int = 50) -> list[JobRecord]:
        res = self._jobs_table().select("*").order("created_at", desc=True).limit(limit).execute()
        return [JobRecord(**row) for row in (res.data or [])]

    def update_job(self, job_id: UUID, changes: dict[str, Any]) -> JobRecord:
        payload = {**changes, "updated_at": (changes.get("updated_at") or utcnow()).isoformat()}
        if "started_at" in payload and payload["started_at"] is not None and hasattr(payload["started_at"], "isoformat"):
            payload["started_at"] = payload["started_at"].isoformat()
        if "finished_at" in payload and payload["finished_at"] is not None and hasattr(payload["finished_at"], "isoformat"):
            payload["finished_at"] = payload["finished_at"].isoformat()
        res = self._jobs_table().update(payload).eq("id", str(job_id)).execute()
        data = (res.data or [None])[0]
        if not data:
            # Some Supabase/PostgREST setups return no updated rows unless select() is chained.
            got = self.get_job(job_id)
            if not got:
                raise RuntimeError(f"Job not found after update: {job_id}")
            return got
        return JobRecord(**data)

    def append_event(self, payload: dict[str, Any]) -> JobEventRecord:
        row = {
            **payload,
            "job_id": str(payload["job_id"]),
            "created_at": (payload.get("created_at") or utcnow()).isoformat(),
        }
        res = self._events_table().insert(row).execute()
        data = (res.data or [None])[0]
        if not data:
            raise RuntimeError("Supabase insert returned no event row")
        return JobEventRecord(**data)

    def list_events(self, job_id: UUID) -> list[JobEventRecord]:
        res = (
            self._events_table()
            .select("*")
            .eq("job_id", str(job_id))
            .order("created_at", desc=False)
            .execute()
        )
        return [JobEventRecord(**row) for row in (res.data or [])]
