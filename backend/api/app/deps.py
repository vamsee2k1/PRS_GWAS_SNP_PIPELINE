from __future__ import annotations

from functools import lru_cache

from app.config import Settings
from app.repositories.memory import MemoryJobRepository
from app.repositories.supabase_repo import SupabaseJobRepository
from app.services.aws_ssm_runner import AwsSsmRunner, MockRunner
from app.services.job_orchestrator import JobOrchestrator


@lru_cache(maxsize=1)
def get_settings() -> Settings:
    settings = Settings()
    settings.validate()
    return settings


@lru_cache(maxsize=1)
def get_repository():
    settings = get_settings()
    if settings.repository_backend == "memory":
        return MemoryJobRepository()
    return SupabaseJobRepository(settings)


@lru_cache(maxsize=1)
def get_runner():
    settings = get_settings()
    if settings.runner_backend == "mock":
        return MockRunner(settings)
    return AwsSsmRunner(settings)


@lru_cache(maxsize=1)
def get_orchestrator() -> JobOrchestrator:
    return JobOrchestrator(get_settings(), get_repository(), get_runner())
