from __future__ import annotations

import os
from dataclasses import dataclass

try:
    from dotenv import load_dotenv
except Exception:  # pragma: no cover - optional in some runtime contexts
    load_dotenv = None

if load_dotenv:
    load_dotenv()


def _get_bool(name: str, default: bool) -> bool:
    raw = os.getenv(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def _get_int(name: str, default: int) -> int:
    raw = os.getenv(name)
    if raw is None or not raw.strip():
        return default
    return int(raw)


@dataclass(frozen=True)
class Settings:
    app_name: str = os.getenv("APP_NAME", "PRS GWAS SNP Pipeline API")
    app_env: str = os.getenv("APP_ENV", "dev")
    app_host: str = os.getenv("APP_HOST", "0.0.0.0")
    app_port: int = _get_int("APP_PORT", 8000)
    repository_backend: str = os.getenv("APP_REPOSITORY_BACKEND", "memory").strip().lower()
    runner_backend: str = os.getenv("APP_RUNNER_BACKEND", "mock").strip().lower()

    supabase_url: str = os.getenv("SUPABASE_URL", "").strip()
    supabase_service_role_key: str = os.getenv("SUPABASE_SERVICE_ROLE_KEY", "").strip()
    supabase_jobs_table: str = os.getenv("SUPABASE_JOBS_TABLE", "pipeline_jobs").strip()
    supabase_job_events_table: str = os.getenv("SUPABASE_JOB_EVENTS_TABLE", "pipeline_job_events").strip()

    aws_region: str = os.getenv("AWS_REGION", "us-east-1").strip()
    aws_ssm_instance_id: str = os.getenv("AWS_SSM_INSTANCE_ID", "").strip()
    aws_ssm_document_name: str = os.getenv("AWS_SSM_DOCUMENT_NAME", "AWS-RunShellScript").strip()
    aws_vm_pipeline_root: str = os.getenv("AWS_VM_PIPELINE_ROOT", "/opt/PRS_GWAS_SNP_PIPELINE").strip()
    aws_vm_python: str = os.getenv("AWS_VM_PYTHON", "python3").strip()
    aws_vm_worker_script: str = os.getenv("AWS_VM_WORKER_SCRIPT", "backend/aws/worker_runner.py").strip()

    api_public_base_url: str = os.getenv("API_PUBLIC_BASE_URL", "http://localhost:8000").rstrip("/")
    pipeline_repo_ref: str = os.getenv("PIPELINE_REPO_REF", "main").strip()

    allow_mock_dispatch: bool = _get_bool("ALLOW_MOCK_DISPATCH", True)

    def validate(self) -> None:
        if self.repository_backend not in {"memory", "supabase"}:
            raise ValueError("APP_REPOSITORY_BACKEND must be one of: memory, supabase")
        if self.runner_backend not in {"mock", "aws_ssm"}:
            raise ValueError("APP_RUNNER_BACKEND must be one of: mock, aws_ssm")
        if self.repository_backend == "supabase":
            if not self.supabase_url or not self.supabase_service_role_key:
                raise ValueError("SUPABASE_URL and SUPABASE_SERVICE_ROLE_KEY are required for Supabase backend")
        if self.runner_backend == "aws_ssm":
            if not self.aws_ssm_instance_id:
                raise ValueError("AWS_SSM_INSTANCE_ID is required for aws_ssm runner backend")
