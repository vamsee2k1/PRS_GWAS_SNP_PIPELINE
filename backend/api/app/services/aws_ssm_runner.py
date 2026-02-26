from __future__ import annotations

import shlex
from typing import Any, Protocol

from app.config import Settings
from app.schemas import JobRecord


class JobRunner(Protocol):
    def dispatch(self, job: JobRecord) -> dict[str, Any]:
        ...


class MockRunner:
    def __init__(self, settings: Settings) -> None:
        self.settings = settings

    def dispatch(self, job: JobRecord) -> dict[str, Any]:
        return {
            "runner": "mock",
            "message": "Mock dispatch completed. No AWS command was sent.",
            "job_id": str(job.id),
        }


class AwsSsmRunner:
    def __init__(self, settings: Settings) -> None:
        self.settings = settings
        try:
            import boto3
        except Exception as exc:  # pragma: no cover - depends on runtime env
            raise RuntimeError("boto3 is required for AWS SSM runner") from exc
        self._ssm = boto3.client("ssm", region_name=settings.aws_region)

    def _command(self, job: JobRecord) -> str:
        pipeline_root = self.settings.aws_vm_pipeline_root
        worker_script = self.settings.aws_vm_worker_script
        py = self.settings.aws_vm_python

        exports = {
            "API_BASE_URL": self.settings.api_public_base_url,
            "JOB_ID": str(job.id),
            "JOB_TOKEN": job.callback_token,
            "PIPELINE_ROOT": pipeline_root,
        }
        export_cmd = " ".join(f"{k}={shlex.quote(v)}" for k, v in exports.items())
        return (
            "set -euo pipefail; "
            f"cd {shlex.quote(pipeline_root)}; "
            f"{export_cmd} {shlex.quote(py)} {shlex.quote(worker_script)}"
        )

    def dispatch(self, job: JobRecord) -> dict[str, Any]:
        command = self._command(job)
        resp = self._ssm.send_command(
            InstanceIds=[self.settings.aws_ssm_instance_id],
            DocumentName=self.settings.aws_ssm_document_name,
            Parameters={"commands": [command]},
            Comment=f"Dispatch pipeline job {job.id}",
        )
        cmd = resp["Command"]
        return {
            "runner": "aws_ssm",
            "aws_region": self.settings.aws_region,
            "instance_id": self.settings.aws_ssm_instance_id,
            "ssm_command_id": cmd["CommandId"],
            "ssm_status": cmd["Status"],
            "commands_preview": [command],
        }
