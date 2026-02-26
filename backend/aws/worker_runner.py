#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import traceback
from pathlib import Path
from urllib import error, parse, request


def env(name: str) -> str:
    val = os.getenv(name, "").strip()
    if not val:
        raise RuntimeError(f"Missing required environment variable: {name}")
    return val


def http_json(method: str, url: str, payload: dict | None = None) -> dict:
    data = None
    headers = {}
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    req = request.Request(url, data=data, headers=headers, method=method)
    with request.urlopen(req, timeout=60) as resp:
        body = resp.read().decode("utf-8")
        return json.loads(body) if body else {}


def post_event(base: str, job_id: str, token: str, event_type: str, message: str, level: str = "info", payload=None):
    url = f"{base}/api/v1/worker/jobs/{job_id}/events"
    body = {
        "token": token,
        "level": level,
        "event_type": event_type,
        "message": message,
        "payload": payload or {},
    }
    return http_json("POST", url, body)


def post_callback(base: str, job_id: str, token: str, status: str, message: str, **kwargs):
    url = f"{base}/api/v1/worker/jobs/{job_id}/callback"
    body = {"token": token, "status": status, "message": message}
    body.update(kwargs)
    return http_json("POST", url, body)


def get_spec(base: str, job_id: str, token: str) -> dict:
    q = parse.urlencode({"token": token})
    url = f"{base}/api/v1/worker/jobs/{job_id}/spec?{q}"
    return http_json("GET", url)


def main() -> int:
    api_base = env("API_BASE_URL").rstrip("/")
    job_id = env("JOB_ID")
    token = env("JOB_TOKEN")
    pipeline_root = Path(env("PIPELINE_ROOT")).resolve()

    backend_runs_dir = pipeline_root / "backend_job_runs"
    backend_runs_dir.mkdir(parents=True, exist_ok=True)
    combined_log = backend_runs_dir / f"{job_id}.combined.log"

    spec = get_spec(api_base, job_id, token)
    post_event(api_base, job_id, token, "worker_started", "Worker started and fetched job spec.")
    post_callback(api_base, job_id, token, "running", "Pipeline run started.", pipeline_log_path=str(combined_log))

    configfile = spec["configfile"]
    output_dir = spec["output_dir"]
    cores = str(spec["cores"])
    run_until = spec.get("run_until")
    snakemake_args = list(spec.get("snakemake_args") or [])

    cmd = [
        str(pipeline_root / "run_pipeline.sh"),
        "--use-conda",
        "--conda-frontend",
        "mamba",
        "--cores",
        cores,
        "--configfile",
        configfile,
    ]
    if run_until:
        cmd.extend(["--until", run_until])
    cmd.extend(snakemake_args)

    post_event(api_base, job_id, token, "pipeline_command", "Launching pipeline command.", payload={"cmd": cmd})

    rc = 0
    try:
        with combined_log.open("w", encoding="utf-8") as logh:
            proc = subprocess.Popen(
                cmd,
                cwd=str(pipeline_root),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
            for line in proc.stdout or []:
                logh.write(line)
            rc = proc.wait()
    except Exception as exc:
        tb = traceback.format_exc()
        try:
            post_event(api_base, job_id, token, "worker_exception", str(exc), level="error", payload={"traceback": tb})
            post_callback(
                api_base,
                job_id,
                token,
                "failed",
                f"Worker exception: {exc}",
                pipeline_log_path=str(combined_log),
                payload={"traceback": tb},
            )
        except Exception:
            pass
        return 1

    output_path = (pipeline_root / output_dir).resolve()
    manifest_path = output_path / "docs" / "run_manifest.txt"
    zip_path = None

    if rc == 0:
        if output_path.exists():
            archive_base = backend_runs_dir / f"{job_id}_results"
            zip_path = shutil.make_archive(str(archive_base), "zip", root_dir=str(output_path))
        msg = "Pipeline run completed successfully."
        post_event(api_base, job_id, token, "pipeline_finished", msg, payload={"return_code": rc})
        post_callback(
            api_base,
            job_id,
            token,
            "succeeded",
            msg,
            pipeline_log_path=str(combined_log),
            run_manifest_path=str(manifest_path) if manifest_path.exists() else None,
            results_zip_path=zip_path,
            payload={"return_code": rc},
        )
        return 0

    msg = f"Pipeline run failed with exit code {rc}."
    post_event(api_base, job_id, token, "pipeline_failed", msg, level="error", payload={"return_code": rc})
    post_callback(
        api_base,
        job_id,
        token,
        "failed",
        msg,
        pipeline_log_path=str(combined_log),
        run_manifest_path=str(manifest_path) if manifest_path.exists() else None,
        payload={"return_code": rc},
    )
    return rc


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except error.HTTPError as exc:
        sys.stderr.write(f"HTTP error: {exc.code} {exc.reason}\n")
        sys.stderr.write(exc.read().decode("utf-8", errors="replace"))
        sys.stderr.write("\n")
        raise
