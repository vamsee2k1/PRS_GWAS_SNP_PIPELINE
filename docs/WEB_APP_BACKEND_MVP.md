# Web App Backend MVP (Supabase + AWS VM)

This document describes the first backend/API layer for turning the pipeline into a web app.

## Architecture (MVP)

- Frontend (future): uploads, job forms, status dashboard
- Backend API (`backend/api`): FastAPI + job orchestration
- Supabase: auth, database (job metadata/events), optional storage later
- AWS EC2 VM: runs Snakemake jobs
- AWS SSM: dispatches commands to the VM without SSH dependency
- Results: zipped on VM and reported back to API (`results_zip_path`) in MVP

## What is implemented now

- FastAPI job API
- Job status/event tracking schema for Supabase
- AWS SSM dispatch integration
- Worker script that runs pipeline and zips results
- In-memory local dev mode (`mock` runner, no Supabase required)

## API endpoints (MVP)

- `POST /api/v1/jobs` create job
- `GET /api/v1/jobs` list jobs
- `GET /api/v1/jobs/{job_id}` get job
- `GET /api/v1/jobs/{job_id}/events` list events
- `POST /api/v1/jobs/{job_id}/dispatch` dispatch to runner
- `GET /api/v1/worker/jobs/{job_id}/spec?token=...` worker fetches job spec
- `POST /api/v1/worker/jobs/{job_id}/events` worker posts events
- `POST /api/v1/worker/jobs/{job_id}/callback` worker updates final status

## Supabase setup

1. Create a Supabase project.
2. Open SQL editor.
3. Run `backend/supabase/schema.sql`.
4. Get:
   - `SUPABASE_URL`
   - `SUPABASE_SERVICE_ROLE_KEY`

Use the service-role key only in the backend API service (not the browser).

## AWS VM setup (MVP)

1. Launch EC2 instance (Ubuntu recommended).
2. Install:
   - Git
   - Miniforge/Mambaforge
   - pipeline dependencies (or bootstrap env)
3. Install and enable SSM Agent (Amazon Linux often has it by default).
4. Attach IAM role with `AmazonSSMManagedInstanceCore`.
5. Clone this repo to a stable path, e.g. `/opt/PRS_GWAS_SNP_PIPELINE`.
6. Ensure the VM can reach your API URL.

## Local backend start (mock mode)

```bash
cd backend/api
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cp .env.example .env
# leave APP_REPOSITORY_BACKEND=memory and APP_RUNNER_BACKEND=mock
uvicorn app.main:app --reload
```

## Example job creation (variant-only)

```bash
curl -s http://localhost:8000/api/v1/jobs \
  -H 'content-type: application/json' \
  -d '{
    "name":"mini-vcf-test",
    "mode":"variant_only",
    "variant_data_mode":"vcf_interpretation",
    "configfile":"config/run.1000g_chr1_variant_only.yaml",
    "output_dir":"results_1000g_chr1_variant_only_test",
    "cores":4,
    "snakemake_args":["--rerun-incomplete"]
  }' | jq
```

Dispatch the returned job id:

```bash
curl -X POST http://localhost:8000/api/v1/jobs/<job-id>/dispatch | jq
```

In mock mode, dispatch is simulated.

## Next production upgrades (recommended)

1. Upload artifacts to object storage instead of only local zip-path reporting:
   - S3 (simplest with AWS VM)
   - or Supabase Storage signed upload URL
2. Add auth (Supabase JWT verification in the API)
3. Add input upload flow (signed URLs)
4. Add per-user access control / RLS
5. Add job cancellation (`ssm cancel-command` or process control on VM)
6. Add timeline parsing endpoint from Snakemake logs
