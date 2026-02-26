# Backend API MVP (Supabase + AWS VM Runner)

This service provides an API layer for submitting pipeline jobs, tracking status in Supabase, and dispatching runs to an AWS VM using SSM.

## What this MVP covers

- `FastAPI` endpoints for job creation, dispatch, status, events, and worker callbacks
- Supabase-backed job metadata storage (with in-memory fallback for local dev)
- AWS SSM command dispatch to a pre-provisioned VM
- Worker script to run the Snakemake pipeline and produce a zipped results folder

## Local development (mock mode)

1. Create a virtualenv and install dependencies:

```bash
cd backend/api
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

2. Copy `.env.example` to `.env` and keep:

- `APP_REPOSITORY_BACKEND=memory`
- `APP_RUNNER_BACKEND=mock`

3. Run the API:

```bash
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

4. Open docs:

- `http://localhost:8000/docs`

## Moving to Supabase + AWS

- Apply SQL in `backend/supabase/schema.sql`
- Set `APP_REPOSITORY_BACKEND=supabase`
- Set `APP_RUNNER_BACKEND=aws_ssm`
- Configure AWS instance and SSM agent
- Copy the repo (including `backend/aws/worker_runner.py`) to the VM

Detailed steps are in `docs/WEB_APP_BACKEND_MVP.md`.
