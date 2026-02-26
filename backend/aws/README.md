# AWS VM Worker (SSM-dispatched)

This directory contains the worker script executed on an AWS VM via SSM.

## VM prerequisites

- AWS SSM Agent installed and online
- IAM role on the EC2 instance with `AmazonSSMManagedInstanceCore`
- Pipeline repo checked out at the path configured in `AWS_VM_PIPELINE_ROOT`
- Conda/Mamba + pipeline dependencies available on the VM
- Network access from VM to the API service URL (`API_PUBLIC_BASE_URL`)

## What `worker_runner.py` does

1. Fetches job spec from the API using a per-job callback token
2. Runs `./run_pipeline.sh` with the configured `--configfile`
3. Writes a combined log to `backend_job_runs/<job_id>.combined.log`
4. Zips the configured `output_dir`
5. Sends final status + zip path back to the API

## Next step (production)

Replace local zip-path reporting with artifact upload (S3 or Supabase Storage signed upload URL), then store the resulting URI in `result_payload`.
