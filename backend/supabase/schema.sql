-- Supabase schema for pipeline job orchestration MVP
-- Run in Supabase SQL editor (project database).
-- This script also creates private Storage buckets used by the web app backend.

create extension if not exists pgcrypto;

create table if not exists public.pipeline_jobs (
    id uuid primary key default gen_random_uuid(),
    name text,
    mode text not null,
    variant_data_mode text not null,
    status text not null default 'queued',
    configfile text not null,
    output_dir text not null,
    cores integer not null,
    run_until text,
    snakemake_args jsonb not null default '[]'::jsonb,
    repo_ref text,
    input_payload jsonb not null default '{}'::jsonb,
    labels jsonb not null default '{}'::jsonb,
    dispatch_payload jsonb not null default '{}'::jsonb,
    result_payload jsonb not null default '{}'::jsonb,
    error_message text,
    callback_token text not null,
    created_at timestamptz not null default now(),
    updated_at timestamptz not null default now(),
    started_at timestamptz,
    finished_at timestamptz
);

create index if not exists pipeline_jobs_status_idx on public.pipeline_jobs (status);
create index if not exists pipeline_jobs_created_at_idx on public.pipeline_jobs (created_at desc);

create table if not exists public.pipeline_job_events (
    id bigint generated always as identity primary key,
    job_id uuid not null references public.pipeline_jobs(id) on delete cascade,
    level text not null default 'info',
    event_type text not null,
    message text not null,
    payload jsonb not null default '{}'::jsonb,
    created_at timestamptz not null default now()
);

create index if not exists pipeline_job_events_job_id_idx on public.pipeline_job_events (job_id, created_at);

create or replace function public.set_updated_at_timestamp()
returns trigger
language plpgsql
as $$
begin
  new.updated_at = now();
  return new;
end;
$$;

drop trigger if exists trg_pipeline_jobs_updated_at on public.pipeline_jobs;
create trigger trg_pipeline_jobs_updated_at
before update on public.pipeline_jobs
for each row
execute procedure public.set_updated_at_timestamp();

-- Optional: add RLS later for user-facing apps.
-- MVP backend uses Supabase service-role key and can operate without RLS policies.

-- Storage buckets (idempotent)
-- Keep private; the backend can use signed URLs later for uploads/downloads.
insert into storage.buckets (id, name, public)
values ('uploads', 'uploads', false)
on conflict (id) do nothing;

insert into storage.buckets (id, name, public)
values ('results', 'results', false)
on conflict (id) do nothing;
