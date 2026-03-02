#!/usr/bin/env bash
set -euo pipefail

# Wrapper for better UX during long runs: shows spinner during quiet phases.
# Keep caches inside the project when running in restricted/sandboxed environments.
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$PWD/.cache}"
mkdir -p "$XDG_CACHE_HOME"

runtime_source_cache="${PWD}/.cache/snakemake/runtime-source-cache"
mkdir -p "$runtime_source_cache"

# Avoid stale local source snapshots in Snakemake runtime cache.
# This has caused old Snakefile content to be reused after updates.
local_cache_root="${runtime_source_cache}/file/${PWD#/}"
if [[ -d "$local_cache_root" ]]; then
  rm -rf "$local_cache_root"
fi

extra_args=()
if [[ " $* " != *" --runtime-source-cache-path "* ]]; then
  extra_args+=(--runtime-source-cache-path "$runtime_source_cache")
fi

loading_cmd=(python3 workflow/scripts/run_with_loading.py)
# Default to concise terminal progress. Set PIPELINE_VERBOSE=1 for full streaming output.
if [[ "${PIPELINE_VERBOSE:-0}" != "1" ]]; then
  loading_cmd+=(--clean-snakemake)
fi

if command -v snakemake >/dev/null 2>&1; then
  "${loading_cmd[@]}" -- snakemake "${extra_args[@]}" "$@"
else
  "${loading_cmd[@]}" -- conda run -n bioinfo-workflow snakemake "${extra_args[@]}" "$@"
fi
