#!/usr/bin/env bash
set -euo pipefail

# Wrapper for better UX during long runs: shows spinner during quiet phases.
# Keep caches inside the project when running in restricted/sandboxed environments.
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$PWD/.cache}"
mkdir -p "$XDG_CACHE_HOME"

runtime_source_cache="${PWD}/.cache/snakemake/runtime-source-cache"
mkdir -p "$runtime_source_cache"

extra_args=()
if [[ " $* " != *" --runtime-source-cache-path "* ]]; then
  extra_args+=(--runtime-source-cache-path "$runtime_source_cache")
fi

if command -v snakemake >/dev/null 2>&1; then
  python3 workflow/scripts/run_with_loading.py -- snakemake "${extra_args[@]}" "$@"
else
  python3 workflow/scripts/run_with_loading.py -- conda run -n bioinfo-workflow snakemake "${extra_args[@]}" "$@"
fi
