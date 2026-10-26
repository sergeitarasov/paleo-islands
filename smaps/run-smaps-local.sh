#!/bin/bash
set -e
cd "$(dirname "$0")/.."
mode="${1:-pilot}"
case "$mode" in
  pilot) run_dir="smaps/runs/posterior-aic-saved-asr-pilot" ;;
  full) run_dir="smaps/runs/posterior-aic-saved-asr" ;;
  *) echo "Usage: bash smaps/run-smaps-local.sh pilot|full" >&2; exit 1 ;;
esac
mkdir -p "$run_dir"
# Prevent simultaneous launchers from writing to the same output files.
if ! mkdir "$run_dir/launcher.lock" 2>/dev/null; then
  echo "Another launcher holds $run_dir/launcher.lock; check running jobs before removing a stale lock." >&2
  exit 1
fi
trap 'rmdir "$run_dir/launcher.lock"' EXIT
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
if [ -f "$run_dir/run.log" ]; then
  mv "$run_dir/run.log" "$run_dir/run-$(date +%Y%m%d-%H%M%S).log"
fi
caffeinate -i Rscript smaps/run-smaps-posterior-parallel.R "$mode" > "$run_dir/run.log" 2>&1
