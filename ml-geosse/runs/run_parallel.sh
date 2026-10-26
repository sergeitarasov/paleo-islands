#!/usr/bin/env bash

set -euo pipefail

# Resolve paths relative to this file, so the launcher works when called from
# either the project root (`bash runs/run_parallel.sh`) or `runs/`.
RUNS_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${RUNS_DIR}/scripts"
LOGS_DIR="${RUNS_DIR}/logs"

mkdir -p "${LOGS_DIR}"
cd "${SCRIPTS_DIR}"

# Run every generated scenario. Current prefixes include Afr16,
# AfrSunken16, Oriental16, and OrientalSunken16, but discovery is dynamic so
# newly generated scenarios are included without another edit.
script_count="$(find . -maxdepth 1 -type f -name '*_gr-*_run-*.R' | wc -l | tr -d ' ')"

if (( script_count == 0 )); then
  echo "No generated scenario scripts found in ${SCRIPTS_DIR}." >&2
  exit 1
fi

printf 'Running %d generated scenario scripts. Logs: %s\n' \
  "${script_count}" "${LOGS_DIR}"

find . -maxdepth 1 -type f -name '*_gr-*_run-*.R' -print0 |
  parallel -0 -j 75% --halt soon,fail=1 \
    'Rscript {} > ../logs/{/.}.log 2>&1'

echo "All analyses completed successfully."
