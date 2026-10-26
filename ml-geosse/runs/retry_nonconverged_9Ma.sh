#!/usr/bin/env bash

set -euo pipefail

RUNS_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${RUNS_DIR}/scripts"

scripts=(
  "Afr16_gr-4_run-4.R"
  "Oriental16_gr-4_run-2.R"
  "Oriental16_gr-6_run-3.R"
  "OrientalSunken16_gr-6_run-2.R"
)

printf '%s\n' "${scripts[@]}" |
  parallel -j 4 --halt never \
    'Rscript {} > ../logs/{.}.retry.log 2>&1'

echo "Retry batch finished at $(date '+%Y-%m-%d %H:%M:%S')"
