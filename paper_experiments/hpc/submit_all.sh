#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"${SCRIPT_DIR}/submit_01_correctness.sh"
"${SCRIPT_DIR}/submit_02_oracle_ablation.sh"
"${SCRIPT_DIR}/submit_03_environment_ablation.sh"
"${SCRIPT_DIR}/submit_04_parallel_scflp.sh"
"${SCRIPT_DIR}/submit_05_parallel_scflp_lp_relax.sh"
