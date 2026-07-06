#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$ROOT/paper_experiments/scripts/run_01_correctness.sh"
"$ROOT/paper_experiments/scripts/run_02_oracle_ablation.sh"
"$ROOT/paper_experiments/scripts/run_03_environment_ablation.sh"
"$ROOT/paper_experiments/scripts/run_04_parallel_scflp.sh"

