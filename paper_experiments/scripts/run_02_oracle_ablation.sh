#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
JULIA_BIN="${JULIA_BIN:-julia}"
PROJECT="$ROOT/paper_experiments"
OUTPUT_DIR="${OUTPUT_DIR:-$PROJECT/results/raw}"
TIME_LIMIT="${TIME_LIMIT:-1800}"
GAP_TOLERANCE="${GAP_TOLERANCE:-1e-4}"
SOLVER_THREADS="${SOLVER_THREADS:-1}"
REPEATS="${REPEATS:-3}"

"$JULIA_BIN" --project="$PROJECT" "$PROJECT/scripts/02_oracle_ablation.jl" \
  --output_dir "$OUTPUT_DIR" \
  --time_limit "$TIME_LIMIT" \
  --gap_tolerance "$GAP_TOLERANCE" \
  --solver_threads "$SOLVER_THREADS" \
  --repeats "$REPEATS"

