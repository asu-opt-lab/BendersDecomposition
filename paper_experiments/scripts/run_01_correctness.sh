#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
JULIA_BIN="${JULIA_BIN:-julia}"
PROJECT="$ROOT/paper_experiments"
OUTPUT_DIR="${OUTPUT_DIR:-$PROJECT/results/raw}"
TIME_LIMIT="${TIME_LIMIT:-600}"
GAP_TOLERANCE="${GAP_TOLERANCE:-1e-6}"
SOLVER_THREADS="${SOLVER_THREADS:-1}"
SOLVER="${SOLVER:-gurobi}"

"$JULIA_BIN" --project="$PROJECT" "$PROJECT/scripts/01_correctness.jl" \
  --output_dir "$OUTPUT_DIR" \
  --time_limit "$TIME_LIMIT" \
  --gap_tolerance "$GAP_TOLERANCE" \
  --solver_threads "$SOLVER_THREADS" \
  --solver "$SOLVER"
