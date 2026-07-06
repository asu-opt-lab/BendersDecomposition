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
ORACLE="${ORACLE:-cfl_knapsack}"
ENV_NAME="${ENV_NAME:-callback}"
THREADS_LIST="${THREADS_LIST:-1 2 4 8}"

for JULIA_THREADS in $THREADS_LIST; do
  "$JULIA_BIN" --threads="$JULIA_THREADS" --project="$PROJECT" "$PROJECT/scripts/04_parallel_scflp.jl" \
    --output_dir "$OUTPUT_DIR" \
    --time_limit "$TIME_LIMIT" \
    --gap_tolerance "$GAP_TOLERANCE" \
    --solver_threads "$SOLVER_THREADS" \
    --repeats "$REPEATS" \
    --oracle "$ORACLE" \
    --env "$ENV_NAME"
done
