#!/bin/sh

set -eu

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT_DIR"

julia --project=. polardcglp/scripts/cflp_milp_fix_scan.jl \
  --instance=T100x100_5_2 \
  --time_limit=60 \
  --threads=7 \
  --values=both
