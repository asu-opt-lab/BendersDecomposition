#!/usr/bin/env bash
set -euo pipefail

# Defaults (override via env if needed)
JULIA_BIN="${JULIA_BIN:-$(command -v julia || true)}"
JULIA_PROJECT="${JULIA_PROJECT:---project=.}"

INSTANCE=""
OUTPUT_DIR=""

# Parse flags: --instance X --output_dir OUT
while [[ $# -gt 0 ]]; do
  case "$1" in
    --instance)   INSTANCE="$2"; shift 2 ;;
    --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "$INSTANCE" || -z "$OUTPUT_DIR" ]]; then
  echo "Usage: $0 --instance <name> --output_dir <dir>" >&2
  exit 2
fi

cp -n -p scripts/MIPLIB_typical.jl "$OUTPUT_DIR/"

# kill job if it runs longer
exec timeout --signal=TERM 7200 \
  "$JULIA_BIN" $JULIA_PROJECT \
  scripts/MIPLIB_typical.jl \
  --instance "$INSTANCE"
