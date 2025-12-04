#!/usr/bin/env bash
# File: run_local_array.sh
# Generate per-instance job scripts and execute them in parallel,
# auto-assigning each job to the first free GPU. Each job logs to its own file.

set -euo pipefail

# ── EDIT ME: experiment tags ──────────────────────────────────────────────────
ROUND_VERSION="Disjunctive"
EXPERIMENT_VERSION="GPU"
# EXPERIMENT_VERSION="CPU"
ROUND_DESCRIPTION="1124_knapsack_PDHG_only_typicaloracle"
EXPERIMENT_DESCRIPTION="Crossover(1), Threads(7), LPWarm(0), PDHGGPU(1)"

# How many jobs to run at once (upper bound)  ← ensure single-task execution
CONCURRENCY="${CONCURRENCY:-1}"

# Pin everything to this GPU (minimal change vs. rewriting dispatcher)
PIN_TO_GPU="${PIN_TO_GPU:-0}"

# Does scripts/cflp_typical_knapsack.sh accept flags --instance/--output_dir ?
#   "flags"     →  ./cflp_typical_knapsack.sh --instance X --output_dir OUT
#   "positional"→  OUTPUT_DIR=OUT ./cflp_typical_knapsack.sh X
CALL_STYLE="${CALL_STYLE:-flags}"

# ── paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT="${REPO_ROOT:-$PWD}"
SCRIPT_PATH="${REPO_ROOT}/scripts/cflp_disjunctive_knapsack.sh"

OUTPUT_DIR="${REPO_ROOT}/experiments/${ROUND_VERSION}/${EXPERIMENT_VERSION}/${ROUND_DESCRIPTION}"
LOG_DIR="${OUTPUT_DIR}/results"
JOB_DIR="${REPO_ROOT}/scripts/local_jobs_${ROUND_DESCRIPTION}"

# Per-run lock dir so old/stale locks never block this run
LOCK_DIR="/tmp/${USER:-user}_gpu_locks_${ROUND_DESCRIPTION}"
mkdir -p "${LOCK_DIR}"

# ── hard stop if experiment directory already exists (like sbatch script) ─────
if [[ -d "${OUTPUT_DIR}" ]]; then
  echo "Error: Experiment directory ${OUTPUT_DIR} already exists. Use a different ROUND_DESCRIPTION/EXPERIMENT_VERSION or remove the existing directory."
  exit 1
fi

# ── instances (edit this list, or load from a file) ───────────────────────────
# To load from file: mapfile -t instances < instances.txt
instances=(
    "T1500x600_10_1" #"T1500x600_10_2""T1500x600_10_3" "T1500x600_10_4" "T1500x600_10_5"
    # "T1500x600_15_1" "T1500x600_15_2" "T1500x600_15_3" "T1500x600_15_4" "T1500x600_15_5"
)

# ── prep ──────────────────────────────────────────────────────────────────────
mkdir -p "${LOG_DIR}"
rm -rf "${JOB_DIR}"
mkdir -p "${JOB_DIR}"
chmod +x "${SCRIPT_PATH}"

# metadata (optional)
mkdir -p "${OUTPUT_DIR}"
cat > "${OUTPUT_DIR}/experiment_metadata.md" << EOF
# Experiment Metadata
- **Round Version**: ${ROUND_VERSION}
- **Round Description**: ${ROUND_DESCRIPTION}
- **Experiment Version**: ${EXPERIMENT_VERSION}
- **Experiment Description**: ${EXPERIMENT_DESCRIPTION}
- **Date**: $(date "+%Y-%m-%d %H:%M:%S")
- **Host**: $(hostname)
- **Concurrency**: ${CONCURRENCY}
EOF

# ── detect GPUs (kept for completeness; we’ll still pin to GPU 0) ──────────────
detect_gpu_count() {
  local n=0
  if command -v nvidia-smi >/dev/null 2>&1; then
    n=$(nvidia-smi -L 2>/dev/null | wc -l | awk '{print $1}')
  fi
  if [[ -z "${n}" || "${n}" -eq 0 ]]; then
    n=2   # fallback
  fi
  echo "${n}"
}
GPU_COUNT="$(detect_gpu_count)"

# Cap concurrency to 1 explicitly (no concurrent jobs)
CONCURRENCY=1

# ── generate per-instance job scripts (no START/END lines here) ───────────────
for inst in "${instances[@]}"; do
  job="${JOB_DIR}/${inst}.sh"

  if [[ "${CALL_STYLE}" == "flags" ]]; then
    cat > "${job}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
cd "__REPO_ROOT__"
stdbuf -oL -eL "__SCRIPT_PATH__" --instance "__INST__" --output_dir "__OUTPUT_DIR__"
EOF
  else
    cat > "${job}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
cd "__REPO_ROOT__"
OUTPUT_DIR="__OUTPUT_DIR__" stdbuf -oL -eL "__SCRIPT_PATH__" "__INST__"
EOF
  fi

  sed -i \
    -e "s|__REPO_ROOT__|${REPO_ROOT}|g" \
    -e "s|__SCRIPT_PATH__|${SCRIPT_PATH}|g" \
    -e "s|__OUTPUT_DIR__|${OUTPUT_DIR}|g" \
    -e "s|__INST__|${inst}|g" \
    "${job}"

  chmod +x "${job}"
done

# ── dispatch jobs (pin to GPU 0; no round-robin) ──────────────────────────────
run_on_any_gpu() {
  local inst="$1"
  local job="$2"
  local gpu="${PIN_TO_GPU}"     # ← always GPU 0 unless overridden

  echo "[$(date '+%F %T')] QUEUED ${inst} (will run on GPU ${gpu})"

  # Single-GPU lock (harmless even with CONCURRENCY=1)
  local lock_path="${LOCK_DIR}/gpu${gpu}.lock"
  while true; do
    exec {lk}>"${lock_path}" || true
    if flock -n "${lk}"; then
      (
        export CUDA_VISIBLE_DEVICES="${gpu}"
        echo "[$(date '+%F %T')] START ${inst} (GPU ${gpu})"
        command -v nvidia-smi >/dev/null 2>&1 && \
          nvidia-smi --query-gpu=index,name,memory.total,memory.used,utilization.gpu --format=csv,noheader -i "${gpu}" || true
        stdbuf -oL -eL "${job}"
        rc=$?
        echo "[$(date '+%F %T')] END   ${inst} (GPU ${gpu}, rc=${rc})"
        exit 0
      )
      flock -u "${lk}" || true
      exec {lk}>&- || true
      return 0
    fi
    exec {lk}>&- || true
    echo "[$(date '+%F %T')] WAIT  ${inst} (GPU ${gpu} busy)"
    sleep 2
  done
}
export -f run_on_any_gpu
export GPU_COUNT LOCK_DIR PIN_TO_GPU

# ── launch ────────────────────────────────────────────────────────────────────
printf "%s\n" "${instances[@]}" | \
  xargs -P "${CONCURRENCY}" -I{} bash -lc '
    set +e
    job="'"${JOB_DIR}"'/{}.sh"
    log="'"${LOG_DIR}"'/{}.log"
    : > "${log}"

    run_on_any_gpu "{}" "${job}" |& tee -a "${log}"
    rc=${PIPESTATUS[0]}
    echo "[${rc}] recorded for {}" >> "${log}"
    exit 0
  '

echo
echo "Done. Logs in: ${LOG_DIR}"
