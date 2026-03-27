#!/usr/bin/env bash
# File: run_local_array_gpu.sh
# Run per-instance jobs in parallel on GPU, auto-assigning each job to the first free GPU.
# Modified: allow MULTIPLE concurrent jobs per GPU via "slots".
set -euo pipefail

# ── EDIT ME: experiment tags ──────────────────────────────────────────────────
PROBLEM_TYPE="MIPLIB_result"
ROUND_VERSION="BD"
EXPERIMENT_VERSION="GPU"
ROUND_DESCRIPTION="260313_warmstart_npresolve_see_divergence"
EXPERIMENT_DESCRIPTION="integrality gap > 10%, CPLEX time >500s, not cut TL, exclude too less NC"

# Upper bound on how many jobs to run at once (across ALL GPUs)
CONCURRENCY="${CONCURRENCY:-8}"

# NEW: how many concurrent jobs you allow PER GPU
# Example: 2 GPUs * SLOTS_PER_GPU=2 => up to 4 simultaneous jobs (subject to CONCURRENCY)
SLOTS_PER_GPU="${SLOTS_PER_GPU:-4}"

# call style for the worker script (flags or positional)
CALL_STYLE="${CALL_STYLE:-flags}"

# ── paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT="${REPO_ROOT:-$PWD}"
SCRIPT_PATH="${REPO_ROOT}/scripts/MIPLIB_typical.sh"

OUTPUT_DIR="${REPO_ROOT}/scripts/${PROBLEM_TYPE}/${ROUND_VERSION}/${EXPERIMENT_VERSION}/${ROUND_DESCRIPTION}"
LOG_DIR="${OUTPUT_DIR}/results"
JOB_DIR="${REPO_ROOT}/scripts/local_jobs_${ROUND_DESCRIPTION}"

# Per-run lock dir so old/stale locks never block this run
LOCK_DIR="/tmp/${USER:-user}_gpu_locks_${ROUND_DESCRIPTION}"
rm -rf "${LOCK_DIR}" && mkdir -p "${LOCK_DIR}"

# fail if experiment dir already exists
if [[ -d "${OUTPUT_DIR}" ]]; then
  echo "Error: Experiment directory ${OUTPUT_DIR} already exists. Use a different ROUND_DESCRIPTION/EXPERIMENT_VERSION or remove the existing directory."
  exit 1
fi

# ── instances (66 - Cut_TL & ETC instances - too less NC = 44) ────────────────────────────────────────────────────────────────
instances=("neos-5188808-nattai" "momentum1" "neos-4321076-ruwer" "supportcase41" "ns1849932" "ns1430538" "neos-2629914-sudost" "neos-4533806-waima" "van" "neos-4555749-wards" "neos-5093327-huahum" "neos-4562542-watut" "zeil" "neos-4408804-prosna" "ns1111636" "neos-4760493-puerua" "neos-4763324-toguru" "rocII-8-11" "roi5alpha10n8" "stockholm" "neos-4355351-swalm" "satellites2-25" "satellites2-40" "neos-5266653-tugela" "neos-5013590-toitoi" "dws008-03" "neos-3025225-shelon" "dws012-01" "dws012-02" "dws012-03" "neos-4300652-rahue" "in" "neos-2991472-kalu" "neos-3209462-rhin" "neos-4409277-trave" "neos-5221106-oparau" "ns1856153" "ns1904248" "ns930473" "snip10x10-35r1budget17" "fhnw-schedule-pairb200" "fhnw-schedule-paira200" "eva1aprime6x6opt" "neos-3695882-vesdre")

# ── prep ──────────────────────────────────────────────────────────────────────
mkdir -p "${LOG_DIR}"
rm -rf "${JOB_DIR}"
mkdir -p "${JOB_DIR}"
chmod +x "${SCRIPT_PATH}"

mkdir -p "${OUTPUT_DIR}"
cat > "${OUTPUT_DIR}/experiment_metadata.md" << EOF
# Experiment Metadata
- Round Version: ${ROUND_VERSION}
- Round Description: ${ROUND_DESCRIPTION}
- Experiment Version: ${EXPERIMENT_VERSION}
- Experiment Description: ${EXPERIMENT_DESCRIPTION}
- Date: $(date "+%Y-%m-%d %H:%M:%S")
- Host: $(hostname)
- Concurrency (upper bound): ${CONCURRENCY}
- Slots per GPU: ${SLOTS_PER_GPU}
EOF

# ── detect GPUs ───────────────────────────────────────────────────────────────
detect_gpu_count() {
  local n=0
  if command -v nvidia-smi >/dev/null 2>&1; then
    n=$(nvidia-smi -L 2>/dev/null | wc -l | awk '{print $1}')
  fi
  if [[ -z "${n}" || "${n}" -eq 0 ]]; then
    n=1   # fallback (change if you prefer)
  fi
  echo "${n}"
}
GPU_COUNT="$(detect_gpu_count)"

# Cap concurrency to total "GPU slots" so you don't launch more wrappers than you can schedule
TOTAL_SLOTS=$(( GPU_COUNT * SLOTS_PER_GPU ))
if (( CONCURRENCY > TOTAL_SLOTS )); then
  CONCURRENCY="${TOTAL_SLOTS}"
fi

echo "Detected GPUs: ${GPU_COUNT}. SLOTS_PER_GPU=${SLOTS_PER_GPU} => TOTAL_SLOTS=${TOTAL_SLOTS}. Running CONCURRENCY=${CONCURRENCY} jobs in parallel."

# ── generate per-instance job scripts ─────────────────────────────────────────
for inst in "${instances[@]}"; do
  job="${JOB_DIR}/${inst}.sh"

  if [[ "${CALL_STYLE}" == "flags" ]]; then
    cat > "${job}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
cd "__REPO_ROOT__"

echo "JOB wrapper PID $$ (instance=__INST__)"
echo "JOB sees CUDA_DEVICE_ORDER=${CUDA_DEVICE_ORDER:-<unset>}"
echo "JOB sees CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>}"

# run the actual worker in foreground
stdbuf -oL -eL "__SCRIPT_PATH__" --instance "__INST__" --output_dir "__OUTPUT_DIR__"
EOF
  else
    cat > "${job}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
cd "__REPO_ROOT__"

echo "JOB wrapper PID $$ (instance=__INST__)"
echo "JOB sees CUDA_DEVICE_ORDER=${CUDA_DEVICE_ORDER:-<unset>}"
echo "JOB sees CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>}"

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

# ── dispatch jobs in parallel (auto-assign any free GPU "slot") ───────────────
run_on_any_gpu() {
  local inst="$1"
  local job="$2"

  echo "[$(date '+%F %T')] QUEUED ${inst} (waiting for free GPU slot; GPUs 0..$((GPU_COUNT-1)), slots 0..$((SLOTS_PER_GPU-1)))"

  while true; do
    for ((gpu=0; gpu<GPU_COUNT; gpu++)); do
      for ((slot=0; slot<SLOTS_PER_GPU; slot++)); do
        local lock_path="${LOCK_DIR}/gpu${gpu}.slot${slot}.lock"

        (
          flock -n 9 || exit 1

          export CUDA_DEVICE_ORDER=PCI_BUS_ID
          export CUDA_VISIBLE_DEVICES="${gpu}"

          echo "[$(date '+%F %T')] START ${inst} (GPU ${gpu}, slot ${slot})"
          command -v nvidia-smi >/dev/null 2>&1 && \
            nvidia-smi --query-gpu=index,name,memory.total,memory.used,utilization.gpu --format=csv,noheader -i "${gpu}" || true

          stdbuf -oL -eL "${job}"
          rc=$?

          echo "[$(date '+%F %T')] END   ${inst} (GPU ${gpu}, slot ${slot}, rc=${rc})"
        ) 9>"${lock_path}" && return 0
      done
    done

    echo "[$(date '+%F %T')] WAIT  ${inst} (all slots busy: ${TOTAL_SLOTS} total)"
    sleep 2
  done
}
export -f run_on_any_gpu
export GPU_COUNT LOCK_DIR SLOTS_PER_GPU TOTAL_SLOTS

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