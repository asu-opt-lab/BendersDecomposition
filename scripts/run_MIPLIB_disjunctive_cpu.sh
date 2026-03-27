#!/usr/bin/env bash
# File: run_local_array_cpu.sh
# Run per-instance jobs in parallel on CPU only.
# Intended: 15 jobs in parallel, CPLEX threads controlled inside Julia.
set -euo pipefail

# ── EDIT ME: experiment tags ──────────────────────────────────────────────────
PROBLEM_TYPE="MIPLIB_result"
ROUND_VERSION="DBD"
EXPERIMENT_VERSION="CPU"
ROUND_DESCRIPTION=""
EXPERIMENT_DESCRIPTION=""

# CPLEX thread budget per job (used only for concurrency calculation)
CPUS_PER_JOB="${CPUS_PER_JOB:-7}"

# desired concurrency (how many jobs to run in parallel)
DESIRED_CONCURRENCY="${DESIRED_CONCURRENCY:-15}"

# call style for the worker script (flags or positional)
CALL_STYLE="${CALL_STYLE:-flags}"

# ── paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT="${REPO_ROOT:-$PWD}"
SCRIPT_PATH="${REPO_ROOT}/scripts/MIPLIB_disjunctive.sh"

OUTPUT_DIR="${REPO_ROOT}/scripts/${PROBLEM_TYPE}/${ROUND_VERSION}/${EXPERIMENT_VERSION}/${ROUND_DESCRIPTION}"
LOG_DIR="${OUTPUT_DIR}/results"
JOB_DIR="${REPO_ROOT}/scripts/local_jobs_${ROUND_DESCRIPTION}"

# fail if experiment dir already exists
if [[ -d "${OUTPUT_DIR}" ]]; then
  echo "Error: Experiment directory ${OUTPUT_DIR} already exists. Use a different ROUND_DESCRIPTION/EXPERIMENT_VERSION or remove the existing directory."
  exit 1
fi

# ── integrality gap > 10%, CPLEX time < 500s ────────────────────────────────────────────────────────────────
instances=("neos-5251015-ogosta" "neos-5188808-nattai" "momentum1" "neos-4321076-ruwer" "supportcase41" "ns1849932" "ns1430538" "hgms62" "hgms30" "neos-2629914-sudost" "neos-4533806-waima" "van" "neos-4555749-wards" "neos-5093327-huahum" "neos-4562542-watut" "ns2124243" "rmatr200-p10" "zeil" "neos-4408804-prosna" "ns1111636" "rmatr200-p5" "neos-4760493-puerua" "neos-4763324-toguru" "rocII-8-11" "roi5alpha10n8" "neos-872648" "stockholm" "satellites3-25" "neos-4355351-swalm" "satellites2-25" "satellites2-40" "neos-5266653-tugela" "satellites4-25" "shipsched" "neos-5013590-toitoi" "dws008-03" "neos-3025225-shelon" "dws012-01" "dws012-02" "dws012-03" "fastxgemm-n3r22s4t6" "fastxgemm-n3r21s3t6" "neos-4300652-rahue" "fastxgemm-n3r23s5t6" "in" "neos-2991472-kalu" "neos-3209462-rhin" "neos-4292145-piako" "neos-4409277-trave" "neos-4535459-waipa" "neos-4545615-waita" "neos-5221106-oparau" "ns1856153" "ns1904248" "ns930473" "snip10x10-35r1budget17" "fhnw-schedule-pairb200" "fhnw-schedule-pairb400" "fhnw-schedule-paira200" "fhnw-schedule-paira400" "eva1aprime6x6opt" "neos-3695882-vesdre" "neos-3208254-reiu")


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
- Desired Concurrency: ${DESIRED_CONCURRENCY}
- CPLEX Threads per job (CPXPARAM_Threads): ${CPUS_PER_JOB}
EOF

# ── detect CPU threads ────────────────────────────────────────────────────────
detect_cpu_count() {
  if command -v nproc >/dev/null 2>&1; then
    nproc
  else
    getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1
  fi
}
CPU_COUNT="$(detect_cpu_count)"

MAX_CONCURRENCY_BY_CPU=$(( CPU_COUNT / CPUS_PER_JOB ))
if (( MAX_CONCURRENCY_BY_CPU < 1 )); then
  echo "Error: CPUS_PER_JOB (${CPUS_PER_JOB}) is larger than available CPU threads (${CPU_COUNT})."
  exit 1
fi

if (( DESIRED_CONCURRENCY > MAX_CONCURRENCY_BY_CPU )); then
  CONCURRENCY="${MAX_CONCURRENCY_BY_CPU}"
else
  CONCURRENCY="${DESIRED_CONCURRENCY}"
fi

echo "Detected CPU threads: ${CPU_COUNT}. Running CONCURRENCY=${CONCURRENCY} jobs in parallel."

# ── generate per-instance job scripts ─────────────────────────────────────────
for inst in "${instances[@]}"; do
  job="${JOB_DIR}/${inst}.sh"

  if [[ "${CALL_STYLE}" == "flags" ]]; then
    cat > "${job}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
cd "__REPO_ROOT__"

echo "JOB wrapper PID $$ (instance=__INST__)"

stdbuf -oL -eL "__SCRIPT_PATH__" --instance "__INST__" --output_dir "__OUTPUT_DIR__"
EOF
  else
    cat > "${job}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
cd "__REPO_ROOT__"

echo "JOB wrapper PID $$ (instance=__INST__)"

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

# ── launch ────────────────────────────────────────────────────────────────────
printf "%s\n" "${instances[@]}" | \
  xargs -P "${CONCURRENCY}" -I{} bash -lc '
    set +e
    job="'"${JOB_DIR}"'/{}.sh"
    log="'"${LOG_DIR}"'/{}.log"
    : > "${log}"

    "${job}" |& tee -a "${log}"
    rc=${PIPESTATUS[0]}
    echo "[${rc}] recorded for {}" >> "${log}"
    exit 0
  '

echo
echo "Done. Logs in: ${LOG_DIR}"
