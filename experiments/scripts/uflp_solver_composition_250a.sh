#!/bin/bash

set -euo pipefail

EXPERIMENT="experiment1"
INSTANCES=(
    ga250a-1 ga250a-2 ga250a-3 ga250a-4 ga250a-5
    gs250a-1 gs250a-2 gs250a-3 gs250a-4 gs250a-5
)
SOLVER_CONFIGS=(
    "cplex cplex"
    "cplex gurobi"
    "gurobi cplex"
    "gurobi gurobi"
)
SEED=1
THREADS=1
TIME_LIMIT=3600

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_DIR="${REPO_ROOT}/experiments/uflp_solver_composition_250a/${EXPERIMENT}"
ERR_OUT_DIR="${OUTPUT_DIR}/results"

[[ ! -d "${OUTPUT_DIR}" ]] || {
    echo "Error: ${OUTPUT_DIR} exists"
    exit 1
}

mkdir -p "${ERR_OUT_DIR}" "${OUTPUT_DIR}/results_csv"
cp "${SCRIPT_DIR}/uflp_solver_composition.jl" "${OUTPUT_DIR}/"
cp "${SCRIPT_DIR}/solver_composition_factories.jl" "${OUTPUT_DIR}/"
cp "${BASH_SOURCE[0]}" "${OUTPUT_DIR}/$(basename "${BASH_SOURCE[0]}")"

cat > "${OUTPUT_DIR}/experiment_metadata.md" <<EOF
# UFLP Solver Composition Experiment

- Instances: ga250a-1 through ga250a-5; gs250a-1 through gs250a-5
- Solver configurations: CPLEX/CPLEX, CPLEX/Gurobi, Gurobi/CPLEX, Gurobi/Gurobi
- Threads: ${THREADS}
- Seed: ${SEED}
- Time limit: ${TIME_LIMIT} seconds
- Date: $(date "+%Y-%m-%d %H:%M:%S")
EOF

for instance in "${INSTANCES[@]}"; do
    for solver_config in "${SOLVER_CONFIGS[@]}"; do
        read -r master_solver sub_solver <<< "${solver_config}"
        job_name="${instance}_master-${master_solver}_sub-${sub_solver}"
        job_file="${OUTPUT_DIR}/${job_name}.sh"

        cat > "${job_file}" <<EOF
#!/bin/bash
#SBATCH -p general
#SBATCH -q grp_gbyeon
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c ${THREADS}
#SBATCH --mem=60G
#SBATCH -t 0-01:10:00
#SBATCH -J ${job_name}
#SBATCH -o ${ERR_OUT_DIR}/${job_name}.out.%j
#SBATCH -e ${ERR_OUT_DIR}/${job_name}.err.%j

set -euo pipefail

module purge
module load julia
module load cplex
module load gurobi

cd "${REPO_ROOT}"
julia --startup-file=no --project="${REPO_ROOT}/experiments/scripts" \
    "${OUTPUT_DIR}/uflp_solver_composition.jl" \
    --instance "${instance}" \
    --master_solver "${master_solver}" \
    --sub_solver "${sub_solver}" \
    --seed ${SEED} \
    --threads ${THREADS} \
    --time_limit ${TIME_LIMIT} \
    --output_dir "${OUTPUT_DIR}"
EOF

        sbatch "${job_file}"
        rm "${job_file}"
    done
done
