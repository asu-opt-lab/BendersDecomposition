#!/bin/bash
#SBATCH -t 0-01:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

# Define variables to make the script more readable and maintainable

ROUND_VERSION="04_parallel_scflp"
ROUND_DESCRIPTION="Parallel SCFLP, generated 100x200 with original stochastic generator, separable oracle"
EXPERIMENT_VERSION="1"
HOUR="02"
TIME_LIMIT="1800"
GAP_TOLERANCE="1e-4"
REPEATS=3
ORACLE="cfl_knapsack"
ENV_NAME="callback"
EXPERIMENT_DESCRIPTION="${HOUR} hr, repeats = ${REPEATS}, env = ${ENV_NAME}, oracle = ${ORACLE}"

FILE_NAME="04_parallel_scflp.jl"
SHELL_FILE_NAME="submit_04_parallel_scflp.sh"
SOLVER_THREADS=1
MEM="100G"
PARTITION="general"
QOS="grp_gbyeon"

# Define variables to make the script more readable and maintainable
OUTPUT_DIR="paper_experiments/results/raw/${ROUND_VERSION}/${EXPERIMENT_VERSION}"
ERR_OUT_DIR="${OUTPUT_DIR}/results"

if [ -d "${OUTPUT_DIR}" ]; then
    echo "Error: Experiment directory ${OUTPUT_DIR} already exists. Please use a different EXPERIMENT_VERSION or remove the existing directory."
    exit 1
fi

# Create necessary directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${ERR_OUT_DIR}"

# Copy experiment scripts to output directory
cp "paper_experiments/scripts/${FILE_NAME}" "${OUTPUT_DIR}/${FILE_NAME}"
cp "paper_experiments/scripts/common.jl" "${OUTPUT_DIR}/common.jl"
cp "paper_experiments/hpc/${SHELL_FILE_NAME}" "${OUTPUT_DIR}/${SHELL_FILE_NAME}"

# Create experiment metadata markdown file
cat > "${OUTPUT_DIR}/experiment_metadata.md" << EOF
# Experiment Metadata

- **Round Version**: ${ROUND_VERSION}
- **Round Description**: ${ROUND_DESCRIPTION}
- **Experiment Version**: ${EXPERIMENT_VERSION}
- **Experiment Description**: ${EXPERIMENT_DESCRIPTION}
- **Date**: $(date "+%Y-%m-%d %H:%M:%S")
- **Time Limit**: ${TIME_LIMIT}
- **Benders Gap Tolerance**: ${GAP_TOLERANCE}
- **Repeats**: ${REPEATS}
- **Environment**: ${ENV_NAME}
- **Oracle**: ${ORACLE}
- **Solver Threads**: ${SOLVER_THREADS}
EOF

scflp_instances=(
    "f100-c200-s256-r5-1" "f100-c200-s256-r5-2" "f100-c200-s256-r5-3" "f100-c200-s256-r5-4" "f100-c200-s256-r5-5"
    "f100-c200-s512-r5-1" "f100-c200-s512-r5-2" "f100-c200-s512-r5-3" "f100-c200-s512-r5-4" "f100-c200-s512-r5-5"
    "f100-c200-s1024-r5-1" "f100-c200-s1024-r5-2" "f100-c200-s1024-r5-3" "f100-c200-s1024-r5-4" "f100-c200-s1024-r5-5"
)

julia_threads_list=(
    "1" "2" "4" "8" "16"
)

for repeat in $(seq 1 "${REPEATS}"); do
    for instance in "${scflp_instances[@]}"; do
        for JULIA_THREADS in "${julia_threads_list[@]}"; do
            JOB_NAME="scflp_${instance}_${ENV_NAME}_${ORACLE}_jt${JULIA_THREADS}_r${repeat}"
            JOBSCRIPT_FILE="${OUTPUT_DIR}/${JOB_NAME}.sh"
            JOB_OUTPUT_DIR="${OUTPUT_DIR}/${JOB_NAME}"
            mkdir -p "${JOB_OUTPUT_DIR}"

            echo "#!/bin/bash" > "${JOBSCRIPT_FILE}"
            echo "#SBATCH -p ${PARTITION}" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -q ${QOS}" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -N 1" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -n 1" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -c 56" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH --nodelist=pcc036,pcc037" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH --mem=${MEM}" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -t 0-${HOUR}:00:00" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -o ${ERR_OUT_DIR}/${JOB_NAME}.out%j" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -e ${ERR_OUT_DIR}/${JOB_NAME}.err%j" >> "${JOBSCRIPT_FILE}"
            echo "module purge" >> "${JOBSCRIPT_FILE}"
            echo "module load julia" >> "${JOBSCRIPT_FILE}"
            echo "module load cplex" >> "${JOBSCRIPT_FILE}"
            echo "module load gurobi" >> "${JOBSCRIPT_FILE}"
            echo "cd ${REPO_ROOT}" >> "${JOBSCRIPT_FILE}"
            echo "julia --threads=${JULIA_THREADS} --project=paper_experiments paper_experiments/scripts/${FILE_NAME} --output_dir=${JOB_OUTPUT_DIR} --time_limit=${TIME_LIMIT} --gap_tolerance=${GAP_TOLERANCE} --solver_threads=${SOLVER_THREADS} --repeats=${REPEATS} --repeat_index=${repeat} --instance=${instance} --oracle=${ORACLE} --env=${ENV_NAME}" >> "${JOBSCRIPT_FILE}"
            sbatch "${JOBSCRIPT_FILE}"
            rm "${JOBSCRIPT_FILE}"
        done
    done
done
