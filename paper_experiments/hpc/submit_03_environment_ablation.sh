#!/bin/bash
#SBATCH -t 0-01:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

# Define variables to make the script more readable and maintainable

ROUND_VERSION="03_environment_ablation"
ROUND_DESCRIPTION="Environment ablation, UFLP, seq vs seq_inout vs callback"
EXPERIMENT_VERSION="1"
HOUR="02"
TIME_LIMIT="1800"
GAP_TOLERANCE="1e-4"
REPEATS=3
EXPERIMENT_DESCRIPTION="${HOUR} hr, repeats = ${REPEATS}, UFLP seq vs seq_inout vs callback"

FILE_NAME="03_environment_ablation.jl"
SHELL_FILE_NAME="submit_03_environment_ablation.sh"
THREADS=7
SOLVER="gurobi"
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
- **Problem**: UFLP
- **Environments**: seq, seq_inout, callback
- **Solver**: ${SOLVER}
- **Threads**: ${THREADS}
EOF

uflp_instances=(
    "ga250c-1" "ga250c-2" "ga250c-3" "ga250c-4" "ga250c-5"
    "gs250c-1" "gs250c-2" "gs250c-3" "gs250c-4" "gs250c-5"
)

for repeat in $(seq 1 "${REPEATS}"); do
    for instance in "${uflp_instances[@]}"; do
        for env_name in "seq" "seq_inout" "callback"; do
            JOB_NAME="uflp_${instance}_${env_name}_r${repeat}"
            JOBSCRIPT_FILE="${OUTPUT_DIR}/${JOB_NAME}.sh"
            JOB_OUTPUT_DIR="${OUTPUT_DIR}/${JOB_NAME}"
            mkdir -p "${JOB_OUTPUT_DIR}"

            echo "#!/bin/bash" > "${JOBSCRIPT_FILE}"
            echo "#SBATCH -p ${PARTITION}" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -q ${QOS}" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -N 1" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -n 1" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -c ${THREADS}" >> "${JOBSCRIPT_FILE}"
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
            echo "julia --project=paper_experiments paper_experiments/scripts/${FILE_NAME} --output_dir=${JOB_OUTPUT_DIR} --time_limit=${TIME_LIMIT} --gap_tolerance=${GAP_TOLERANCE} --solver_threads=${THREADS} --solver=${SOLVER} --repeats=${REPEATS} --repeat_index=${repeat} --problem=uflp --instance=${instance} --env=${env_name}" >> "${JOBSCRIPT_FILE}"
            sbatch "${JOBSCRIPT_FILE}"
            rm "${JOBSCRIPT_FILE}"
        done
    done
done
