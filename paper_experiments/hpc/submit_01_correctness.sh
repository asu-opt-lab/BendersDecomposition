#!/bin/bash
#SBATCH -t 0-01:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

# Define variables to make the script more readable and maintainable

ROUND_VERSION="01_correctness"
ROUND_DESCRIPTION="Correctness check, UFLP/CFLP p1-p71, fresh extensive-form baseline"
EXPERIMENT_VERSION="1"
HOUR="01"
TIME_LIMIT="600"
GAP_TOLERANCE="1e-6"
EXPERIMENT_DESCRIPTION="${HOUR} hr, p1-p71, baseline + Benders variants"

FILE_NAME="01_correctness.jl"
SHELL_FILE_NAME="submit_01_correctness.sh"
THREADS=7
MEM="40G"
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
- **Threads**: ${THREADS}
EOF

# Define arrays of problem and instance names
problems=(
    "uflp_correctness"
    "cflp_correctness"
)

instances=()
for i in $(seq 1 71); do
    instances+=("p${i}")
done

# Loop through cases and create one job script for each problem-instance pair
for problem in "${problems[@]}"; do
    for instance in "${instances[@]}"; do
        JOB_NAME="${problem}_${instance}"
        JOBSCRIPT_FILE="${OUTPUT_DIR}/${JOB_NAME}.sh"
        JOB_OUTPUT_DIR="${OUTPUT_DIR}/${JOB_NAME}"

        mkdir -p "${JOB_OUTPUT_DIR}"

        # Create job script file
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

        # Load necessary modules
        echo "module purge" >> "${JOBSCRIPT_FILE}"
        echo "module load julia" >> "${JOBSCRIPT_FILE}"
        echo "module load cplex" >> "${JOBSCRIPT_FILE}"
        echo "module load gurobi" >> "${JOBSCRIPT_FILE}"

        # Run Julia script with algorithm parameters
        echo "cd ${REPO_ROOT}" >> "${JOBSCRIPT_FILE}"
        echo "julia --project=paper_experiments paper_experiments/scripts/${FILE_NAME} --output_dir=${JOB_OUTPUT_DIR} --time_limit=${TIME_LIMIT} --gap_tolerance=${GAP_TOLERANCE} --solver_threads=${THREADS} --problem=${problem} --instance=${instance}" >> "${JOBSCRIPT_FILE}"

        # Submit job
        sbatch "${JOBSCRIPT_FILE}"
        rm "${JOBSCRIPT_FILE}"
    done
done

