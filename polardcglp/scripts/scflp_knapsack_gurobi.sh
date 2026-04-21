#!/bin/sh
#SBATCH -t 0-01:00:00

# Define variables to make the script more readable and maintainable

ROUND_VERSION="scflp_knapsack_gurobi"
ROUND_DESCRIPTION="SCFLP Benders BnB, knapsack oracle, 7 private cores, Gurobi"
EXPERIMENT_VERSION="1"
SEED="1"
HOUR="04"
EXPERIMENT_DESCRIPTION="${HOUR} hr, seed = ${SEED}"

FILE_NAME="scflp_knapsack_gurobi.jl"
THREADS=7

# Define variables to make the script more readable and maintainable
OUTPUT_DIR="polardcglp/experiment/${ROUND_VERSION}/${EXPERIMENT_VERSION}"
ERR_OUT_DIR="${OUTPUT_DIR}/results"

if [ -d "${OUTPUT_DIR}" ]; then
    echo "Error: Experiment directory ${OUTPUT_DIR} already exists. Please use a different EXPERIMENT_VERSION or remove the existing directory."
    exit 1
fi

# Create necessary directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${ERR_OUT_DIR}"

# Copy src directory to output directory
cp -r polardcglp/scripts/${FILE_NAME} "${OUTPUT_DIR}/${FILE_NAME}"

# Create experiment metadata markdown file
cat > "${OUTPUT_DIR}/experiment_metadata.md" << EOF
# Experiment Metadata

- **Round Version**: ${ROUND_VERSION}
- **Round Description**: ${ROUND_DESCRIPTION}
- **Experiment Version**: ${EXPERIMENT_VERSION}
- **Experiment Description**: ${EXPERIMENT_DESCRIPTION}
- **Date**: $(date "+%Y-%m-%d %H:%M:%S")
EOF

# Define an array of instance names
instances=(

    # f10-c10 (tiny)
    # "f10-c10-s5-r3-1"
    # "f10-c10-s25-r3-1"

    # f25-c50
    "f25-c50-s64-r10-1" "f25-c50-s64-r10-2" "f25-c50-s64-r10-3" "f25-c50-s64-r10-4" "f25-c50-s64-r10-5"
    # "f25-c50-s64-r3-1" "f25-c50-s64-r3-2" "f25-c50-s64-r3-3" "f25-c50-s64-r3-4" "f25-c50-s64-r3-5"
    # "f25-c50-s64-r5-1" "f25-c50-s64-r5-2" "f25-c50-s64-r5-3" "f25-c50-s64-r5-4" "f25-c50-s64-r5-5"

    # "f25-c50-s128-r3-1" "f25-c50-s128-r3-2" "f25-c50-s128-r3-3" "f25-c50-s128-r3-4" "f25-c50-s128-r3-5"
    # "f25-c50-s128-r5-1" "f25-c50-s128-r5-2" "f25-c50-s128-r5-3" "f25-c50-s128-r5-4" "f25-c50-s128-r5-5"
    # "f25-c50-s128-r10-1" "f25-c50-s128-r10-2" "f25-c50-s128-r10-3" "f25-c50-s128-r10-4" "f25-c50-s128-r10-5"

    # "f25-c50-s254-r3-1" "f25-c50-s254-r3-2" "f25-c50-s254-r3-3" "f25-c50-s254-r3-4" "f25-c50-s254-r3-5"
    # "f25-c50-s254-r5-1" "f25-c50-s254-r5-2" "f25-c50-s254-r5-3" "f25-c50-s254-r5-4" "f25-c50-s254-r5-5"
    # "f25-c50-s254-r10-1" "f25-c50-s254-r10-2" "f25-c50-s254-r10-3" "f25-c50-s254-r10-4" "f25-c50-s254-r10-5"

    # "f25-c50-s512-r3-1" "f25-c50-s512-r3-2" "f25-c50-s512-r3-3" "f25-c50-s512-r3-4" "f25-c50-s512-r3-5"
    # "f25-c50-s512-r5-1" "f25-c50-s512-r5-2" "f25-c50-s512-r5-3" "f25-c50-s512-r5-4" "f25-c50-s512-r5-5"
    # "f25-c50-s512-r10-1" "f25-c50-s512-r10-2" "f25-c50-s512-r10-3" "f25-c50-s512-r10-4" "f25-c50-s512-r10-5"

    # "f25-c50-s1024-r3-1" "f25-c50-s1024-r3-2" "f25-c50-s1024-r3-3" "f25-c50-s1024-r3-4" "f25-c50-s1024-r3-5"
    # "f25-c50-s1024-r5-1" "f25-c50-s1024-r5-2" "f25-c50-s1024-r5-3" "f25-c50-s1024-r5-4" "f25-c50-s1024-r5-5"
    # "f25-c50-s1024-r10-1" "f25-c50-s1024-r10-2" "f25-c50-s1024-r10-3" "f25-c50-s1024-r10-4" "f25-c50-s1024-r10-5"

    # f50-c50
    # "f50-c50-s256-r3-1" "f50-c50-s256-r3-2" "f50-c50-s256-r3-3" "f50-c50-s256-r3-4" "f50-c50-s256-r3-5"
    # "f50-c50-s256-r5-1" "f50-c50-s256-r5-2" "f50-c50-s256-r5-3" "f50-c50-s256-r5-4" "f50-c50-s256-r5-5"
    # "f50-c50-s256-r10-1" "f50-c50-s256-r10-2" "f50-c50-s256-r10-3" "f50-c50-s256-r10-4" "f50-c50-s256-r10-5"

    # "f50-c50-s512-r3-1" "f50-c50-s512-r3-2" "f50-c50-s512-r3-3" "f50-c50-s512-r3-4" "f50-c50-s512-r3-5"
    # "f50-c50-s512-r5-1" "f50-c50-s512-r5-2" "f50-c50-s512-r5-3" "f50-c50-s512-r5-4" "f50-c50-s512-r5-5"
    # "f50-c50-s512-r10-1" "f50-c50-s512-r10-2" "f50-c50-s512-r10-3" "f50-c50-s512-r10-4" "f50-c50-s512-r10-5"

    # f50-c100
    # "f50-c100-s256-r3-1" "f50-c100-s256-r3-2" "f50-c100-s256-r3-3" "f50-c100-s256-r3-4" "f50-c100-s256-r3-5"
    # "f50-c100-s256-r5-1" "f50-c100-s256-r5-2" "f50-c100-s256-r5-3" "f50-c100-s256-r5-4" "f50-c100-s256-r5-5"
    # "f50-c100-s256-r10-1" "f50-c100-s256-r10-2" "f50-c100-s256-r10-3" "f50-c100-s256-r10-4" "f50-c100-s256-r10-5"

    # "f50-c100-s512-r3-1" "f50-c100-s512-r3-2" "f50-c100-s512-r3-3" "f50-c100-s512-r3-4" "f50-c100-s512-r3-5"
    # "f50-c100-s512-r5-1" "f50-c100-s512-r5-2" "f50-c100-s512-r5-3" "f50-c100-s512-r5-4" "f50-c100-s512-r5-5"
    # "f50-c100-s512-r10-1" "f50-c100-s512-r10-2" "f50-c100-s512-r10-3" "f50-c100-s512-r10-4" "f50-c100-s512-r10-5"

    # f100-c100
    # "f100-c100-s256-r3-1" "f100-c100-s256-r3-2" "f100-c100-s256-r3-3" "f100-c100-s256-r3-4" "f100-c100-s256-r3-5"
    # "f100-c100-s256-r5-1" "f100-c100-s256-r5-2" "f100-c100-s256-r5-3" "f100-c100-s256-r5-4" "f100-c100-s256-r5-5"
    # "f100-c100-s256-r10-1" "f100-c100-s256-r10-2" "f100-c100-s256-r10-3" "f100-c100-s256-r10-4" "f100-c100-s256-r10-5"

    # "f100-c100-s512-r3-1" "f100-c100-s512-r3-2" "f100-c100-s512-r3-3" "f100-c100-s512-r3-4" "f100-c100-s512-r3-5"
    # "f100-c100-s512-r5-1" "f100-c100-s512-r5-2" "f100-c100-s512-r5-3" "f100-c100-s512-r5-4" "f100-c100-s512-r5-5"
    # "f100-c100-s512-r10-1" "f100-c100-s512-r10-2" "f100-c100-s512-r10-3" "f100-c100-s512-r10-4" "f100-c100-s512-r10-5"

)

# Loop through the instances and create a job script for each
for instance in "${instances[@]}"; do
    JOBSCRIPT_FILE="${OUTPUT_DIR}/${instance}.sh"

    # Create job script file
    echo "#!/bin/bash" > "${JOBSCRIPT_FILE}"

    echo "#SBATCH -p htc" >> "${JOBSCRIPT_FILE}"
    # echo "#SBATCH -q grp_gbyeon" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -N 1" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -n 1" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -c ${THREADS}" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH --nodelist=pcc036,pcc037" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH --mem=60G" >> "${JOBSCRIPT_FILE}"

    echo "#SBATCH -t 0-${HOUR}:00:00" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -o ${ERR_OUT_DIR}/${instance}.out%j" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -e ${ERR_OUT_DIR}/${instance}.err%j" >> "${JOBSCRIPT_FILE}"

    # Load necessary modules
    echo "module purge" >> "${JOBSCRIPT_FILE}"
    echo "module load julia" >> "${JOBSCRIPT_FILE}"
    echo "module load gurobi" >> "${JOBSCRIPT_FILE}"

    # Run Julia script with algorithm parameters
    echo "julia --project=polardcglp polardcglp/scripts/${FILE_NAME} --instance=${instance} --output_dir=${OUTPUT_DIR} --seed=${SEED}" >> "${JOBSCRIPT_FILE}"

    # Submit job
    sbatch "${JOBSCRIPT_FILE}"
    rm "${JOBSCRIPT_FILE}"
done
