#!/bin/sh
#SBATCH -t 0-01:00:00

# Define variables to make the script more readable and maintainable

ROUND_VERSION="vertical_reverse_polar/snip_classical"
ROUND_DESCRIPTION="vertical_reverse_polar SNIP, classical separable oracle, freq500, 56 cores SLURM, julia --threads=8, CPLEX"
EXPERIMENT_VERSION="new_true_parallel"
SEED="1"
HOUR="01"
TIME_LIMIT="14400"
DCGLP_TIME_LIMIT="100"
DCGLP_ITER_LIMIT="250"
DCGLP_HALT_LIMIT="3"
FREQUENCY="500"
REUSE_DCGLP="false"
STRENGTHENED="true"
LIFT="false"
EXPERIMENT_DESCRIPTION="${HOUR} hr, seed = ${SEED}, vertical_reverse_polar"

FILE_NAME="snip_classical.jl"
SHELL_FILE_NAME="snip_classical.sh"
THREADS=7
JULIA_THREADS=7

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

# Copy julia and shell script to output directory
cp -r polardcglp/scripts/vertical_reverse_polar/${FILE_NAME} "${OUTPUT_DIR}/${FILE_NAME}"
cp -r polardcglp/scripts/vertical_reverse_polar/${SHELL_FILE_NAME} "${OUTPUT_DIR}/${SHELL_FILE_NAME}"

# Create experiment metadata markdown file
cat > "${OUTPUT_DIR}/experiment_metadata.md" << EOF
# Experiment Metadata

- **Round Version**: ${ROUND_VERSION}
- **Round Description**: ${ROUND_DESCRIPTION}
- **Experiment Version**: ${EXPERIMENT_VERSION}
- **Experiment Description**: ${EXPERIMENT_DESCRIPTION}
- **Date**: $(date "+%Y-%m-%d %H:%M:%S")
EOF

# Define instance parameter arrays
instance_nos=(0 1 2 3 4)
snip_nos=(3 4)
budgets=(30.0 40.0 50.0 60.0 70.0 80.0 90.0)
# budgets=(100.0 110.0 120.0 130.0 140.0 150.0)

# Loop through the instances and create a job script for each combination
for instance_no in "${instance_nos[@]}"; do
    for snip_no in "${snip_nos[@]}"; do
        for budget in "${budgets[@]}"; do
            TAG="inst${instance_no}_snip${snip_no}_bud${budget}"
            JOBSCRIPT_FILE="${OUTPUT_DIR}/${TAG}.sh"

            # Create job script file
            echo "#!/bin/bash" > "${JOBSCRIPT_FILE}"

            echo "#SBATCH -q grp_gbyeon" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -p general" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -N 1" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -n 1" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -c ${THREADS}" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH --nodelist=pcc036,pcc037" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH --mem=100G" >> "${JOBSCRIPT_FILE}"

            echo "#SBATCH -t 0-${HOUR}:30:00" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -o ${ERR_OUT_DIR}/${TAG}.out%j" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -e ${ERR_OUT_DIR}/${TAG}.err%j" >> "${JOBSCRIPT_FILE}"

            # Load necessary modules
            echo "module purge" >> "${JOBSCRIPT_FILE}"
            echo "module load julia" >> "${JOBSCRIPT_FILE}"
            echo "module load cplex" >> "${JOBSCRIPT_FILE}"

            # Run Julia script with algorithm parameters
            echo "julia --threads=${JULIA_THREADS} --project=polardcglp polardcglp/scripts/vertical_reverse_polar/${FILE_NAME} --instance_no=${instance_no} --snip_no=${snip_no} --budget=${budget} --seed=${SEED} --time_limit=${TIME_LIMIT} --dcglp_time_limit=${DCGLP_TIME_LIMIT} --dcglp_iter_limit=${DCGLP_ITER_LIMIT} --dcglp_halt_limit=${DCGLP_HALT_LIMIT} --frequency=${FREQUENCY} --reuse_dcglp=${REUSE_DCGLP} --strengthened=${STRENGTHENED} --lift=${LIFT}" >> "${JOBSCRIPT_FILE}"

            # Submit job
            sbatch "${JOBSCRIPT_FILE}"
            rm "${JOBSCRIPT_FILE}"
        done
    done
done
