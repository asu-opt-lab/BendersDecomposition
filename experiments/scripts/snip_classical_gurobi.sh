#!/bin/sh
#SBATCH -t 0-01:00:00

# Define variables to make the script more readable and maintainable

ROUND_VERSION="snip_classical_gurobi"
ROUND_DESCRIPTION="snip classical BendersBnB + LazyCallback + BendersSeq root, 56 cores SLURM, julia --threads=8, Gurobi"
EXPERIMENT_VERSION="1"
SEED="1"
HOUR="01"
EXPERIMENT_DESCRIPTION="${HOUR} hr, seed = ${SEED}"

FILE_NAME="snip_classical_gurobi.jl"
THREADS=56
JULIA_THREADS=8

# Define variables to make the script more readable and maintainable
OUTPUT_DIR="experiments/${ROUND_VERSION}/${EXPERIMENT_VERSION}"
ERR_OUT_DIR="${OUTPUT_DIR}/results"

if [ -d "${OUTPUT_DIR}" ]; then
    echo "Error: Experiment directory ${OUTPUT_DIR} already exists. Please use a different EXPERIMENT_VERSION or remove the existing directory."
    exit 1
fi

# Create necessary directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${ERR_OUT_DIR}"

# Copy julia script to output directory
cp -r scripts/${FILE_NAME} "${OUTPUT_DIR}/${FILE_NAME}"

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

# Loop through the instances and create a job script for each combination
for instance_no in "${instance_nos[@]}"; do
    for snip_no in "${snip_nos[@]}"; do
        for budget in "${budgets[@]}"; do
            TAG="inst${instance_no}_snip${snip_no}_bud${budget}"
            JOBSCRIPT_FILE="${OUTPUT_DIR}/${TAG}.sh"

            # Create job script file
            echo "#!/bin/bash" > "${JOBSCRIPT_FILE}"

            echo "#SBATCH -q grp_gbyeon" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -N 1" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -n 1" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -c ${THREADS}" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH --nodelist=pcc036,pcc037" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH --mem=200G" >> "${JOBSCRIPT_FILE}"

            echo "#SBATCH -t 0-${HOUR}:00:00" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -o ${ERR_OUT_DIR}/${TAG}.out%j" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -e ${ERR_OUT_DIR}/${TAG}.err%j" >> "${JOBSCRIPT_FILE}"

            # Load necessary modules
            echo "module purge" >> "${JOBSCRIPT_FILE}"
            echo "module load julia" >> "${JOBSCRIPT_FILE}"
            echo "module load gurobi" >> "${JOBSCRIPT_FILE}"

            # Run Julia script with algorithm parameters
            echo "julia --threads=${JULIA_THREADS} --project=experiments/scripts experiments/scripts/${FILE_NAME} --instance_no ${instance_no} --snip_no ${snip_no} --budget ${budget} --output_dir ${OUTPUT_DIR} --seed ${SEED}" >> "${JOBSCRIPT_FILE}"

            # Submit job
            sbatch "${JOBSCRIPT_FILE}"
            rm "${JOBSCRIPT_FILE}"
        done
    done
done
