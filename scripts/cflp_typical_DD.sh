#!/bin/sh
#SBATCH -t 0-01:00:00

ROUND_VERSION="Dual_Decomposition"
ROUND_DESCRIPTION="DD_w_strict_params"
EXPERIMENT_VERSION="1006_DD_strict_typicalKnapsack_isapprox_50threads_T1000100051_T1500300102"
EXPERIMENT_DESCRIPTION=" strict params are applied, noViolation(DD) runs knapsack, check isapprox between dual value and t value, 50 threads"

# Define variables to make the script more readable and maintainable
OUTPUT_DIR="experiments/${ROUND_VERSION}/${EXPERIMENT_VERSION}"
ERR_OUT_DIR="${OUTPUT_DIR}/results"

# Check if experiment directory already exists
if [ -d "${OUTPUT_DIR}" ]; then
    echo "Error: Experiment directory ${OUTPUT_DIR} already exists. Please use a different EXPERIMENT_VERSION or remove the existing directory."
    exit 1
fi

# Create necessary directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${ERR_OUT_DIR}"

# Define job script directory
JOBSCRIPT_DIR="scripts"
# mkdir -p "${JOBSCRIPT_DIR}"

# Copy src directory to output directory
cp -r scripts/cflp_typical_DD.jl "${OUTPUT_DIR}/cflp_typical_DD.jl"

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
    "T1000x1000_5_1" "T1500x300_10_2"
    # "T100x100_3_1" "T100x100_3_2" "T100x100_3_3" "T100x100_3_4" "T100x100_3_5"
    # "T100x100_5_1" "T100x100_5_2" "T100x100_5_3" "T100x100_5_4" "T100x100_5_5"
    # "T100x100_10_1" "T100x100_10_2" "T100x100_10_3" "T100x100_10_4" "T100x100_10_5"

    # "T200x100_3_1" "T200x100_3_2" "T200x100_3_3" "T200x100_3_4" "T200x100_3_5"
    # "T200x100_5_1" "T200x100_5_2" "T200x100_5_3" "T200x100_5_4" "T200x100_5_5"
    # "T200x100_10_1" "T200x100_10_2" "T200x100_10_3" "T200x100_10_4" "T200x100_10_5"

    # "T500x100_3_1" "T500x100_3_2" "T500x100_3_3" "T500x100_3_4" "T500x100_3_5"
    # "T500x100_5_1" "T500x100_5_2" "T500x100_5_3" "T500x100_5_4" "T500x100_5_5"
    # "T500x100_10_1" "T500x100_10_2" "T500x100_10_3" "T500x100_10_4" "T500x100_10_5"

    # "T200x200_3_1" "T200x200_3_2" "T200x200_3_3" "T200x200_3_4" "T200x200_3_5"
    # "T200x200_5_1" "T200x200_5_2" "T200x200_5_3" "T200x200_5_4" "T200x200_5_5"
    # "T200x200_10_1" "T200x200_10_2" "T200x200_10_3" "T200x200_10_4" "T200x200_10_5"

    # "T500x200_3_1" "T500x200_3_2" "T500x200_3_3" "T500x200_3_4" "T500x200_3_5"
    # "T500x200_5_1" "T500x200_5_2" "T500x200_5_3" "T500x200_5_4" "T500x200_5_5"
    # "T500x200_10_1" "T500x200_10_2" "T500x200_10_3" "T500x200_10_4" "T500x200_10_5"

    # "T300x300_5_1" "T300x300_5_2" "T300x300_5_3" "T300x300_5_4" "T300x300_5_5"
    # "T300x300_10_1" "T300x300_10_2" "T300x300_10_3" "T300x300_10_4" "T300x300_10_5"
    # "T300x300_15_1" "T300x300_15_2" "T300x300_15_3" "T300x300_15_4" "T300x300_15_5"
    # "T300x300_20_1" "T300x300_20_2" "T300x300_20_3" "T300x300_20_4" "T300x300_20_5"

    # "T1500x300_5_1" "T1500x300_5_2" "T1500x300_5_3" "T1500x300_5_4" "T1500x300_5_5"
    # "T1500x300_10_1" "T1500x300_10_2" "T1500x300_10_3" "T1500x300_10_4" "T1500x300_10_5"
    # "T1500x300_15_1" "T1500x300_15_2" "T1500x300_15_3" "T1500x300_15_4" "T1500x300_15_5"
    # "T1500x300_20_1" "T1500x300_20_2" "T1500x300_20_3" "T1500x300_20_4" "T1500x300_20_5"

    # "T500x500_5_1" "T500x500_5_2" "T500x500_5_3" "T500x500_5_4" "T500x500_5_5"
    # "T500x500_10_1" "T500x500_10_2" "T500x500_10_3" "T500x500_10_4" "T500x500_10_5"
    # "T500x500_15_1" "T500x500_15_2" "T500x500_15_3" "T500x500_15_4" "T500x500_15_5"
    # "T500x500_20_1" "T500x500_20_2" "T500x500_20_3" "T500x500_20_4" "T500x500_20_5"

    # "T1500x600_5_1" "T1500x600_5_2" "T1500x600_5_3" "T1500x600_5_4" "T1500x600_5_5"
    # "T1500x600_10_1" "T1500x600_10_2" "T1500x600_10_3" "T1500x600_10_4" "T1500x600_10_5"
    # "T1500x600_15_1" "T1500x600_15_2" "T1500x600_15_3" "T1500x600_15_4" "T1500x600_15_5"
    # "T1500x600_20_1" "T1500x600_20_2" "T1500x600_20_3" "T1500x600_20_4" "T1500x600_20_5"

    # "T700x700_5_1" "T700x700_5_2" "T700x700_5_3" "T700x700_5_4" "T700x700_5_5"
    # "T700x700_10_1" "T700x700_10_2" "T700x700_10_3" "T700x700_10_4" "T700x700_10_5"
    # "T700x700_15_1" "T700x700_15_2" "T700x700_15_3" "T700x700_15_4" "T700x700_15_5"
    # "T700x700_20_1" "T700x700_20_2" "T700x700_20_3" "T700x700_20_4" "T700x700_20_5"

    # "T1000x1000_5_1" "T1000x1000_5_2" "T1000x1000_5_3" "T1000x1000_5_4" "T1000x1000_5_5"
    # "T1000x1000_10_1" "T1000x1000_10_2" "T1000x1000_10_3" "T1000x1000_10_4" "T1000x1000_10_5"
    # "T1000x1000_15_1" "T1000x1000_15_2" "T1000x1000_15_3" "T1000x1000_15_4" "T1000x1000_15_5"
    # "T1000x1000_20_1" "T1000x1000_20_2" "T1000x1000_20_3" "T1000x1000_20_4" "T1000x1000_20_5"
)

# Loop through the instances and create a job script for each
for instance in "${instances[@]}"; do
    JOBSCRIPT_FILE="${JOBSCRIPT_DIR}/${instance}.sh"
    
    # Create job script file
    echo "#!/bin/bash" > "${JOBSCRIPT_FILE}"

    echo "#SBATCH -p general" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -q grp_gbyeon" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -N 1" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -n 1" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -c 50" >> "${JOBSCRIPT_FILE}"
    # echo "#SBATCH --nodelist=pcc036" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH --mem=60G" >> "${JOBSCRIPT_FILE}"

    echo "#SBATCH -t 0-01:30:00" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -o ${ERR_OUT_DIR}/${instance}.out%j" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -e ${ERR_OUT_DIR}/${instance}.err%j" >> "${JOBSCRIPT_FILE}"

    # Load necessary modules
    echo "module purge" >> "${JOBSCRIPT_FILE}"
    echo "module load julia" >> "${JOBSCRIPT_FILE}"
    echo "module load cplex" >> "${JOBSCRIPT_FILE}"
    echo "module load gurobi" >> "${JOBSCRIPT_FILE}"

    # Run Julia script with algorithm parameters
    echo "JULIA_NUM_THREADS=50 julia --project=. scripts/cflp_typical_DD.jl --instance ${instance} --output_dir ${OUTPUT_DIR}" >> "${JOBSCRIPT_FILE}"

    # Submit job
    sbatch "${JOBSCRIPT_FILE}"
    rm "${JOBSCRIPT_FILE}"
done