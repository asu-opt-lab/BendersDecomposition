#!/bin/sh
#SBATCH -t 0-01:00:00

# Define variables to make the script more readable and maintainable

ROUND_VERSION="milp/cflp_gurobi"
ROUND_DESCRIPTION="CFLP MILP, 7 private cores, Gurobi"
EXPERIMENT_VERSION="1"
SEED="1"
HOUR="04"
EXPERIMENT_DESCRIPTION="${HOUR} hr, seed = ${SEED}"

FILE_NAME="cflp_gurobi.jl"
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
cp -r polardcglp/scripts/milp/${FILE_NAME} "${OUTPUT_DIR}/${FILE_NAME}"

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

    # "T100x100_3_1" "T100x100_3_2" "T100x100_3_3" "T100x100_3_4" "T100x100_3_5"
    "T100x100_5_1" "T100x100_5_2" "T100x100_5_3" "T100x100_5_4" "T100x100_5_5"
    # "T100x100_10_1" "T100x100_10_2" "T100x100_10_3" "T100x100_10_4" "T100x100_10_5"

    # "T200x100_3_1" "T200x100_3_2" "T200x100_3_3" "T200x100_3_4" "T200x100_3_5"
    # "T200x100_5_1" "T200x100_5_2" "T200x100_5_3" "T200x100_5_4" "T200x100_5_5"
    # "T200x100_10_1" "T200x100_10_2" "T200x100_10_3" "T200x100_10_4" "T200x100_10_5"

    # "T200x200_3_1" "T200x200_3_2" "T200x200_3_3" "T200x200_3_4" "T200x200_3_5"
    # "T200x200_5_1" "T200x200_5_2" "T200x200_5_3" "T200x200_5_4" "T200x200_5_5"
    # "T200x200_10_1" "T200x200_10_2" "T200x200_10_3" "T200x200_10_4" "T200x200_10_5"

    # "T300x300_5_1" "T300x300_5_2" "T300x300_5_3" "T300x300_5_4" "T300x300_5_5"
    # "T300x300_10_1" "T300x300_10_2" "T300x300_10_3" "T300x300_10_4" "T300x300_10_5"
    # "T300x300_15_1" "T300x300_15_2" "T300x300_15_3" "T300x300_15_4" "T300x300_15_5"
    # "T300x300_20_1" "T300x300_20_2" "T300x300_20_3" "T300x300_20_4" "T300x300_20_5"

    # "T500x100_3_1" "T500x100_3_2" "T500x100_3_3" "T500x100_3_4" "T500x100_3_5"
    # "T500x100_5_1" "T500x100_5_2" "T500x100_5_3" "T500x100_5_4" "T500x100_5_5"
    # "T500x100_10_1" "T500x100_10_2" "T500x100_10_3" "T500x100_10_4" "T500x100_10_5"

    # "T500x200_3_1" "T500x200_3_2" "T500x200_3_3" "T500x200_3_4" "T500x200_3_5"
    # "T500x200_5_1" "T500x200_5_2" "T500x200_5_3" "T500x200_5_4" "T500x200_5_5"
    # "T500x200_10_1" "T500x200_10_2" "T500x200_10_3" "T500x200_10_4" "T500x200_10_5"

    # "T500x500_5_1" "T500x500_5_2" "T500x500_5_3" "T500x500_5_4" "T500x500_5_5"
    # "T500x500_10_1" "T500x500_10_2" "T500x500_10_3" "T500x500_10_4" "T500x500_10_5"
    # "T500x500_15_1" "T500x500_15_2" "T500x500_15_3" "T500x500_15_4" "T500x500_15_5"
    # "T500x500_20_1" "T500x500_20_2" "T500x500_20_3" "T500x500_20_4" "T500x500_20_5"

    # "T700x700_5_1" "T700x700_5_2" "T700x700_5_3" "T700x700_5_4" "T700x700_5_5"
    # "T700x700_10_1" "T700x700_10_2" "T700x700_10_3" "T700x700_10_4" "T700x700_10_5"
    # "T700x700_15_1" "T700x700_15_2" "T700x700_15_3" "T700x700_15_4" "T700x700_15_5"
    # "T700x700_20_1" "T700x700_20_2" "T700x700_20_3" "T700x700_20_4" "T700x700_20_5"

    # "T1000x1000_5_1" "T1000x1000_5_2" "T1000x1000_5_3" "T1000x1000_5_4" "T1000x1000_5_5"
    # "T1000x1000_10_1" "T1000x1000_10_2" "T1000x1000_10_3" "T1000x1000_10_4" "T1000x1000_10_5"
    # "T1000x1000_15_1" "T1000x1000_15_2" "T1000x1000_15_3" "T1000x1000_15_4" "T1000x1000_15_5"
    # "T1000x1000_20_1" "T1000x1000_20_2" "T1000x1000_20_3" "T1000x1000_20_4" "T1000x1000_20_5"

    # "T1500x300_5_1" "T1500x300_5_2" "T1500x300_5_3" "T1500x300_5_4" "T1500x300_5_5"
    # "T1500x300_10_1" "T1500x300_10_2" "T1500x300_10_3" "T1500x300_10_4" "T1500x300_10_5"
    # "T1500x300_15_1" "T1500x300_15_2" "T1500x300_15_3" "T1500x300_15_4" "T1500x300_15_5"
    # "T1500x300_20_1" "T1500x300_20_2" "T1500x300_20_3" "T1500x300_20_4" "T1500x300_20_5"

    # "T1500x600_5_1" "T1500x600_5_2" "T1500x600_5_3" "T1500x600_5_4" "T1500x600_5_5"
    # "T1500x600_10_1" "T1500x600_10_2" "T1500x600_10_3" "T1500x600_10_4" "T1500x600_10_5"
    # "T1500x600_15_1" "T1500x600_15_2" "T1500x600_15_3" "T1500x600_15_4" "T1500x600_15_5"
    # "T1500x600_20_1" "T1500x600_20_2" "T1500x600_20_3" "T1500x600_20_4" "T1500x600_20_5"
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
    echo "julia --project=polardcglp polardcglp/scripts/milp/${FILE_NAME} --instance=${instance} --output_dir=${OUTPUT_DIR} --seed=${SEED}" >> "${JOBSCRIPT_FILE}"

    # Submit job
    sbatch "${JOBSCRIPT_FILE}"
    rm "${JOBSCRIPT_FILE}"
done
