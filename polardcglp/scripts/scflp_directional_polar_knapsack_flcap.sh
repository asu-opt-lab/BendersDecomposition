#!/bin/sh
#SBATCH -t 0-01:00:00

# Define variables to make the script more readable and maintainable

ROUND_VERSION="scflp_directional_polar_knapsack_flcap"
ROUND_DESCRIPTION="DirectionalPolarDCGLP SCFLP, knapsack typical oracle, FLCAP-based stochastic data, 7 private cores, Gurobi"
EXPERIMENT_VERSION="1"
SEED="1"
HOUR="04"
EXPERIMENT_DESCRIPTION="${HOUR} hr, seed = ${SEED}"

FILE_NAME="scflp_directional_polar_knapsack_flcap.jl"
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

# Define an array of instance names (FLCAP-based stochastic dataset)
instances=(

    # cap10x (25 facilities)
    "cap101-s256" "cap102-s256" "cap103-s256" "cap104-s256"
    # "cap101-s512" "cap102-s512" "cap103-s512" "cap104-s512"
    # "cap101-s1024" "cap102-s1024" "cap103-s1024" "cap104-s1024"

    # cap12x (50 facilities)
    # "cap121-s256" "cap122-s256" "cap123-s256" "cap124-s256"
    # "cap121-s512" "cap122-s512" "cap123-s512" "cap124-s512"
    # "cap121-s1024" "cap122-s1024" "cap123-s1024" "cap124-s1024"

    # cap13x (50 facilities)
    # "cap131-s256" "cap132-s256" "cap133-s256" "cap134-s256"
    # "cap131-s512" "cap132-s512" "cap133-s512" "cap134-s512"
    # "cap131-s1024" "cap132-s1024" "cap133-s1024" "cap134-s1024"

    # cap6x (16 facilities)
    # "cap61-s256" "cap62-s256" "cap63-s256" "cap64-s256"
    # "cap61-s512" "cap62-s512" "cap63-s512" "cap64-s512"
    # "cap61-s1024" "cap62-s1024" "cap63-s1024" "cap64-s1024"

    # cap7x (16 facilities)
    # "cap71-s256" "cap72-s256" "cap73-s256" "cap74-s256"
    # "cap71-s512" "cap72-s512" "cap73-s512" "cap74-s512"
    # "cap71-s1024" "cap72-s1024" "cap73-s1024" "cap74-s1024"

    # cap9x (16 facilities)
    # "cap91-s256" "cap92-s256" "cap93-s256" "cap94-s256"
    # "cap91-s512" "cap92-s512" "cap93-s512" "cap94-s512"
    # "cap91-s1024" "cap92-s1024" "cap93-s1024" "cap94-s1024"

    # capax (16 facilities)
    # "capa1-s256" "capa2-s256" "capa3-s256" "capa4-s256"
    # "capa1-s512" "capa2-s512" "capa3-s512" "capa4-s512"
    # "capa1-s1024" "capa2-s1024" "capa3-s1024" "capa4-s1024"

    # capbx (16 facilities)
    # "capb1-s256" "capb2-s256" "capb3-s256" "capb4-s256"
    # "capb1-s512" "capb2-s512" "capb3-s512" "capb4-s512"
    # "capb1-s1024" "capb2-s1024" "capb3-s1024" "capb4-s1024"

    # capcx (16 facilities)
    # "capc1-s256" "capc2-s256" "capc3-s256" "capc4-s256"
    # "capc1-s512" "capc2-s512" "capc3-s512" "capc4-s512"
    # "capc1-s1024" "capc2-s1024" "capc3-s1024" "capc4-s1024"

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
    echo "julia --project=polardcglp polardcglp/scripts/${FILE_NAME} --instance=${instance} --seed=${SEED}" >> "${JOBSCRIPT_FILE}"

    # Submit job
    sbatch "${JOBSCRIPT_FILE}"
    rm "${JOBSCRIPT_FILE}"
done
