#!/bin/sh
#SBATCH -t 0-01:00:00

ROUND_VERSION="round7"
ROUND_DESCRIPTION="For new instances"
EXPERIMENT_VERSION="uflp_milp"
EXPERIMENT_DESCRIPTION="benchmark"

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
JOBSCRIPT_DIR="./job_scripts"
# mkdir -p "${JOBSCRIPT_DIR}"

# Copy src directory to output directory
cp -r scripts/uflp_milp.jl "${OUTPUT_DIR}/uflp_milp.jl"

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

    "ga250a-1" "ga250a-2" "ga250a-3" "ga250a-4" "ga250a-5"
    "ga250b-1" "ga250b-2" "ga250b-3" "ga250b-4" "ga250b-5"
    "ga250c-1" "ga250c-2" "ga250c-3" "ga250c-4" "ga250c-5"

    "ga500a-1" "ga500a-2" "ga500a-3" "ga500a-4" "ga500a-5"
    "ga500b-1" "ga500b-2" "ga500b-3" "ga500b-4" "ga500b-5"
    "ga500c-1" "ga500c-2" "ga500c-3" "ga500c-4" "ga500c-5"

    "ga750a-1" "ga750a-2" "ga750a-3" "ga750a-4" "ga750a-5"
    "ga750b-1" "ga750b-2" "ga750b-3" "ga750b-4" "ga750b-5"
    "ga750c-1" "ga750c-2" "ga750c-3" "ga750c-4" "ga750c-5"

    "gs250a-1" "gs250a-2" "gs250a-3" "gs250a-4" "gs250a-5"
    "gs250b-1" "gs250b-2" "gs250b-3" "gs250b-4" "gs250b-5"
    "gs250c-1" "gs250c-2" "gs250c-3" "gs250c-4" "gs250c-5"

    "gs500a-1" "gs500a-2" "gs500a-3" "gs500a-4" "gs500a-5"
    "gs500b-1" "gs500b-2" "gs500b-3" "gs500b-4" "gs500b-5"
    "gs500c-1" "gs500c-2" "gs500c-3" "gs500c-4" "gs500c-5"

    "gs750a-1" "gs750a-2" "gs750a-3" "gs750a-4" "gs750a-5"
    "gs750b-1" "gs750b-2" "gs750b-3" "gs750b-4" "gs750b-5"
)

# Loop through the instances and create a job script for each
for instance in "${instances[@]}"; do
    JOBSCRIPT_FILE="${JOBSCRIPT_DIR}/${instance}.sh"
    
    # Create job script file
    echo "#!/bin/bash" > "${JOBSCRIPT_FILE}"

    echo "#SBATCH -p htc" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -q grp_gbyeon" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -N 1" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -n 1" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -c 7" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH --nodelist=pcc037" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH --mem=60G" >> "${JOBSCRIPT_FILE}"

    echo "#SBATCH -t 0-04:00:00" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -o ${OUTPUT_DIR}/${instance}.out%j" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -e ${OUTPUT_DIR}/${instance}.err%j" >> "${JOBSCRIPT_FILE}"

    # Load necessary modules
    echo "module purge" >> "${JOBSCRIPT_FILE}"
    echo "module load julia" >> "${JOBSCRIPT_FILE}"
    echo "module load cplex" >> "${JOBSCRIPT_FILE}"
    echo "module load gurobi" >> "${JOBSCRIPT_FILE}"

    # Run Julia script with algorithm parameters
    echo "julia --project=. scripts/uflp_milp.jl --instance ${instance} --output_dir ${OUTPUT_DIR}" >> "${JOBSCRIPT_FILE}"

    # Submit job
    sbatch "${JOBSCRIPT_FILE}"
    rm "${JOBSCRIPT_FILE}"
done
