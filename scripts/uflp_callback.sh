#!/bin/sh
#SBATCH -t 0-01:00:00

# Define variables to make the script more readable and maintainable

# config_path="scripts/config/config_uflp_disjunctive.toml"
config_path="scripts/config/config_uflp_knapsack.toml"

JOBSCRIPT_DIR="scripts"
OUTPUT_DIR="knapsack_cplex_HF1_250107"
# OUTPUT_DIR="knapsack_cplex_LBH1_250107"

# Create necessary directories
mkdir -p "${OUTPUT_DIR}"

# Copy config file to output directory
cp "${config_path}" "${OUTPUT_DIR}/config.toml"

# Define an array of instance names
instances=(
    "ga500a-1" "ga500a-2" "ga500a-3" "ga500a-4" "ga500a-5"
    "ga500b-1" "ga500b-2" "ga500b-3" "ga500b-4" "ga500b-5"
    "ga750a-1" "ga750a-2" "ga750a-3" "ga750a-4" "ga750a-5"
    "ga750b-1" "ga750b-2" "ga750b-3" "ga750b-4" "ga750b-5"
    "ga750c-1" "ga750c-2" "ga750c-3" "ga750c-4" "ga750c-5"
)

# instances=(
#     "ga250a-3"
# )

# Loop through the instances and create a job script for each
for instance in "${instances[@]}"; do
    JOBSCRIPT_FILE="${JOBSCRIPT_DIR}/${OUTPUT_DIR}${instance}.sh"
    
    # Create job script file
    echo "#!/bin/bash" > "${JOBSCRIPT_FILE}"
    echo "#SBATCH -N 1" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -n 20" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -t 0-03:00:00" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -o ${OUTPUT_DIR}/${instance}.out%j" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -e ${OUTPUT_DIR}/${instance}.err%j" >> "${JOBSCRIPT_FILE}"

    # Load necessary modules
    echo "module purge" >> "${JOBSCRIPT_FILE}"
    echo "module load julia" >> "${JOBSCRIPT_FILE}"
    echo "module load cplex" >> "${JOBSCRIPT_FILE}"
    echo "module load gurobi" >> "${JOBSCRIPT_FILE}"

    # Run Julia script with algorithm parameters
    echo "julia --project=. scripts/uflp_callback.jl --instance ${instance} --output_dir ${OUTPUT_DIR} --config ${config_path}" >> "${JOBSCRIPT_FILE}"

    # Submit job
    sbatch "${JOBSCRIPT_FILE}"
    rm "${JOBSCRIPT_FILE}"
done
