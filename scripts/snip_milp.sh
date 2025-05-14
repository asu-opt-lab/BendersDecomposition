#!/bin/sh
#SBATCH -t 0-01:00:00

# Define variables to make the script more readable and maintainable

OUTPUT_DIR="experiments_snip/milp"

# Define job script directory
JOBSCRIPT_DIR="scripts"

# Create necessary directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${JOBSCRIPT_DIR}"

# Copy src directory to output directory
cp -r src "${OUTPUT_DIR}/src"
cp scripts/snip_milp.jl "${OUTPUT_DIR}/src"

# Define an array of instance names
instances=(0 1 2 3 4)
snipnos=(1 2 3 4)
budgets=(30.0 40.0 50.0 60.0 70.0 80.0 90.0)

# Loop through the instances and create a job script for each
for ((i=0; i<${#instances[@]}; i++)); do
    instance="${instances[$i]}"
    snipno="${snipnos[$i]}"
    budget="${budgets[$i]}"

    budget_filename=$(echo "${budget}" | tr '.' '_')
    JOBSCRIPT_FILE="${JOBSCRIPT_DIR}/${instance}-${snipno}-${budget_filename}.sh"
    
    # Create job script file
    echo "#!/bin/bash" > "${JOBSCRIPT_FILE}"

    echo "#SBATCH -p htc" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -q grp_gbyeon" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -N 1" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -n 1" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -c 14" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH --nodelist=pcc037" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH --mem=100G" >> "${JOBSCRIPT_FILE}"

    echo "#SBATCH -t 0-04:00:00" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -o ${OUTPUT_DIR}/${instance}-${snipno}-${budget_filename}.out%j" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -e ${OUTPUT_DIR}/${instance}-${snipno}-${budget_filename}.err%j" >> "${JOBSCRIPT_FILE}"

    # Load necessary modules
    echo "module purge" >> "${JOBSCRIPT_FILE}"
    echo "module load julia" >> "${JOBSCRIPT_FILE}"
    echo "module load cplex" >> "${JOBSCRIPT_FILE}"
    echo "module load gurobi" >> "${JOBSCRIPT_FILE}"

    # Run Julia script with algorithm parameters
    echo "julia --project=. scripts/cflp_milp.jl --snip_instance ${instance} --snip_no ${snipno} --snip_budget ${budget} --output_dir ${OUTPUT_DIR}" >> "${JOBSCRIPT_FILE}"

    # Submit job
    sbatch "${JOBSCRIPT_FILE}"
    rm "${JOBSCRIPT_FILE}"
done
