#!/bin/sh
#SBATCH -t 0-01:00:00
#SBATCH -N 1
#SBATCH -n 10

# Define variables to make the script more readable and maintainable

config_path="scripts/config/milp.toml"
# config_path="scripts/config/config_cflp_disjunctive.toml"

JOBSCRIPT_DIR="scripts"
# OUTPUT_DIR="experiments/cflp_callback_benchmark_hard_cplex"
OUTPUT_DIR="experiments/snip_milp_sequential"
# OUTPUT_DIR="experiments/snip_milp_callback_cplex_turnoff"

# Create necessary directories
mkdir -p "${OUTPUT_DIR}"

# Copy config file to output directory
cp "${config_path}" "${OUTPUT_DIR}/config.toml"

# Define arrays of parameters
instances=(0 1 2 3 4)
snipno=(1 2 3 4)
budgets=(30.0 40.0 50.0 60.0 70.0 80.0 90.0)

# snipno=(1 2)
# budgets=(70.0 80.0 90.0)

# instances=(0)
# snipno=(1)
# budgets=(30.0)


# Loop through all combinations
for instance in "${instances[@]}"; do
    for snip in "${snipno[@]}"; do
        for budget in "${budgets[@]}"; do
            # Create unique job name
            job_name="s${snip}_b${budget}_i${instance}"
            JOBSCRIPT_FILE="${JOBSCRIPT_DIR}/${job_name}.sh"
            
            # Create job script file
            echo "#!/bin/bash" > "${JOBSCRIPT_FILE}"
            echo "#SBATCH -N 1" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -n 20" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -t 0-03:00:00" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH --mem=48G" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -o ${OUTPUT_DIR}/${job_name}.out%j" >> "${JOBSCRIPT_FILE}"
            echo "#SBATCH -e ${OUTPUT_DIR}/${job_name}.err%j" >> "${JOBSCRIPT_FILE}"

            # Load necessary modules
            echo "module purge" >> "${JOBSCRIPT_FILE}"
            echo "module load julia" >> "${JOBSCRIPT_FILE}"
            echo "module load cplex" >> "${JOBSCRIPT_FILE}"
            echo "module load gurobi" >> "${JOBSCRIPT_FILE}"

            # Run Julia script with all parameters
            echo "julia --project=. scripts/snip_callback.jl --instance ${instance} --snip_no ${snip} --budget ${budget} --output_dir ${OUTPUT_DIR} --config ${config_path}" >> "${JOBSCRIPT_FILE}"

            # Submit job
            sbatch "${JOBSCRIPT_FILE}"
            rm "${JOBSCRIPT_FILE}"
        done
    done
done
