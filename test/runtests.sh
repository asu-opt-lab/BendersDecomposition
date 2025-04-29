#!/bin/sh
#SBATCH -t 0-01:00:00

# Define variables to make the script more readable and maintainable

OUTPUT_DIR="test/results/round2"
# Define job script directory
JOBSCRIPT_DIR="./job_scripts"

# Create necessary directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${JOBSCRIPT_DIR}"

instances=(
    # "1_test_sequential_typical"
    # "2_test_sequential_in_out_typical"
    # "3_test_sequential_disjunctive"
    "4_test_callback_typical"
    "5_test_callback_disjunctive"
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
    echo "#SBATCH -c 14" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH --nodelist=pcc036" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH --mem=100G" >> "${JOBSCRIPT_FILE}"

    echo "#SBATCH -t 0-04:00:00" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -o ${OUTPUT_DIR}/${instance}.out%j" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -e ${OUTPUT_DIR}/${instance}.err%j" >> "${JOBSCRIPT_FILE}"

    # Load necessary modules
    echo "module purge" >> "${JOBSCRIPT_FILE}"
    echo "module load julia" >> "${JOBSCRIPT_FILE}"
    echo "module load cplex" >> "${JOBSCRIPT_FILE}"
    echo "module load gurobi" >> "${JOBSCRIPT_FILE}"

    # Run Julia script with algorithm parameters
    echo "julia --project=. test/${instance}/runtest.jl" >> "${JOBSCRIPT_FILE}"

    # Submit job
    sbatch "${JOBSCRIPT_FILE}"
    rm "${JOBSCRIPT_FILE}"
done
