#!/bin/sh
#SBATCH -t 0-01:00:00

# Define variables to make the script more readable and maintainable

# config_path="scripts/config/config_cflp_benchmark.toml"
config_path="scripts/config/config_scflp_disjunctive_1.toml"

JOBSCRIPT_DIR="scripts"
# OUTPUT_DIR="experiments/cflp_callback_benchmark_hard_cplex"
# OUTPUT_DIR="experiments/scflp_sequential_benchmark_cplex"
OUTPUT_DIR="experiments/scflp_sequential_disjunctive_gurobi_3"

# Create necessary directories
mkdir -p "${OUTPUT_DIR}"

# Copy config file to output directory
cp "${config_path}" "${OUTPUT_DIR}/config.toml"

# Define an array of instance names
# instances=(
    
#     # # 64 scenarios
#     # "f25-c50-s64-r3-1" "f25-c50-s64-r3-2" "f25-c50-s64-r3-3" "f25-c50-s64-r3-4" "f25-c50-s64-r3-5"
#     # "f25-c50-s64-r5-1" "f25-c50-s64-r5-2" "f25-c50-s64-r5-3" "f25-c50-s64-r5-4" "f25-c50-s64-r5-5"
#     # "f25-c50-s64-r10-1" "f25-c50-s64-r10-2" "f25-c50-s64-r10-3" "f25-c50-s64-r10-4" "f25-c50-s64-r10-5"

#     # 128 scenarios
#     "f25-c50-s128-r3-1" "f25-c50-s128-r3-2" "f25-c50-s128-r3-3" "f25-c50-s128-r3-4" "f25-c50-s128-r3-5"
#     "f25-c50-s128-r5-1" "f25-c50-s128-r5-2" "f25-c50-s128-r5-3" "f25-c50-s128-r5-4" "f25-c50-s128-r5-5"
#     "f25-c50-s128-r10-1" "f25-c50-s128-r10-2" "f25-c50-s128-r10-3" "f25-c50-s128-r10-4" "f25-c50-s128-r10-5"

#     # 254 scenarios
#     "f25-c50-s254-r3-1" "f25-c50-s254-r3-2" "f25-c50-s254-r3-3" "f25-c50-s254-r3-4" "f25-c50-s254-r3-5"
#     "f25-c50-s254-r5-1" "f25-c50-s254-r5-2" "f25-c50-s254-r5-3" "f25-c50-s254-r5-4" "f25-c50-s254-r5-5"
#     "f25-c50-s254-r10-1" "f25-c50-s254-r10-2" "f25-c50-s254-r10-3" "f25-c50-s254-r10-4" "f25-c50-s254-r10-5"

#     # 512 scenarios
#     # "f25-c50-s512-r3-1" "f25-c50-s512-r3-2" "f25-c50-s512-r3-3" "f25-c50-s512-r3-4" "f25-c50-s512-r3-5"
#     "f25-c50-s512-r5-1" "f25-c50-s512-r5-2" "f25-c50-s512-r5-3" "f25-c50-s512-r5-4" "f25-c50-s512-r5-5"
#     # "f25-c50-s512-r10-1" "f25-c50-s512-r10-2" "f25-c50-s512-r10-3" "f25-c50-s512-r10-4" "f25-c50-s512-r10-5"

#     # 1024 scenarios
#     "f25-c50-s1024-r3-1" "f25-c50-s1024-r3-2" "f25-c50-s1024-r3-3" "f25-c50-s1024-r3-4" "f25-c50-s1024-r3-5"
#     # "f25-c50-s1024-r5-1" "f25-c50-s1024-r5-2" "f25-c50-s1024-r5-3" "f25-c50-s1024-r5-4" "f25-c50-s1024-r5-5"
#     # "f25-c50-s1024-r10-1" "f25-c50-s1024-r10-2" "f25-c50-s1024-r10-3" "f25-c50-s1024-r10-4" "f25-c50-s1024-r10-5"
# )

instances=(
    # 256 scenarios
    "f50-c50-s256-r3-1" "f50-c50-s256-r3-2" "f50-c50-s256-r3-3" "f50-c50-s256-r3-4" "f50-c50-s256-r3-5"
    "f50-c50-s256-r5-1" "f50-c50-s256-r5-2" "f50-c50-s256-r5-3" "f50-c50-s256-r5-4" "f50-c50-s256-r5-5"
    "f50-c50-s256-r10-1" "f50-c50-s256-r10-2" "f50-c50-s256-r10-3" "f50-c50-s256-r10-4" "f50-c50-s256-r10-5"

    # 512 scenarios
    "f50-c50-s512-r3-1" "f50-c50-s512-r3-2" "f50-c50-s512-r3-3" "f50-c50-s512-r3-4" "f50-c50-s512-r3-5"
    "f50-c50-s512-r5-1" "f50-c50-s512-r5-2" "f50-c50-s512-r5-3" "f50-c50-s512-r5-4" "f50-c50-s512-r5-5"
    "f50-c50-s512-r10-1" "f50-c50-s512-r10-2" "f50-c50-s512-r10-3" "f50-c50-s512-r10-4" "f50-c50-s512-r10-5"

    # 256 scenarios
    "f50-c100-s256-r3-1" "f50-c100-s256-r3-2" "f50-c100-s256-r3-3" "f50-c100-s256-r3-4" "f50-c100-s256-r3-5"
    "f50-c100-s256-r5-1" "f50-c100-s256-r5-2" "f50-c100-s256-r5-3" "f50-c100-s256-r5-4" "f50-c100-s256-r5-5"
    "f50-c100-s256-r10-1" "f50-c100-s256-r10-2" "f50-c100-s256-r10-3" "f50-c100-s256-r10-4" "f50-c100-s256-r10-5"

    # 512 scenarios
    "f50-c100-s512-r3-1" "f50-c100-s512-r3-2" "f50-c100-s512-r3-3" "f50-c100-s512-r3-4" "f50-c100-s512-r3-5"
    "f50-c100-s512-r5-1" "f50-c100-s512-r5-2" "f50-c100-s512-r5-3" "f50-c100-s512-r5-4" "f50-c100-s512-r5-5"
    "f50-c100-s512-r10-1" "f50-c100-s512-r10-2" "f50-c100-s512-r10-3" "f50-c100-s512-r10-4" "f50-c100-s512-r10-5"

    # 256 scenarios
    "f100-c100-s256-r3-1" "f100-c100-s256-r3-2" "f100-c100-s256-r3-3" "f100-c100-s256-r3-4" "f100-c100-s256-r3-5"
    "f100-c100-s256-r5-1" "f100-c100-s256-r5-2" "f100-c100-s256-r5-3" "f100-c100-s256-r5-4" "f100-c100-s256-r5-5"
    "f100-c100-s256-r10-1" "f100-c100-s256-r10-2" "f100-c100-s256-r10-3" "f100-c100-s256-r10-4" "f100-c100-s256-r10-5"
)


# Loop through the instances and create a job script for each
for instance in "${instances[@]}"; do
    JOBSCRIPT_FILE="${JOBSCRIPT_DIR}/${instance}.sh"
    
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
    echo "julia --project=. scripts/scflp_sequential.jl --instance ${instance} --output_dir ${OUTPUT_DIR} --config ${config_path}" >> "${JOBSCRIPT_FILE}"

    # Submit job
    sbatch "${JOBSCRIPT_FILE}"
    rm "${JOBSCRIPT_FILE}"
done
