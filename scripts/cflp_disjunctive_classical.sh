#!/bin/sh
#SBATCH -t 0-01:00:00

ROUND_VERSION="round5"
ROUND_DESCRIPTION="For instances 700 and 1000"
EXPERIMENT_VERSION="cflp_disjunctive_classical_01"
EXPERIMENT_DESCRIPTION="try it first"

# Define variables to make the script more readable and maintainable
OUTPUT_DIR="experiments/${ROUND_VERSION}/${EXPERIMENT_VERSION}"
ERR_OUT_DIR="${OUTPUT_DIR}/results"

# Create necessary directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${ERR_OUT_DIR}"

# Define job script directory
JOBSCRIPT_DIR="./job_scripts"
# mkdir -p "${JOBSCRIPT_DIR}"

# Copy src directory to output directory
cp -r scripts/cflp_disjunctive_classical.jl "${OUTPUT_DIR}/cflp_disjunctive_classical.jl"

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
#     # # 10 facilities, 10 customers
#     # "f10-c10-r3-1" "f10-c10-r3-2"
#     # "f10-c10-r5-1" "f10-c10-r5-2"
#     # "f10-c10-r10-1" "f10-c10-r10-2"

#     # # 25 facilities, 25 customers
#     # "f25-c25-r3-1" "f25-c25-r3-2"
#     # "f25-c25-r5-1" "f25-c25-r5-2"
#     # "f25-c25-r10-1" "f25-c25-r10-2"

#     # # 50 facilities, 50 customers
#     # "f50-c50-r3-1" "f50-c50-r3-2"
#     # "f50-c50-r5-1" "f50-c50-r5-2"
#     # "f50-c50-r10-1" "f50-c50-r10-2"
    
#     # 100 facilities, 100 customers
    # "f100-c100-r3-1" "f100-c100-r3-2" "f100-c100-r3-3" "f100-c100-r3-4" "f100-c100-r3-5"
    # "f100-c100-r5-1" "f100-c100-r5-2" "f100-c100-r5-3" "f100-c100-r5-4" "f100-c100-r5-5"
    # "f100-c100-r10-1" "f100-c100-r10-2" "f100-c100-r10-3" "f100-c100-r10-4" "f100-c100-r10-5"
    
# #     # 100 facilities, 200 customers
#     "f100-c200-r3-1" "f100-c200-r3-2" "f100-c200-r3-3" "f100-c200-r3-4" "f100-c200-r3-5"
#     "f100-c200-r5-1" "f100-c200-r5-2" "f100-c200-r5-3" "f100-c200-r5-4" "f100-c200-r5-5"
#     "f100-c200-r10-1" "f100-c200-r10-2" "f100-c200-r10-3" "f100-c200-r10-4" "f100-c200-r10-5"

# # #     # 100 facilities, 500 customers
#     "f100-c500-r3-1" "f100-c500-r3-2" "f100-c500-r3-3" "f100-c500-r3-4" "f100-c500-r3-5"
#     "f100-c500-r5-1" "f100-c500-r5-2" "f100-c500-r5-3" "f100-c500-r5-4" "f100-c500-r5-5"
#     "f100-c500-r10-1" "f100-c500-r10-2" "f100-c500-r10-3" "f100-c500-r10-4" "f100-c500-r10-5"

# # # #     # 200 facilities, 200 customers
#     "f200-c200-r3-1" "f200-c200-r3-2" "f200-c200-r3-3" "f200-c200-r3-4" "f200-c200-r3-5"
#     "f200-c200-r5-1" "f200-c200-r5-2" "f200-c200-r5-3" "f200-c200-r5-4" "f200-c200-r5-5"
#     "f200-c200-r10-1" "f200-c200-r10-2" "f200-c200-r10-3" "f200-c200-r10-4" "f200-c200-r10-5"

# # #     # 200 facilities, 500 customers
#     "f200-c500-r3-1" "f200-c500-r3-2" "f200-c500-r3-3" "f200-c500-r3-4" "f200-c500-r3-5"
#     "f200-c500-r5-1" "f200-c500-r5-2" "f200-c500-r5-3" "f200-c500-r5-4" "f200-c500-r5-5"
#     "f200-c500-r10-1" "f200-c500-r10-2" "f200-c500-r10-3" "f200-c500-r10-4" "f200-c500-r10-5"
    
#     # 300 facilities, 300 customers
#     "f300-c300-r3-1" "f300-c300-r3-2" "f300-c300-r3-3" "f300-c300-r3-4" "f300-c300-r3-5"
#     "f300-c300-r5-1" "f300-c300-r5-2" "f300-c300-r5-3" "f300-c300-r5-4" "f300-c300-r5-5"
#     "f300-c300-r10-1" "f300-c300-r10-2" "f300-c300-r10-3" "f300-c300-r10-4" "f300-c300-r10-5"

#     # 300 facilities, 1500 customers
#     "f300-c1500-r3-1" "f300-c1500-r3-2" "f300-c1500-r3-3" "f300-c1500-r3-4" "f300-c1500-r3-5"
#     "f300-c1500-r5-1" "f300-c1500-r5-2" "f300-c1500-r5-3" "f300-c1500-r5-4" "f300-c1500-r5-5"
#     "f300-c1500-r10-1" "f300-c1500-r10-2" "f300-c1500-r10-3" "f300-c1500-r10-4" "f300-c1500-r10-5"

    # 500 facilities, 500 customers
    # "f500-c500-r3-1" "f500-c500-r3-2" "f500-c500-r3-3" "f500-c500-r3-4" "f500-c500-r3-5"
    # "f500-c500-r5-1" "f500-c500-r5-2" "f500-c500-r5-3" "f500-c500-r5-4" "f500-c500-r5-5"
    # "f500-c500-r10-1" "f500-c500-r10-2" "f500-c500-r10-3" "f500-c500-r10-4" "f500-c500-r10-5"

    # # 600 facilities, 1500 customers
    # "f600-c1500-r3-1" "f600-c1500-r3-2" "f600-c1500-r3-3" "f600-c1500-r3-4" "f600-c1500-r3-5"
    # "f600-c1500-r5-1" "f600-c1500-r5-2" "f600-c1500-r5-3" "f600-c1500-r5-4" "f600-c1500-r5-5"
    # "f600-c1500-r10-1" "f600-c1500-r10-2" "f600-c1500-r10-3" "f600-c1500-r10-4" "f600-c1500-r10-5"

    # 700 facilities, 700 customers
    "f700-c700-r3-1" "f700-c700-r3-2" "f700-c700-r3-3" "f700-c700-r3-4" "f700-c700-r3-5"
    "f700-c700-r5-1" "f700-c700-r5-2" "f700-c700-r5-3" "f700-c700-r5-4" "f700-c700-r5-5"
    "f700-c700-r10-1" "f700-c700-r10-2" "f700-c700-r10-3" "f700-c700-r10-4" "f700-c700-r10-5"

    # # 1000 facilities, 1000 customers
    # "f1000-c1000-r3-1" "f1000-c1000-r3-2" "f1000-c1000-r3-3" "f1000-c1000-r3-4" "f1000-c1000-r3-5"
    # "f1000-c1000-r5-1" "f1000-c1000-r5-2" "f1000-c1000-r5-3" "f1000-c1000-r5-4" "f1000-c1000-r5-5"
    # "f1000-c1000-r10-1" "f1000-c1000-r10-2" "f1000-c1000-r10-3" "f1000-c1000-r10-4" "f1000-c1000-r10-5"

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
    echo "#SBATCH --nodelist=pcc036" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH --mem=60G" >> "${JOBSCRIPT_FILE}"

    echo "#SBATCH -t 0-04:00:00" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -o ${ERR_OUT_DIR}/${instance}.out%j" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -e ${ERR_OUT_DIR}/${instance}.err%j" >> "${JOBSCRIPT_FILE}"

    # Load necessary modules
    echo "module purge" >> "${JOBSCRIPT_FILE}"
    echo "module load julia" >> "${JOBSCRIPT_FILE}"
    echo "module load cplex" >> "${JOBSCRIPT_FILE}"
    echo "module load gurobi" >> "${JOBSCRIPT_FILE}"

    # Run Julia script with algorithm parameters
    echo "julia --project=. scripts/cflp_callback_disjunctive.jl --instance ${instance} --output_dir ${ERR_OUT_DIR}" >> "${JOBSCRIPT_FILE}"

    # Submit job
    sbatch "${JOBSCRIPT_FILE}"
    rm "${JOBSCRIPT_FILE}"
done

