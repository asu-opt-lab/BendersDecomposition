#!/bin/sh
#SBATCH -t 0-01:00:00

ROUND_VERSION="simplex_norm_test/uflp_knapsack"
ROUND_DESCRIPTION="SimplexNormTestDCGLP hard UFL, freq500, knapsack, 7 private cores"
EXPERIMENT_VERSION="1"
SEED="1"
HOUR="04"
EXPERIMENT_DESCRIPTION="${HOUR} hr, seed = ${SEED}"

FILE_NAME="uflp_knapsack.jl"
THREADS=7

OUTPUT_DIR="polardcglp/experiment/${ROUND_VERSION}/${EXPERIMENT_VERSION}"
ERR_OUT_DIR="${OUTPUT_DIR}/results"

if [ -d "${OUTPUT_DIR}" ]; then
    echo "Error: Experiment directory ${OUTPUT_DIR} already exists. Please use a different EXPERIMENT_VERSION or remove the existing directory."
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${ERR_OUT_DIR}"

cp -r polardcglp/scripts/simplex_norm_test/${FILE_NAME} "${OUTPUT_DIR}/${FILE_NAME}"

cat > "${OUTPUT_DIR}/experiment_metadata.md" << EOF
# Experiment Metadata

- **Round Version**: ${ROUND_VERSION}
- **Round Description**: ${ROUND_DESCRIPTION}
- **Experiment Version**: ${EXPERIMENT_VERSION}
- **Experiment Description**: ${EXPERIMENT_DESCRIPTION}
- **Date**: $(date "+%Y-%m-%d %H:%M:%S")
EOF

instances=(

    # "ga250a-1" "ga250a-2" "ga250a-3" "ga250a-4" "ga250a-5"
    # "ga250b-1" "ga250b-2" "ga250b-3" "ga250b-4" "ga250b-5"
    # "ga250c-1" "ga250c-2" "ga250c-3" "ga250c-4" "ga250c-5"

    # "gs250a-1" "gs250a-2" "gs250a-3" "gs250a-4" "gs250a-5"
    # "gs250b-1" "gs250b-2" "gs250b-3" "gs250b-4" "gs250b-5"
    # "gs250c-1" "gs250c-2" "gs250c-3" "gs250c-4" "gs250c-5"

    "ga500a-1" "ga500a-2" "ga500a-3" "ga500a-4" "ga500a-5"
    "ga500b-1" "ga500b-2" "ga500b-3" "ga500b-4" "ga500b-5"
    "ga500c-1" "ga500c-2" "ga500c-3" "ga500c-4" "ga500c-5"

    "gs500a-1" "gs500a-2" "gs500a-3" "gs500a-4" "gs500a-5"
    "gs500b-1" "gs500b-2" "gs500b-3" "gs500b-4" "gs500b-5"
    "gs500c-1" "gs500c-2" "gs500c-3" "gs500c-4" "gs500c-5"
)

for instance in "${instances[@]}"; do
    JOBSCRIPT_FILE="${OUTPUT_DIR}/${instance}.sh"

    echo "#!/bin/bash" > "${JOBSCRIPT_FILE}"

    echo "#SBATCH -p htc" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -N 1" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -n 1" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -c ${THREADS}" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH --mem=60G" >> "${JOBSCRIPT_FILE}"

    echo "#SBATCH -t 0-${HOUR}:00:00" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -o ${ERR_OUT_DIR}/${instance}.out%j" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -e ${ERR_OUT_DIR}/${instance}.err%j" >> "${JOBSCRIPT_FILE}"

    echo "module purge" >> "${JOBSCRIPT_FILE}"
    echo "module load julia" >> "${JOBSCRIPT_FILE}"
    echo "module load cplex" >> "${JOBSCRIPT_FILE}"
    echo "module load gurobi" >> "${JOBSCRIPT_FILE}"

    echo "julia --project=polardcglp polardcglp/scripts/simplex_norm_test/${FILE_NAME} --instance=${instance} --seed=${SEED}" >> "${JOBSCRIPT_FILE}"

    sbatch "${JOBSCRIPT_FILE}"
    rm "${JOBSCRIPT_FILE}"
done
