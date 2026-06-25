#!/bin/sh
#SBATCH -t 0-01:00:00

ROUND_VERSION="vertical_reverse_polar/scflp_bodur"
ROUND_DESCRIPTION="vertical_reverse_polar SCFLP, Bodur data, freq250, knapsack/classical oracle, 7 private cores, CPLEX"
EXPERIMENT_VERSION="0_noroot"
SEED="1"
HOUR="04"
TIME_LIMIT="14400"
DCGLP_TIME_LIMIT="1000"
DCGLP_ITER_LIMIT="250"
DCGLP_HALT_LIMIT="3"
FREQUENCY="250"
REUSE_DCGLP="false"
STRENGTHENED="true"
LIFT="false"
ORACLE="knapsack"
EXPERIMENT_DESCRIPTION="${HOUR} hr, seed = ${SEED}, vertical_reverse_polar, oracle = ${ORACLE}"

FILE_NAME="scflp_bodur.jl"
SHELL_FILE_NAME="scflp_bodur.sh"
THREADS=7

OUTPUT_DIR="polardcglp/experiment/${ROUND_VERSION}/${EXPERIMENT_VERSION}"
ERR_OUT_DIR="${OUTPUT_DIR}/results"

if [ -d "${OUTPUT_DIR}" ]; then
    echo "Error: Experiment directory ${OUTPUT_DIR} already exists. Please use a different EXPERIMENT_VERSION or remove the existing directory."
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${ERR_OUT_DIR}"

cp -r "polardcglp/scripts/vertical_reverse_polar/${FILE_NAME}" "${OUTPUT_DIR}/${FILE_NAME}"
cp -r "polardcglp/scripts/vertical_reverse_polar/${SHELL_FILE_NAME}" "${OUTPUT_DIR}/${SHELL_FILE_NAME}"

cat > "${OUTPUT_DIR}/experiment_metadata.md" << EOF
# Experiment Metadata

- **Round Version**: ${ROUND_VERSION}
- **Round Description**: ${ROUND_DESCRIPTION}
- **Experiment Version**: ${EXPERIMENT_VERSION}
- **Experiment Description**: ${EXPERIMENT_DESCRIPTION}
- **Date**: $(date "+%Y-%m-%d %H:%M:%S")
EOF

BASE_INSTANCES=(
    "cap101" "cap102" "cap103" "cap104"
    "cap121" "cap122" "cap123" "cap124"
    "cap131" "cap132" "cap133" "cap134"
    # "cap61" "cap62" "cap63" "cap64"
    # "cap71" "cap72" "cap73" "cap74"
    # "cap91" "cap92" "cap93" "cap94"
    # "cap123" "cap124"
)

SCENARIO_SIZES=(
    "250" "500" "1500"
    # "500" "1500"
)

for base_instance in "${BASE_INSTANCES[@]}"; do
for scenario_size in "${SCENARIO_SIZES[@]}"; do
    instance="${base_instance}-s${scenario_size}"
    JOBSCRIPT_FILE="${OUTPUT_DIR}/${instance}.sh"

    echo "#!/bin/bash" > "${JOBSCRIPT_FILE}"
    echo "#SBATCH -p general" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -q grp_gbyeon" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -N 1" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -n 1" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -c ${THREADS}" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH --nodelist=pcc036,pcc037" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH --mem=60G" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -t 0-${HOUR}:30:00" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -o ${ERR_OUT_DIR}/${instance}.out%j" >> "${JOBSCRIPT_FILE}"
    echo "#SBATCH -e ${ERR_OUT_DIR}/${instance}.err%j" >> "${JOBSCRIPT_FILE}"
    echo "module purge" >> "${JOBSCRIPT_FILE}"
    echo "module load julia" >> "${JOBSCRIPT_FILE}"
    echo "module load cplex" >> "${JOBSCRIPT_FILE}"
    echo "julia --project=polardcglp polardcglp/scripts/vertical_reverse_polar/${FILE_NAME} --instance=${instance} --seed=${SEED} --time_limit=${TIME_LIMIT} --dcglp_time_limit=${DCGLP_TIME_LIMIT} --dcglp_iter_limit=${DCGLP_ITER_LIMIT} --dcglp_halt_limit=${DCGLP_HALT_LIMIT} --frequency=${FREQUENCY} --threads=${THREADS} --reuse_dcglp=${REUSE_DCGLP} --strengthened=${STRENGTHENED} --lift=${LIFT} --oracle=${ORACLE}" >> "${JOBSCRIPT_FILE}"

    sbatch "${JOBSCRIPT_FILE}"
    rm "${JOBSCRIPT_FILE}"
done
done
