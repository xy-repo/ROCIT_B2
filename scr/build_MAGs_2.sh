#!/usr/bin/env bash
#SBATCH --job-name=sqmeta-restart
#SBATCH --output=/group/sbms004/yxia/GUT/logs/MAGs_restart_%A_%a.out
#SBATCH --error=/group/sbms004/yxia/GUT/logs/MAGs_restart_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --array=1-30

set -euo pipefail

# -------- PATHS --------
SAMPLE_LIST_DIR="/group/sbms004/yxia/GUT/sample_list"
PROJECT_ROOT="/group/sbms004/yxia/GUT/SqueezeMeta_projects"
GTDBTK_DB="/group/sbms004/yxia/GUT/db/release232/release232"
LOG_DIR="/group/sbms004/yxia/GUT/logs"
# -----------------------

mkdir -p "${LOG_DIR}"

# -------- Activate conda --------
source ~/miniconda3/etc/profile.d/conda.sh

set +u
conda activate SqueezeMeta18
set -u

export GTDBTK_DATA_PATH="${GTDBTK_DB}"

LIST="${SLURM_ARRAY_TASK_ID}"
SAMPLE_LIST="${SAMPLE_LIST_DIR}/sample_list_${LIST}.tsv"
PROJECT_GROUP_DIR="${PROJECT_ROOT}/sample_list_${LIST}"

echo "=========================================="
echo "[+] Starting SqueezeMeta restart-only job"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Sample list: ${SAMPLE_LIST}"
echo "Project group dir: ${PROJECT_GROUP_DIR}"
echo "GTDBTK_DATA_PATH: ${GTDBTK_DATA_PATH}"
echo "Start time: $(date)"
echo "=========================================="

echo
echo "[1/3] Checking sample list..."
if [[ ! -s "${SAMPLE_LIST}" ]]; then
    echo "ERROR: sample list does not exist or is empty: ${SAMPLE_LIST}"
    exit 1
fi

echo "OK: sample list exists"
wc -l "${SAMPLE_LIST}"
head "${SAMPLE_LIST}"

echo
echo "[2/3] Checking project group directory..."
if [[ ! -d "${PROJECT_GROUP_DIR}" ]]; then
    echo "ERROR: project group directory does not exist: ${PROJECT_GROUP_DIR}"
    exit 1
fi
echo "OK: project group directory exists"

echo
echo "[3/3] Collecting sample names..."

mapfile -t SAMPLES < <(
    awk -F'\t' 'NF >= 2 && $1 != "" {print $1}' "${SAMPLE_LIST}" | sort -u
)

if (( ${#SAMPLES[@]} == 0 )); then
    echo "ERROR: no sample names found in ${SAMPLE_LIST}"
    exit 1
fi

echo
echo "Samples found:"
printf '  %s\n' "${SAMPLES[@]}"

restarted=0
skipped=0
failed=0

for sample in "${SAMPLES[@]}"; do
    echo
    echo "------------------------------------------"
    echo "[+] Checking sample: ${sample}"
    echo "------------------------------------------"

    # Your real structure:
    # /group/sbms004/yxia/GUT/SqueezeMeta_projects/sample_list_1/RFMT001/RFMT001
    PROJECT_PARENT_DIR="${PROJECT_GROUP_DIR}/${sample}"
    PROJECT_DIR="${PROJECT_PARENT_DIR}/${sample}"
    RESULTS_DIR="${PROJECT_DIR}/results"

    echo "Project parent dir: ${PROJECT_PARENT_DIR}"
    echo "Project dir:        ${PROJECT_DIR}"
    echo "Results dir:        ${RESULTS_DIR}"

    if [[ ! -d "${PROJECT_PARENT_DIR}" ]]; then
        echo "[SKIP] Project parent folder does not exist:"
        echo "       ${PROJECT_PARENT_DIR}"
        skipped=$((skipped + 1))
        continue
    fi

    if [[ ! -d "${PROJECT_DIR}" ]]; then
        echo "[SKIP] Nested SqueezeMeta project folder does not exist:"
        echo "       ${PROJECT_DIR}"
        skipped=$((skipped + 1))
        continue
    fi

    if [[ ! -d "${RESULTS_DIR}" ]]; then
        echo "[SKIP] Results folder does not exist:"
        echo "       ${RESULTS_DIR}"
        echo "       Not starting a fresh SqueezeMeta run."
        skipped=$((skipped + 1))
        continue
    fi

    echo "[RUN] Results folder exists:"
    echo "      ${RESULTS_DIR}"
    echo "[RUN] Restarting SqueezeMeta project:"
    echo "      ${PROJECT_DIR}"
    echo
    echo "Command:"
    echo "SqueezeMeta.pl -p \"${PROJECT_DIR}\" --restart"

    if SqueezeMeta.pl \
        -p "${PROJECT_DIR}" \
        --restart
    then
        echo "[OK] Restart finished for sample: ${sample}"
        restarted=$((restarted + 1))
    else
        echo "[ERROR] Restart failed for sample: ${sample}"
        failed=$((failed + 1))
    fi

    echo "Time: $(date)"
done

echo
echo "=========================================="
echo "[+] Restart-only job finished"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Project group dir: ${PROJECT_GROUP_DIR}"
echo "Restarted samples: ${restarted}"
echo "Skipped samples: ${skipped}"
echo "Failed samples: ${failed}"
echo "End time: $(date)"
echo "=========================================="

if (( failed > 0 )); then
    echo "ERROR: one or more samples failed during restart."
    exit 1
fi