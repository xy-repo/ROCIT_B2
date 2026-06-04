#!/usr/bin/env bash
#SBATCH --job-name=sqmeta-array
#SBATCH --output=/group/sbms004/yxia/GUT/logs/MAGs_annotation_%A_%a.out
#SBATCH --error=/group/sbms004/yxia/GUT/logs/MAGs_annotation_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --array=1-300

set -eo pipefail

# -------- PATHS --------
INPUT_DIR="/group/sbms004/yxia/GUT/cleaned_fastq"
SAMPLE_LIST_DIR="/group/sbms004/yxia/GUT/sample_list"
PROJECT_ROOT="/group/sbms004/yxia/GUT/MAG_per_sample"
GTDBTK_DB="/group/sbms004/yxia/GUT/db/release232/release232"
LOG_DIR="/group/sbms004/yxia/GUT/logs"
# -----------------------

mkdir -p "${LOG_DIR}"
mkdir -p "${PROJECT_ROOT}"

# -------- Activate conda --------
source ~/miniconda3/etc/profile.d/conda.sh

# Avoid PERL5LIB unbound variable during conda activation
set +u
conda activate SqueezeMeta18
set -u

# Important for GTDB-Tk
export GTDBTK_DATA_PATH="${GTDBTK_DB}"

LIST="${SLURM_ARRAY_TASK_ID}"
THREADS="${SLURM_CPUS_PER_TASK}"

SAMPLE_LIST="${SAMPLE_LIST_DIR}/sample_list_${LIST}.tsv"
PROJECT_GROUP_DIR="${PROJECT_ROOT}/sample_list_${LIST}"

mkdir -p "${PROJECT_GROUP_DIR}"

echo "=========================================="
echo "[+] Starting SqueezeMeta array task"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Input reads dir: ${INPUT_DIR}"
echo "Sample list: ${SAMPLE_LIST}"
echo "Project group dir: ${PROJECT_GROUP_DIR}"
echo "GTDB-Tk DB: ${GTDBTK_DB}"
echo "GTDBTK_DATA_PATH: ${GTDBTK_DATA_PATH}"
echo "Threads: ${THREADS}"
echo "Start time: $(date)"
echo "=========================================="

echo
echo "[1/5] Checking input directory..."
if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "ERROR: input directory does not exist: ${INPUT_DIR}"
    exit 1
fi
echo "OK: input directory exists"

echo
echo "[2/5] Checking sample list..."
if [[ ! -s "${SAMPLE_LIST}" ]]; then
    echo "ERROR: sample list does not exist or is empty: ${SAMPLE_LIST}"
    exit 1
fi

echo "OK: sample list exists"
echo "Sample list line count:"
wc -l "${SAMPLE_LIST}"

echo
echo "Sample list preview:"
head "${SAMPLE_LIST}"

echo
echo "[3/5] Checking GTDB-Tk database..."
if [[ ! -d "${GTDBTK_DB}" ]]; then
    echo "ERROR: GTDB-Tk database directory does not exist: ${GTDBTK_DB}"
    exit 1
fi
echo "OK: GTDB-Tk database directory exists"

echo
echo "[4/5] Checking FASTQ files listed in sample list..."

missing=0
checked=0

while IFS=$'\t' read -r sample fastq pair; do
    [[ -z "${sample:-}" ]] && continue
    [[ -z "${fastq:-}" ]] && continue

    fq_path="${INPUT_DIR}/${fastq}"
    checked=$((checked + 1))

    if [[ -s "${fq_path}" ]]; then
        echo "OK: ${sample} ${pair:-NA} -> ${fq_path}"
    else
        echo "MISSING/EMPTY: ${sample} ${pair:-NA} -> ${fq_path}"
        missing=$((missing + 1))
    fi
done < "${SAMPLE_LIST}"

echo
echo "FASTQ files checked: ${checked}"
echo "Missing/empty FASTQ files: ${missing}"

if (( missing > 0 )); then
    echo "ERROR: Some FASTQ files are missing or empty."
    exit 1
fi

echo
echo "[5/5] Running or restarting SqueezeMeta sample by sample..."

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

for sample in "${SAMPLES[@]}"; do
    echo
    echo "=========================================="
    echo "[+] Processing sample: ${sample}"
    echo "=========================================="

    PROJECT_DIR="${PROJECT_GROUP_DIR}/${sample}"
    RESULTS_DIR="${PROJECT_DIR}/results"
    SAMPLE_ONLY_LIST="${PROJECT_GROUP_DIR}/${sample}.sample_list.tsv"

    echo "Project dir: ${PROJECT_DIR}"
    echo "Results dir: ${RESULTS_DIR}"
    echo "Sample-only list: ${SAMPLE_ONLY_LIST}"

    awk -F'\t' -v s="${sample}" 'BEGIN{OFS="\t"} $1 == s {print $1, $2, $3}' "${SAMPLE_LIST}" > "${SAMPLE_ONLY_LIST}"

    if [[ ! -s "${SAMPLE_ONLY_LIST}" ]]; then
        echo "ERROR: failed to create sample-only list for ${sample}"
        exit 1
    fi

    echo
    echo "Sample-only list content:"
    cat "${SAMPLE_ONLY_LIST}"

    echo
    echo "Checking FASTQ files for ${sample}..."

    sample_missing=0

    while IFS=$'\t' read -r s fastq pair; do
        [[ -z "${s:-}" ]] && continue
        [[ -z "${fastq:-}" ]] && continue

        fq_path="${INPUT_DIR}/${fastq}"

        if [[ -s "${fq_path}" ]]; then
            echo "OK: ${s} ${pair:-NA} -> ${fq_path}"
        else
            echo "MISSING/EMPTY: ${s} ${pair:-NA} -> ${fq_path}"
            sample_missing=$((sample_missing + 1))
        fi
    done < "${SAMPLE_ONLY_LIST}"

    if (( sample_missing > 0 )); then
        echo "ERROR: sample ${sample} has missing FASTQ files."
        exit 1
    fi

    if [[ -d "${RESULTS_DIR}" ]]; then
        echo
        echo "[+] Results folder exists for ${sample}"
        echo "[+] Restarting SqueezeMeta"
        echo "Command:"
        echo "SqueezeMeta.pl -p \"${PROJECT_DIR}\" --restart"

        SqueezeMeta.pl \
            -p "${PROJECT_DIR}" \
            --restart

    else
        echo
        echo "[+] Results folder does not exist for ${sample}"
        echo "[+] Starting fresh SqueezeMeta project"
        echo "Command:"
        cat <<EOF
SqueezeMeta.pl \\
    -f "${INPUT_DIR}" \\
    -s "${SAMPLE_ONLY_LIST}" \\
    -p "${PROJECT_DIR}" \\
    -m sequential \\
    -t "${THREADS}" \\
    -c 2500 \\
    --onlybins \\
    --gtdbtk \\
    -gtdbtk_data_path "${GTDBTK_DB}" \\
    -binners metabat2
EOF

        SqueezeMeta.pl \
            -f "${INPUT_DIR}" \
            -s "${SAMPLE_ONLY_LIST}" \
            -p "${PROJECT_DIR}" \
            -m sequential \
            -t "${THREADS}" \
            -c 2500 \
            --onlybins \
            --gtdbtk \
            -gtdbtk_data_path "${GTDBTK_DB}" \
            -binners metabat2
    fi

    echo
    echo "[+] Finished sample: ${sample}"
    echo "Project dir: ${PROJECT_DIR}"
    echo "Results dir: ${RESULTS_DIR}"
    echo "Time: $(date)"
done

echo
echo "=========================================="
echo "[+] Finished SqueezeMeta array task"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Project group dir: ${PROJECT_GROUP_DIR}"
echo "End time: $(date)"
echo "=========================================="