#!/usr/bin/env bash
#SBATCH --job-name=sqmeta-list7-bigmem
#SBATCH --output=/group/sbms004/yxia/GUT/logs/MAGs_annotation_list7_bigmem_%j.out
#SBATCH --error=/group/sbms004/yxia/GUT/logs/MAGs_annotation_list7_bigmem_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G

set -eo pipefail

# -------- PATHS --------
INPUT_DIR="/group/sbms004/yxia/GUT/cleaned_fastq"
SAMPLE_LIST="/group/sbms004/yxia/GUT/sample_list/sample_list_6.tsv"
PROJECT_ROOT="/group/sbms004/yxia/GUT/SqueezeMeta_projects"
PROJECT_DIR="${PROJECT_ROOT}/sample_list_6_bigmem"
GTDBTK_DB="/group/sbms004/yxia/GUT/db/release226"
LOG_DIR="/group/sbms004/yxia/GUT/logs"
# -----------------------

mkdir -p "${LOG_DIR}"
mkdir -p "${PROJECT_ROOT}"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate SqueezeMeta18

THREADS="${SLURM_CPUS_PER_TASK}"

echo "=========================================="
echo "[+] Starting SqueezeMeta rerun for sample_list_6"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Input reads dir: ${INPUT_DIR}"
echo "Sample list: ${SAMPLE_LIST}"
echo "Project dir: ${PROJECT_DIR}"
echo "GTDB-Tk DB: ${GTDBTK_DB}"
echo "Threads: ${THREADS}"
echo "Memory requested: 256G"
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
wc -l "${SAMPLE_LIST}"
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
        echo "OK: ${sample} ${pair} -> ${fq_path}"
    else
        echo "MISSING/EMPTY: ${sample} ${pair} -> ${fq_path}"
        missing=$((missing + 1))
    fi
done < "${SAMPLE_LIST}"

echo "FASTQ files checked: ${checked}"
echo "Missing/empty FASTQ files: ${missing}"

if (( missing > 0 )); then
    echo "ERROR: Some FASTQ files are missing or empty."
    exit 1
fi

echo
echo "[5/5] Checking project directory..."
if [[ -d "${PROJECT_DIR}" ]]; then
    echo "ERROR: project directory already exists:"
    echo "${PROJECT_DIR}"
    echo "Remove it or change PROJECT_DIR before rerunning."
    exit 1
fi

echo
echo "[+] Running SqueezeMeta..."
echo "Command:"
cat <<EOF
SqueezeMeta.pl \\
    -f "${INPUT_DIR}" \\
    -s "${SAMPLE_LIST}" \\
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
    -s "${SAMPLE_LIST}" \
    -p "${PROJECT_DIR}" \
    -m sequential \
    -t "${THREADS}" \
    -c 2500 \
    --onlybins \
    --gtdbtk \
    -gtdbtk_data_path "${GTDBTK_DB}" \
    -binners metabat2

echo
echo "=========================================="
echo "[+] Finished SqueezeMeta rerun for sample_list_7"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Project dir: ${PROJECT_DIR}"
echo "End time: $(date)"
echo "=========================================="