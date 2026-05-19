#!/usr/bin/env bash

set -eo pipefail

# -------- DRY RUN SETTINGS --------
LIST="${SLURM_ARRAY_TASK_ID:-1}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

INPUT_DIR="/group/sbms004/yxia/GUT/cleaned_fastq"
SAMPLE_LIST="/group/sbms004/yxia/GUT/sample_list/sample_list_${LIST}.tsv"
PROJECT_DIR="/group/sbms004/yxia/GUT/SqueezeMeta_projects/sample_list_${LIST}"
GTDBTK_DB="/group/sbms004/yxia/GUT/db/release226"
LOG_DIR="/group/sbms004/yxia/GUT/logs"
# ----------------------------------

echo "=========================================="
echo "[DRY RUN] SqueezeMeta precheck"
echo "Array task/list ID: ${LIST}"
echo "Threads: ${THREADS}"
echo "Input reads dir: ${INPUT_DIR}"
echo "Sample list: ${SAMPLE_LIST}"
echo "Project dir: ${PROJECT_DIR}"
echo "GTDB-Tk DB: ${GTDBTK_DB}"
echo "Log dir: ${LOG_DIR}"
echo "Date: $(date)"
echo "=========================================="

mkdir -p "${LOG_DIR}"
mkdir -p "$(dirname "${PROJECT_DIR}")"

echo
echo "[1/6] Checking conda environment..."
if command -v conda >/dev/null 2>&1; then
    echo "OK: conda found: $(command -v conda)"
else
    echo "WARNING: conda not found in current shell."
    echo "Try: source ~/miniconda3/etc/profile.d/conda.sh"
fi

echo
echo "[2/6] Checking SqueezeMeta.pl..."
if command -v SqueezeMeta.pl >/dev/null 2>&1; then
    echo "OK: SqueezeMeta.pl found: $(command -v SqueezeMeta.pl)"
else
    echo "WARNING: SqueezeMeta.pl not found in PATH."
    echo "If needed, run:"
    echo "  source ~/miniconda3/etc/profile.d/conda.sh"
    echo "  conda activate SqueezeMeta"
fi

echo
echo "[3/6] Checking directories..."
if [[ -d "${INPUT_DIR}" ]]; then
    echo "OK: input directory exists"
else
    echo "ERROR: input directory does not exist: ${INPUT_DIR}"
    exit 1
fi

if [[ -d "${GTDBTK_DB}" ]]; then
    echo "OK: GTDB-Tk DB directory exists"
else
    echo "ERROR: GTDB-Tk DB directory does not exist: ${GTDBTK_DB}"
    exit 1
fi

echo
echo "[4/6] Checking sample list..."
if [[ ! -s "${SAMPLE_LIST}" ]]; then
    echo "ERROR: sample list does not exist or is empty: ${SAMPLE_LIST}"
    exit 1
fi

echo "OK: sample list exists and is not empty"
echo "Lines in sample list:"
wc -l "${SAMPLE_LIST}"

echo
echo "Preview:"
head "${SAMPLE_LIST}"

echo
echo "[5/6] Checking FASTQ files listed in sample list..."
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

echo
echo "FASTQ files checked: ${checked}"
echo "Missing/empty FASTQ files: ${missing}"

if (( missing > 0 )); then
    echo "ERROR: Some FASTQ files are missing or empty."
    exit 1
fi

echo
echo "[6/6] Command that would run:"
echo
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

echo
echo "=========================================="
echo "[DRY RUN PASSED]"
echo "No SqueezeMeta job was started."
echo "Project would be written to:"
echo "${PROJECT_DIR}"
echo "=========================================="