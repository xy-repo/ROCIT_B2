#!/usr/bin/env bash
#SBATCH --job-name=clean
#SBATCH --array=1-30%15
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=/group/sbms004/yxia/GUT/logs/clean_reads_%A_%a.out
#SBATCH --error=/group/sbms004/yxia/GUT/logs/clean_reads_%A_%a.err

source activate py3
set -euo pipefail
shopt -s nullglob

INPUT_DIR="/group/sbms004/yxia/GUT/merged_fastq"
HOST_INDEX="/group/sbms004/yxia/GUT/db/GRCh38_noalt_as/GRCh38_noalt_as"
OUTPUT_DIR="/group/sbms004/yxia/GUT/cleaned_fastq"
THREADS="${SLURM_CPUS_PER_TASK}"

mkdir -p "${OUTPUT_DIR}"
mkdir -p /group/sbms004/yxia/GUT/logs

# Build list of R1 files
R1_LIST=("${INPUT_DIR}"/*_R1.fastq.gz)
TOTAL=${#R1_LIST[@]}

if (( TOTAL == 0 )); then
    echo "ERROR: No R1 files found in ${INPUT_DIR}"
    exit 1
fi

# Ensure task ID does not exceed sample count
if (( SLURM_ARRAY_TASK_ID > TOTAL )); then
    echo "Task ID ${SLURM_ARRAY_TASK_ID} exceeds number of samples (${TOTAL}). Exiting."
    exit 0
fi

# Select files for this task
R1="${R1_LIST[$((SLURM_ARRAY_TASK_ID - 1))]}"
R2="${R1%_R1.fastq.gz}_R2.fastq.gz"
BASENAME="$(basename "${R1}" _R1.fastq.gz)"

FASTP_R1="${OUTPUT_DIR}/${BASENAME}_R1.fastp.fastq.gz"
FASTP_R2="${OUTPUT_DIR}/${BASENAME}_R2.fastp.fastq.gz"

CLEAN_TMP_PREFIX="${OUTPUT_DIR}/${BASENAME}_clean.fastq.gz"
CLEAN_TMP_R1="${OUTPUT_DIR}/${BASENAME}_clean.fastq.1.gz"
CLEAN_TMP_R2="${OUTPUT_DIR}/${BASENAME}_clean.fastq.2.gz"

CLEAN_R1="${OUTPUT_DIR}/${BASENAME}_R1_clean.fastq.gz"
CLEAN_R2="${OUTPUT_DIR}/${BASENAME}_R2_clean.fastq.gz"

echo "=========================================="
echo "[+] Processing sample (${SLURM_ARRAY_TASK_ID}/${TOTAL}): ${BASENAME}"
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Output R1: ${CLEAN_R1}"
echo "Output R2: ${CLEAN_R2}"
echo "=========================================="

# Skip if final R2 clean output already exists
if [[ -s "${CLEAN_R2}" ]]; then
    echo "[+] ${CLEAN_R2} already exists and is not empty."
    echo "[+] Skipping sample: ${BASENAME}"
    exit 0
fi

# Optional stronger skip: if both final clean files exist, skip
if [[ -s "${CLEAN_R1}" && -s "${CLEAN_R2}" ]]; then
    echo "[+] Clean R1 and R2 already exist."
    echo "[+] Skipping sample: ${BASENAME}"
    exit 0
fi

# Check R2 exists
if [[ ! -s "${R2}" ]]; then
    echo "ERROR: R2 file does not exist or is empty:"
    echo "${R2}"
    exit 1
fi

# Remove incomplete intermediate files from previous failed runs
rm -f "${FASTP_R1}" "${FASTP_R2}"
rm -f "${CLEAN_TMP_R1}" "${CLEAN_TMP_R2}"

# Step 1: Quality trimming & adapter removal with fastp
fastp \
    -i "${R1}" \
    -I "${R2}" \
    -o "${FASTP_R1}" \
    -O "${FASTP_R2}" \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --cut_right \
    --cut_mean_quality 20 \
    --length_required 50 \
    --thread "${THREADS}" \
    --html "${OUTPUT_DIR}/${BASENAME}_fastp.html" \
    --json "${OUTPUT_DIR}/${BASENAME}_fastp.json"

# Check fastp outputs
if [[ ! -s "${FASTP_R1}" || ! -s "${FASTP_R2}" ]]; then
    echo "ERROR: fastp output missing or empty for sample: ${BASENAME}"
    exit 1
fi

# Step 2: Remove host contamination with Bowtie2
bowtie2 \
    -x "${HOST_INDEX}" \
    -1 "${FASTP_R1}" \
    -2 "${FASTP_R2}" \
    --very-sensitive \
    --threads "${THREADS}" \
    --un-conc-gz "${CLEAN_TMP_PREFIX}" \
    -S /dev/null

# Check bowtie2 outputs before renaming
if [[ ! -s "${CLEAN_TMP_R1}" || ! -s "${CLEAN_TMP_R2}" ]]; then
    echo "ERROR: Bowtie2 clean output missing or empty for sample: ${BASENAME}"
    echo "Expected:"
    echo "${CLEAN_TMP_R1}"
    echo "${CLEAN_TMP_R2}"
    exit 1
fi

# Rename output files to *_R1_clean.fastq.gz and *_R2_clean.fastq.gz
mv "${CLEAN_TMP_R1}" "${CLEAN_R1}"
mv "${CLEAN_TMP_R2}" "${CLEAN_R2}"

# Check final outputs before deleting intermediate fastp files
if [[ ! -s "${CLEAN_R1}" || ! -s "${CLEAN_R2}" ]]; then
    echo "ERROR: Final clean output missing or empty for sample: ${BASENAME}"
    echo "fastp intermediate files will NOT be deleted."
    exit 1
fi

# Remove intermediate fastp output files
rm -f "${FASTP_R1}" "${FASTP_R2}"

echo "[+] Finished sample: ${BASENAME}"
echo "Final clean R1: ${CLEAN_R1}"
echo "Final clean R2: ${CLEAN_R2}"