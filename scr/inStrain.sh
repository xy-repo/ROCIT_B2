#!/usr/bin/env bash
#SBATCH --job-name=inS
#SBATCH --output=/group/sbms004/yxia/GUT/logs/inStrain_%A_%a.out
#SBATCH --error=/group/sbms004/yxia/GUT/logs/inStrain_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --array=1-31

set -eo pipefail
shopt -s nullglob

# -------------------------
# 1. Activate environment
# -------------------------
source activate mgs

# -------------------------
# 2. Get individual ID: RFMT001 to RFMT031
# -------------------------
IND=$(printf "RFMT%03d" "${SLURM_ARRAY_TASK_ID}")

# -------------------------
# 3. Paths
# -------------------------
THREADS="${SLURM_CPUS_PER_TASK}"

MAG_ROOT="/group/sbms004/yxia/GUT/MAG_per_sample"
CLEAN_DIR="/group/sbms004/yxia/GUT/cleaned_fastq"
WORK_ROOT="/group/sbms004/yxia/GUT/inStrain"
LOG_DIR="/group/sbms004/yxia/GUT/logs"

IND_WORK="${WORK_ROOT}/${IND}"
OUT_DIR="${IND_WORK}/results"
COMBINED_MAGS="${IND_WORK}/${IND}_all_MAGs.fasta"
INDEX_PREFIX="${IND_WORK}/${IND}_all_MAGs_index"
SAMPLE_LIST="${IND_WORK}/${IND}_samples.txt"
VALID_SAMPLE_LIST="${IND_WORK}/${IND}_valid_samples.txt"
MAG_LIST="${IND_WORK}/${IND}_MAG_files.txt"

mkdir -p "${LOG_DIR}" "${IND_WORK}" "${OUT_DIR}"

echo "=========================================="
echo "[+] Starting inStrain pipeline"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Individual: ${IND}"
echo "Node: $(hostname)"
echo "MAG root: ${MAG_ROOT}"
echo "Clean reads dir: ${CLEAN_DIR}"
echo "Work dir: ${IND_WORK}"
echo "Output dir: ${OUT_DIR}"
echo "Threads: ${THREADS}"
echo "Start time: $(date)"
echo "=========================================="

# -------------------------
# 4. Find samples for this individual
# -------------------------
echo
echo "[1/6] Finding samples for individual ${IND}..."

: > "${SAMPLE_LIST}"
: > "${VALID_SAMPLE_LIST}"
: > "${MAG_LIST}"

for bins_dir in "${MAG_ROOT}"/sample_list_*/"${IND}"_*/"${IND}"_*/results/bins; do
    [[ -d "${bins_dir}" ]] || continue

    sample="$(basename "$(dirname "$(dirname "$(dirname "${bins_dir}")")")")"
    echo "${sample}" >> "${SAMPLE_LIST}"
done

sort -u "${SAMPLE_LIST}" -o "${SAMPLE_LIST}"

if [[ ! -s "${SAMPLE_LIST}" ]]; then
    echo "WARNING: no sample folders found for ${IND}; skipping this individual."
    echo "Expected pattern:"
    echo "  ${MAG_ROOT}/sample_list_*/${IND}_*/${IND}_*/results/bins"
    echo "End time: $(date)"
    exit 0
fi

echo "All candidate samples found:"
cat "${SAMPLE_LIST}"

# -------------------------
# 5. Check MAG folders and collect MAG files
# -------------------------
echo
echo "[2/6] Checking MAG folders and collecting MAG files..."

while read -r sample; do
    [[ -z "${sample}" ]] && continue

    echo
    echo "Checking sample: ${sample}"

    matching_bins_dirs=( "${MAG_ROOT}"/sample_list_*/"${sample}"/"${sample}"/results/bins )

    found_valid_bins=0

    for bins_dir in "${matching_bins_dirs[@]}"; do
        [[ -d "${bins_dir}" ]] || continue

        echo "MAG folder found:"
        echo "  ${bins_dir}"

        mag_files=( "${bins_dir}"/*.fa "${bins_dir}"/*.fna "${bins_dir}"/*.fasta )

        if (( ${#mag_files[@]} == 0 )); then
            echo "WARNING: MAG folder exists but is empty:"
            echo "  ${bins_dir}"
            continue
        fi

        found_valid_bins=1
        echo "${sample}" >> "${VALID_SAMPLE_LIST}"

        for mag in "${mag_files[@]}"; do
            [[ -s "${mag}" ]] || continue
            echo "${sample}"$'\t'"${mag}" >> "${MAG_LIST}"
        done
    done

    if (( found_valid_bins == 0 )); then
        echo "WARNING: no non-empty MAG folder found for ${sample}; skipping this sample."
    fi

done < "${SAMPLE_LIST}"

sort -u "${VALID_SAMPLE_LIST}" -o "${VALID_SAMPLE_LIST}"
sort -u "${MAG_LIST}" -o "${MAG_LIST}"

if [[ ! -s "${VALID_SAMPLE_LIST}" ]]; then
    echo "WARNING: no valid samples with non-empty MAG folders found for ${IND}; skipping this individual."
    echo "End time: $(date)"
    exit 0
fi

if [[ ! -s "${MAG_LIST}" ]]; then
    echo "WARNING: no MAG fasta files found for ${IND}; skipping this individual."
    echo "End time: $(date)"
    exit 0
fi

echo
echo "Valid samples with non-empty MAG folders:"
cat "${VALID_SAMPLE_LIST}"

echo
echo "MAG files collected:"
cat "${MAG_LIST}"

MAG_COUNT=$(wc -l < "${MAG_LIST}")
echo
echo "Total MAG files: ${MAG_COUNT}"

# -------------------------
# 6. Combine MAGs
# -------------------------
echo
echo "[3/6] Combining MAGs into:"
echo "  ${COMBINED_MAGS}"

: > "${COMBINED_MAGS}"

while IFS=$'\t' read -r sample mag; do
    [[ -z "${sample}" ]] && continue
    [[ -z "${mag}" ]] && continue

    if [[ ! -s "${mag}" ]]; then
        echo "WARNING: MAG file missing or empty, skipping:"
        echo "  ${mag}"
        continue
    fi

    mag_base="$(basename "${mag}")"
    mag_id="${mag_base%.*}"

    echo "Adding MAG:"
    echo "  sample: ${sample}"
    echo "  mag: ${mag}"

    sed "s/^>/>${sample}_${mag_id}_/" "${mag}" >> "${COMBINED_MAGS}"

done < "${MAG_LIST}"

if [[ ! -s "${COMBINED_MAGS}" ]]; then
    echo "WARNING: combined MAG fasta is empty for ${IND}; skipping this individual."
    echo "End time: $(date)"
    exit 0
fi

echo
echo "Combined MAG fasta created:"
ls -lh "${COMBINED_MAGS}"

# -------------------------
# 7. Build Bowtie2 index
# -------------------------
echo
echo "[4/6] Building Bowtie2 index..."

if [[ -s "${INDEX_PREFIX}.1.bt2" || -s "${INDEX_PREFIX}.1.bt2l" ]]; then
    echo "Existing Bowtie2 index found. Skipping index build."
else
    bowtie2-build \
        "${COMBINED_MAGS}" \
        "${INDEX_PREFIX}"
fi

# -------------------------
# 8. Map reads and run inStrain
# -------------------------
echo
echo "[5/6] Mapping reads and running inStrain profile..."

while read -r sample; do
    [[ -z "${sample}" ]] && continue

    echo
    echo "=========================================="
    echo "[+] Processing sample: ${sample}"
    echo "=========================================="

    R1="${CLEAN_DIR}/${sample}_R1_clean.fastq.gz"
    R2="${CLEAN_DIR}/${sample}_R2_clean.fastq.gz"

    BAM="${OUT_DIR}/${sample}.bam"
    IS_OUT="${OUT_DIR}/${sample}_IS"

    if [[ ! -s "${R1}" ]]; then
        echo "WARNING: R1 not found or empty:"
        echo "  ${R1}"
        echo "Skipping ${sample}"
        continue
    fi

    if [[ ! -s "${R2}" ]]; then
        echo "WARNING: R2 not found or empty:"
        echo "  ${R2}"
        echo "Skipping ${sample}"
        continue
    fi

    echo "R1: ${R1}"
    echo "R2: ${R2}"
    echo "BAM: ${BAM}"
    echo "inStrain output: ${IS_OUT}"

    if [[ -s "${BAM}" && -s "${BAM}.bai" ]]; then
        echo "Existing BAM and index found. Skipping bowtie2 mapping."
    else
        echo "Running bowtie2 mapping..."

        bowtie2 \
            -x "${INDEX_PREFIX}" \
            -1 "${R1}" \
            -2 "${R2}" \
            -p "${THREADS}" \
            | samtools view -@ "${THREADS}" -bS - \
            | samtools sort -@ "${THREADS}" -o "${BAM}"

        samtools index "${BAM}"
    fi

    if [[ -d "${IS_OUT}" ]]; then
        echo "Existing inStrain output folder found:"
        echo "  ${IS_OUT}"
        echo "Skipping inStrain profile for ${sample}"
    else
        echo "Running inStrain profile for ${sample}..."

        inStrain profile \
            "${BAM}" \
            "${COMBINED_MAGS}" \
            -p "${THREADS}" \
            -o "${IS_OUT}"
    fi

    echo "[+] Finished sample: ${sample}"

done < "${VALID_SAMPLE_LIST}"

echo
echo "[6/6] Done."
echo "Individual: ${IND}"
echo "Results:"
echo "  ${OUT_DIR}"
echo "End time: $(date)"
echo "=========================================="