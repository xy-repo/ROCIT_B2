#!/usr/bin/env bash
#SBATCH --job-name=inS
#SBATCH --output=/group/sbms004/yxia/GUT/logs/inStrain_%A_%a.out
#SBATCH --error=/group/sbms004/yxia/GUT/logs/inStrain_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=512G
#SBATCH --array=20-31

set -eo pipefail
shopt -s nullglob

# -------------------------
# 1. Activate environment
# -------------------------
# If this does not work on Kaya, replace with:
# source activate mgs
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mgs

# -------------------------
# 2. RFMT group ID
# Example:
# 1  -> RFMT001
# 31 -> RFMT031
#
# Each RFMT group should contain donor + recipient samples.
# -------------------------
IND=$(printf "RFMT%03d" "${SLURM_ARRAY_TASK_ID}")

# -------------------------
# 3. Paths
# -------------------------
THREADS="${SLURM_CPUS_PER_TASK}"

MAG_ROOT="/group/sbms004/yxia/GUT/MAG_per_sample"
CLEAN_DIR="/group/sbms004/yxia/GUT/cleaned_fastq"
WORK_ROOT="/group/sbms004/yxia/GUT/inStrain_final"
LOG_DIR="/group/sbms004/yxia/GUT/logs"

IND_WORK="${WORK_ROOT}/${IND}"
OUT_DIR="${IND_WORK}/results"
COMPARE_OUT="${IND_WORK}/compare_out"

COMBINED_MAGS="${IND_WORK}/${IND}_all_MAGs.fasta"
STB_FILE="${IND_WORK}/${IND}_all_MAGs.stb"
INDEX_PREFIX="${IND_WORK}/${IND}_all_MAGs_index"

SAMPLE_LIST="${IND_WORK}/${IND}_samples.txt"
VALID_SAMPLE_LIST="${IND_WORK}/${IND}_valid_samples.txt"
MAG_LIST="${IND_WORK}/${IND}_MAG_files.txt"
COMPLETED_PROFILE_LIST="${IND_WORK}/${IND}_completed_profiles.txt"

mkdir -p "${LOG_DIR}" "${IND_WORK}" "${OUT_DIR}"

echo "=========================================="
echo "[+] Starting inStrain pipeline"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "RFMT group: ${IND}"
echo "Node: $(hostname)"
echo "MAG root: ${MAG_ROOT}"
echo "Clean reads dir: ${CLEAN_DIR}"
echo "Work dir: ${IND_WORK}"
echo "Output dir: ${OUT_DIR}"
echo "Compare output dir: ${COMPARE_OUT}"
echo "Threads: ${THREADS}"
echo "Start time: $(date)"
echo "=========================================="

# -------------------------
# 4. Find samples for this RFMT group
# Expected structure:
# ${MAG_ROOT}/sample_list_*/RFMT001_3537/RFMT001_3537/results/bins
# -------------------------
echo
echo "[1/7] Finding samples for RFMT group ${IND}..."

: > "${SAMPLE_LIST}"
: > "${VALID_SAMPLE_LIST}"
: > "${MAG_LIST}"
: > "${COMPLETED_PROFILE_LIST}"

for bins_dir in "${MAG_ROOT}"/sample_list_*/"${IND}"_*/"${IND}"_*/results/bins; do
    [[ -d "${bins_dir}" ]] || continue

    sample="$(basename "$(dirname "$(dirname "$(dirname "${bins_dir}")")")")"
    echo "${sample}" >> "${SAMPLE_LIST}"
done

sort -u "${SAMPLE_LIST}" -o "${SAMPLE_LIST}"

if [[ ! -s "${SAMPLE_LIST}" ]]; then
    echo "WARNING: no sample folders found for ${IND}; skipping this RFMT group."
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
echo "[2/7] Checking MAG folders and collecting MAG files..."

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
    echo "WARNING: no valid samples with non-empty MAG folders found for ${IND}; skipping this RFMT group."
    echo "End time: $(date)"
    exit 0
fi

VALID_SAMPLE_COUNT=$(wc -l < "${VALID_SAMPLE_LIST}")

if (( VALID_SAMPLE_COUNT < 2 )); then
    echo "WARNING: only ${VALID_SAMPLE_COUNT} valid sample found for ${IND}; skipping this RFMT group."
    echo "Reason: inStrain compare needs at least 2 samples."
    echo "Valid samples:"
    cat "${VALID_SAMPLE_LIST}"
    echo "End time: $(date)"
    exit 0
fi

if [[ ! -s "${MAG_LIST}" ]]; then
    echo "WARNING: no MAG fasta files found for ${IND}; skipping this RFMT group."
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
# 6. Combine MAGs and generate STB file
# STB format:
# scaffold_name<TAB>MAG_ID
# -------------------------
echo
echo "[3/7] Combining MAGs and generating STB..."
echo "Combined FASTA:"
echo "  ${COMBINED_MAGS}"
echo "STB file:"
echo "  ${STB_FILE}"

: > "${COMBINED_MAGS}"
: > "${STB_FILE}"

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

    # MAG/bin ID used in STB file.
    # Example:
    # sample = RFMT001_3537
    # mag_id = bin.1
    # bin_id = RFMT001_3537_bin.1
    bin_id="${sample}_${mag_id}"

    echo "Adding MAG:"
    echo "  sample: ${sample}"
    echo "  mag: ${mag}"
    echo "  bin_id: ${bin_id}"

    awk -v bin_id="${bin_id}" -v stb="${STB_FILE}" '
        /^>/ {
            header=$0
            sub(/^>/, "", header)
            split(header, a, /[ \t]/)
            old_contig=a[1]

            new_contig=bin_id "_" old_contig

            print ">" new_contig
            print new_contig "\t" bin_id >> stb

            next
        }
        {
            print
        }
    ' "${mag}" >> "${COMBINED_MAGS}"

done < "${MAG_LIST}"

if [[ ! -s "${COMBINED_MAGS}" ]]; then
    echo "WARNING: combined MAG FASTA is empty for ${IND}; skipping this RFMT group."
    echo "End time: $(date)"
    exit 0
fi

if [[ ! -s "${STB_FILE}" ]]; then
    echo "WARNING: STB file is empty for ${IND}; skipping this RFMT group."
    echo "End time: $(date)"
    exit 0
fi

echo
echo "Combined MAG FASTA created:"
ls -lh "${COMBINED_MAGS}"

echo
echo "STB file created:"
ls -lh "${STB_FILE}"

echo
echo "Preview FASTA headers:"
grep '^>' "${COMBINED_MAGS}" | head || true

echo
echo "Preview STB:"
head "${STB_FILE}" || true

# -------------------------
# 7. Build Bowtie2 index
# Rebuild if FASTA has changed
# -------------------------
echo
echo "[4/7] Building Bowtie2 index..."

if [[ -s "${INDEX_PREFIX}.1.bt2" && "${INDEX_PREFIX}.1.bt2" -nt "${COMBINED_MAGS}" ]]; then
    echo "Existing Bowtie2 index is newer than FASTA. Skipping index build."
elif [[ -s "${INDEX_PREFIX}.1.bt2l" && "${INDEX_PREFIX}.1.bt2l" -nt "${COMBINED_MAGS}" ]]; then
    echo "Existing large Bowtie2 index is newer than FASTA. Skipping index build."
else
    echo "Building/rebuilding Bowtie2 index..."
    rm -f "${INDEX_PREFIX}".*.bt2 "${INDEX_PREFIX}".*.bt2l

    bowtie2-build \
        "${COMBINED_MAGS}" \
        "${INDEX_PREFIX}"
fi

# -------------------------
# 8. Map reads and run inStrain profile
# -------------------------
echo
echo "[5/7] Mapping reads and running inStrain profile..."

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

    if [[ -s "${BAM}" && -s "${BAM}.bai" && "${BAM}" -nt "${COMBINED_MAGS}" ]]; then
        echo "Existing BAM and BAI found, and BAM is newer than reference. Skipping mapping."
    else
        echo "Running bowtie2 mapping..."

        rm -f "${BAM}" "${BAM}.bai"

        bowtie2 \
            -x "${INDEX_PREFIX}" \
            -1 "${R1}" \
            -2 "${R2}" \
            -p "${THREADS}" \
            | samtools view -@ "${THREADS}" -bS - \
            | samtools sort -@ "${THREADS}" -o "${BAM}"

        samtools index "${BAM}"
    fi

    # Check whether inStrain profile is already complete.
    # The exact output filenames can vary by inStrain version, so this checks common table names.
    PROFILE_COMPLETE=0

    if [[ -d "${IS_OUT}" ]]; then
        if find "${IS_OUT}" -type f \( \
            -name "*genome_info.tsv" -o \
            -name "*scaffold_info.tsv" -o \
            -name "*SNVs.tsv" -o \
            -name "genome_info.tsv" -o \
            -name "scaffold_info.tsv" \
        \) -size +0c | grep -q .; then
            PROFILE_COMPLETE=1
        fi
    fi

    if (( PROFILE_COMPLETE == 1 )); then
        echo "Existing completed inStrain profile found. Skipping profile:"
        echo "  ${IS_OUT}"
    else
        if [[ -d "${IS_OUT}" ]]; then
            echo "WARNING: existing inStrain folder appears incomplete. Removing:"
            echo "  ${IS_OUT}"
            rm -rf "${IS_OUT}"
        fi

        echo "Running inStrain profile for ${sample}..."

        inStrain profile \
            "${BAM}" \
            "${COMBINED_MAGS}" \
            -p "${THREADS}" \
            -o "${IS_OUT}" \
            -s "${STB_FILE}"
    fi

    # Record only completed profile folders for compare.
    PROFILE_COMPLETE=0

    if [[ -d "${IS_OUT}" ]]; then
        if find "${IS_OUT}" -type f \( \
            -name "*genome_info.tsv" -o \
            -name "*scaffold_info.tsv" -o \
            -name "*SNVs.tsv" -o \
            -name "genome_info.tsv" -o \
            -name "scaffold_info.tsv" \
        \) -size +0c | grep -q .; then
            PROFILE_COMPLETE=1
        fi
    fi

    if (( PROFILE_COMPLETE == 1 )); then
        echo "${IS_OUT}" >> "${COMPLETED_PROFILE_LIST}"
        echo "[+] Completed profile recorded: ${IS_OUT}"
    else
        echo "WARNING: no completed inStrain profile detected for ${sample}; not adding to compare list."
    fi

    echo "[+] Finished sample: ${sample}"

done < "${VALID_SAMPLE_LIST}"

sort -u "${COMPLETED_PROFILE_LIST}" -o "${COMPLETED_PROFILE_LIST}"

echo
echo "Completed inStrain profiles:"
cat "${COMPLETED_PROFILE_LIST}" || true

PROFILE_COUNT=$(wc -l < "${COMPLETED_PROFILE_LIST}" || echo 0)

# -------------------------
# 9. Run inStrain compare
# -------------------------
echo
echo "[6/7] Running inStrain compare..."

if (( PROFILE_COUNT < 2 )); then
    echo "WARNING: fewer than 2 completed inStrain profiles for ${IND}; skipping inStrain compare."
    echo "Completed profile count: ${PROFILE_COUNT}"
else
    echo "Completed profile count: ${PROFILE_COUNT}"
    echo "Compare output:"
    echo "  ${COMPARE_OUT}"

    COMPARE_COMPLETE=0

    if [[ -d "${COMPARE_OUT}" ]]; then
        if find "${COMPARE_OUT}" -type f \( \
            -name "*comparisons*.tsv" -o \
            -name "*compare*.tsv" -o \
            -name "*genome*.tsv" -o \
            -name "comparisonsTable.tsv" \
        \) -size +0c | grep -q .; then
            COMPARE_COMPLETE=1
        fi
    fi

    if (( COMPARE_COMPLETE == 1 )); then
        echo "Existing completed compare output found. Skipping compare:"
        echo "  ${COMPARE_OUT}"
    else
        if [[ -d "${COMPARE_OUT}" ]]; then
            echo "WARNING: existing compare folder appears incomplete. Removing:"
            echo "  ${COMPARE_OUT}"
            rm -rf "${COMPARE_OUT}"
        fi

        mapfile -t PROFILE_DIRS < "${COMPLETED_PROFILE_LIST}"

        echo "Running inStrain compare with these profiles:"
        printf '%s\n' "${PROFILE_DIRS[@]}"

        inStrain compare \
            -i "${PROFILE_DIRS[@]}" \
            -o "${COMPARE_OUT}" \
            -p "${THREADS}"
    fi
fi

# -------------------------
# 10. Final summary
# -------------------------
echo
echo "[7/7] Done."
echo "RFMT group: ${IND}"
echo "Combined MAG FASTA:"
echo "  ${COMBINED_MAGS}"
echo "STB file:"
echo "  ${STB_FILE}"
echo "Profile results:"
echo "  ${OUT_DIR}"
echo "Compare results:"
echo "  ${COMPARE_OUT}"
echo "End time: $(date)"
echo "=========================================="