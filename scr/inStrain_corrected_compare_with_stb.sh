#!/usr/bin/env bash
#SBATCH --job-name=inS_cmp
#SBATCH --output=/group/sbms004/yxia/GUT/logs/inStrain_compare_%A_%a.out
#SBATCH --error=/group/sbms004/yxia/GUT/logs/inStrain_compare_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --array=1-20

set -eo pipefail
shopt -s nullglob

# -------------------------
# 1. Activate environment
# -------------------------
# Use the one that works on your cluster.

source ~/miniconda3/etc/profile.d/conda.sh
conda activate mgs

# If the above fails, comment it out and use:
# source activate mgs

# -------------------------
# 2. RFMT group ID
# -------------------------
IND=$(printf "RFMT%03d" "${SLURM_ARRAY_TASK_ID}")

THREADS="${SLURM_CPUS_PER_TASK}"

# -------------------------
# 3. Paths
# -------------------------
WORK_ROOT="/group/sbms004/yxia/GUT/inStrain_final"
LOG_DIR="/group/sbms004/yxia/GUT/logs"

IND_WORK="${WORK_ROOT}/${IND}"
PROFILE_ROOT="${IND_WORK}/results"

STB_FILE="${IND_WORK}/${IND}_all_MAGs.stb"

# Keep your old compare_out untouched.
# New output goes here:
COMPARE_OUT="${IND_WORK}/compare_out_with_stb"

COMPLETED_PROFILE_LIST="${IND_WORK}/${IND}_completed_profiles_for_compare_with_stb.txt"

mkdir -p "${LOG_DIR}" "${IND_WORK}"

echo "=========================================="
echo "[+] Starting inStrain compare-only step"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "RFMT group: ${IND}"
echo "Profile root: ${PROFILE_ROOT}"
echo "STB file: ${STB_FILE}"
echo "Compare output: ${COMPARE_OUT}"
echo "Threads: ${THREADS}"
echo "Start time: $(date)"
echo "=========================================="

# -------------------------
# 4. Check required inputs
# -------------------------
if [[ ! -d "${PROFILE_ROOT}" ]]; then
    echo "WARNING: profile root folder does not exist:"
    echo "  ${PROFILE_ROOT}"
    echo "Skipping ${IND}."
    echo "End time: $(date)"
    exit 0
fi

if [[ ! -s "${STB_FILE}" ]]; then
    echo "WARNING: STB file missing or empty:"
    echo "  ${STB_FILE}"
    echo "Skipping ${IND}."
    echo "End time: $(date)"
    exit 0
fi

# -------------------------
# 5. Collect completed inStrain profile folders
# -------------------------
echo
echo "[1/3] Collecting completed profile folders..."

: > "${COMPLETED_PROFILE_LIST}"

for IS_OUT in "${PROFILE_ROOT}"/*_IS; do
    [[ -d "${IS_OUT}" ]] || continue

    PROFILE_COMPLETE=0

    # Check for common inStrain profile output files.
    # Different versions may use slightly different names.
    if find "${IS_OUT}" -type f \( \
        -name "*genome_info.tsv" -o \
        -name "genome_info.tsv" -o \
        -name "*scaffold_info.tsv" -o \
        -name "scaffold_info.tsv" -o \
        -name "*SNVs.tsv" -o \
        -name "SNVs.tsv" \
    \) -size +0c | grep -q .; then
        PROFILE_COMPLETE=1
    fi

    if (( PROFILE_COMPLETE == 1 )); then
        echo "${IS_OUT}" >> "${COMPLETED_PROFILE_LIST}"
        echo "Completed profile found:"
        echo "  ${IS_OUT}"
    else
        echo "WARNING: profile folder found but no expected output table detected. Skipping:"
        echo "  ${IS_OUT}"
    fi
done

sort -u "${COMPLETED_PROFILE_LIST}" -o "${COMPLETED_PROFILE_LIST}"

PROFILE_COUNT=$(wc -l < "${COMPLETED_PROFILE_LIST}" || echo 0)

echo
echo "Completed profile count: ${PROFILE_COUNT}"

# Need at least 2 samples/profiles for inStrain compare.
if (( PROFILE_COUNT < 2 )); then
    echo "WARNING: fewer than 2 completed profiles for ${IND}; skipping inStrain compare."
    echo "Profiles found:"
    cat "${COMPLETED_PROFILE_LIST}" || true
    echo "End time: $(date)"
    exit 0
fi

echo
echo "Profiles that will be compared:"
cat "${COMPLETED_PROFILE_LIST}"

# -------------------------
# 6. Prepare compare output
# -------------------------
echo
echo "[2/3] Preparing compare output folder..."

if [[ -d "${COMPARE_OUT}" ]]; then
    BACKUP="${COMPARE_OUT}_old_$(date +%Y%m%d_%H%M%S)"
    echo "Existing compare output folder found. Moving it to:"
    echo "  ${BACKUP}"
    mv "${COMPARE_OUT}" "${BACKUP}"
fi

# -------------------------
# 7. Run inStrain compare with STB
# -------------------------
echo
echo "[3/3] Running inStrain compare with STB..."

mapfile -t PROFILE_DIRS < "${COMPLETED_PROFILE_LIST}"

inStrain compare \
    -i "${PROFILE_DIRS[@]}" \
    -o "${COMPARE_OUT}" \
    -p "${THREADS}" \
    -s "${STB_FILE}"

echo
echo "[+] inStrain compare finished."

# -------------------------
# 8. Report output files
# -------------------------
echo
echo "All compare output files:"
find "${COMPARE_OUT}" -type f | sort || true

echo
echo "Possible genome-level / MAG-level output files:"
find "${COMPARE_OUT}" -type f | grep -Ei "genome|genomeWide|wide|compare|comparisons|tsv" || true

echo
echo "Preview TSV headers:"
for f in "${COMPARE_OUT}"/output/*.tsv "${COMPARE_OUT}"/*.tsv; do
    [[ -s "${f}" ]] || continue
    echo "==== ${f} ===="
    head -n 1 "${f}"
done

echo
echo "Done."
echo "RFMT group: ${IND}"
echo "Compare output:"
echo "  ${COMPARE_OUT}"
echo "End time: $(date)"
echo "=========================================="