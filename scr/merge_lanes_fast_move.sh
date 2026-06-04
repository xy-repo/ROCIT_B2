#!/usr/bin/env bash
#SBATCH --job-name=merge_lanes_array
#SBATCH --output=/group/sbms004/yxia/GUT/logs/merge_lanes_%A_%a.out
#SBATCH --error=/group/sbms004/yxia/GUT/logs/merge_lanes_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=1-300

source activate gut

set -euo pipefail
shopt -s nullglob

# ─── CONFIG ─────────────────────────────────────────────────────────────
INPUT_DIR="/group/sbms004/yxia/GUT/AGRF_download"
OUTPUT_DIR="/group/sbms004/yxia/GUT/merged_fastq"
LOG_DIR="/group/sbms004/yxia/GUT/logs"

THREADS="${SLURM_CPUS_PER_TASK:-4}"
# ───────────────────────────────────────────────────────────────────────

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

echo "=========================================="
echo "Merge lanes array job"
echo "Input dir: $INPUT_DIR"
echo "Output dir: $OUTPUT_DIR"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Threads per sample: $THREADS"
echo "Original raw files will NOT be moved or deleted."
echo "=========================================="

# ─── Build sample list ─────────────────────────────────────────────────
# Example filename:
#   RMFT001_3521_xxx_L001_R1.fastq.gz
#
# Sample ID:
#   RMFT001_3521

mapfile -t ALL_SAMPLES < <(
    for f in "$INPUT_DIR"/*_R1.fastq.gz; do
        basename "$f" | awk -F'_' '{print $1"_"$2}'
    done | sort -u
)

TOTAL_SAMPLES=${#ALL_SAMPLES[@]}

echo "Total samples found: $TOTAL_SAMPLES"

if (( TOTAL_SAMPLES == 0 )); then
    echo "No R1 files found in $INPUT_DIR"
    exit 0
fi

TASK_ID="${SLURM_ARRAY_TASK_ID}"

if (( TASK_ID > TOTAL_SAMPLES )); then
    echo "Array task $TASK_ID is greater than total samples $TOTAL_SAMPLES"
    echo "Nothing to do."
    exit 0
fi

sample="${ALL_SAMPLES[$((TASK_ID - 1))]}"

echo "Selected sample: $sample"

check_sample_id_match() {
    local sample="$1"
    shift

    local f
    local base
    local detected_sample

    for f in "$@"; do
        base="$(basename "$f")"
        detected_sample="$(echo "$base" | awk -F'_' '{print $1"_"$2}')"

        if [[ "$detected_sample" != "$sample" ]]; then
            echo "ERROR: Sample ID mismatch!"
            echo "Expected sample: $sample"
            echo "Detected sample: $detected_sample"
            echo "File: $f"
            return 1
        fi

        if [[ "$base" != "${sample}_"* ]]; then
            echo "ERROR: File does not start with exact sample prefix."
            echo "Expected prefix: ${sample}_"
            echo "File: $f"
            return 1
        fi
    done
}

merge_sample() {
    local sample="$1"

    echo "[$(date '+%F %T')] Starting sample: $sample"

    local r1_lanes=()
    local r2_lanes=()

    r1_lanes=( "$INPUT_DIR/${sample}"_*_R1.fastq.gz )
    r2_lanes=( "$INPUT_DIR/${sample}"_*_R2.fastq.gz )

    local r1_out="$OUTPUT_DIR/${sample}_R1.fastq.gz"
    local r2_out="$OUTPUT_DIR/${sample}_R2.fastq.gz"

    local r1_tmp="${r1_out}.tmp"
    local r2_tmp="${r2_out}.tmp"

    if (( ${#r1_lanes[@]} == 0 )); then
        echo "ERROR: No R1 files found for sample: $sample"
        return 1
    fi

    if (( ${#r2_lanes[@]} == 0 )); then
        echo "ERROR: No R2 files found for sample: $sample"
        return 1
    fi

    if (( ${#r1_lanes[@]} != ${#r2_lanes[@]} )); then
        echo "ERROR: R1/R2 lane count mismatch for sample: $sample"
        echo "R1 files: ${#r1_lanes[@]}"
        printf '  %s\n' "${r1_lanes[@]}"
        echo "R2 files: ${#r2_lanes[@]}"
        printf '  %s\n' "${r2_lanes[@]}"
        return 1
    fi

    check_sample_id_match "$sample" "${r1_lanes[@]}" "${r2_lanes[@]}"

    if [[ -s "$r1_out" && -s "$r2_out" ]]; then
        echo "[$(date '+%F %T')] Output already exists for sample: $sample"
        echo "Existing R1: $r1_out"
        echo "Existing R2: $r2_out"
        echo "Skipping sample."
        return 0
    fi

    if [[ -e "$r1_out" || -e "$r2_out" ]]; then
        echo "ERROR: Only one output exists or one output is empty for sample: $sample"
        [[ -e "$r1_out" ]] && echo "Existing R1: $r1_out"
        [[ -e "$r2_out" ]] && echo "Existing R2: $r2_out"
        echo "Please check manually before rerunning."
        return 1
    fi

    rm -f "$r1_tmp" "$r2_tmp"

    echo "[$(date '+%F %T')] R1 input files:"
    printf '  %s\n' "${r1_lanes[@]}"

    echo "[$(date '+%F %T')] R2 input files:"
    printf '  %s\n' "${r2_lanes[@]}"

    echo "[$(date '+%F %T')] Merging R1 for sample: $sample"

    if command -v pigz >/dev/null 2>&1; then
        if ! zcat "${r1_lanes[@]}" | pigz -p "$THREADS" > "$r1_tmp"; then
            echo "ERROR: R1 merge failed for sample: $sample"
            rm -f "$r1_tmp" "$r2_tmp"
            return 1
        fi
    else
        if ! zcat "${r1_lanes[@]}" | gzip > "$r1_tmp"; then
            echo "ERROR: R1 merge failed for sample: $sample"
            rm -f "$r1_tmp" "$r2_tmp"
            return 1
        fi
    fi

    echo "[$(date '+%F %T')] Merging R2 for sample: $sample"

    if command -v pigz >/dev/null 2>&1; then
        if ! zcat "${r2_lanes[@]}" | pigz -p "$THREADS" > "$r2_tmp"; then
            echo "ERROR: R2 merge failed for sample: $sample"
            rm -f "$r1_tmp" "$r2_tmp"
            return 1
        fi
    else
        if ! zcat "${r2_lanes[@]}" | gzip > "$r2_tmp"; then
            echo "ERROR: R2 merge failed for sample: $sample"
            rm -f "$r1_tmp" "$r2_tmp"
            return 1
        fi
    fi

    if [[ ! -s "$r1_tmp" ]]; then
        echo "ERROR: R1 tmp output is empty for sample: $sample"
        rm -f "$r1_tmp" "$r2_tmp"
        return 1
    fi

    if [[ ! -s "$r2_tmp" ]]; then
        echo "ERROR: R2 tmp output is empty for sample: $sample"
        rm -f "$r1_tmp" "$r2_tmp"
        return 1
    fi

    echo "[$(date '+%F %T')] Checking gzip integrity"

    if ! gzip -t "$r1_tmp"; then
        echo "ERROR: R1 tmp gzip check failed for sample: $sample"
        rm -f "$r1_tmp" "$r2_tmp"
        return 1
    fi

    if ! gzip -t "$r2_tmp"; then
        echo "ERROR: R2 tmp gzip check failed for sample: $sample"
        rm -f "$r1_tmp" "$r2_tmp"
        return 1
    fi

    mv "$r1_tmp" "$r1_out"
    mv "$r2_tmp" "$r2_out"

    if [[ ! -s "$r1_out" || ! -s "$r2_out" ]]; then
        echo "ERROR: Final output is empty for sample: $sample"
        return 1
    fi

    echo "[$(date '+%F %T')] Merge successful for sample: $sample"
    echo "Final R1: $r1_out"
    echo "Final R2: $r2_out"
    echo "Original raw files kept in: $INPUT_DIR"
}

merge_sample "$sample"

echo "=========================================="
echo "Done."
echo "=========================================="