#!/usr/bin/env bash
#SBATCH --job-name=merge_lanes_fast_move
#SBATCH --output=/group/sbms004/yxia/GUT/logs/merge_lanes_fast_move_%j.out
#SBATCH --error=/group/sbms004/yxia/GUT/logs/merge_lanes_fast_move_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00

source activate gut

set -euo pipefail
shopt -s nullglob

# ─── CONFIG ─────────────────────────────────────────────────────────────
INPUT_DIR="/group/sbms004/yxia/GUT/AGRF_download"
OUTPUT_DIR="/group/sbms004/yxia/GUT/merged_fastq"
RAW_BACKUP_DIR="/group/sbms004/yxia/GUT/raw_moved_after_successful_merge"

SAMPLES_PER_RUN=19
PARALLEL_JOBS=19

# 1 = only print what would be moved
# 0 = actually move raw files to backup
DRY_RUN="${DRY_RUN:-0}"
# ───────────────────────────────────────────────────────────────────────

mkdir -p "$OUTPUT_DIR"
mkdir -p "$RAW_BACKUP_DIR"
mkdir -p /group/sbms004/yxia/GUT/logs

echo "=========================================="
echo "Fast merge + move script"
echo "Input dir: $INPUT_DIR"
echo "Output dir: $OUTPUT_DIR"
echo "Raw backup dir: $RAW_BACKUP_DIR"
echo "Samples per run: $SAMPLES_PER_RUN"
echo "Parallel jobs: $PARALLEL_JOBS"
echo "DRY_RUN: $DRY_RUN"
echo "=========================================="

if [[ "$DRY_RUN" == "1" ]]; then
    echo "DRY_RUN=1: raw files will NOT be moved."
else
    echo "DRY_RUN=0: raw files will be moved to backup after successful merge."
fi

# ─── Get current remaining sample names from INPUT_DIR ─────────────────
# Important:
# After first batch, INPUT_DIR has changed.
# This script always takes the first 19 samples currently remaining.
mapfile -t ALL_SAMPLES < <(
    for f in "$INPUT_DIR"/*_R1.fastq.gz; do
        basename "$f" | cut -d'_' -f1
    done | sort -u
)

TOTAL_SAMPLES=${#ALL_SAMPLES[@]}

echo "Remaining samples found in INPUT_DIR: $TOTAL_SAMPLES"

if (( TOTAL_SAMPLES == 0 )); then
    echo "No remaining R1 files found in $INPUT_DIR"
    echo "Nothing to do."
    exit 0
fi

mapfile -t SAMPLES < <(
    printf "%s\n" "${ALL_SAMPLES[@]}" |
    head -n "$SAMPLES_PER_RUN"
)

echo "Samples selected for this run:"
printf "%s\n" "${SAMPLES[@]}"

export INPUT_DIR OUTPUT_DIR RAW_BACKUP_DIR DRY_RUN

check_sample_id_match() {
    local sample="$1"
    shift

    local f
    local base
    local detected_sample

    for f in "$@"; do
        base="$(basename "$f")"
        detected_sample="$(echo "$base" | cut -d'_' -f1)"

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

export -f check_sample_id_match

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

    # ─── Check 1: R1/R2 must both exist ────────────────────────────────
    if (( ${#r1_lanes[@]} == 0 )); then
        echo "ERROR: No R1 files found for sample: $sample"
        return 1
    fi

    if (( ${#r2_lanes[@]} == 0 )); then
        echo "ERROR: No R2 files found for sample: $sample"
        return 1
    fi

    # ─── Check 2: R1/R2 lane count must match ─────────────────────────
    if (( ${#r1_lanes[@]} != ${#r2_lanes[@]} )); then
        echo "ERROR: R1/R2 lane count mismatch for sample: $sample"
        echo "R1 files: ${#r1_lanes[@]}"
        printf '  %s\n' "${r1_lanes[@]}"
        echo "R2 files: ${#r2_lanes[@]}"
        printf '  %s\n' "${r2_lanes[@]}"
        return 1
    fi

    # ─── Check 3: sample ID must match exactly before merge ────────────
    check_sample_id_match "$sample" "${r1_lanes[@]}" "${r2_lanes[@]}"

    # ─── Check 4: do not overwrite existing outputs ───────────────────
    if [[ -e "$r1_out" || -e "$r2_out" ]]; then
        echo "ERROR: Output already exists for sample: $sample"
        [[ -e "$r1_out" ]] && echo "Existing: $r1_out"
        [[ -e "$r2_out" ]] && echo "Existing: $r2_out"
        echo "Raw files will NOT be moved."
        return 1
    fi

    rm -f "$r1_tmp" "$r2_tmp"

    echo "[$(date '+%F %T')] R1 input files for $sample:"
    printf '  %s\n' "${r1_lanes[@]}"

    echo "[$(date '+%F %T')] R2 input files for $sample:"
    printf '  %s\n' "${r2_lanes[@]}"

    # ─── Merge R1 ─────────────────────────────────────────────────────
    echo "[$(date '+%F %T')] Merging R1 for sample: $sample"

    if command -v pigz >/dev/null 2>&1; then
        if ! zcat "${r1_lanes[@]}" | pigz -p 1 > "$r1_tmp"; then
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

    # ─── Merge R2 ─────────────────────────────────────────────────────
    echo "[$(date '+%F %T')] Merging R2 for sample: $sample"

    if command -v pigz >/dev/null 2>&1; then
        if ! zcat "${r2_lanes[@]}" | pigz -p 1 > "$r2_tmp"; then
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

    # ─── Check 5: tmp outputs must not be empty ───────────────────────
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

    # ─── Move tmp files to final outputs ──────────────────────────────
    mv "$r1_tmp" "$r1_out"
    mv "$r2_tmp" "$r2_out"

    if [[ ! -s "$r1_out" || ! -s "$r2_out" ]]; then
        echo "ERROR: Final output is empty for sample: $sample"
        return 1
    fi

    echo "[$(date '+%F %T')] Merge successful for sample: $sample"
    echo "Final R1: $r1_out"
    echo "Final R2: $r2_out"

    # ─── Final sample ID check before moving raw files ────────────────
    check_sample_id_match "$sample" "${r1_lanes[@]}" "${r2_lanes[@]}"

    # ─── Move raw files only after BOTH R1 and R2 succeeded ───────────
    local sample_backup_dir="$RAW_BACKUP_DIR/$sample"
    mkdir -p "$sample_backup_dir"

    if [[ "$DRY_RUN" == "1" ]]; then
        echo "[$(date '+%F %T')] DRY_RUN=1: raw files will NOT be moved."
        echo "Files that would be moved:"
        printf '  %s\n' "${r1_lanes[@]}" "${r2_lanes[@]}"
    else
        echo "[$(date '+%F %T')] Moving raw files to:"
        echo "$sample_backup_dir"

        mv "${r1_lanes[@]}" "${r2_lanes[@]}" "$sample_backup_dir/"

        echo "[$(date '+%F %T')] Raw files moved safely for sample: $sample"
    fi

    echo "[$(date '+%F %T')] Finished sample: $sample"
}

export -f merge_sample

# ─── Run selected 19 remaining samples in parallel ─────────────────────
printf "%s\n" "${SAMPLES[@]}" | xargs -P "$PARALLEL_JOBS" -I {} bash -c 'merge_sample "$@"' _ {}

echo "=========================================="
echo "✅ Done! This batch completed."
echo "Merged files are in: $OUTPUT_DIR"
echo "Raw backup folder: $RAW_BACKUP_DIR"
echo "=========================================="