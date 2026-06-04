#!/usr/bin/env bash

set -eo pipefail
shopt -s nullglob

SAMPLE_LIST_DIR="/group/sbms004/yxia/GUT/sample_list"
PROJECT_ROOT="/group/sbms004/yxia/GUT/MAG_per_sample"

REPORT="${PROJECT_ROOT}/squeezemeta_syslog_check.tsv"

echo -e "sample_list\tsample\tproject_dir\tsyslog\tbins_dir\tstatus\tnote" > "${REPORT}"

echo "Checking SqueezeMeta syslog completion..."
echo "Report:"
echo "${REPORT}"
echo

for SAMPLE_LIST in "${SAMPLE_LIST_DIR}"/sample_list_*.tsv; do
    [[ -s "${SAMPLE_LIST}" ]] || continue

    list_base="$(basename "${SAMPLE_LIST}" .tsv)"
    PROJECT_GROUP_DIR="${PROJECT_ROOT}/${list_base}"

    echo "Checking ${list_base}..."

    mapfile -t SAMPLES < <(
        awk -F'\t' 'NF >= 2 && $1 != "" {print $1}' "${SAMPLE_LIST}" | sort -u
    )

    for sample in "${SAMPLES[@]}"; do
        PROJECT_DIR="${PROJECT_GROUP_DIR}/${sample}/${sample}"
        SYSLOG="${PROJECT_DIR}/syslog"
        BINS_DIR="${PROJECT_DIR}/results/bins"

        status="UNKNOWN"
        note=""

        if [[ ! -d "${PROJECT_DIR}" ]]; then
            status="NOT_RUN"
            note="project directory missing"

        elif [[ ! -s "${SYSLOG}" ]]; then
            status="NO_SYSLOG"
            note="syslog missing or empty"

        elif grep -q "INFO: Done." "${SYSLOG}" && \
             grep -q "INFO: Removing intermediate files." "${SYSLOG}" && \
             grep -q "INFO: Intermediate files removed." "${SYSLOG}"; then

            bins=( "${BINS_DIR}"/*.fa "${BINS_DIR}"/*.fna "${BINS_DIR}"/*.fasta )

            if [[ ! -d "${BINS_DIR}" ]]; then
                status="DONE_NO_BINS_DIR"
                note="syslog done, but bins directory missing"

            elif (( ${#bins[@]} == 0 )); then
                status="DONE_EMPTY_BINS"
                note="syslog done, but bins directory has no MAG fasta files"

            else
                status="DONE_OK"
                note="syslog done and MAG fasta files found: ${#bins[@]}"
            fi

        else
            status="NOT_DONE"
            note="syslog does not contain all completion lines"
        fi

        echo -e "${list_base}\t${sample}\t${PROJECT_DIR}\t${SYSLOG}\t${BINS_DIR}\t${status}\t${note}" >> "${REPORT}"
    done
done

echo
echo "Summary:"
cut -f6 "${REPORT}" | tail -n +2 | sort | uniq -c

echo
echo "Incomplete / problematic samples:"
awk -F'\t' 'NR == 1 || $6 != "DONE_OK"' "${REPORT}"

echo
echo "Done."
echo "Full report:"
echo "${REPORT}"