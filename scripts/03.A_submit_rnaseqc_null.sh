#!/usr/bin/env bash

set -o nounset -o errexit

# Source the config file
CONFIG_FILE="${1:-}"
if [[ -z "${CONFIG_FILE}" ]]; then
    echo "Usage: $(basename "$0") <config_file>"
    exit 1
fi
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Build optional args
extra_args=()
if [[ -n "${reference_fasta:-}" ]]; then
    extra_args+=(--reference_fasta "${reference_fasta}")
fi
if [[ -n "${chr_sizes:-}" ]]; then
    extra_args+=(--chr_sizes "${chr_sizes}")
fi
if [[ -n "${intervals_bed:-}" ]]; then
    extra_args+=(--intervals_bed "${intervals_bed}")
fi

# Call the shared batch submission utility
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    "${CONFIG_FILE}" \
    "${code_dir}/utils/batch_process_files.sh" \
    "rnaseqc_null" \
    "${code_dir}/03.A_quantify_rnaseqc_null.sh" \
    --local_reference_dir "${reference_dir}" \
    --output_dir "${output_dir}" \
    --code_dir "${code_dir}" \
    --genes_gtf "${genes_gtf}" \
    "${extra_args[@]}" \
    --bam_file_end "${bam_file_end}"
