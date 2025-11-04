#!/usr/bin/env bash

set -o nounset -o errexit

# Source the config file
CONFIG_FILE="${1:-}"
if [[ -z "${CONFIG_FILE}" ]]; then
    echo "Usage: $(basename "$0") <config_file>"
    exit 1
fi
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Call the shared batch submission utility
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    "${CONFIG_FILE}" \
    "${code_dir}/utils/batch_process_files.sh" \
    "quantify_bam" \
    "${code_dir}/03_quantify_bam.sh" \
    --reference_dir "${reference_dir}" \
    --vcf_dir "${vcf_dir}" \
    --output_dir "${output_dir}" \
    --code_dir "${code_dir}" \
    --reference_fasta "${reference_fasta}" \
    --genes_gtf "${genes_gtf}" \
    --chr_sizes "${chr_sizes}" \
    --intervals_bed "${intervals_bed}" \
    --ipa_annotation "${ipa_annotation}" \
    --editing_bed "${editing_bed}"
