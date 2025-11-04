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
    "call_mutect" \
    "${code_dir}/03.B_call_mutect.sh" \
    --reference_dir "${reference_dir}" \
    --output_dir "${output_dir}" \
    --code_dir "${code_dir}" \
    --vcf_dir "${vcf_dir}" \
    --reference_fasta "${reference_fasta}" \
    --dbsnp "${dbsnp}" \
    --indels_mills "${indels_mills}" \
    --indels_decoy "${indels_decoy}" \
    --gene_intervals_bed "${gene_intervals_bed}" \
    --full_vcf_file "${full_vcf_file}" \
    --exac_reference "${exac_reference}"
