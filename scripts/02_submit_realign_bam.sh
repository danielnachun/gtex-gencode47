#!/usr/bin/env bash

set -o nounset -o errexit

# Source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/02_realign_truncated.sh"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Call the shared batch submission utility
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    "${CONFIG_FILE}" \
    "${code_dir}/utils/batch_process_files.sh" \
    "realign_bam" \
    "${code_dir}/02_realign_bam.sh" \
    --reference_dir "${reference_dir}" \
    --vcf_dir "${vcf_dir}" \
    --output_dir "${output_dir}" \
    --code_dir "${code_dir}" \
    --reference_fasta "${reference_fasta}" \
    --rsem_ref_dir "${rsem_ref_dir}" \
    --star_index "${star_index}"
