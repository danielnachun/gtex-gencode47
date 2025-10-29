#!/usr/bin/env bash

set -o nounset -o errexit

# Source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/03.A_quantify_null_v11.sh"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Call the shared batch submission utility
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    "${CONFIG_FILE}" \
    "${code_dir}/utils/batch_process_files.sh" \
    "rnaseqc_null" \
    "${code_dir}/03.A_quantify_rnaseqc_null.sh" \
    --reference_dir "${reference_dir}" \
    --output_dir "${output_dir}" \
    --code_dir "${code_dir}" \
    --reference_fasta "${reference_fasta}" \
    --chr_sizes "${chr_sizes}" \
    --genes_gtf "${genes_gtf}" \
    --intervals_bed "${intervals_bed}" \
    --bam_file_end "${bam_file_end}"
