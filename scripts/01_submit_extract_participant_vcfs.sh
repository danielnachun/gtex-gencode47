#!/usr/bin/env bash

set -o nounset -o errexit

# Source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/01_extract_vcfs_all.sh"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Call the shared batch submission utility
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    "${CONFIG_FILE}" \
    "${code_dir}/utils/batch_process_files.sh" \
    "extract_participant_vcfs" \
    "${code_dir}/01_extract_participant_vcfs.sh" \
    --participant_id_list "${output_dir}/participants_to_process.txt" \
    --full_vcf "${full_vcf}" \
    --output_dir "${output_dir}" \
    --code_dir "${code_dir}"
