#!/usr/bin/env bash

set -o nounset -o errexit

# Source the config file
CONFIG_FILE="${1:-}"
if [[ -z "${CONFIG_FILE}" ]]; then
    echo "Usage: $(basename "$0") <config_file>"
    exit 1
fi
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Create completion directory
mkdir -p "${output_dir}/completed/extract_participant_vcfs"

# Submit batch jobs: helper will filter completed items and divide into batches
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    --config_file "${CONFIG_FILE}" \
    --batch_script "${code_dir}/utils/batch_process_files.sh" \
    --job_name "extract_participant_vcfs" \
    --items_list_file "${participant_id_list}" \
    --processing_script "${code_dir}/01_extract_participant_vcfs.sh" \
    --file_param "--participant_id" \
    --completion_dir "${output_dir}/completed/extract_participant_vcfs" \
    --regenerate_all "${regenerate_all:-false}" \
    -- \
    --full_vcf "${full_vcf}" \
    --output_dir "${output_dir}" \
    --code_dir "${code_dir}"
