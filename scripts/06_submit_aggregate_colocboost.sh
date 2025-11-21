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
mkdir -p "${coloc_base_dir}/completed/aggregate_colocboost"

# Submit batch jobs: helper will filter completed items and divide into batches
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    --config_file "${CONFIG_FILE}" \
    --batch_script "${code_dir}/utils/batch_process_files.sh" \
    --job_name "aggregate_colocboost" \
    --items_list_file "${tissue_id_list}" \
    --processing_script "${code_dir}/06_aggregate_colocboosts.sh" \
    --file_param "--tissue_id" \
    --completion_dir "${coloc_base_dir}/completed/aggregate_colocboost" \
    --regenerate_all "${regenerate:-false}" \
    -- \
    --tissue_id_list "${tissue_id_list}" \
    --coloc_base_dir "${coloc_base_dir}" \
    --code_dir "${code_dir}" \
    --regenerate "${regenerate}" \
    --skip_robust "${skip_robust}"
