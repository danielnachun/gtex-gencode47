#!/usr/bin/env bash

set -o nounset -o errexit

# Source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/06_aggregate_colocboost_gwas.sh"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Call the shared batch submission utility
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    "${CONFIG_FILE}" \
    "${code_dir}/utils/batch_process_files.sh" \
    "aggregate_colocboost" \
    "${code_dir}/06_aggregate_colocboosts.sh" \
    --tissue_id_list "${tissue_id_list}" \
    --coloc_base_dir "${coloc_base_dir}" \
    --code_dir "${code_dir}" \
    --regenerate "${regenerate}" \
    --skip_robust "${skip_robust}"
