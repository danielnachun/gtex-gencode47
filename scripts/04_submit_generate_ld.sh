#!/usr/bin/env bash

set -o nounset -o errexit

# Source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/04_ld_gtex_eur.sh"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Call the shared batch submission utility
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    "${CONFIG_FILE}" \
    "${code_dir}/utils/batch_process_files.sh" \
    "generate_ld" \
    "${code_dir}/04_generate_ld.sh" \
    --genotype_prefix "${genotype_prefix}" \
    --sample_ids "${out_dir}/european_subject_ids.keep.txt" \
    --output_dir "${out_dir}" \
    --code_dir "${code_dir}" \
    --regenerate "${regenerate_all}"