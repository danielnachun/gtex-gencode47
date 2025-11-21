#!/usr/bin/env bash

set -o nounset -o errexit

# Source the config file
CONFIG_FILE="${1:-}"
if [[ -z "${CONFIG_FILE}" ]]; then
    echo "Usage: $(basename "$0") <config_file>"
    exit 1
fi
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Generate list of all LD regions to process
# Convert BED file (0-based, half-open) to region strings (chr:start-end, 1-based inclusive)
ld_region_list_file="${out_dir}/file_lists/ld_region_list.txt"
mkdir -p "$(dirname "${ld_region_list_file}")"
if [ ! -f "${ld_region_list_file}" ]; then
    awk 'NR>1 {printf "%s:%d-%d\n", $1, $2+1, $3}' "${ld_blocks_bed}" > "${ld_region_list_file}"
    echo "Created LD region list file: ${ld_region_list_file}"
fi

# Create completion directory
mkdir -p "${out_dir}/completed/generate_ld"

# Submit batch jobs: helper will normalize regions for completion checking, filter completed, and divide into batches
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    --config_file "${CONFIG_FILE}" \
    --batch_script "${code_dir}/utils/batch_process_files.sh" \
    --job_name "generate_ld" \
    --items_list_file "${ld_region_list_file}" \
    --processing_script "${code_dir}/04_generate_ld.sh" \
    --file_param "--region" \
    --completion_dir "${out_dir}/completed/generate_ld" \
    --regenerate_all "${regenerate_all:-false}" \
    -- \
    --genotype_prefix "${genotype_prefix}" \
    --sample_ids "${out_dir}/european_subject_ids.keep.txt" \
    --output_dir "${out_dir}" \
    --code_dir "${code_dir}" \
    --regenerate "${regenerate_all}"