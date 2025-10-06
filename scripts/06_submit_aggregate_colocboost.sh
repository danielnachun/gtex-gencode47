#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/aggregate_colocboost.sh"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

num_tissues=$(wc -l < "${tissue_id_list}")
echo "Number of tissues: ${num_tissues}"

if [ "${num_tissues}" -eq 0 ]; then
    echo "No tissues to aggregate"
    exit 0
fi

sbatch_params=(
    --output "${output_dir}/logs/aggregate_colocboost/%A_%a.log"
    --error "${output_dir}/logs/aggregate_colocboost/%A_%a.log"
    --array "1-${num_tissues}%250"
    --time 12:00:00
    --cpus-per-task 8
    --mem 128G
    --job-name aggregate_colocboost
    ${code_dir}/06_aggregate_colocboost.sh \
        --tissue_id_list ${tissue_id_list} \
        --coloc_base_dir ${coloc_base_dir} \
        --code_dir ${code_dir}
)


# Submit on either sherlock or scg
if [ "${submit_on}" = 'sherlock' ]; then
    # Additional parameters for sherlock
    sbatch \
        --partition normal,owners \
        --tmp 200G \
        "${sbatch_params[@]}" 
elif [ "${submit_on}" = 'scg' ]; then
    # Additional parameters for scg
    sbatch \
        --account smontgom \
        --partition batch \
        --constraint="nvme" \
        "${sbatch_params[@]}" 
else
    echo "must submit on either 'sherlock' or 'scg'"
fi
