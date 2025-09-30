#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/coloc_single_tissue.sh"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }



# List all LD regions that do not have a completion file

completion_dir="${output_dir}/${tissue_id}/completed/single_tissue_coloc"
mkdir -p "${completion_dir}"

# Extract all LD regions from the BED file (skip header), convert to 1-based region string
ld_regions=$(awk 'NR>1 {printf "%s:%d-%d\n", $1, $2+1, $3}' "${ld_region_list}")
total_count=$(wc -l < "${ld_region_list}")
regenerate_all=${regenerate_all:-false}
if [ "${regenerate_all}" = true ]; then
    # run all the bams in the input folder
    missing_ld_regions=${ld_regions}
else
    # only coloc a block if the completion marker file does not already exist
    all_completion_files=$(printf "%s\n" ${ld_regions} | sed "s|^|${completion_dir}/${tissue_id}.|;s|\$|.completed|")
    missing_ld_regions=$(paste <(printf "%s\n" ${ld_regions}) <(printf "%s\n" ${all_completion_files}) | awk '{if(system("[ -f \""$2"\" ]")==0) next; print $1}')
fi

# Write missing LD regions to a file
missing_ld_regions_file="${output_dir}/file_lists/file_list_coloc_${tissue_id}_ld_blocks.txt"
mkdir -p $(dirname "${missing_ld_regions_file}")
printf "%s\n" ${missing_ld_regions} > "${missing_ld_regions_file}"
num_ld_blocks=$(wc -l < "${missing_ld_regions_file}")

echo "Total LD regions: ${total_count}"
echo "To be processed: ${num_ld_blocks}"
# Check if num_ld_blocks is greater than max_array_size
if [ "${num_ld_blocks}" -gt "${max_array_size}" ]; then
    num_ld_blocks="${max_array_size}" 
fi
echo "Submitting: ${num_ld_blocks}"


sbatch_params=(
    --output "${output_dir}/logs/coloc/%A_%a.log"
    --error "${output_dir}/logs/coloc/%A_%a.log"
    --array "1-${num_ld_blocks}%250"
    --time 2:00:00
    --cpus-per-task 4
    --mem 64G
    --job-name coloc_${tissue_id}_ld_blocks
    ${code_dir}/colocalize_regions_single_tissue.sh
        --ld_region_list ${missing_ld_regions_file}
        --tissue_id ${tissue_id}
        --genotype_stem ${genotype_stem}
        --covariate_path ${covariate_dir_path}/${tissue_id}.v11.covariates.txt
        --expression_path ${expression_dir_path}/${tissue_id}.v11.normalized_expression.bed.gz
        --all_v39_genes_path ${all_v39_genes_path}
        --region_padding ${region_padding}
        --association_padding ${association_padding}
        --output_dir ${output_dir}/${tissue_id}
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