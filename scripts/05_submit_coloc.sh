#!/usr/bin/env bash

set -o nounset -o errexit

# Unified submission script for all co-localization analyses
# Usage: 05_submit_coloc.sh <config_file>

# Source config file
CONFIG_FILE="${1:-/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/05_coloc_all.sh}"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Determine output directory based on boolean flags
if [ "${multi_tissue:-FALSE}" = "TRUE" ] || [ "${multi_tissue:-FALSE}" = "true" ]; then
    if [ "${run_separate_gwas:-FALSE}" = "TRUE" ] || [ "${run_separate_gwas:-FALSE}" = "true" ]; then
        output_dir="${output_dir}/multi_tissue_gwas"
    elif [ "${run_joint_gwas:-FALSE}" = "TRUE" ] || [ "${run_joint_gwas:-FALSE}" = "true" ]; then
        output_dir="${output_dir}/multi_tissue_gwas"
    elif [ "${run_xqtl_only:-FALSE}" = "TRUE" ] || [ "${run_xqtl_only:-FALSE}" = "true" ]; then
        output_dir="${output_dir}/multi_tissue_qtl_only"
    else
        echo "Error: Multi-tissue mode requires either run_separate_gwas=TRUE, run_joint_gwas=TRUE, or run_xqtl_only=TRUE"
        exit 1
    fi
else
    if [ "${run_separate_gwas:-FALSE}" = "TRUE" ] || [ "${run_separate_gwas:-FALSE}" = "true" ]; then
        if [ "${strong_only:-FALSE}" = "TRUE" ] || [ "${strong_only:-FALSE}" = "true" ]; then
            output_dir="${output_dir}/single_tissue_gwas_strong_only"
        else
            output_dir="${output_dir}/single_tissue_gwas"
        fi
    elif [ "${run_joint_gwas:-FALSE}" = "TRUE" ] || [ "${run_joint_gwas:-FALSE}" = "true" ]; then
        if [ "${strong_only:-FALSE}" = "TRUE" ] || [ "${strong_only:-FALSE}" = "true" ]; then
            output_dir="${output_dir}/single_tissue_gwas_strong_only"
        else
            output_dir="${output_dir}/single_tissue_gwas"
        fi
    elif [ "${run_xqtl_only:-FALSE}" = "TRUE" ] || [ "${run_xqtl_only:-FALSE}" = "true" ]; then
        output_dir="${output_dir}/single_tissue_qtl_only"
    else
        echo "Error: Single-tissue mode requires either run_separate_gwas=TRUE, run_joint_gwas=TRUE, or run_xqtl_only=TRUE"
        exit 1
    fi
fi

# Generate list of all LD regions to process
# Convert BED file (0-based, half-open) to region strings (chr:start-end, 1-based inclusive)
ld_region_list_file="${output_dir}/file_lists/file_list_coloc_ld_blocks.txt"
mkdir -p "$(dirname "${ld_region_list_file}")"
if [ ! -f "${ld_region_list_file}" ]; then
    awk 'NR>1 {printf "%s:%d-%d\n", $1, $2+1, $3}' "${ld_region_list}" > "${ld_region_list_file}"
    echo "Created LD region list file: ${ld_region_list_file}"
fi

# Create completion directory (using job name)
job_name="coloc_$(basename "${output_dir}")"
mkdir -p "${output_dir}/completed/${job_name}"

# Submit batch jobs: helper will normalize regions for completion checking, filter completed, and divide into batches
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    --config_file "${CONFIG_FILE}" \
    --batch_script "${code_dir}/utils/batch_process_files.sh" \
    --job_name "${job_name}" \
    --items_list_file "${ld_region_list_file}" \
    --processing_script "${code_dir}/05_colocalize_regions.sh" \
    --file_param "--ld_region" \
    --completion_dir "${output_dir}/completed/${job_name}" \
    --regenerate_all "${regenerate_all:-false}" \
    -- \
    --ld_region_list "${ld_region_list_file}" \
    --tissue_id_list "${tissue_id_list}" \
    --genotype_stem "${genotype_stem}" \
    --covariate_dir "${covariate_dir}" \
    --expression_dir "${expression_dir}" \
    --all_v39_genes_path "${all_v39_genes_path}" \
    --region_padding "${region_padding}" \
    --output_dir "${output_dir}" \
    --code_dir "${code_dir}" \
    --multi_tissue "${multi_tissue:-FALSE}" \
    --run_v39 "${run_v39:-FALSE}" \
    --run_individual "${run_individual:-FALSE}" \
    --run_xqtl_only "${run_xqtl_only:-FALSE}" \
    --run_separate_gwas "${run_separate_gwas:-FALSE}" \
    --run_joint_gwas "${run_joint_gwas:-FALSE}" \
    ${gwas_params-}
