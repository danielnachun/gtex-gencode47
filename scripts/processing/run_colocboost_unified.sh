#!/usr/bin/env bash

# Unified colocboost script for all analysis modes
# Usage: run_colocboost.sh [parameters...]

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o nounset -o pipefail -o errexit

# Helper functions
check_for_file() {
    argument_name="${1}"
    file_path="${2}"
    if [[ ${file_path} != "none" ]] && [[ ! -f ${file_path} ]]; then
        echo "Error: file ${file_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

check_for_directory() {
    argument_name="${1}"
    directory_path="${2}"
    if [[ ${directory_path} != "none" ]] && [[ ! -d ${directory_path} ]]; then
        echo "Error: directory ${directory_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

# Define all possible parameters
options_array=(
    # Core parameters
    tissue_id
    ld_region
    association_region
    gene_region
    gene_bed_list
    covariate_list
    covariate_path
    genotype_stem
    v39_gene_id_path
    output_dir
    code_dir
    # GWAS-specific parameters (optional)
    gwas_id_list
    gwas_phenotype_list
    gwas_meta
    ld_meta
    gwas_column_matching
    # Analysis mode flags
    multi_tissue
    run_v39
    run_individual
    run_xqtl_only
    run_separate_gwas
    run_joint_gwas
    # Quality control parameters
    maf_cutoff
    mac_cutoff
    xvar_cutoff
    imiss_cutoff
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'colocboost' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        # Core parameters
        --tissue_id ) tissue_id="${2}"; shift 2 ;;
        --ld_region ) ld_region="${2}"; shift 2 ;;
        --association_region ) association_region="${2}"; shift 2 ;;
        --gene_region ) gene_region="${2}"; shift 2 ;;
        --gene_bed_list ) gene_bed_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --covariate_list ) covariate_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --covariate_path ) covariate_path="${2}"; shift 2 ;;
        --genotype_stem ) genotype_stem="${2}"; shift 2 ;;
        --v39_gene_id_path ) v39_gene_id_path="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_dir ) output_dir="${2}"; shift 2 ;;
        --code_dir ) code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        # GWAS-specific parameters
        --gwas_id_list ) gwas_id_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --gwas_phenotype_list ) gwas_phenotype_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --gwas_meta ) gwas_meta="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --ld_meta ) ld_meta="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --gwas_column_matching ) gwas_column_matching="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        # Analysis mode flags
        --multi_tissue ) multi_tissue="${2}"; shift 2 ;;
        --run_v39 ) run_v39="${2}"; shift 2 ;;
        --run_individual ) run_individual="${2}"; shift 2 ;;
        --run_xqtl_only ) run_xqtl_only="${2}"; shift 2 ;;
        --run_separate_gwas ) run_separate_gwas="${2}"; shift 2 ;;
        --run_joint_gwas ) run_joint_gwas="${2}"; shift 2 ;;
        # Quality control parameters
        --maf_cutoff ) maf_cutoff="${2}"; shift 2 ;;
        --mac_cutoff ) mac_cutoff="${2}"; shift 2 ;;
        --xvar_cutoff ) xvar_cutoff="${2}"; shift 2 ;;
        --imiss_cutoff ) imiss_cutoff="${2}"; shift 2 ;;
        --) shift; break;;
        * ) echo "Invalid argument ${1} ${2}" >&2; exit 1
    esac
done

# Set default values
multi_tissue=${multi_tissue:-FALSE}
run_v39=${run_v39:-FALSE}
run_individual=${run_individual:-FALSE}
run_xqtl_only=${run_xqtl_only:-FALSE}
run_separate_gwas=${run_separate_gwas:-FALSE}
run_joint_gwas=${run_joint_gwas:-FALSE}
maf_cutoff=${maf_cutoff:-0.01}
mac_cutoff=${mac_cutoff:-10}
xvar_cutoff=${xvar_cutoff:-0}
imiss_cutoff=${imiss_cutoff:-0.9}

# Create output directory
mkdir -p "${output_dir}"

echo $(date +"[%b %d %H:%M:%S] Running colocboost on LD region ${ld_region} in tissue ${tissue_id} on the gene bed list ${gene_bed_list}")

# Determine which R script to call based on analysis mode
if [ "${run_xqtl_only}" = "TRUE" ] || [ "${run_xqtl_only}" = "true" ] || [ "${run_separate_gwas}" = "TRUE" ] || [ "${run_separate_gwas}" = "true" ] || [ "${run_joint_gwas}" = "TRUE" ] || [ "${run_joint_gwas}" = "true" ]; then
    # Use unified compare script for QTL and GWAS analyses
    r_script="${code_dir}/colocboost.compare.R"
else
    # Default to basic colocboost
    r_script="${code_dir}/colocboost.R"
fi

# Build the R script command
r_command="${r_script} \
    --tissue_id ${tissue_id} \
    --gene_region ${gene_region} \
    --variant_region ${association_region} \
    --ld_region ${ld_region} \
    --genotype_stem ${genotype_stem} \
    --phenotype_list ${gene_bed_list} \
    --covariate_list ${covariate_list} \
    --covariate_path ${covariate_path} \
    --output_dir ${output_dir} \
    --maf_cutoff ${maf_cutoff} \
    --mac_cutoff ${mac_cutoff} \
    --xvar_cutoff ${xvar_cutoff} \
    --imiss_cutoff ${imiss_cutoff} \
    --run_single_gene ${run_individual} \
    --run_v39_genes ${run_v39} \
    --v39_gene_id_path ${v39_gene_id_path}"

# Add GWAS-specific parameters if running GWAS analysis
if [ "${run_separate_gwas}" = "TRUE" ] || [ "${run_separate_gwas}" = "true" ] || [ "${run_joint_gwas}" = "TRUE" ] || [ "${run_joint_gwas}" = "true" ]; then
    r_command="${r_command} \
        --gwas_id_list ${gwas_id_list} \
        --gwas_phenotype_list ${gwas_phenotype_list} \
        --gwas_meta ${gwas_meta} \
        --ld_meta ${ld_meta} \
        --gwas_column_matching ${gwas_column_matching}"
fi

# Add analysis mode flags to R script
r_command="${r_command} \
    --multi_tissue ${multi_tissue} \
    --run_v39 ${run_v39} \
    --run_individual ${run_individual} \
    --run_xqtl_only ${run_xqtl_only} \
    --run_separate_gwas ${run_separate_gwas} \
    --run_joint_gwas ${run_joint_gwas}"

# Execute the R script
eval "${r_command}"

echo $(date +"[%b %d %H:%M:%S] Done with gene bed list ${gene_bed_list}")
