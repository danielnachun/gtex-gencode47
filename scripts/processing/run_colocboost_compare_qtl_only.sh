#!/usr/bin/env bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o nounset -o pipefail -o errexit

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

options_array=(
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
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'fraser' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --tissue_id )
            tissue_id="${2}"; shift 2 ;;
        --ld_region )
            ld_region="${2}"; shift 2 ;;
        --association_region )
            association_region="${2}"; shift 2 ;;
        --gene_region )
            gene_region="${2}"; shift 2 ;;
        --gene_bed_list )
            gene_bed_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --covariate_list )
            covariate_list="${2}"; shift 2 ;;
        --covariate_path )
            covariate_path="${2}"; shift 2 ;;
        --genotype_stem )
            genotype_stem="${2}"; shift 2 ;;
        --v39_gene_id_path )
            v39_gene_id_path="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

# Normalize optional covariate inputs and validate that at least one is provided
covariate_list="${covariate_list:-none}"
covariate_path="${covariate_path:-none}"

if [[ "${covariate_list}" == "none" && "${covariate_path}" == "none" ]]; then
    echo "Error: you must provide either --covariate_list or --covariate_path." >&2
    exit 1
fi

if [[ "${covariate_list}" != "none" ]]; then
    check_for_file "--covariate_list" "${covariate_list}"
fi

if [[ "${covariate_path}" != "none" ]]; then
    check_for_file "--covariate_path" "${covariate_path}"
fi

mkdir -p ${output_dir}

echo $(date +"[%b %d %H:%M:%S] Running colocboost on LD region ${ld_region} in tissue ${tissue_id} on the gene bed list ${gene_bed_list}")

${code_dir}/colocboost_compare_qtl_only.R \
    --tissue_id ${tissue_id} \
    --gene_region ${gene_region} \
    --variant_region ${association_region} \
    --ld_region ${ld_region} \
    --genotype_stem ${genotype_stem} \
    --phenotype_list ${gene_bed_list} \
    --covariate_list ${covariate_list} \
    --covariate_path ${covariate_path} \
    --output_dir ${output_dir} \
    --maf_cutoff 0.01 \
    --mac_cutoff 10 \
    --xvar_cutoff 0 \
    --imiss_cutoff 0.9 \
    --run_single_gene TRUE \
    --run_v39_genes TRUE \
    --v39_gene_id_path ${v39_gene_id_path}

echo $(date +"[%b %d %H:%M:%S] Done with gene bed list ${gene_bed_list}")

