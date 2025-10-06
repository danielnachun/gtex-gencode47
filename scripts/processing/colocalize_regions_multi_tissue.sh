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
    ld_region_list
    tissue_id_list
    gwas_id_list
    genotype_stem
    covariate_dir
    expression_dir
    gwas_dir
    gwas_meta
    ld_meta
    gwas_column_matching
    all_v39_genes_path
    region_padding
    association_padding
    output_dir
    code_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'fraser' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --ld_region_list )
            ld_region_list="${2}"; shift 2 ;;
        --tissue_id_list )
            tissue_id_list="${2}"; shift 2 ;;
        --gwas_id_list )
            gwas_id_list="${2}"; shift 2 ;;
        --genotype_stem )
            genotype_stem="${2}"; shift 2 ;;
        --covariate_dir )
            covariate_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --expression_dir )
            expression_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --gwas_dir )
            gwas_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --gwas_meta )
            gwas_meta="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --ld_meta )
            ld_meta="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --gwas_column_matching )
            gwas_column_matching="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --all_v39_genes_path )
            all_v39_genes_path="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --region_padding )
            region_padding="${2}"; shift 2 ;;
        --association_padding )
            association_padding="${2}"; shift 2 ;;
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

mkdir -p ${output_dir}

line_number="${SLURM_ARRAY_TASK_ID}"
working_dir="${TMPDIR}/coloc_${line_number}"
mkdir -p "${working_dir}"
rsync -PrhLtv "${genotype_stem}"* "${working_dir}/"
local_genotype_stem="${working_dir}/$(basename "${genotype_stem}")"


# activate the pixi enviroment
source <(pixi shell-hook --environment pecotmr-dev --manifest-path ${code_dir}/pixi.toml)

# map job id to line number and then to participant id
line_number=${SLURM_ARRAY_TASK_ID}
ld_region_bed_line=$(sed "${line_number}q; d" "${ld_region_list}")

# Parse the region from colon-separated format (e.g., "chr1:16104-1170341")
region_chr=$(echo "${ld_region_bed_line}" | cut -d':' -f1)
region_coords=$(echo "${ld_region_bed_line}" | cut -d':' -f2)
region_start=$(echo "${region_coords}" | cut -d'-' -f1)
region_end=$(echo "${region_coords}" | cut -d'-' -f2)
ld_region="${region_chr}:${region_start}-${region_end}"
echo $(date +"[%b %d %H:%M:%S] Running colocboost with and without v39 genes on region ${ld_region}")

# get region for variants
association_region_start=$((region_start - association_padding))
association_region_end=$((region_end + association_padding))
# Ensure padded region start is not less than 0
if [ "${association_region_start}" -lt 0 ]; then
    association_region_start=0
fi
association_region="${region_chr}:${association_region_start}-${association_region_end}"

# get regions for genes
region_start_padded=$((region_start - region_padding))
region_end_padded=$((region_end + region_padding))
# Ensure padded region start is not less than 0
if [ "${region_start_padded}" -lt 0 ]; then
    region_start_padded=0
fi
gene_region="${region_chr}:${region_start_padded}-${region_end_padded}"


# Create a list of GWAS phenotype file paths for this job's region
gwas_phenotype_list="${working_dir}/gwas_phenotype_list.txt"
rm -f "${gwas_phenotype_list}"
touch "${gwas_phenotype_list}"

while read -r gwas_id; do
    gwas_file="${gwas_dir}/imputed_${gwas_id}.txt.gz"
    echo "${gwas_file}" >> "${gwas_phenotype_list}"
done < "${gwas_id_list}"


completion_dir="${output_dir}/completed"
mkdir -p "${completion_dir}"

# Build aggregated gene_bed_list and covariate_list across all tissues
gene_bed_list_all="${working_dir}/all_tissues.${ld_region}.gene_bed_list.txt"
covariate_list_all="${working_dir}/all_tissues.${ld_region}.covariate_list.txt"
rm -f "${gene_bed_list_all}" "${covariate_list_all}"
touch "${gene_bed_list_all}" "${covariate_list_all}"

mkdir -p "${working_dir}/gene_beds"

while read -r tissue_id; do
    expression_path="${expression_dir}/${tissue_id}.v11.normalized_expression.bed.gz"
    covariate_path="${covariate_dir}/${tissue_id}.v11.covariates.txt"

    expression_unzipped="${working_dir}/$(basename "${expression_path}" .gz)"
    if [ ! -f "${expression_unzipped}" ]; then
        gunzip -c "${expression_path}" > "${expression_unzipped}"
    fi
    expression_header=$(head -n 1 "${expression_unzipped}")

    awk -v chr="${region_chr}" -v start="${region_start_padded}" -v end="${region_end_padded}" '($1 == chr && $2 <= end && $3 >= start){print NR}' "${expression_unzipped}" | while read -r line_num; do
        gene_line=$(sed -n "${line_num}p" "${expression_unzipped}")
        gene_id=$(echo "${gene_line}" | awk '{print $4}')
        echo "line number ${line_num} for gene ${gene_id} in region ${ld_region}"
        gene_bed_path="${working_dir}/gene_beds/${tissue_id}.${gene_id}.bed.gz"
        (echo "${expression_header}"; echo "${gene_line}") | bgzip > "${gene_bed_path}"
        tabix -s 1 -b 2 -e 3 "${gene_bed_path}"
        echo "${gene_bed_path}" >> "${gene_bed_list_all}"
        echo "${covariate_path}" >> "${covariate_list_all}"
    done
done < "${tissue_id_list}"

echo "$(wc -l < "${gene_bed_list_all}") tissue×gene phenotypes selected for analysis."

completion_file="${completion_dir}/${ld_region}.completed"
if [ ! -f "${completion_file}" ]; then
    echo "Running colocboost aggregated on LD region ${ld_region} across all tissues"

    ${code_dir}/run_colocboost_compare.sh \
        --tissue_id all_tissues \
        --ld_region ${ld_region} \
        --association_region ${association_region} \
        --gene_region ${gene_region} \
        --gene_bed_list ${gene_bed_list_all} \
        --covariate_list ${covariate_list_all} \
        --covariate_path none \
        --genotype_stem ${local_genotype_stem} \
        --gwas_id_list ${gwas_id_list} \
        --gwas_phenotype_list ${gwas_phenotype_list} \
        --gwas_meta ${gwas_meta} \
        --ld_meta ${ld_meta} \
        --gwas_column_matching ${gwas_column_matching} \
        --v39_gene_id_path ${all_v39_genes_path} \
        --output_dir ${output_dir}/all_tissues \
        --code_dir ${code_dir}

    echo $(date +"[%b %d %H:%M:%S] Done with all tissues for region ${ld_region}")
    echo "Processing completed successfully for ld region ${ld_region} across all tissues" > "${completion_file}"
    echo "Completion marker created: ${completion_file}"
else
    echo "Skipping region ${ld_region} because ${completion_file} already exists."
fi