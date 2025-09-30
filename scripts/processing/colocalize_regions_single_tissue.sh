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
    tissue_id
    genotype_stem
    covariate_path
    expression_path
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
        --tissue_id )
            tissue_id="${2}"; shift 2 ;;
        --genotype_stem )
            genotype_stem="${2}"; shift 2 ;;
        --covariate_path )
            covariate_path="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --expression_path )
            expression_path="${2}"; shift 2 ;;
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
rsync -PrhLtv "${covariate_path}" "${working_dir}/"
local_covariate_path="${working_dir}/$(basename "${covariate_path}")"


# activate the pixi enviroment
source <(pixi shell-hook --environment pecotmr --manifest-path ${code_dir}/pixi.toml)

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

# Unzip the expression file to a temporary location if not already uncompressed
expression_unzipped="${working_dir}/$(basename "${expression_path}" .gz)"
if [ ! -f "${expression_unzipped}" ]; then
    gunzip -c "${expression_path}" > "${expression_unzipped}"
fi
expression_header=$(head -n 1 "${expression_unzipped}")

# Initialize gene bed list
gene_bed_list="${working_dir}/${tissue_id}.${ld_region}.gene_bed_list.txt"
rm -f "${gene_bed_list}" 
touch "${gene_bed_list}" 

# Process expression file to get one phenotype file per gene
mkdir -p "${working_dir}/gene_beds"
awk -v chr="${region_chr}" -v start="${region_start_padded}" -v end="${region_end_padded}" '($1 == chr && $2 <= end && $3 >= start){print NR}' "${expression_unzipped}" | while read -r line_num; do
    # Extract the gene line
    gene_line=$(sed -n "${line_num}p" "${expression_unzipped}")
    # Get gene_id (assumed to be column 4)
    gene_id=$(echo "${gene_line}" | awk '{print $4}')
    echo "line number ${line_num} for gene ${gene_id} in region ${ld_region}"
    # Output file path
    gene_bed_path="${working_dir}/gene_beds/${tissue_id}.${gene_id}.bed.gz"
    # Write header and gene line, compress, and save
    (echo "${expression_header}"; echo "${gene_line}") | bgzip > "${gene_bed_path}"
    tabix -s 1 -b 2 -e 3 "${gene_bed_path}"
    # Add gene bed filename to all genes list
    echo "${gene_bed_path}" >> "${gene_bed_list}"
done

echo "$(wc -l < "${gene_bed_list}") genes selected for all analysis."

echo "Running colocboost on LD region ${ld_region} in tissue ${tissue_id}"

${code_dir}/run_colocboost_compare.sh \
    --tissue_id ${tissue_id} \
    --ld_region ${ld_region} \
    --association_region ${association_region} \
    --gene_region ${gene_region} \
    --gene_bed_list ${gene_bed_list} \
    --covariate_path ${local_covariate_path} \
    --genotype_stem ${local_genotype_stem} \
    --v39_gene_id_path ${all_v39_genes_path} \
    --output_dir ${output_dir} \
    --code_dir ${code_dir} 

echo $(date +"[%b %d %H:%M:%S] Done with all genes")


# Create completion marker file to indicate successful processing
completion_dir="${output_dir}/completed/single_tissue_coloc"
mkdir -p "${completion_dir}"
completion_file="${completion_dir}/${tissue_id}.${ld_region}.completed"
echo "Processing completed successfully for tissue ${tissue_id}, ld region ${ld_region}" > "${completion_file}"
echo "Completion marker created: ${completion_file}"



# # hardcoded test
# ld_region="chr1:1-1190341"
# region_padding=1000000
# working_dir="/home/klawren/oak/gtex/output/test/coloc_test"
# expression_path="/home/klawren/oak/gtex/output/test/coloc_test/Brain_Caudate_basal_ganglia.v11.normalized_expression.bed.gz"
# tissue_id="Brain_Caudate_basal_ganglia"
# code_dir="/home/klawren/oak/gtex/scripts/processing"
# genotype_stem="/home/klawren/oak/gtex/output/test/coloc_test/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.MAF01"
# covariate_path="/home/klawren/oak/gtex/output/test/coloc_test/Brain_Caudate_basal_ganglia.v11.covariates.txt"
# output_dir="/home/klawren/oak/gtex/output/test/coloc_test/results"
# all_v39_genes_path="/home/klawren/oak/gtex/data/other_references/gencode/v39_genes.txt"