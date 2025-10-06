#!/usr/bin/env bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o nounset -o pipefail -o errexit

# Helper function to check for file existence
check_for_file() {
    argument_name="${1}"
    file_path="${2}"
    if [[ ${file_path} != "none" ]] && [[ ! -f ${file_path} ]]; then
        echo "Error: file ${file_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

# Helper function to check for directory existence
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
    genotype_stem
    covariate_dir
    expression_dir
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
        --genotype_stem )
            genotype_stem="${2}"; shift 2 ;;
        --covariate_dir )
            covariate_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --expression_dir )
            expression_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
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

mkdir -p "${output_dir}"

# Ensure TMPDIR has a sane default to avoid unbound variable exit
TMPDIR="${TMPDIR:-/tmp}"

line_number="${SLURM_ARRAY_TASK_ID}"
working_dir="${TMPDIR}/coloc_${line_number}"
mkdir -p "${working_dir}"
rsync -PrhLtv "${genotype_stem}"* "${working_dir}/"
local_genotype_stem="${working_dir}/$(basename "${genotype_stem}")"

# Activate the pixi environment
source <(pixi shell-hook --environment pecotmr-dev --manifest-path ${code_dir}/pixi.toml)

# Map job id to line number and then to participant id
line_number=${SLURM_ARRAY_TASK_ID}
ld_region_bed_line=$(sed "${line_number}q; d" "${ld_region_list}")

# Parse the region from colon-separated format (e.g., "chr1:16104-1170341")
region_chr=$(echo "${ld_region_bed_line}" | cut -d':' -f1)
region_coords=$(echo "${ld_region_bed_line}" | cut -d':' -f2)
region_start=$(echo "${region_coords}" | cut -d'-' -f1)
region_end=$(echo "${region_coords}" | cut -d'-' -f2)
ld_region="${region_chr}:${region_start}-${region_end}"
echo $(date +"[%b %d %H:%M:%S] Running colocboost with and without v39 genes on region ${ld_region}")

# Get region for variants
association_region_start=$((region_start - association_padding))
association_region_end=$((region_end + association_padding))
if [ "${association_region_start}" -lt 0 ]; then
    association_region_start=0
fi
association_region="${region_chr}:${association_region_start}-${association_region_end}"

# Get regions for genes
region_start_padded=$((region_start - region_padding))
region_end_padded=$((region_end + region_padding))
if [ "${region_start_padded}" -lt 0 ]; then
    region_start_padded=0
fi
gene_region="${region_chr}:${region_start_padded}-${region_end_padded}"

completion_dir="${output_dir}/completed"
mkdir -p "${completion_dir}"

# Function to process a single tissue_id
process_tissue() {
    tissue_id="$1"
    completion_file="${completion_dir}/${tissue_id}_${ld_region}.completed"
    all_genes_file="${output_dir}/${tissue_id}/${tissue_id}.${ld_region}.all_genes.colocboost.rds"
    if [ -f "${completion_file}" ] || [ -f "${all_genes_file}" ]; then
        if [ -f "${completion_file}" ]; then
            echo "Skipping tissue ${tissue_id} for region ${ld_region} because ${completion_file} already exists."
        fi
        if [ -f "${all_genes_file}" ]; then
            echo "Skipping tissue ${tissue_id} for region ${ld_region} because ${all_genes_file} already exists."
        fi
        return
    fi

    expression_path="${expression_dir}/${tissue_id}.v11.normalized_expression.bed.gz"
    covariate_path="${covariate_dir}/${tissue_id}.v11.covariates.txt"

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
        gene_line=$(sed -n "${line_num}p" "${expression_unzipped}")
        gene_id=$(echo "${gene_line}" | awk '{print $4}')
        echo "line number ${line_num} for gene ${gene_id} in region ${ld_region}"
        gene_bed_path="${working_dir}/gene_beds/${tissue_id}.${gene_id}.bed.gz"
        (echo "${expression_header}"; echo "${gene_line}") | bgzip > "${gene_bed_path}"
        tabix -s 1 -b 2 -e 3 "${gene_bed_path}"
        echo "${gene_bed_path}" >> "${gene_bed_list}"
    done

    echo "$(wc -l < "${gene_bed_list}") genes selected for all analysis."
    echo "Running colocboost on LD region ${ld_region} in tissue ${tissue_id}"

    ${code_dir}/run_colocboost_compare_qtl_only.sh \
        --tissue_id ${tissue_id} \
        --ld_region ${ld_region} \
        --association_region ${association_region} \
        --gene_region ${gene_region} \
        --gene_bed_list ${gene_bed_list} \
        --covariate_path ${covariate_path} \
        --genotype_stem ${local_genotype_stem} \
        --v39_gene_id_path ${all_v39_genes_path} \
        --output_dir ${output_dir}/${tissue_id} \
        --code_dir ${code_dir}

    echo $(date +"[%b %d %H:%M:%S] Done with all genes for tissue ${tissue_id}")
    echo "Processing completed successfully for ld region ${ld_region} in tissue ${tissue_id}" > "${completion_file}"
    echo "Completion marker created: ${completion_file}"
}

export -f process_tissue
export working_dir region_chr region_start_padded region_end_padded ld_region \
    association_region gene_region completion_dir output_dir code_dir \
    all_v39_genes_path local_genotype_stem expression_dir covariate_dir \
    association_padding region_padding

# Read tissue IDs into an array
mapfile -t tissue_ids < "${tissue_id_list}"

# Use number of threads (if set) for parallel jobs, otherwise default to 1
if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
    parallel_jobs=$((SLURM_CPUS_PER_TASK - 1))
elif [[ -n "${OMP_NUM_THREADS:-}" ]]; then
    parallel_jobs=$((OMP_NUM_THREADS - 1))
else
    parallel_jobs=1
fi

echo "Running process_tissue in parallel for ${parallel_jobs} jobs, across ${#tissue_ids[@]} tissues"

# Set R to use only 1 thread per job to avoid exceeding CPU limits
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Run process_tissue in parallel for each tissue_id
printf "%s\n" "${tissue_ids[@]}" | parallel --tag --line-buffer -j "${parallel_jobs}" process_tissue {} 2>&1

# Create completion marker file to indicate successful processing
overall_completion_file="${completion_dir}/${ld_region}.completed"
echo "Processing completed successfully for ld region ${ld_region}" > "${overall_completion_file}"
echo "Completion marker created: ${overall_completion_file}"