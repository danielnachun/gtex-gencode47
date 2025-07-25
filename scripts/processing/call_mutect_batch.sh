#!/usr/bin/env bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace  -o nounset -o pipefail -o errexit

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
    reference_dir
    output_dir
    code_dir
    vcf_dir
    bam_list_paths
    reference_fasta
    dbsnp
    indels_mills
    indels_decoy
    gene_intervals_bed
    full_vcf_file
    exac_reference
    step_size
)


longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'call_edsites_batch' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --reference_dir )
            reference_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --vcf_dir )
            vcf_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --bam_list_paths )
            bam_list_paths="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --reference_fasta )
            reference_fasta="${2}"; shift 2 ;;
        --dbsnp )
            dbsnp="${2}"; shift 2 ;;
        --indels_mills )
            indels_mills="${2}"; shift 2 ;;
        --indels_decoy )
            indels_decoy="${2}"; shift 2 ;;
        --gene_intervals_bed )
            gene_intervals_bed="${2}"; shift 2 ;;
        --full_vcf_file )
            full_vcf_file="${2}"; shift 2 ;;   
        --exac_reference )
            exac_reference="${2}"; shift 2 ;;            
        --step_size )
            step_size="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

# get the paths to the bam files for this batch
line_number="${SLURM_ARRAY_TASK_ID}"
bam_list="$(sed "${line_number}q; d" "${bam_list_paths}")"

# copy over the references to this node
# folder specific to this batch
reference_dir_prefix="${TMPDIR}/references_${line_number}"
mkdir -p "${reference_dir_prefix}"
rsync -PrhLtv "${reference_dir}"/* "${reference_dir_prefix}/"

echo "Step size: ${step_size}"
echo "Requested CPUs: ${SLURM_CPUS_PER_TASK}"

# run the batch
cat "${bam_list}" | parallel -j"${step_size}" --ungroup --verbose \
        "${code_dir}/call_mutect.sh" \
        --local_reference_dir "${reference_dir_prefix}/" \
        --output_dir "${output_dir}" \
        --code_dir "${code_dir}" \
        --vcf_dir "${vcf_dir}" \
        --bam_file {} \
        --reference_fasta "${reference_fasta}" \
        --dbsnp "${dbsnp}" \
        --indels_mills "${indels_mills}" \
        --indels_decoy "${indels_decoy} " \
        --gene_intervals_bed "${gene_intervals_bed}" \
        --full_vcf_file "${full_vcf_file}" \
        --exac_reference "${exac_reference}"

echo "Batch finished"
