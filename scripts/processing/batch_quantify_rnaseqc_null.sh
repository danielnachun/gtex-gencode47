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
    bam_list_paths
    reference_fasta
    chr_sizes
    genes_gtf
    intervals_bed
    step_size
    bam_file_end
)


longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'bam_process' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --reference_dir )
            reference_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --bam_list_paths )
            bam_list_paths="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --reference_fasta )
            reference_fasta="${2}"; shift 2 ;;
        --chr_sizes )
            chr_sizes="${2}"; shift 2 ;;
        --genes_gtf )
            genes_gtf="${2}"; shift 2 ;;
        --intervals_bed )
            intervals_bed="${2}"; shift 2 ;;
        --step_size )
            step_size="${2}"; shift 2 ;;
        --bam_file_end )
            bam_file_end="${2}"; shift 2 ;;
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

# I'm getting errors requesting too many jobs?
echo "Step size: ${step_size}"
echo "Requested CPUs: ${SLURM_CPUS_PER_TASK}"

# run the batch
cat "${bam_list}" | parallel -j"${step_size}" --eta --ungroup \
        "${code_dir}/quantify_rnaseqc_null.sh" \
        --local_reference_dir "${reference_dir_prefix}/" \
        --output_dir "${output_dir}" \
        --code_dir "${code_dir}" \
        --bam_file {} \
        --reference_fasta "${reference_fasta}" \
        --chr_sizes "${chr_sizes}" \
        --genes_gtf "${genes_gtf}" \
        --intervals_bed "${intervals_bed}" \
        --bam_file_end "${bam_file_end}"

echo "Batch finished"
