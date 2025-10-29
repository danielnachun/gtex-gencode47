#!/usr/bin/env bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o pipefail -o errexit

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
    initial_bam_file
    star_aligned_bam
    sample_id
    tmp_dir
    output_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'fastq_to_star' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --initial_bam_file )
            initial_bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --star_aligned_bam )
            star_aligned_bam="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --sample_id )
            sample_id="${2}"; shift 2 ;;
        --tmp_dir )
            tmp_dir="${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

mkdir -p ${tmp_dir}
mkdir -p ${output_dir}

echo $(date +"[%b %d %H:%M:%S] Running bamsync")
run_bamsync.sh ${initial_bam_file} ${star_aligned_bam} ${tmp_dir}/${sample_id}
echo $(date +"[%b %d %H:%M:%S] Running samtools flagstat")
samtools flagstat ${tmp_dir}/${sample_id}.Aligned.sortedByCoord.out.patched.bam > ${output_dir}/${sample_id}.Aligned.sortedByCoord.out.patched.bam.flagstat.txt
echo $(date +"[%b %d %H:%M:%S] Done")
