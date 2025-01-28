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
    sample_id
    dir_prefix
    duplicate_marked_bam
    output_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'fastq_to_star' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --sample_id )
            sample_id="${2}"; shift 2 ;;
        --dir_prefix )
            dir_prefix="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --duplicate_marked_bam )
            duplicate_marked_bam="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

pixi shell

mkdir -p ${output_dir}
mkdir -p ${dir_prefix}/tmp/filtered_leafcutter_bam/

echo $(date +"[%b %d %H:%M:%S] Extracting junctions for sample ${sample_id}")
# select uniquely mapped reads that pass WASP filters
samtools view -h -q 255 ${duplicate_marked_bam} | grep -v "vW:i:[2-7]" | samtools view -b > ${dir_prefix}/tmp/filtered_leafcutter_bam/${sample_id}_filtered.bam
samtools index ${dir_prefix}/tmp/filtered_leafcutter_bam/${sample_id}_filtered.bam
# run leafcutter
regtools junctions extract -a 8 -m 50 -M 500000 -s 0 ${dir_prefix}/tmp/filtered_leafcutter_bam/${sample_id}_filtered.bam | gzip -c > ${output_dir}/${sample_id}.regtools_junc.txt.gz
echo $(date +"[%b %d %H:%M:%S] Done")