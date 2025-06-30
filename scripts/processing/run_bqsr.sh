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
    duplicate_marked_bam
    sample_id
    reference_fasta
    dbsnp
    indels_mills
    indels_decoy
    tmp_dir
    output_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'bsqr' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --split_bam )
            split_bam="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --sample_id )
            sample_id="${2}"; shift 2 ;;
        --reference_fasta )
            reference_fasta="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --dbsnp )
            dbsnp="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --indels_mills )
            indels_mills="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --indels_decoy )
            indels_decoy="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
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


echo $(date +"[%b %d %H:%M:%S] Running BaseRecalibrator for ${sample_id}")
gatk BaseRecalibrator \
    -R ${reference_fasta} \
    -I ${split_bam} \
    --use-original-qualities \
    -O ${tmp_dir}/${sample_id}.recal_data.csv \
    --known-sites ${dbsnp} \
    --known-sites ${indels_mills} \
    --known-sites ${indels_decoy}


echo $(date +"[%b %d %H:%M:%S] Running ApplyBQSR for ${sample_id}")
gatk ApplyBQSR \
    --add-output-sam-program-record \
    -R ${reference_fasta} \
    -I ${split_bam} \
    --use-original-qualities \
    -O ${output_dir}/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.split.recalibrated.bam \
    --bqsr-recal-file ${tmp_dir}/${sample_id}.recal_data.csv 

echo $(date +"[%b %d %H:%M:%S] Done")