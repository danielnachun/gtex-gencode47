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
    fastq_1
    fastq_2
    sample_id
    tmp_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'ipa_finder' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --fastq_1 )
            fastq_1="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --fastq_2 )
            fastq_2="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --sample_id )
            sample_id="${2}"; shift 2 ;;
        --tmp_dir )
            tmp_dir="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

mkdir -p ${tmp_dir}

clean_1=${tmp_dir}/${sample_id}.clean_1.fastq.gz
clean_2=${tmp_dir}/${sample_id}.clean_2.fastq.gz

echo $(date +"[%b %d %H:%M:%S] Running fastp on fastqs")
# qc on fastqs with fastp
fastp -i ${fastq_1} \
    -I ${fastq_2} \
    -o ${clean_1} \
    -O ${clean_2} \
    -t 1 \
    -c -p -5 -3 \
    -j ${tmp_dir}/${sample_id}_fasstp_qc.json \
    -h ${tmp_dir}/${sample_id}_fastp_qc.html

# unzip fastqs for sprint
gunzip ${clean_1}
gunzip ${clean_2}

echo $(date +"[%b %d %H:%M:%S] Done")


