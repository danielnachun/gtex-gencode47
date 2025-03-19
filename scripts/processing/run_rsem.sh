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
    rsem_ref_dir
    transcriptome_bam
    sample_id
    output_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'fastq_to_star' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --rsem_ref_dir )
            rsem_ref_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --transcriptome_bam )
            transcriptome_bam="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --sample_id )
            sample_id="${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

mkdir -p ${output_dir}
echo $(date +"[%b %d %H:%M:%S] Running RSEM for ${sample_id}")
# RSEM transcript quantification
run_RSEM.py \
    ${rsem_ref_dir} \
    ${transcriptome_bam} \
    ${sample_id} \
    --output_dir ${output_dir} \
    --threads 4
echo $(date +"[%b %d %H:%M:%S] Done")
# # question: where is rsem ref/ when in pipeline should we make this?
# run_RSEM.py \
#         /data/rsem_reference \
#         /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/tmp/star/GTEX-1A3MV-0005-SM-7PC1O.Aligned.toTranscriptome.out.bam \
#         GTEX-1A3MV-0005-SM-7PC1O \
#         --output_dir /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/output/rsem
#         --threads 4
