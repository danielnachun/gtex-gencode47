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
    bam_file
    sample_id
    reference_fasta
    output_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'bam_process' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --bam_file )
            bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --sample_id )
            sample_id="${2}"; shift 2 ;;
        --reference_fasta )
            reference_fasta="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
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

# question: do we have to make the output dir first? 
# or is that handled by the script which will call this script?

mkdir -p ${output_dir}
echo $(date +"[%b %d %H:%M:%S] getting fastqs for ${sample_id}")
run_SamToFastq.py ${bam_file} \
    --jar $(pwd)/.pixi/envs/default/share/picard-2.27.1-0/picard.jar \
    -p ${sample_id} \
    --reference_fasta ${reference_fasta} \
    --output_dir ${output_dir} \
    --memory 64
echo $(date +"[%b %d %H:%M:%S] Done")

# run_SamToFastq.py /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/data/test_bams/GTEX-1A3MV-0005-SM-7PC1O.Aligned.sortedByCoord.out.patched.md.bam \
#     --jar $(pwd)/.pixi/envs/default/share/picard-2.27.1-0/picard.jar \
#     -p GTEX-1A3MV-0005-SM-7PC1O \
#     --reference_fasta /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
#     --output_dir /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/tmp/fastq \
#     --memory 64
