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
    genome_bam_file
    output_prefix
    output_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'fastq_to_star' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --genome_bam_file )
            genome_bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_prefix )
            output_prefix="${2}"; shift 2 ;;
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
echo $(date +"[%b %d %H:%M:%S] Marking duplicated for sample ${sample_id}")
# mark duplicates with Picard
run_MarkDuplicates.py ${genome_bam_file} \
    --jar $(pwd)/.pixi/envs/default/share/picard-2.27.1-0/picard.jar \
    ${output_prefix} \
    --output_dir ${output_dir}

echo $(date +"[%b %d %H:%M:%S] Indexing for ${sample_id}")
samtools index ${output_dir}/${output_prefix}.bam
echo $(date +"[%b %d %H:%M:%S] Done")

# run_MarkDuplicates.py \
#        --jar $(pwd)/.pixi/envs/default/share/picard-2.27.1-0/picard.jar \
#         /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/tmp/star/GTEX-1A3MV-0005-SM-7PC1O.Aligned.sortedByCoord.out.patched.bam \
#         GTEX-1A3MV-0005-SM-7PC1O.Aligned.sortedByCoord.out.patched.v11md \
#         --output_dir oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/output/genomebam

