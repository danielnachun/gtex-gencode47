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
    repeat_bed
    reference_fasta
    output_dir
    code_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'ipa_finder' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --duplicate_marked_bam )
            duplicate_marked_bam="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --repeat_bed )
            repeat_bed="${2}"; shift 2 ;;
        --reference_fasta )
            reference_fasta="${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

mkdir -p ${output_dir}
samtools_path=`which samtools`

echo $(date +"[%b %d %H:%M:%S] Running SPRINT from aligned BAM")

# run sprint on aligned bam to get regular editing
# the expected argument order does not match the manual
# python /home/klawren/oak/gtex/scripts/processing/.pixi/envs/calledsites/lib/python2.7/site-packages/sprint/sprint_from_bam.py \

python ${code_dir}/sprint_from_bam.py \
    -rp ${repeat_bed} \
    ${duplicate_marked_bam} \
    ${reference_fasta} \
    ${output_dir} \
    ${samtools_path}

echo $(date +"[%b %d %H:%M:%S] Done")

# code_dir='/home/klawren/oak/gtex/scripts/processing'
# output_dir='output/test_bams/output_kate/sprint/GTEX-1C4CL-2126-SM-7IGQC/aligned'
# duplicate_marked_bam='output/test_bams/output_kate/genome_bam/GTEX-1C4CL-2126-SM-7IGQC.Aligned.sortedByCoord.out.patched.v11md.bam'
# repeat_bed="data/edsite_references/hg38_repeat.bed"
# reference_fasta="data/edsite_references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta"
