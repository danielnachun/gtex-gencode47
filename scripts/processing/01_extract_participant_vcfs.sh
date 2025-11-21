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
    full_vcf
    output_dir
    code_dir
    participant_id
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'fastq_to_star' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --full_vcf )
            full_vcf="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --participant_id )
            participant_id="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

# activate the pixi enviroment
source <(pixi shell-hook --environment extractvcf --manifest-path ${code_dir}/pixi.toml)

echo $(date +"[%b %d %H:%M:%S] Generating participant VCF (SNPs only)")
# select SNPs, filter out missing sites
bcftools view --no-update -s ${participant_id} -v snps ${full_vcf} -Ou | bcftools view --no-update -e 'GT=".|."' -Oz -o ${output_dir}/${participant_id}.snps.vcf.gz
tabix ${output_dir}/${participant_id}.snps.vcf.gz

echo $(date +"[%b %d %H:%M:%S] Subsetting biallelic het sites for ASE")
bcftools view --no-update -i 'GT="het"' ${output_dir}/${participant_id}.snps.vcf.gz -Ou | bcftools norm -m+ -Ou | bcftools view -m2 -M2 -Oz -o ${output_dir}/${participant_id}.snps.het.vcf.gz
tabix ${output_dir}/${participant_id}.snps.het.vcf.gz

# Create completion marker
completion_dir="${output_dir}/completed/extract_participant_vcfs"
mkdir -p "${completion_dir}"
completion_file="${completion_dir}/${participant_id}.completed"
echo "Processing completed successfully for participant: ${participant_id}" > "${completion_file}"
echo "Completion marker created: ${completion_file}"

echo $(date +"[%b %d %H:%M:%S] Done")
