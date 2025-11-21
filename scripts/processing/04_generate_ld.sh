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
    genotype_prefix
    sample_ids
    output_dir
    code_dir
    regenerate
    region
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'fastq_to_star' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --genotype_prefix )
            genotype_prefix="${2}"; shift 2 ;;
        --sample_ids )
            sample_ids="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --regenerate )
            regenerate="${2}"; shift 2 ;;
        --region )
            region="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

# activate the pixi enviroment
source <(pixi shell-hook --environment plink --manifest-path ${code_dir}/pixi.toml)


echo "Processing region: ${region}"
# Parse region string format: chr:start-end (1-based inclusive)
# Example: chr1:1000-2000 or 1:1000-2000
chr_id=$(echo "${region}" | cut -d':' -f 1)
region_range=$(echo "${region}" | cut -d':' -f 2)
from_bp=$(echo "${region_range}" | cut -d'-' -f 1)
to_bp=$(echo "${region_range}" | cut -d'-' -f 2)

# Remove "chr" prefix if present for PLINK
chr_id=$(echo "$chr_id" | sed 's/^chr//')

out_file="${output_dir}/LD_chr${chr_id}_${from_bp}_${to_bp}.ld.gz"
echo "Checking for ${out_file}"

if [[ ${regenerate} = false ]] && [[ -f ${out_file} ]]; then
    echo $(date +"[%b %d %H:%M:%S] LD file ${out_file} already exists, skipping")
else
    echo $(date +"[%b %d %H:%M:%S] Running PLINK LD for chr${chr_id}:${from_bp}-${to_bp}")
    ${code_dir}/run_04_plink_ld.sh \
        --genotype_prefix "$genotype_prefix" \
        --sample_ids "$sample_ids" \
        --chr_id "$chr_id" \
        --from_bp "$from_bp" \
        --to_bp "$to_bp" \
        --output_dir "$output_dir"
fi

# Create completion marker using normalized region string for filesystem safety
# Normalize the region string to handle colons and other special characters
completion_dir="${output_dir}/completed"
mkdir -p "${completion_dir}"
normalized_region=$(echo "$region" | tr '\t' '_' | tr '\n' ' ' | sed 's/[:/]/_/g' | tr -d '\r' | sed 's/  */_/g' | sed 's/^_\|_$//g')
completion_file="${completion_dir}/${normalized_region}.completed"
echo "Processing completed successfully for region: ${region}" > "${completion_file}"
echo "Completion marker created: ${completion_file}"

echo $(date +"[%b %d %H:%M:%S] Done")
