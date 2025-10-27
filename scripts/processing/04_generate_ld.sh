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
    regions_bed
    regenerate
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
        --regions_bed )
            regions_bed="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --regenerate )
            regenerate="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

# activate the pixi enviroment
source <(pixi shell-hook --environment plink --manifest-path ${code_dir}/pixi.toml)

# map job id to line number and then to region
line_number=${SLURM_ARRAY_TASK_ID}

# pull chr, from_bp, to_bp from the line of the regions_bed
region="$(sed "${line_number}q; d" "${regions_bed}")" 
echo "Processing region: ${region}"
chr_id=$(echo "${region}" | cut -f 1)
from_bp=$(echo "${region}" | cut -f 2)
to_bp=$(echo "${region}" | cut -f 3)

out_file="${output_dir}/LD_chr${chr_id}_${from_bp}_${to_bp}.ld.gz"
echo "Checking for ${out_file}"

if [[ ${regenerate} = false ]] && [[ -f ${out_file} ]]; then
    echo $(date +"[%b %d %H:%M:%S] LD file ${out_file} already exists, skipping")
else
    echo $(date +"[%b %d %H:%M:%S] Running PLINK LD for chr${chr_id}:${from_bp}-${to_bp}")
    ${code_dir}/run_plink_ld.sh \
        --genotype_prefix "$genotype_prefix" \
        --sample_ids "$sample_ids" \
        --chr_id "$chr_id" \
        --from_bp "$from_bp" \
        --to_bp "$to_bp" \
        --output_dir "$output_dir"
fi

echo $(date +"[%b %d %H:%M:%S] Done")
