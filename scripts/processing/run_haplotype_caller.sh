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
    bqsr_bam
    sample_id
    reference_fasta
    gene_intervals_bed
    dbsnp
    output_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'haplotype_caller' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --bqsr_bam )
            bqsr_bam="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --sample_id )
            sample_id="${2}"; shift 2 ;;
        --reference_fasta )
            reference_fasta="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --gene_intervals_bed )
            gene_intervals_bed="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --dbsnp )
            dbsnp="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
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

echo $(date +"[%b %d %H:%M:%S] Running HaplotypeCaller for ${sample_id}")

gatk HaplotypeCaller \
    -R ${reference_fasta} \
    -I ${bqsr_bam} \
    -L ${gene_intervals_bed} \
    -O ${output_dir}/${sample_id}.hc.vcf.gz \
    -A ReadPosition \
    -A FragmentLength \
    -dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 0 \
    --min-base-quality-score 20 \
    --dbsnp ${dbsnp}


echo $(date +"[%b %d %H:%M:%S] Filtering HaplotypeCaller")
bcftools view ${output_dir}/${sample_id}.hc.vcf.gz \
    --include "INFO/DP>=10 & MQ>=40 & MQRankSum>=-12.5 & QD>=2 & ReadPosRankSum>=-8 & ReadPosRankSum<8 & AD > 3 & MPOS > 6" \
    --min-alleles 2 \
    --max-alleles 2 \
    --output ${output_dir}/${sample_id}.hc_filtered.vcf.gz

echo $(date +"[%b %d %H:%M:%S] Done")


