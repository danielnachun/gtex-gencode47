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
    vcf_file
    exac_reference
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
        --vcf_file )
            vcf_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --exac_reference )
            exac_reference="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
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

echo $(date +"[%b %d %H:%M:%S] Running Mutect2 for ${sample_id}")

gatk Mutect2 \
    --input "${bqsr_bam}" \
    --output "${output_dir}/${sample_id}.mutect2.vcf" \
    --reference "${reference_fasta}" \
    --intervals "${gene_intervals_bed}" \
    --dont-use-soft-clipped-bases \
    --panel-of-normals "${vcf_file}" \
    --f1r2-tar-gz "${output_dir}/${sample_id}.f1r2.tar.gz" \
    -A OrientationBiasReadCounts \
    -A ReadPosRankSumTest \
    -A MappingQualityRankSumTest \
    -A RMSMappingQuality \
    -A QualByDepth \
    -A FragmentLength


echo $(date +"[%b %d %H:%M:%S] Getting pileup summaries.")
gatk GetPileupSummaries \
    --input "${bqsr_bam}" \
    --variant "${exac_reference}" \
    --intervals "${exac_reference}" \
    --reference "${reference_fasta}" \
    --output ${output_dir}/${sample_id}.pileups.table

echo $(date +"[%b %d %H:%M:%S] Calculating cross sample contamination.")
gatk CalculateContamination \
        --input "${output_dir}/${sample_id}.pileups.table" \
        --output "${output_dir}/${sample_id}.contamination.table"


echo $(date +"[%b %d %H:%M:%S] Learning read orinetation bias model.")
gatk LearnReadOrientationModel \
        --input ${output_dir}/${sample_id}.f1r2.tar.gz \
        --output ${output_dir}/${sample_id}.artifact_prior.tar.gz

    
echo $(date +"[%b %d %H:%M:%S] Filtering somatic variants with FilterMutectCalls.")
gatk FilterMutectCalls \
        --variant "${output_dir}/${sample_id}.mutect2.vcf" \
        --output "${output_dir}/${sample_id}.mutect2.filtered.vcf" \
        --contamination-table "${output_dir}/${sample_id}.contamination.table" \
        --ob-priors  "${output_dir}/${sample_id}.artifact_prior.tar.gz" \
        --reference "${reference_fasta}" \
        --min-median-base-quality 20 \
        --max-alt-allele-count 1 \
        --min-median-read-position 6 \
        --unique-alt-read-count 3 \
        --read-filter NotSupplementaryAlignmentReadFilter 


echo $(date +"[%b %d %H:%M:%S] Done")

