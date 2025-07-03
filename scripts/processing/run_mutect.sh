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
    --panel-of-normals "${vcf_file}" 

echo $(date +"[%b %d %H:%M:%S] Done")

# bash scripts/processing/run_mutect.sh \
#     --bqsr_bam /home/klawren/oak/gtex/output/test_bams/output_kate/bqsr/GTEX-1C4CL-2126-SM-7IGQC.Aligned.sortedByCoord.out.patched.v11md.recalibrated.bam \
#     --sample_id GTEX-1C4CL-2126-SM-7IGQC \
#     --reference_fasta /home/klawren/oak/gtex/data/edsite_references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
#     --gene_intervals_bed /home/klawren/oak/gtex/data/edsite_references/gencode.v47.merged.bed \
#     --vcf_file /home/klawren/oak/gtex/data/processed/vcfs/GTEX-1C4CL.snps.het.vcf.gz \
#     --output_dir /home/klawren/oak/gtex/output/test_bams/output_kate/mutect

