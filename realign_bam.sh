#!/usr/bin/env bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o pipefail -o errexit

# example args
# dir_prefix = /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow
# sample_id = GTEX-1A3MV-0005-SM-7PC1O
# reference_fasta = /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/references/reference.fasta
# rsem_ref_dir = ??
# star_index = ??


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
    dir_prefix
    sample_id
    reference_fasta
    rsem_ref_dir
    star_index
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'bam_process' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --dir_prefix )
            dir_prefix="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --sample_id )
            sample_id="${2}"; shift 2 ;;
        --reference_fasta )
            reference_fasta="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --rsem_ref_dir )
            rsem_ref_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --star_index )
            star_index="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

# process bams to fastqs
run_bam_to_fastq.sh --bam_file ${dir_prefix}/data/old_bams/${sample_id}.Aligned.sortedByCoord.out.patched.md.bam \
    --sample_id ${sample_id} \
    --reference_fasta ${reference_fasta} \
    --output_dir ${dir_prefix}/tmp/fastq

# align with star
run_fastq_to_star.sh --star_index ${star_index} \
    --fastq_1 ${dir_prefix}/tmp/fastq/${sample_id}.fastq.gz \
    --fastq_2 ${dir_prefix}/tmp/fastq/${sample_id}.fastq.gz \
    --sample_id ${sample_id} \
    --output_dir ${dir_prefix}/tmp/star

# mark duplicates, get the genome bam that we save
run_mark_duplicates.sh --genome_bam_file ${dir_prefix}/tmp/star/${sample_id}.Aligned.sortedByCoord.out.patched.bam \
    --output_prefix ${sample_id}.Aligned.sortedByCoord.out.patched.v11md \
    --output_dir ${dir_prefix}/output/genome_bam

# run rsem, save isoform quantification
run_rsem.sh --rsem_ref_dir ${rsem_ref_dir} \
    --transcriptome_bam ${dir_prefix}/tmp/star/${sample_id}.Aligned.toTranscriptome.out.bam\
    --sample_id ${sample_id} \
    --output_dir ${dir_prefix}/output/rsem

