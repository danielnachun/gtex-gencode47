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
    dir_prefix
    sample_id
    reference_fasta
    chr_sizes
    genes_gtf
    intervals_bed
    vcf_dir
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
        --chr_sizes )
            chr_sizes="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --genes_gtf )
            genes_gtf="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --intervals_bed )
            intervals_bed="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --vcf_dir )
            vcf_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done



# run rnaseq qc
run_rnaseq_qc.sh \
    --duplicate_marked_bam ${dir_prefix}/output/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
     --genes_gtf ${genes_gtf} \
     --genome_fasta ${reference_fasta} \
     --sample_id ${sample_id} \
     --output_dir ${dir_prefix}/output/rnaseq_qc

# run coverage
run_bam_to_coverage.sh \
    --duplicate_marked_bam ${dir_prefix}/output/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
    --chr_sizes ${chr_sizes} \
    --sample_id ${sample_id} \
    --intervals_bed ${intervals_bed} \
    --output_dir ${dir_prefix}/output/coverage

# run gatk
run_gatk.sh \
    --sample_id ${sample_id} \
    --dir_prefix ${dir_prefix} \
    --genome_fasta ${reference_fasta} \
    --vcf_dir ${vcf_dir} \
    --duplicate_marked_bam ${dir_prefix}/output/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
    --output_dir ${dir_prefix}/output/gatk

# run regtools
run_regtools.sh \
    --sample_id ${sample_id} \
    --dir_prefix ${dir_prefix} \
    --duplicate_marked_bam ${dir_prefix}/output/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
    --output_dir ${dir_prefix}/output/leafcutter
