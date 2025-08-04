
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
    bam_files
    participant_id
    chromosome
    genotype_vcf
    gene_model_bed
    output_dir
    working_dir
    num_threads
    # haplotype_blacklist_bed
    # phaser_blacklist_bed
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'phaser' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --bam_files )
            bam_files="${2}"; shift 2 ;;
        --participant_id )
            participant_id="${2}"; shift 2 ;;
        --chromosome )
            chromosome="${2}"; shift 2 ;;
        --genotype_vcf )
            genotype_vcf="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --gene_model_bed )
            gene_model_bed="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --working_dir )
            working_dir="${2}"; shift 2 ;;
        --num_threads )
            num_threads="${2}"; shift 2 ;;
        # --haplotype_blacklist_bed )
        #     haplotype_blacklist_bed="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        # --phaser_blacklist_bed )
        #     phaser_blacklist_bed="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

echo $(date +"[%b %d %H:%M:%S] Subsetting chromosome from bam files")
bam_stagin_dir=${working_dir}/bam_staging_${chromosome}
# mkdir -p ${bam_stagin_dir}
# for bam_file in ${bam_files}; do
#     echo ${bam_file}
#     samtools view -h $bam_file ${chromosome} | \
#     grep -v "vW:i:[2-7]" | \
#     samtools view -h1 | samtools sort > ${bam_stagin_dir}/$(basename $bam_file)
#     samtools index -@ ${num_threads} ${bam_stagin_dir}/$(basename $bam_file)
# done

echo $(date +"[%b %d %H:%M:%S] Running phASER")
phaser_wrapper phase ${participant_id} ${bam_stagin_dir}/*.bam \
    ${genotype_vcf} ${gene_model_bed} ${output_dir}  \
    --chr ${chromosome} \
    # --haplo-count-blacklist + haplotype_blacklist_bed} \
    # --blacklist + phaser_blacklist_bed} \
