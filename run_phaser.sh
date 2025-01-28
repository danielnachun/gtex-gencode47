
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
    sample_id
    chromosome
    haplotype_blacklist_bed
    phaser_blacklist_bed
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'bam_process' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --sample_id )
            sample_id="${2}"; shift 2 ;;
        --chromosome )
            chromosome="${2}"; shift 2 ;;
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

# echo $(date +"[%b %d %H:%M:%S] Preparing indices")
# for index in ${sep=" " bam_indices}; do
#     touch $index
# done

echo $(date +"[%b %d %H:%M:%S] Preparing bam files")
mkdir ./bam_staging
for bam_file in ${sep=" " bam_files}; do
    samtools view -h $bam_file ${chromosome} | \
    grep -v "vW:i:[2-7]" | \
    samtools view -h1 | samtools sort > ./bam_staging/$(basename $bam_file)
    samtools index -@ ${num_threads} ./bam_staging/$(basename $bam_file)
done

echo $(date +"[%b %d %H:%M:%S] Running phASER")

# this was being called with python 2.7

python2.7 phaser/wrapper.py phase ${sample_id} ./bam_staging/*.bam \
    ${genotype_vcf} ${gene_model_bed} .  \
    ${"--haplo-count-blacklist=" + haplotype_blacklist_bed} \
    ${"--blacklist=" + phaser_blacklist_bed} \
    --chr ${chromosome}