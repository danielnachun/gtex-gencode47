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
    dir_prefix
    genome_fasta
    het_vcf
    duplicate_marked_bam
    output_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'fastq_to_star' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --sample_id )
            sample_id="${2}"; shift 2 ;;
        --dir_prefix )
            dir_prefix="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --genome_fasta )
            genome_fasta="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --het_vcf )
            het_vcf="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --duplicate_marked_bam )
            duplicate_marked_bam="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

mkdir -p ${dir_prefix}/tmp/filtered_gatk_bam/
mkdir -p ${output_dir}

# filter_wasp is true by default
echo $(date +"[%b %d %H:%M:%S] Filtering out reads with allelic mapping bias")
samtools view -h ${duplicate_marked_bam} | grep -v "vW:i:[2-7]" | samtools view -1 > ${dir_prefix}/tmp/filtered_gatk_bam/${sample_id}_filtered.bam
samtools index ${dir_prefix}/tmp/filtered_gatk_bam/${sample_id}_filtered.bam
echo $(date +"[%b %d %H:%M:%S] Running GATK on ${sample_id}")
# dan will patch conda so gatk_jar is not needed
run_GATK_ASEReadCounter.py ${genome_fasta} \
    ${het_vcf} \
    ${dir_prefix}/tmp/filtered_gatk_bam/${sample_id}_filtered.bam \
    ${output_dir}/${sample_id}

# filter out chrX
mv ${output_dir}/${sample_id}.readcounts.txt.gz ${output_dir}/${sample_id}.readcounts.all.txt.gz
zcat ${output_dir}/${sample_id}.readcounts.all.txt.gz | awk '$1!="chrX" && $1!="X" {print $0}' | gzip -c > ${output_dir}/${sample_id}.readcounts.txt.gz
zcat ${output_dir}/${sample_id}.readcounts.all.txt.gz | awk '$1=="contig" || $1=="chrX" || $1=="X" {print $0}' | gzip -c > ${output_dir}/${sample_id}.readcounts.chrX.txt.gz
echo $(date +"[%b %d %H:%M:%S] Done")
