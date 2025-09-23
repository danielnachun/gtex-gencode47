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
    clean_fastq_1
    clean_fastq_2
    sample_id
    repeat_bed
    reference_fasta
    output_dir
    code_dir
    ss
    get_hyperedited_reads
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'ipa_finder' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --clean_fastq_1 )
            clean_fastq_1="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --clean_fastq_2 )
            clean_fastq_2="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --sample_id )
            sample_id="${2}"; shift 2 ;;
        --repeat_bed )
            repeat_bed="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --reference_fasta )
            reference_fasta="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; shift 2 ;;
        --ss )
            ss="${2}"; shift 2 ;;
        --get_hyperedited_reads )
            get_hyperedited_reads="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

mkdir -p ${output_dir}

bwa_path=`which bwa`
samtools_path=`which samtools`

echo $(date +"[%b %d %H:%M:%S] Running SPRINT from fastqs")

# run sprint to indentify RES (RNA editing sites)
python ${code_dir}/sprint_main.py \
    -1 ${clean_fastq_1} \
    -2 ${clean_fastq_2} \
    -ss ${ss} \
    -rp ${repeat_bed} \
    -c 6 \
    -p 6 \
    ${reference_fasta} \
    ${output_dir} \
    ${bwa_path} \
    ${samtools_path}
echo $(date +"[%b %d %H:%M:%S] Done")

if [ "${get_hyperedited_reads}" == "true" ]; then
    echo $(date +"[%b %d %H:%M:%S] Getting the hypereditied reads into BAM format")
    # make a bam file of the hypereditied reads
    python ${code_dir}/zz2sam.py ${output_dir}/tmp/all_combined.zz
    samtools view -H ${output_dir}/tmp/genome/all.bam > ${output_dir}/tmp/SAMheader.txt
    cat ${output_dir}/tmp/SAMheader.txt ${output_dir}/tmp/all_combined.zz.sam > ${output_dir}/tmp/all_combined.zz.sam.header
    samtools view -bS ${output_dir}/tmp/all_combined.zz.sam.header > ${output_dir}/tmp/all_combined.zz.bam
    samtools sort ${output_dir}/tmp/all_combined.zz.bam -f ${output_dir}/${sample_id}.hyperedited_reads.sorted.bam
    samtools index ${output_dir}/${sample_id}.hyperedited_reads.sorted.bam
    echo $(date +"[%b %d %H:%M:%S] Done")
fi


# sample_id='GTEX-1C4CL-2126-SM-7IGQC'
# code_dir='/home/klawren/oak/gtex/scripts/processing'
# output_dir='output/test_bams/output_kate/sprint/GTEX-1C4CL-2126-SM-7IGQC/unaligned'
# fastq_1="output/test_bams/output_kate/sprint/GTEX-1C4CL-2126-SM-7IGQC/GTEX-1C4CL-2126-SM-7IGQC.unaligned_1.fastq"
# fastq_2="output/test_bams/output_kate/sprint/GTEX-1C4CL-2126-SM-7IGQC/GTEX-1C4CL-2126-SM-7IGQC.unaligned_2.fastq"
# repeat_bed="data/sprint_references/hg38_repeat.bed"
# reference_fasta="data/sprint_references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta"

