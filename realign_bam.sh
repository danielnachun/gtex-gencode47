#!/usr/bin/env bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace  -o nounset -o pipefail -o errexit

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
    reference_dir
    vcf_dir
    output_dir
    code_dir
    bam_list
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
        --reference_dir )
            reference_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --vcf_dir )
            vcf_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --bam_list )
            bam_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --reference_fasta )
            reference_fasta="${2}"; shift 2 ;;
        --rsem_ref_dir )
            rsem_ref_dir="${2}"; shift 2 ;;
        --star_index )
            star_index="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

# activate the pixi enviroment
source <(pixi shell-hook --environment realignbam --manifest-path ${code_dir}/pixi.toml)


# map job id to line number and then to sample id
line_number=${SLURM_ARRAY_TASK_ID}
bam_file="$(sed "${line_number}q; d" "${bam_list}")"
sample_id=$(basename $(echo ${bam_file} | sed 's/\.Aligned\.sortedByCoord\.out\.patched\.md\.bam//'))
participant_id=$(echo ${sample_id} | cut -d '-' -f1,2)
vcf_file=${participant_id}.snps.vcf.gz
vcf_index=${participant_id}.snps.vcf.gz.tbi

# check for vcf file and vcf file index
check_for_file "vcf_file" "${vcf_dir}/${vcf_file}"
check_for_file "vcf_index" "${vcf_dir}/${vcf_index}"

# check for bam file
check_for_file "bam_file" "${bam_file}"
check_for_file "bam_file_index" "${bam_file}.bai"

# make tmp dir
dir_prefix=${TMPDIR}/${sample_id}
vcf_dir_tmp=${dir_prefix}/vcfs
mkdir -p ${vcf_dir_tmp}
mkdir -p ${dir_prefix}/raw
mkdir -p ${dir_prefix}/tmp

# copy references and data to temop direcotry in compute node
rsync -PrhLtv ${reference_dir} ${dir_prefix}
rsync -PrhLtv ${vcf_dir}/${vcf_file} ${vcf_dir_tmp}
rsync -PrhLtv ${vcf_dir}/${vcf_file}.tbi ${vcf_dir_tmp}
rsync -PrhLtv ${bam_file} ${dir_prefix}/raw
rsync -PrhLtv ${bam_file}.bai ${dir_prefix}/raw

# process bams to fastqs
bash ${code_dir}/run_bam_to_fastq.sh --bam_file ${dir_prefix}/raw/${sample_id}.Aligned.sortedByCoord.out.patched.md.bam \
    --sample_id ${sample_id} \
    --reference_fasta ${dir_prefix}/references/${reference_fasta} \
    --tmp_dir ${dir_prefix}/tmp/fastq

# align with star
bash ${code_dir}/run_fastq_to_star.sh \
    --star_index ${dir_prefix}/references/${star_index} \
    --fastq_1 ${dir_prefix}/tmp/fastq/${sample_id}_1.fastq.gz \
    --fastq_2 ${dir_prefix}/tmp/fastq/${sample_id}_2.fastq.gz \
    --sample_id ${sample_id} \
    --vcf_file ${vcf_dir_tmp}/${vcf_file} \
    --tmp_dir ${dir_prefix}/tmp/star

# sync bams
bash ${code_dir}/run_bam_sync.sh \
    --initial_bam_file ${dir_prefix}/raw/${sample_id}.Aligned.sortedByCoord.out.patched.md.bam \
    --star_aligned_bam ${dir_prefix}/tmp/star/${sample_id}.Aligned.sortedByCoord.out.bam \
    --sample_id ${sample_id} \
    --tmp_dir ${dir_prefix}/tmp/bamsync \
    --output_dir ${dir_prefix}/output/flagstat

# mark duplicates, get the genome bam that we save
bash ${code_dir}/run_mark_duplicates.sh \
    --genome_bam_file ${dir_prefix}/tmp/bamsync/${sample_id}.Aligned.sortedByCoord.out.patched.bam \
    --output_prefix ${sample_id}.Aligned.sortedByCoord.out.patched.v11md \
    --output_dir ${dir_prefix}/output/genome_bam

# run rsem, save isoform quantification
bash ./run_rsem.sh \
    --rsem_ref_dir ${dir_prefix}/references/${rsem_ref_dir} \
    --transcriptome_bam ${dir_prefix}/tmp/star/${sample_id}.Aligned.toTranscriptome.out.bam \
    --sample_id ${sample_id} \
    --output_dir ${dir_prefix}/output/rsem


rsync -Prhltv ${dir_prefix}/output/ ${output_dir}
