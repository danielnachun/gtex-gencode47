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
    reference_dir
    vcf_dir
    output_dir
    code_dir
    bam_list
    reference_fasta
    chr_sizes
    genes_gtf
    intervals_bed
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
            output_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --bam_list )
            bam_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --reference_fasta )
            reference_fasta="${2}"; shift 2 ;;
        --chr_sizes )
            chr_sizes="${2}"; shift 2 ;;
        --genes_gtf )
            genes_gtf="${2}"; shift 2 ;;
        --intervals_bed )
            intervals_bed="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

# activate the pixi enviroment
source <(pixi shell-hook --environment quantifybam --manifest-path ${code_dir}/pixi.toml)


# map job id to line number and then to sample id
line_number=${SLURM_ARRAY_TASK_ID}
bam_file="$(sed "${line_number}q; d" "${bam_list}")"
sample_id=$(basename $(echo ${bam_file} | sed 's/\.Aligned\.sortedByCoord\.out\.patched\.v11md\.bam//'))
participant_id=$(echo ${sample_id} | cut -d '-' -f1,2)

vcf_file=${participant_id}.snps.vcf.gz
vcf_index=${participant_id}.snps.vcf.gz.tbi

# check for vcf file and vcf file index
check_for_file "vcf_file" "${vcf_dir}/${vcf_file}"
check_for_file "vcf_index" "${vcf_dir}/${vcf_index}"

# make tmp dir
dir_prefix=${TMPDIR}/${sample_id}
vcf_dir_tmp=${dir_prefix}/vcfs
mkdir -p ${vcf_dir_tmp}
mkdir -p ${dir_prefix}/output/genome_bam
mkdir -p ${dir_prefix}/tmp

# copy references and data to temop direcotry in compute node
rsync -PrhLtv ${reference_dir} ${dir_prefix}
rsync -PrhLtv ${vcf_dir}/${vcf_file} ${vcf_dir_tmp}
rsync -PrhLtv ${vcf_dir}/${vcf_file}.tbi ${vcf_dir_tmp}
rsync -PrhLtv ${bam_file} ${dir_prefix}/output/genome_bam
rsync -PrhLtv ${bam_file}.bai ${dir_prefix}/output/genome_bam

# run rnaseq qc
bash ${code_dir}/run_rnaseq_qc.sh \
    --duplicate_marked_bam ${dir_prefix}/output/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
    --genes_gtf ${dir_prefix}/references/${genes_gtf} \
    --genome_fasta ${dir_prefix}/references/${reference_fasta} \
    --intervals_bed ${dir_prefix}/references/${intervals_bed} \
    --sample_id ${sample_id} \
    --output_dir ${dir_prefix}/output/rnaseq_qc

# run coverage
bash ${code_dir}/run_bam_to_coverage.sh \
    --duplicate_marked_bam ${dir_prefix}/output/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
    --chr_sizes ${dir_prefix}/references/${chr_sizes} \
    --sample_id ${sample_id} \
    --output_dir ${dir_prefix}/output/coverage

# run gatk
bash ${code_dir}/run_gatk.sh \
    --sample_id ${sample_id} \
    --dir_prefix ${dir_prefix} \
    --genome_fasta ${dir_prefix}/references/${reference_fasta} \
    --vcf_dir ${vcf_dir} \
    --duplicate_marked_bam ${dir_prefix}/output/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
    --output_dir ${dir_prefix}/output/gatk

# run regtools
bash ${code_dir}/run_regtools.sh \
    --sample_id ${sample_id} \
    --dir_prefix ${dir_prefix} \
    --duplicate_marked_bam ${dir_prefix}/output/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
    --output_dir ${dir_prefix}/output/leafcutter

rsync -Prhltv ${dir_prefix}/output/ ${output_dir}
