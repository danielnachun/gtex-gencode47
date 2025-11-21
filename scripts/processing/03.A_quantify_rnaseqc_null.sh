#!/usr/bin/env bash

# useful for running though the null

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
    local_reference_dir
    output_dir
    code_dir
    sample_id
    bam_dir
    reference_fasta
    chr_sizes
    genes_gtf
    intervals_bed
    bam_file_end
    completion_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'bam_process' -- "$@")
eval set -- "${arguments}"

# Initialize optional vars to avoid nounset when not provided
intervals_bed=${intervals_bed-}
reference_fasta=${reference_fasta-}
chr_sizes=${chr_sizes-}

while true; do
    case "${1}" in
        --local_reference_dir )
            local_reference_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --sample_id )
            sample_id="${2}"; shift 2 ;;
        --bam_dir )
            bam_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --reference_fasta )
            reference_fasta="${2}"; shift 2 ;;
        --chr_sizes )
            chr_sizes="${2}"; shift 2 ;;
        --genes_gtf )
            genes_gtf="${2}"; shift 2 ;;
        --intervals_bed )
            intervals_bed="${2}"; shift 2 ;;
        --bam_file_end )
            bam_file_end="${2}"; shift 2 ;;
        --completion_dir )
            completion_dir="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

# activate the pixi enviroment
source <(pixi shell-hook --environment quantifybam --manifest-path ${code_dir}/pixi.toml)

# Construct BAM file path from bam_dir, sample_id, and bam_file_end
bam_file="${bam_dir}/${sample_id}.${bam_file_end}"

# Check that BAM file exists
if [[ ! -f ${bam_file} ]]; then
    echo "Error: BAM file ${bam_file} does not exist."
    exit 1
fi

participant_id=$(echo ${sample_id} | cut -d '-' -f1,2)

# make tmp dir
dir_prefix=${TMPDIR}/${sample_id}
mkdir -p ${dir_prefix}/output/genome_bam
mkdir -p ${dir_prefix}/tmp

# copy data to temp direcotry in compute node
rsync -PrhLtv ${bam_file} ${dir_prefix}/output/genome_bam
rsync -PrhLtv ${bam_file}.bai ${dir_prefix}/output/genome_bam

# Build rnaseqc args
rnaseqc_args=(
    --duplicate_marked_bam ${dir_prefix}/output/genome_bam/${sample_id}.${bam_file_end}
    --genes_gtf ${local_reference_dir}/${genes_gtf}
    --sample_id ${sample_id}
    --output_dir ${dir_prefix}/output/rnaseqc_null
)
if [[ -n "${reference_fasta:-}" ]]; then
    rnaseqc_args+=(--genome_fasta ${local_reference_dir}/${reference_fasta})
fi
if [[ -n "${intervals_bed:-}" ]]; then
    rnaseqc_args+=(--intervals_bed ${local_reference_dir}/${intervals_bed})
fi

# run rnaseq qc
bash ${code_dir}/run_03_rnaseqc.sh "${rnaseqc_args[@]}"

mkdir -p ${output_dir}/rnaseqc_null/
rsync -Prhltv ${dir_prefix}/output/rnaseqc_null/* ${output_dir}/rnaseqc_null/

# Create completion marker
completion_dir="${completion_dir:-${output_dir}/completed/rnaseqc_null}"
mkdir -p "${completion_dir}"
completion_file="${completion_dir}/${sample_id}.completed"
echo "Processing completed successfully for sample: ${sample_id}" > "${completion_file}"
echo "Completion marker created: ${completion_file}"
