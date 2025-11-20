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

# check for file, warn instead of exiting
check_vcf_file() {
    file_path="${1}"
    file_desc="${2}"
    if [[ ! -f ${file_path} ]]; then
        echo "Warning: ${file_desc} file ${file_path} does not exist."
        return 1
    fi
    return 0
}

options_array=(
    local_reference_dir
    vcf_dir
    output_dir
    code_dir
    bam_file
    reference_fasta
    chr_sizes
    genes_gtf
    intervals_bed
    ipa_annotation
    editing_bed
    bam_file_end
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'bam_process' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --local_reference_dir )
            local_reference_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --vcf_dir )
            vcf_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --bam_file )
            bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --reference_fasta )
            reference_fasta="${2}"; shift 2 ;;
        --chr_sizes )
            chr_sizes="${2}"; shift 2 ;;
        --genes_gtf )
            genes_gtf="${2}"; shift 2 ;;
        --intervals_bed )
            intervals_bed="${2}"; shift 2 ;;
        --ipa_annotation )
            ipa_annotation="${2}"; shift 2 ;;
        --editing_bed )
            editing_bed="${2}"; shift 2 ;;
        --bam_file_end )
            bam_file_end="${2}"; shift 2 ;;
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
sample_id=$(basename $(echo ${bam_file} | sed "s/.${bam_file_end}//"))
participant_id=$(echo ${sample_id} | cut -d '-' -f1,2)


# make tmp dir
dir_prefix=${TMPDIR}/${sample_id}
mkdir -p ${dir_prefix}/references/genome_bam
mkdir -p ${dir_prefix}/tmp

# copy data to temp direcotry in compute node
rsync -PrhLtv ${bam_file} ${dir_prefix}/references/genome_bam
rsync -PrhLtv ${bam_file}.bai ${dir_prefix}/references/genome_bam

# run rnaseq qc
bash ${code_dir}/run_03_rnaseqc.sh \
    --duplicate_marked_bam ${dir_prefix}/references/genome_bam/${sample_id}.${bam_file_end} \
    --genes_gtf ${local_reference_dir}/${genes_gtf} \
    --genome_fasta ${local_reference_dir}/${reference_fasta} \
    --intervals_bed ${local_reference_dir}/${intervals_bed} \
    --sample_id ${sample_id} \
    --output_dir ${dir_prefix}/output/rnaseqc

# copy out results
rsync -Prhltv ${dir_prefix}/output/ ${output_dir}

# run coverage
bash ${code_dir}/run_03_bam_to_coverage.sh \
    --duplicate_marked_bam ${dir_prefix}/references/genome_bam/${sample_id}.${bam_file_end} \
    --chr_sizes ${local_reference_dir}/${chr_sizes} \
    --sample_id ${sample_id} \
    --output_dir ${dir_prefix}/output/coverage

# copy out results
rsync -Prhltv ${dir_prefix}/output/ ${output_dir}


# Check for VCF files and run GATK only if they exist
vcf_file=${participant_id}.snps.vcf.gz 
vcf_index=${participant_id}.snps.vcf.gz.tbi
vcf_path="${vcf_dir}/${vcf_file}"
vcf_index_path="${vcf_dir}/${vcf_index}"

if check_vcf_file "$vcf_path" "VCF" && check_vcf_file "$vcf_index_path" "VCF index"; then
    echo "VCF files found. Running GATK step..."
    vcf_dir_tmp=${dir_prefix}/vcfs
    mkdir -p ${vcf_dir_tmp}
    rsync -PrhLtv ${vcf_dir}/${vcf_file} ${vcf_dir_tmp}
    rsync -PrhLtv ${vcf_dir}/${vcf_file}.tbi ${vcf_dir_tmp}

    # run gatk ase
    bash ${code_dir}/run_03_gatk_ase.sh \
        --sample_id ${sample_id} \
        --dir_prefix ${dir_prefix} \
        --genome_fasta ${local_reference_dir}/${reference_fasta} \
        --vcf_dir ${vcf_dir} \
        --duplicate_marked_bam ${dir_prefix}/references/genome_bam/${sample_id}.${bam_file_end} \
        --output_dir ${dir_prefix}/output/ase
else
    echo "Warning: Skipping GATK step because the required VCF files are missing."
fi

# copy out results
rsync -Prhltv ${dir_prefix}/output/ ${output_dir}

# run regtools
bash ${code_dir}/run_03_regtools.sh \
    --sample_id ${sample_id} \
    --dir_prefix ${dir_prefix} \
    --duplicate_marked_bam ${dir_prefix}/references/genome_bam/${sample_id}.${bam_file_end} \
    --output_dir ${dir_prefix}/output/leafcutter

# copy out results
rsync -Prhltv ${dir_prefix}/output/ ${output_dir}

# run fraser quantification
bash ${code_dir}/run_03_fraser.sh \
    --duplicate_marked_bam ${dir_prefix}/references/genome_bam/${sample_id}.${bam_file_end} \
    --sample_id ${sample_id} \
    --output_dir ${dir_prefix}/output/fraser \
    --working_dir ${local_reference_dir}/${sample_id} \
    --code_dir ${code_dir} 

# copy out results
rsync -Prhltv ${dir_prefix}/output/ ${output_dir}

# run ipafinder
bash ${code_dir}/run_03_ipafinder.sh \
    --duplicate_marked_bam ${dir_prefix}/references/genome_bam/${sample_id}.${bam_file_end} \
    --sample_id ${sample_id} \
    --ipa_annotation ${local_reference_dir}/${ipa_annotation} \
    --output_dir ${dir_prefix}/output/ipafinder

# copy out results
rsync -Prhltv ${dir_prefix}/output/ ${output_dir}

# # run edsite quantification
# bash ${code_dir}/run_edsite_pileup.sh \
#     --duplicate_marked_bam ${dir_prefix}/references/genome_bam/${sample_id}.${bam_file_end} \
#     --sample_id ${sample_id} \
#     --reference_fasta ${local_reference_dir}/${reference_fasta} \
#     --editing_bed ${local_reference_dir}/${editing_bed} \
#     --output_dir ${dir_prefix}/output/rnaediting

# # copy out results
rsync -Prhltv ${dir_prefix}/output/ ${output_dir}

# Create completion marker file to indicate successful processing
completion_dir="${output_dir}/completed/quantifications"
mkdir -p "${completion_dir}"
completion_file="${completion_dir}/${sample_id}.completed"
echo "Processing completed successfully for sample: ${sample_id}" > "${completion_file}"
echo "Completion marker created: ${completion_file}"
