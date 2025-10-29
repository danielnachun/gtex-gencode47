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
    rsem_ref_dir
    star_index
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
            output_dir="${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --bam_file )
            bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
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

sample_id=$(basename $(echo ${bam_file} | sed 's/\.Aligned\.sortedByCoord\.out\.patched\.v11md\.bam//'))
participant_id=$(echo ${sample_id} | cut -d '-' -f1,2)


# make tmp dir
dir_prefix=${TMPDIR}/${sample_id}
mkdir -p ${dir_prefix}/raw
mkdir -p ${dir_prefix}/tmp

# copy data to temp direcotry in compute node
rsync -PrhLtv ${bam_file} ${dir_prefix}/raw
rsync -PrhLtv ${bam_file}.bai ${dir_prefix}/raw

# process bams to fastqs
bash ${code_dir}/run_02_bam_to_fastq.sh --bam_file ${dir_prefix}/raw/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
    --sample_id ${sample_id} \
    --reference_fasta ${local_reference_dir}/${reference_fasta} \
    --tmp_dir ${dir_prefix}/tmp/fastq


# Check for VCF files and run star with WASP if the exist
vcf_file=${participant_id}.snps.vcf.gz 
vcf_index=${participant_id}.snps.vcf.gz.tbi
vcf_path="${vcf_dir}/${vcf_file}"
vcf_index_path="${vcf_dir}/${vcf_index}"

if check_vcf_file "$vcf_path" "VCF" && check_vcf_file "$vcf_index_path" "VCF index"; then
    echo "VCF files found. Running STAR with WASP..."
    vcf_dir_tmp=${dir_prefix}/vcfs
    mkdir -p ${vcf_dir_tmp}
    rsync -PrhLtv ${vcf_dir}/${vcf_file} ${vcf_dir_tmp}
    rsync -PrhLtv ${vcf_dir}/${vcf_file}.tbi ${vcf_dir_tmp}

    # align with star
    bash ${code_dir}/run_02_star.sh \
        --star_index ${local_reference_dir}/${star_index} \
        --fastq_1 ${dir_prefix}/tmp/fastq/${sample_id}_1.fastq.gz \
        --fastq_2 ${dir_prefix}/tmp/fastq/${sample_id}_2.fastq.gz \
        --sample_id ${sample_id} \
        --vcf_file ${vcf_dir_tmp}/${vcf_file} \
        --tmp_dir ${dir_prefix}/tmp/star
else
    echo "Warning: VCF files not found, running STAR without WASP..."
    # align with star
    bash ${code_dir}/run_02_star_no_vcf.sh \
        --star_index ${local_reference_dir}/${star_index} \
        --fastq_1 ${dir_prefix}/tmp/fastq/${sample_id}_1.fastq.gz \
        --fastq_2 ${dir_prefix}/tmp/fastq/${sample_id}_2.fastq.gz \
        --sample_id ${sample_id} \
        --tmp_dir ${dir_prefix}/tmp/star
fi

# # sync bams
# bash ${code_dir}/run_02_bam_sync.sh \
#     --initial_bam_file ${dir_prefix}/raw/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
#     --star_aligned_bam ${dir_prefix}/tmp/star/${sample_id}.Aligned.sortedByCoord.out.bam \
#     --sample_id ${sample_id} \
#     --tmp_dir ${dir_prefix}/tmp/bamsync \
#     --output_dir ${dir_prefix}/output/flagstat

# # mark duplicates, get the genome bam that we save
# bash ${code_dir}/run_02_mark_duplicates.sh \
#     --genome_bam_file ${dir_prefix}/tmp/bamsync/${sample_id}.Aligned.sortedByCoord.out.patched.bam \
#     --output_prefix ${sample_id}.Aligned.sortedByCoord.out.patched.v11md \
#     --output_dir ${dir_prefix}/output/genome_bam

# run rsem, save isoform quantification
bash ${code_dir}/run_02_rsem.sh \
    --rsem_ref_dir ${local_reference_dir}/${rsem_ref_dir} \
    --transcriptome_bam ${dir_prefix}/tmp/star/${sample_id}.Aligned.toTranscriptome.out.bam \
    --sample_id ${sample_id} \
    --output_dir ${dir_prefix}/output/rsem

echo "Copying out results"
rsync -Prhltv ${dir_prefix}/output/ ${output_dir}

# Copy out all STAR result files except the sorted BAM and its index, handle gzipped & non-gzipped files
echo "Copying out STAR results"
mkdir -p "${output_dir}/star"
for ext in Aligned.toTranscriptome.out.bam \
           Chimeric.out.junction \
           Chimeric.out.junction.gz \
           Log.final.out \
           Log.out \
           Log.progress.out \
           ReadsPerGene.out.tab \
           ReadsPerGene.out.tab.gz \
           SJ.out.tab \
           SJ.out.tab.gz \
           _STARpass1
do
    file_path="${dir_prefix}/tmp/star/${sample_id}.${ext}"
    if [[ -e "${file_path}" ]]; then
        rsync -Prhltv "${file_path}" "${output_dir}/star/"
    else
        echo "Warning: ${file_path} not found. Skipping."
    fi
done

# Create completion marker
completion_dir="${output_dir}/completed/realign"
mkdir -p "${completion_dir}"
completion_file="${completion_dir}/${sample_id}.completed"
echo "Processing completed successfully for sample: ${sample_id}" > "${completion_file}"
echo "Completion marker created: ${completion_file}"

echo "Done"

