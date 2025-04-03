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
#source <(pixi shell-hook --environment realignbam --manifest-path ${code_dir}/pixi.toml)

# Define variables
max_attempts=5
attempt=1
wait_time=30  # seconds to wait between attempts

# Function to attempt pixi environment activation
activate_pixi_env() {
    echo "Attempt $attempt of $max_attempts: Activating pixi environment..."
    
    # Try to activate the environment using shell-hook
    if eval "$(pixi shell-hook --environment realignbam --manifest-path ${code_dir}/pixi.toml)"; then
        # Check if environment was properly activated
        if command -v picard >/dev/null 2>&1; then
            echo "Successfully activated pixi environment"
            return 0
        else
            echo "Environment activated, but picard not available"
            return 1
        fi
    else
        echo "Failed to activate pixi environment"
        return 1
    fi
}

# Main loop to try activation with retries
while [ $attempt -le $max_attempts ]; do
    if activate_pixi_env; then
        # Success - break out of the loop
        break
    else
        # Failed attempt
        if [ $attempt -eq $max_attempts ]; then
            echo "ERROR: Failed to activate pixi environment after $max_attempts attempts."
            echo "Please check your pixi installation and environment configuration."
            exit 1
        fi
        
        # Increment attempt counter and wait before retry (can't use sleep due to sherlock requirements)
        echo "Waiting $wait_time seconds before retry..."
        end_time=$(($(date +%s) + $wait_time))
        while [ $(date +%s) -lt $end_time ]; do
            # Perform trivial computation 
            for i in {1..1000}; do echo $i > /dev/null; done
        done
        attempt=$((attempt+1))
    fi
done

echo "Pixi environment activated successfully after $attempt attempt(s)."



# map job id to line number and then to sample id
line_number=${SLURM_ARRAY_TASK_ID}
bam_file="$(sed "${line_number}q; d" "${bam_list}")"
sample_id=$(basename $(echo ${bam_file} | sed 's/\.Aligned\.sortedByCoord\.out\.patched\.md\.bam//'))
participant_id=$(echo ${sample_id} | cut -d '-' -f1,2)

# check for bam file
check_for_file "bam_file" "${bam_file}"
check_for_file "bam_file_index" "${bam_file}.bai"

# make tmp dir
dir_prefix=${TMPDIR}/${sample_id}
mkdir -p ${dir_prefix}/raw
mkdir -p ${dir_prefix}/tmp

# copy references and data to temop direcotry in compute node
rsync -PrhLtv ${reference_dir}/* ${dir_prefix}/references/
rsync -PrhLtv ${bam_file} ${dir_prefix}/raw
rsync -PrhLtv ${bam_file}.bai ${dir_prefix}/raw

# process bams to fastqs
bash ${code_dir}/run_bam_to_fastq.sh --bam_file ${dir_prefix}/raw/${sample_id}.Aligned.sortedByCoord.out.patched.md.bam \
    --sample_id ${sample_id} \
    --reference_fasta ${dir_prefix}/references/${reference_fasta} \
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
    bash ${code_dir}/run_fastq_to_star.sh \
        --star_index ${dir_prefix}/references/${star_index} \
        --fastq_1 ${dir_prefix}/tmp/fastq/${sample_id}_1.fastq.gz \
        --fastq_2 ${dir_prefix}/tmp/fastq/${sample_id}_2.fastq.gz \
        --sample_id ${sample_id} \
        --vcf_file ${vcf_dir_tmp}/${vcf_file} \
        --tmp_dir ${dir_prefix}/tmp/star
else
    echo "Warning: VCF files not found, running STAR without WASP..."
    # align with star
    bash ${code_dir}/run_fastq_to_star_no_vcf.sh \
        --star_index ${dir_prefix}/references/${star_index} \
        --fastq_1 ${dir_prefix}/tmp/fastq/${sample_id}_1.fastq.gz \
        --fastq_2 ${dir_prefix}/tmp/fastq/${sample_id}_2.fastq.gz \
        --sample_id ${sample_id} \
        --tmp_dir ${dir_prefix}/tmp/star
fi

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
bash ${code_dir}/run_rsem.sh \
    --rsem_ref_dir ${dir_prefix}/references/${rsem_ref_dir} \
    --transcriptome_bam ${dir_prefix}/tmp/star/${sample_id}.Aligned.toTranscriptome.out.bam \
    --sample_id ${sample_id} \
    --output_dir ${dir_prefix}/output/rsem


rsync -Prhltv ${dir_prefix}/output/ ${output_dir}

# delete the corresponding bam in the extra bam folder
extra_bam_dir="/home/klawren/oak/gtex/data/raw/bam_copy"
extra_bam_file="${extra_bam_dir}/${sample_id}.Aligned.sortedByCoord.out.patched.md.bam"
extra_bam_index="${extra_bam_dir}/${sample_id}.Aligned.sortedByCoord.out.patched.md.bam.bai"
rm -f ${extra_bam_file}
rm -f ${$extra_bam_index}
