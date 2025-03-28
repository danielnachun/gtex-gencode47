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
#source <(pixi shell-hook --environment quantifybam --manifest-path ${code_dir}/pixi.toml)
# Define variables
max_attempts=5
attempt=1
wait_time=30  # seconds to wait between attempts

# Function to attempt pixi environment activation
activate_pixi_env() {
    echo "Attempt $attempt of $max_attempts: Activating pixi environment..."
    
    # Try to activate the environment using shell-hook
    if eval "$(pixi shell-hook --environment quantifybam --manifest-path ${code_dir}/pixi.toml)"; then
        # Check if environment was properly activated 
        if command -v rnaseqc >/dev/null 2>&1; then
            echo "Successfully activated pixi environment"
            return 0
        else
            echo "Environment activated, but rnaseqc not available"
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
sample_id=$(basename $(echo ${bam_file} | sed 's/\.Aligned\.sortedByCoord\.out\.patched\.v11md\.bam//'))
participant_id=$(echo ${sample_id} | cut -d '-' -f1,2)


# check for bam file
check_for_file "bam_file" "${bam_file}"
check_for_file "bam_file_index" "${bam_file}.bai"

# make tmp dir
dir_prefix=${TMPDIR}/${sample_id}
mkdir -p ${dir_prefix}/output/genome_bam
mkdir -p ${dir_prefix}/tmp

# copy references and data to temop direcotry in compute node
rsync -PrhLtv ${reference_dir}/* ${dir_prefix}/references/
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

mkdir -p ${output_dir}/rnaseq_qc/
rsync -Prhltv ${dir_prefix}/output/rnaseq_qc/* ${output_dir}/rnaseq_qc/
