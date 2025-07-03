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
    output_dir
    code_dir
    vcf_dir
    bam_file
    reference_fasta
    dbsnp
    indels_mills
    indels_decoy
    gene_intervals_bed
    full_vcf_file
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'call_edsites' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --local_reference_dir )
            local_reference_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --vcf_dir )
            vcf_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --bam_file )
            bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --reference_fasta )
            reference_fasta="${2}"; shift 2 ;;
        --dbsnp )
            dbsnp="${2}"; shift 2 ;;
        --indels_mills )
            indels_mills="${2}"; shift 2 ;;
        --indels_decoy )
            indels_decoy="${2}"; shift 2 ;;
        --gene_intervals_bed )
            gene_intervals_bed="${2}"; shift 2 ;;
        --full_vcf_file )
            full_vcf_file="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

# activate the pixi enviroment
max_attempts=5
attempt=1
wait_time=30  # seconds to wait between attempts

# Function to attempt pixi environment activation
activate_pixi_env() {
    echo "Attempt $attempt of $max_attempts: Activating pixi environment..."
    
    # Try to activate the environment using shell-hook
    # TODO make pixi
    if eval "$(pixi shell-hook --environment calledsites --manifest-path ${code_dir}/pixi.toml)"; then
        # Check if environment was properly activated 
        if command -v gatk >/dev/null 2>&1; then
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
sample_id=$(basename $(echo ${bam_file} | sed 's/\.Aligned\.sortedByCoord\.out\.patched\.v11md\.bam//'))
participant_id=$(echo ${sample_id} | cut -d '-' -f1,2)


# make tmp dir
dir_prefix=${TMPDIR}/${sample_id}
mkdir -p ${dir_prefix}/references/genome_bam
mkdir -p ${dir_prefix}/tmp

# copy data to temp direcotry in compute node
rsync -PrhLtv ${bam_file} ${dir_prefix}/references/genome_bam
rsync -PrhLtv ${bam_file}.bai ${dir_prefix}/references/genome_bam

echo $(date +"[%b %d %H:%M:%S] Callind edSites for ${sample_id}")

# run SplitNCigarReads
# 2 hour run time
bash ${code_dir}/run_split_reads.sh \
    --duplicate_marked_bam ${dir_prefix}/references/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
    --sample_id ${sample_id} \
    --reference_fasta ${local_reference_dir}/${reference_fasta} \
    --output_dir ${dir_prefix}/tmp/split_bam

# run bqsr
# 1 hour 15 minute run time
# then 22 minute run time
bash ${code_dir}/run_bqsr.sh \
    --split_bam ${dir_prefix}/tmp/split_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.split.bam \
    --sample_id ${sample_id} \
    --reference_fasta ${local_reference_dir}/${reference_fasta} \
    --dbsnp ${local_reference_dir}/${dbsnp} \
    --indels_mills ${local_reference_dir}/${indels_mills} \
    --indels_decoy ${local_reference_dir}/${indels_decoy} \
    --tmp_dir ${dir_prefix}/tmp/bqsr_recal_data \
    --output_dir ${dir_prefix}/output/bqsr


# Check for VCF files and run mutect with participant vcf if it exists, otherwise full vcf
vcf_file=${participant_id}.snps.vcf.gz 
vcf_index=${participant_id}.snps.vcf.gz.tbi
vcf_path="${vcf_dir}/${vcf_file}"
vcf_index_path="${vcf_dir}/${vcf_index}"

if check_vcf_file "$vcf_path" "VCF" && check_vcf_file "$vcf_index_path" "VCF index"; then
    echo "VCF files found. Running MUTECT2 with PON from participant VCF..."
    vcf_dir_tmp=${dir_prefix}/vcfs
    mkdir -p ${vcf_dir_tmp}
    rsync -PrhLtv ${vcf_dir}/${vcf_file} ${vcf_dir_tmp}
    rsync -PrhLtv ${vcf_dir}/${vcf_file}.tbi ${vcf_dir_tmp}

    bash ${code_dir}/run_mutect.sh \
        --bqsr_bam ${dir_prefix}/output/bqsr/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.recalibrated.bam \
        --sample_id ${sample_id} \
        --reference_fasta ${local_reference_dir}/${reference_fasta} \
        --gene_intervals_bed ${local_reference_dir}/${gene_intervals_bed} \
        --vcf_file ${vcf_dir_tmp}/${vcf_file} \
        --output_dir ${dir_prefix}/output/mutect

else
    echo "Warning: VCF files not found, running MUTECT with PON from all participant combined VCF ..."
    bash ${code_dir}/run_mutect.sh \
        --bqsr_bam ${dir_prefix}/output/bqsr/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.recalibrated.bam \
        --sample_id ${sample_id} \
        --reference_fasta ${local_reference_dir}/${reference_fasta} \
        --gene_intervals_bed ${local_reference_dir}/${gene_intervals_bed} \
        --vcf_file ${local_reference_dir}/${full_vcf_file} \
        --output_dir ${dir_prefix}/output/mutect
fi


# 44 minute run time
bash ${code_dir}/run_haplotype_caller.sh \
    --bqsr_bam ${dir_prefix}/output/bqsr/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.recalibrated.bam \
    --sample_id ${sample_id} \
    --reference_fasta ${local_reference_dir}/${reference_fasta} \
    --gene_intervals_bed ${local_reference_dir}/${gene_intervals_bed} \
    --dbsnp ${local_reference_dir}/${dbsnp} \
    --output_dir ${dir_prefix}/output/haplotype_caller


rsync -Prhltv ${dir_prefix}/output/ ${output_dir}
