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
    output_dir
    code_dir
    bam_file
    reference_fasta
    masked_genome_dir
    splicesites
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'bam_process' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --local_reference_dir )
            local_reference_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --bam_file )
            bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --reference_fasta )
            reference_fasta="${2}"; shift 2 ;;
        --masked_genome_dir )
            masked_genome_dir="${2}"; shift 2 ;;
        --splicesites )
            splicesites="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

sample_id=$(basename $(echo ${bam_file} | sed 's/\.Aligned\.sortedByCoord\.out\.patched\.v11md\.bam//'))
participant_id=$(echo ${sample_id} | cut -d '-' -f1,2)

# make tmp dir
dir_prefix=${TMPDIR}/${sample_id}
mkdir -p ${dir_prefix}/references/genome_bam
mkdir -p ${dir_prefix}/tmp

# copy bam data to temp direcotry in compute node
rsync -PrhLtv ${bam_file} ${dir_prefix}/references/genome_bam
rsync -PrhLtv ${bam_file}.bai ${dir_prefix}/references/genome_bam



# activate the pixi enviroment for extracting fastqs 
source <(pixi shell-hook --environment realignbam --manifest-path ${code_dir}/pixi.toml)

# running on just unliagned reads as a smaller test case for unstranded version 

# get unaligned reads from bam
bash ${code_dir}/run_get_unaligned_bam.sh \
    --bam_file ${dir_prefix}/references/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
    --sample_id ${sample_id} \
    --tmp_dir ${dir_prefix}/tmp/unaligned_bam

# process unaligned bams to fastqs
bash ${code_dir}/run_bam_to_fastq.sh \
    --bam_file ${dir_prefix}/tmp/unaligned_bam/${sample_id}.unaligned.bam \
    --sample_id ${sample_id} \
    --reference_fasta ${local_reference_dir}/${reference_fasta} \
    --tmp_dir ${dir_prefix}/tmp/unaligned_fastq

# do we get different resutls if we run on all reads (i.e. not just unaligned reads)?
# process all bam reads to fastqs
bash ${code_dir}/run_bam_to_fastq.sh \
    --bam_file ${dir_prefix}/references/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
    --sample_id ${sample_id} \
    --reference_fasta ${local_reference_dir}/${reference_fasta} \
    --tmp_dir ${dir_prefix}/tmp/fastq

# activate the pixi enviroment for hyperediting
source <(pixi shell-hook --environment calledsites --manifest-path ${code_dir}/pixi.toml)

# Run ADAR-bismark 
/home/klawren/oak/bioinformatics-tools/adar_bismark/Bismark/bismark --adar --hisat2 --non_directional --local \
  ${local_reference_dir}/${masked_genome_dir} \
  -1 ${dir_prefix}/tmp/unaligned_fastq/${sample_id}_1.fastq.gz \
  -2 ${dir_prefix}/tmp/unaligned_fastq/${sample_id}_2.fastq.gz \
  -o ${output_dir}/hyperediting/unaligned_adar_masked/ \
  --known-splicesite-infile ${local_reference_dir}/${splicesites} \
  -B ${sample_id}

# Run GA-bismark
/home/klawren/oak/bioinformatics-tools/adar_bismark/Bismark/bismark --custom_ref G --custom_alt A --hisat2 --non_directional --local \
  ${local_reference_dir}/${masked_genome_dir} \
  -1 ${dir_prefix}/tmp/unaligned_fastq/${sample_id}_1.fastq.gz \
  -2 ${dir_prefix}/tmp/unaligned_fastq/${sample_id}_2.fastq.gz \
  -o ${output_dir}/hyperediting/unaligned_GA_masked/ \
  --known-splicesite-infile ${local_reference_dir}/${splicesites} \
  -B ${sample_id}




# Run ADAR-bismark 
/home/klawren/oak/bioinformatics-tools/adar_bismark/Bismark/bismark --adar --hisat2 --non_directional --local \
  ${local_reference_dir}/${masked_genome_dir} \
  -1 ${dir_prefix}/tmp/fastq/${sample_id}_1.fastq.gz \
  -2 ${dir_prefix}/tmp/fastq/${sample_id}_2.fastq.gz \
  -o ${output_dir}/hyperediting/adar_masked/ \
  --known-splicesite-infile ${local_reference_dir}/${splicesites} \
  -B ${sample_id}

# Run GA-bismark
/home/klawren/oak/bioinformatics-tools/adar_bismark/Bismark/bismark --custom_ref G --custom_alt A --hisat2 --non_directional --local \
  ${local_reference_dir}/${masked_genome_dir} \
  -1 ${dir_prefix}/tmp/fastq/${sample_id}_1.fastq.gz \
  -2 ${dir_prefix}/tmp/fastq/${sample_id}_2.fastq.gz \
  -o ${output_dir}/hyperediting/GA_masked/ \
  --known-splicesite-infile ${local_reference_dir}/${splicesites} \
  -B ${sample_id}

echo "Copying out results"
rsync -Prhltv ${dir_prefix}/output/ ${output_dir}