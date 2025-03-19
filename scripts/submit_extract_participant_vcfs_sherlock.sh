#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/extract_vcfs_caudate.sh"
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    echo "Error: Config file $CONFIG_FILE not found!"
    exit 1
fi

mkdir -p ${output_dir}
mkdir -p ${output_dir}/logs

participant_id_length=$(wc -l < ${participant_id_list})
echo "BAM array length: ${participant_id_length}"

# check to see if vcf already exists
snps_vcf="${output_dir}/${participant_id}.snps.vcf.gz"
snps_vcf_index="${output_dir}/${participant_id}.snps.vcf.gz.tbi"
het_vcf="${output_dir}/${participant_id}.snps.het.vcf.gz"
het_vcf_index="${output_dir}/${participant_id}.snps.het.vcf.gz.tbi"

check_file() {
    [ -f "$1" ] && [ -s "$1" ]
}

# Check if all output files exist and are non-empty
if check_file "$snps_vcf" && check_file "$snps_vcf_index" && \
   check_file "$het_vcf" && check_file "$het_vcf_index"; then
    echo "$(date +"[%b %d %H:%M:%S] All output files already exist. Skipping processing.")"
    exit 0
fi




sbatch --output "${output_dir}/logs/%A_%a.log" \
    --error "${output_dir}/logs/%A_%a.log" \
    --array "1-${participant_id_length}%250" \
    --time 4:00:00 \
    --partition normal,owners \
    --cpus-per-task 1 \
    --mem 64G \
    --job-name participant_vcf \
    ${code_dir}/extract_participant_vcfs.sh \
        --participant_id_list ${participant_id_list} \
        --full_vcf ${full_vcf} \
        --output_dir ${output_dir} \
        --code_dir ${code_dir}
