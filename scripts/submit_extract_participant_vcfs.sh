#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/extract_vcfs_all_tissues.sh"
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    echo "Error: Config file $CONFIG_FILE not found!"
    exit 1
fi

# Check for duplicates in participant_id_list, duplicated can cuase errors as files colide
if [ $(sort "${participant_id_list}" | uniq -d | wc -l) -gt 0 ]; then
    echo "Error: Duplicate entries found in ${participant_id_list}:"
    sort "${participant_id_list}" | uniq -d
    exit 1
fi

# if true do all in participant_id_list
# if false, do all in participant_id_list that do not already vcfs and vcf indexs
regenerate_all=false

mkdir -p ${output_dir}
mkdir -p ${output_dir}/logs


# Create new participant list for processing
new_participant_list="${output_dir}/participants_to_process.txt"
> "${new_participant_list}"  # Clear the file

# Function to check if file exists and is non-empty
check_file() {
    [ -f "$1" ] && [ -s "$1" ]
}

# Process each participant
while read -r participant_id; do
    # Define expected output files
    snps_vcf="${output_dir}/${participant_id}.snps.vcf.gz"
    snps_vcf_index="${output_dir}/${participant_id}.snps.vcf.gz.tbi"
    het_vcf="${output_dir}/${participant_id}.snps.het.vcf.gz"
    het_vcf_index="${output_dir}/${participant_id}.snps.het.vcf.gz.tbi"

    if [ "$regenerate_all" = true ]; then
        # Add all participants if regenerating everything
        echo "${participant_id}" >> "${new_participant_list}"
    else
        # Check if any required file is missing or empty
        if ! check_file "$snps_vcf" || \
           ! check_file "$snps_vcf_index" || \
           ! check_file "$het_vcf" || \
           ! check_file "$het_vcf_index"; then
            echo "${participant_id}" >> "${new_participant_list}"
        fi
    fi
done < "${participant_id_list}"

original_count=$(wc -l < "${participant_id_list}")
to_process_count=$(wc -l < "${new_participant_list}")
completed_count=$((original_count - to_process_count))
echo "Original participant count: ${original_count}"
echo "Already completed: ${completed_count}"
echo "To be processed: ${to_process_count}"

# submit on either sherlock or scg
if [ "${submit_on}" = 'sherlock' ]; then
    # submit on sherlock
    sbatch --output "${output_dir}/logs/%A_%a.log" \
        --error "${output_dir}/logs/%A_%a.log" \
        --array "1-${to_process_count}%250" \
        --time 4:00:00 \
        --partition normal,owners \
        --cpus-per-task 1 \
        --mem 64G \
        --job-name participant_vcf \
        ${code_dir}/extract_participant_vcfs.sh \
            --participant_id_list ${new_participant_list} \
            --full_vcf ${full_vcf} \
            --output_dir ${output_dir} \
            --code_dir ${code_dir}
elif [ "${submit_on}" = 'scg' ]; then
    # submit on scg
    sbatch --output "${output_dir}/logs/%A_%a.log" \
        --error "${output_dir}/logs/%A_%a.log" \
        --array "1-${to_process_count}%250" \
        --time 4:00:00 \
        --account smontgom \
        --partition batch \
        --cpus-per-task 1 \
        --mem 64G \
        --job-name participant_vcf \
        ${code_dir}/extract_participant_vcfs.sh \
            --participant_id_list ${new_participant_list} \
            --full_vcf ${full_vcf} \
            --output_dir ${output_dir} \
            --code_dir ${code_dir}
else
    echo "must submit on either 'sherlock' or 'scg'"
fi
