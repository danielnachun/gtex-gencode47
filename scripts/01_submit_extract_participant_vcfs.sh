#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/01_extract_vcfs_all.sh"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }


# Check for duplicates in participant_id_list, duplicated can cuase errors as files colide
if [ $(sort "${participant_id_list}" | uniq -d | wc -l) -gt 0 ]; then
    echo "Error: Duplicate entries found in ${participant_id_list}:"
    sort "${participant_id_list}" | uniq -d
    exit 1
fi

# if true do all in participant_id_list
# if false, do all in participant_id_list that do not already vcfs and vcf indexs
regenerate_all=false

logs_dir="${output_dir}/logs/extract_participant_vcfs"
mkdir -p "${logs_dir}"


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

# Exit early if nothing to process
if [ "${to_process_count}" -eq 0 ]; then
    echo "All participants processed"
    exit 0
fi

# Build shared sbatch parameters
sbatch_params=(
    --output "${logs_dir}/%A_%a.log"
    --error "${logs_dir}/%A_%a.log"
    --array "1-${to_process_count}%250"
    --time 4:00:00
    --cpus-per-task 1
    --mem 64G
    --job-name participant_vcf
    ${code_dir}/01_extract_participant_vcfs.sh \
        --participant_id_list ${new_participant_list} \
        --full_vcf ${full_vcf} \
        --output_dir ${output_dir} \
        --code_dir ${code_dir}
)

# Submit on either sherlock or scg
if [ "${submit_on}" = 'sherlock' ]; then
    sbatch \
        --partition normal,owners \
        --tmp 200G \
        "${sbatch_params[@]}" 
elif [ "${submit_on}" = 'scg' ]; then
    sbatch \
        --account smontgom \
        --partition batch \
        --constraint="nvme" \
        "${sbatch_params[@]}" 
else
    echo "must submit on either 'sherlock' or 'scg'"
fi
