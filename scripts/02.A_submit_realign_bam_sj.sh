#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# Source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/02.A_realign_star_outputs.sh"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Define function to get files to process
get_files_to_process() {
    # Get all the v10 bam paths from bam dir
    bam_dir_bam_list=$(ls "$bam_dir" | grep 'Aligned.sortedByCoord.out.patched.v11md.bam$')
    # Filter to only those in gtex_ids
    full_bam_list=$(grep -F -f "$gtex_ids" <<< "$bam_dir_bam_list")
    original_count=$(echo "$full_bam_list" | wc -l)
    echo "Original sample count: ${original_count}"
    
    # Determine which BAMs to process based on regenerate_all setting
    if [ "${regenerate_all}" = true ]; then
        # Process all BAMs in input folder
        echo "$full_bam_list" | sed "s|^|${bam_dir}/|"
    else
        # Only process BAMs that don't have completion markers
        completion_dir="${output_dir}/completed/star"
        all_completion_files=$(printf "%s\n" ${full_bam_list} | sed "s|^|${completion_dir}/|;s|\.Aligned\.sortedByCoord\.out\.patched\.v11md\.bam$|.completed|")
        printf "%s\n" ${full_bam_list} | paste - <(printf "%s\n" ${all_completion_files}) | awk '{if(system("[ -f \""$2"\" ]")==0) next; print $1}' | sed "s|^|${bam_dir}/|"
    fi
}

# Call the shared batch submission utility
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    "${CONFIG_FILE}" \
    "${code_dir}/utils/batch_process_files.sh" \
    "realign_bam_sj" \
    "${code_dir}/02.A_realign_bam_sj.sh" \
    --reference_dir "${reference_dir}" \
    --vcf_dir "${vcf_dir}" \
    --output_dir "${output_dir}" \
    --code_dir "${code_dir}" \
    --reference_fasta "${reference_fasta}" \
    --star_index "${star_index}"
