#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/test/sprint_test_bams.sh"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# get all the realigned bam paths
bam_dir_bam_list=$(ls "$realign_bam_dir" | grep 'Aligned.sortedByCoord.out.patched.v11md.bam$')
# filter to only those in gtex_ids
full_bam_list=$(grep -F -f "$gtex_ids" <<< "$bam_dir_bam_list")
original_count=$(echo "$full_bam_list" | wc -l)
echo "Original sample count: ${original_count}"


# if true do all in gtex_ids
# if false, do all in gtex ids that do not already have a pileup output in the edsites_pileup folder (that is the last step)
regenerate_all=${regenerate_all:-false}
if [ "${regenerate_all}" = true ]; then
    # run all the bams in the input folder
    bams_to_sprint=$(sed "s|^|${realign_bam_dir}/|" <<< "$full_bam_list")
else
    # only sprint a bam if hyperedited reads output doesn't exist
    bams_to_sprint=$(grep -v -F -f <(ls "${output_dir}/sprint/unaligned/" | sed 's|\.hyperedited_reads\.sorted\.bam$|.Aligned.sortedByCoord.out.patched.v11md.bam|') <<< "$full_bam_list" | sed "s|^|${realign_bam_dir}/|")
fi

# Check if bams_to_sprint is empty
if [ -z "$bams_to_sprint" ]; then
    echo "To be processed: 0"
    echo "All bams realigned"
    exit 0  # Exit the script with a success status
else
    # If not empty, count the number of lines (bams to realign)
    to_process_count=$(echo "$bams_to_sprint" | wc -l)
    echo "To be processed: $to_process_count"
fi

completed_count=$((original_count - to_process_count))

# create a folder with a file per step, with one bam path per line in the file
bam_list_folder="$output_dir/file_lists/file_lists_sprint"
rm -rf "${bam_list_folder}"
mkdir -p "${bam_list_folder}"
split -l "${step_size}" --additional-suffix=".txt" <(echo "${bams_to_sprint}") "${bam_list_folder}/bam_list_" 

# create a file with one folder path per line
bam_list_paths="${output_dir}/file_lists/file_list_paths_sprint.txt"
rm -rf "${bam_list_paths}"
printf "%s\n" "${bam_list_folder}"/* > "${bam_list_paths}"
num_batches=$(wc -l < "${bam_list_paths}")
echo "Batches created: ${num_batches}"

# Check if num_batches is greater than max_array_size
if [ "${num_batches}" -gt "${max_array_size}" ]; then
    num_batches="${max_array_size}" 
fi

echo "Already completed: ${completed_count}"
echo "To be quantified: ${to_process_count}"
echo "Batches needed: $(( (to_process_count + step_size - 1) / step_size ))"
echo "Batches created: ${num_batches}"


sbatch_params=(
    --output "${output_dir}/logs/sprint/%A_%a.log"
    --error "${output_dir}/logs/sprint/%A_%a.log"
    --array "1-${num_batches}%250"
    --time 12:00:00
    --cpus-per-task "${step_size}"
    --mem 128G
    --job-name edsite_bam_batch
    ${code_dir}/batch_sprint_edsites.sh
        --reference_dir ${reference_dir}
        --output_dir ${output_dir}
        --code_dir ${code_dir}
        --bam_list_paths ${bam_list_paths}
        --reference_fasta ${reference_fasta}
        --repeat_bed ${repeat_bed}
        --step_size ${step_size}
)

# Submit on either sherlock or scg
if [ "${submit_on}" = 'sherlock' ]; then
    # Additional parameters for sherlock
    sbatch \
        --partition normal,owners \
        --tmp 200G \
        "${sbatch_params[@]}" 
elif [ "${submit_on}" = 'scg' ]; then
    # Additional parameters for scg
    sbatch \
        --account smontgom \
        --partition batch \
        --constraint="nvme" \
        "${sbatch_params[@]}" 
else
    echo "must submit on either 'sherlock' or 'scg'"
fi
