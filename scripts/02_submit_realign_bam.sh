#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/realign_truncated.sh"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }


# get all the v10 bam paths from bam dir
bam_dir_bam_list=$(ls "$bam_dir" | grep 'Aligned.sortedByCoord.out.patched.v11md.bam$')
# filter to only those in gtex_ids
full_bam_list=$(grep -F -f "$gtex_ids" <<< "$bam_dir_bam_list")
original_count=$(echo "$full_bam_list" | wc -l)
echo "Original sample count: ${original_count}"

# if true do all in gtex_ids
# if false, do all in gtex ids that do not already have a genome_bam in the output folder
regenerate_all=${regenerate_all:-false}
if [ "${regenerate_all}" = true ]; then
    # run all the bams in the input folder
    bams_to_realign=$(sed "s|^|${bam_dir}/|" <<< "$full_bam_list")
else
    # only realign a bam if the v11 genome bam does not already exist
    bams_to_realign=$(grep -v -F -f <(find "${output_dir}/genome_bam/" -name "*.v11md.bam" -exec basename {} \; | sed 's|\.v11md\.bam$|.md.bam|') <<< "$full_bam_list" | sed "s|^|${bam_dir}/|")
fi

# Check if bams_to_realign is empty
if [ -z "$bams_to_realign" ]; then
    echo "To be processed: 0"
    echo "All bams realigned"
    exit 0  # Exit the script with a success status
else
    # If not empty, count the number of lines (bams to realign)
    to_process_count=$(echo "$bams_to_realign" | wc -l)
    echo "To be processed: $to_process_count"
fi

completed_count=$((original_count - to_process_count))

# create a folder with a file per step, with one bam path per line in the file
bam_list_folder="$output_dir/file_lists/file_lists_realign"
rm -rf "${bam_list_folder}"
mkdir -p "${bam_list_folder}"
split -l "${step_size}" --additional-suffix=".txt" <(echo "${bams_to_realign}") "${bam_list_folder}/bam_list_" 

# create a file with one folder path per line
bam_list_paths="${output_dir}/file_lists/file_list_paths_realign.txt"
rm -rf "${bam_list_paths}"
printf "%s\n" "${bam_list_folder}"/* > "${bam_list_paths}"
num_batches=$(wc -l < "${bam_list_paths}")
echo "Batches created: ${num_batches}"

# Check if num_batches is greater than max_array_size
if [ "${num_batches}" -gt "${max_array_size}" ]; then
    num_batches="${max_array_size}" 
fi

echo "Already completed: ${completed_count}"
echo "To be realigned: ${to_process_count}"
echo "Batches needed: $(( (to_process_count + step_size - 1) / step_size ))"
echo "Batches created: ${num_batches}"


sbatch_params=(
    --output "${output_dir}/logs/realign/%A_%a.log"
    --error "${output_dir}/logs/realign/%A_%a.log"
    --array "1-${num_batches}%250"
    --time 12:00:00
    --cpus-per-task "${step_size}"
    --mem 128G
    --job-name realign_bam_batch
    ${code_dir}/batch_realign_bam.sh \
        --reference_dir ${reference_dir} \
        --vcf_dir ${vcf_dir} \
        --output_dir ${output_dir} \
        --code_dir ${code_dir} \
        --bam_list_paths ${bam_list_paths} \
        --reference_fasta ${reference_fasta} \
        --rsem_ref_dir ${rsem_ref_dir} \
        --star_index ${star_index} \
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
