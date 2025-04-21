#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/realign_all_tissues_copied_only.sh"
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    echo "Error: Config file $CONFIG_FILE not found!"
    exit 1
fi


mkdir -p ${output_dir}
mkdir -p ${output_dir}/logs

# get all the v10 bam paths
full_bam_list=$(ls "$bam_dir" | grep 'Aligned.sortedByCoord.out.patched.md.bam$')
original_count=$(echo "$full_bam_list" | wc -l)
echo "Original sample count: ${original_count}"

# if true do all in gtex_ids
# if false, do all in gtex ids that do not already have a genome_bam in the output folder
regenerate_all=${regenerate_all:-false}
if [ "${regenerate_all}" = true ]; then
    # run all the bams in the input folder
    bams_to_realign="${full_bam_list}"
else
    # only realign a bam if the v11 genome bam does not already exist
    bams_to_realign=$(grep -v -F -f <(ls "${output_dir}/genome_bam/") <(sed 's|\.md\.bam$|.v11md.bam|' <<< "$full_bam_list"))
fi
bams_to_realign=$(grep -v -F -f <(ls "${output_dir}/genome_bam/" | sed 's|\.v11md\.bam$|.md.bam|') <<< "$full_bam_list" | sed "s|^|${bam_dir}/|")

to_process_count=$(echo "$bams_to_realign" | wc -l)
echo "To be processed: ${to_process_count}"
completed_count=$((original_count - to_process_count))
echo "Already completed: ${completed_count}"

echo "Batches needed: $(( (to_process_count + step_size - 1) / step_size ))"

# create a folder with a file per step, with one bam path per line in the file
bam_list_folder="$output_dir/bam_file_lists"
rm -rf "${bam_list_folder}"
mkdir -p "${bam_list_folder}"
split -l "${step_size}" --additional-suffix=".txt" <(echo "${bams_to_realign}") "${bam_list_folder}/bam_list_" 

# create a file with one folder path per line
bam_list_paths="${output_dir}/bam_list_paths.txt"
rm -rf "${bam_list_paths}"
printf "%s\n" "${bam_list_folder}"/* > "${bam_list_paths}"
num_batches=$(wc -l < "${bam_list_paths}")
echo "Batches created: ${num_batches}"

# Check if num_batches is greater than max_array_size
if [ "${num_batches}" -gt "${max_array_size}" ]; then
    num_batches="${max_array_size}" 
fi
echo "Batches running: ${num_batches}"

sbatch --output="${output_dir}/logs/%A_%a.log" \
            --error="${output_dir}/logs/%A_%a.log" \
            --array="1-${num_batches}%250" \
            --time=24:00:00 \
            --cpus-per-task="${step_size}" \
            --partition=normal,owners \
            --mem=256G \
            --tmp=250G \
            --job-name=realign_bam_batch \
            ${code_dir}/realign_bam_batch.sh \
                --reference_dir ${reference_dir} \
                --vcf_dir ${vcf_dir} \
                --output_dir ${output_dir} \
                --code_dir ${code_dir} \
                --bam_list_paths ${bam_list_paths} \
                --reference_fasta ${reference_fasta} \
                --rsem_ref_dir ${rsem_ref_dir} \
                --star_index ${star_index} \
                --step_size ${step_size}
