#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/quantify_all_tissues.sh"
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    echo "Error: Config file $CONFIG_FILE not found!"
    exit 1
fi


mkdir -p ${output_dir}
mkdir -p ${output_dir}/logs


# get all the realigned bam paths
bam_dir_bam_list=$(ls "$realign_bam_dir" | grep 'Aligned.sortedByCoord.out.patched.v11md.bam$')
# filter to only those in gtex_ids
full_bam_list=$(grep -F -f "$gtex_ids" <<< "$bam_dir_bam_list")
original_count=$(echo "$full_bam_list" | wc -l)
echo "Original sample count: ${original_count}"


# if true do all in gtex_ids
# if false, do all in gtex ids that do not already have a regtools output in the leafcutter folder (that is the last step)
regenerate_all=${regenerate_all:-false}
if [ "${regenerate_all}" = true ]; then
    # run all the bams in the input folder
    bams_to_quantify="${full_bam_list}"
else
    # only quantify a bam if the regtools output does not already exist
    bams_to_quantify=$(grep -v -F -f <(ls "${output_dir}/leafcutter/" | sed 's|\.regtools_junc\.txt\.gz$|.Aligned.sortedByCoord.out.patched.v11md.bam|') <<< "$full_bam_list" | sed "s|^|${realign_bam_dir}/|")
fi

# Check if bams_to_quantify is empty
if [ -z "$bams_to_quantify" ]; then
    echo "To be processed: 0"
    echo "All bams realigned"
    exit 0  # Exit the script with a success status
else
    # If not empty, count the number of lines (bams to realign)
    to_process_count=$(echo "$bams_to_quantify" | wc -l)
    echo "To be processed: $to_process_count"
fi

completed_count=$((original_count - to_process_count))
echo "Already completed: ${completed_count}"
echo "Batches needed: $(( (to_process_count + step_size - 1) / step_size ))"

# create a folder with a file per step, with one bam path per line in the file
bam_list_folder="$output_dir/quantify_sherlock_bam_file_lists"
rm -rf "${bam_list_folder}"
mkdir -p "${bam_list_folder}"
split -l "${step_size}" --additional-suffix=".txt" <(echo "${bams_to_quantify}") "${bam_list_folder}/bam_list_" 

# create a file with one folder path per line
bam_list_paths="${output_dir}/quantify_sherlock_bam_list_paths.txt"
rm -rf "${bam_list_paths}"
printf "%s\n" "${bam_list_folder}"/* > "${bam_list_paths}"
num_batches=$(wc -l < "${bam_list_paths}")
echo "Batches created: ${num_batches}"


# Check if num_batches is greater than max_array_size
if [ "${num_batches}" -gt "${max_array_size}" ]; then
    num_batches="${max_array_size}" 
fi
echo "Batches running: ${num_batches}"

sbatch --output "${output_dir}/logs/%A_%a.log" \
        --error "${output_dir}/logs/%A_%a.log" \
        --array="1-${num_batches}%250" \
        --time 6:00:00 \
        --cpus-per-task="${step_size}" \
        --partition normal,owners \
        --mem=128G \
        --tmp=200G \
        --job-name qunatify_bam_batch \
        ${code_dir}/quantify_bam_batch.sh \
            --reference_dir ${reference_dir} \
            --vcf_dir ${vcf_dir} \
            --output_dir ${output_dir} \
            --code_dir ${code_dir} \
            --bam_list_paths ${bam_list_paths} \
            --reference_fasta ${reference_fasta} \
            --genes_gtf ${genes_gtf} \
            --chr_sizes ${chr_sizes} \
            --intervals_bed ${intervals_bed} \
            --step_size ${step_size}

