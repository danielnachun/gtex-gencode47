#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/realign_test_bams.sh"
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    echo "Error: Config file $CONFIG_FILE not found!"
    exit 1
fi

# if true do all in gtex_ids
# if false, do all in gtex ids that do not already have a genome_bam in the output folder
regenerate_all=false

# make a bam list
mkdir -p ${output_dir}
mkdir -p ${output_dir}/logs
bam_list="$output_dir/bam_list"
> "$bam_list" # clear the file


if [ "$regenerate_all" = true ]; then
    # all gtex ids have a corresponding bam
    awk -v dir="$bam_dir" '{print dir "/" $0 ".Aligned.sortedByCoord.out.patched.md.bam"}' "$gtex_ids" > "$bam_list"
else
    # only add a gtex id if the genome bam does not already exist
    > "$bam_list"  
    while read -r gtex_id; do
        bam_file="${output_dir}/genome_bam/${gtex_id}.Aligned.sortedByCoord.out.patched.v11md.bam"
        if [ ! -f "$bam_file" ]; then
            echo "$bam_dir/$gtex_id.Aligned.sortedByCoord.out.patched.md.bam" >> "$bam_list"
        fi
    done < "$gtex_ids"
fi

original_count=$(wc -l < "${gtex_ids}")
to_process_count=$(wc -l < "${bam_list}")
completed_count=$((original_count - to_process_count))
echo "Original sample count: ${original_count}"
echo "Already completed: ${completed_count}"
echo "To be processed: ${to_process_count}"


sbatch --output "${output_dir}/logs/%A_%a.log" \
    --error "${output_dir}/logs/%A_%a.log" \
    --array "1-${to_process_count}" \
    --time 24:00:00 \
    --account smontgom \
    --partition batch \
    --cpus-per-task 1 \
    --mem 64G \
    --constraint="nvme" \
    --job-name realign_bam \
    ${code_dir}/realign_bam.sh \
        --reference_dir ${reference_dir} \
        --vcf_dir ${vcf_dir} \
        --output_dir ${output_dir} \
        --code_dir ${code_dir} \
        --bam_list ${bam_list} \
        --reference_fasta ${reference_fasta} \
        --rsem_ref_dir ${rsem_ref_dir} \
        --star_index ${star_index}
