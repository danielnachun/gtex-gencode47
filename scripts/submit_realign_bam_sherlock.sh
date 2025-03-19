#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/realign_caudate
.sh"
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    echo "Error: Config file $CONFIG_FILE not found!"
    exit 1
fi


# if true do all in gtex_ids
# if false, do all in gtex ids that do not already have a genome_bam in the output folder
regenerate_all=flase

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

bam_array_length=$(wc -l < "${bam_list}")
echo "BAM array length: ${bam_array_length}"

sbatch --output "${output_dir}/logs/%A_%a.log" \
    --error "${output_dir}/logs/%A_%a.log" \
    --array "1-${bam_array_length}%250" \
    --time 12:00:00 \
    --cpus-per-task 1 \
    --partition normal,owners \
    --mem 64G \
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
