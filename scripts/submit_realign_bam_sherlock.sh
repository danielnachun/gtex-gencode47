#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/realign_all_tissues.sh"
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

# # sherlock only lets up to 1k?ish jobs at a time
# if [ "$to_process_count" -gt 1000 ]; then
#   to_process_count=1000
# fi


step=10
max_array_size=1000 

# sequential submission of 1k batches

previous_jobid=""
for batch_start in $(seq 0 $((max_array_size * step)) $to_process_count); do
    batch_end=$((batch_start + (max_array_size - 1) * step))
    if [ "$batch_end" -gt "$to_process_count" ]; then
        batch_end=$to_process_count
    fi

    if [ -z "$previous_jobid" ]; then
        # First batch, no dependency
        jobid=$(sbatch --parsable --output="${output_dir}/logs/%A_%a.log" \
            --error="${output_dir}/logs/%A_%a.log" \
            --array="${batch_start}-${batch_end}:${step}%250" \
            --time=24:00:00 \
            --ntasks-per-node=1 \
            --nodes=10 \
            --partition=normal,owners \
            --mem=100G \
            --job-name=realign_bam_batch \
            ${code_dir}/realign_bam_batch.sh \
                --reference_dir ${reference_dir} \
                --vcf_dir ${vcf_dir} \
                --output_dir ${output_dir} \
                --code_dir ${code_dir} \
                --bam_list ${bam_list} \
                --reference_fasta ${reference_fasta} \
                --rsem_ref_dir ${rsem_ref_dir} \
                --star_index ${star_index})
    else
        # Subsequent batches submit only after first batch is done
        jobid=$(sbatch --parsable --dependency=afterany:${previous_jobid} \
            --output="${output_dir}/logs/%A_%a.log" \
            --error="${output_dir}/logs/%A_%a.log" \
            --array="${batch_start}-${batch_end}:${step}%250" \
            --time=24:00:00 \
            --ntasks-per-node=1 \
            --nodes=10 \
            --partition=normal,owners \
            --mem=100G \
            --job-name=realign_bam_batch \
            ${code_dir}/realign_bam_batch.sh \
                --reference_dir ${reference_dir} \
                --vcf_dir ${vcf_dir} \
                --output_dir ${output_dir} \
                --code_dir ${code_dir} \
                --bam_list ${bam_list} \
                --reference_fasta ${reference_fasta} \
                --rsem_ref_dir ${rsem_ref_dir} \
                --star_index ${star_index})
    fi
    previous_jobid=$jobid
done
