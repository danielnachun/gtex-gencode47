#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/quantify_caudate_v10.sh"
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    echo "Error: Config file $CONFIG_FILE not found!"
    exit 1
fi

# if true do all in gtex_ids
# if false, do all in gtex ids that do not already have a regtools output in the leafcutter folder (that is the last step)
regenerate_all=false

# make a bam list
mkdir -p ${output_dir}
mkdir -p ${output_dir}/logs
bam_list="$output_dir/realigned_bam_list"
> "$bam_list" # clear the file


if [ "$regenerate_all" = true ]; then
awk -v dir="$realign_bam_dir" -v fileEnd="$bam_file_end" '{print dir "/" $0 file_end}' "$gtex_ids" > "$bam_list"
else
    # only add a gtex id if the null rnaseqc doesn't exist
    > "$bam_list"  
    while read -r gtex_id; do
        bam_file="${output_dir}/rnaseq_qc/${gtex_id}.gene_tpm.gct.gz"
        if [ ! -f "$bam_file" ]; then
            echo "$realign_bam_dir/$gtex_id.$bam_file_end" >> "$bam_list"
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
    --array "1-${to_process_count}%250" \
    --time 6:00:00 \
    --partition normal,owners \
    --cpus-per-task 1 \
    --mem 64G \
    --job-name quantiy_bam \
    ${code_dir}/quantify_rnaseqc.sh \
        --reference_dir ${reference_dir} \
        --output_dir ${output_dir} \
        --code_dir ${code_dir} \
        --reference_fasta ${reference_fasta} \
        --bam_list=${bam_list} \
        --genes_gtf ${genes_gtf} \
        --chr_sizes ${chr_sizes} \
        --intervals_bed ${intervals_bed} \
        --bam_file_end ${bam_file_end}
