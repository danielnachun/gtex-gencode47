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


bam_list="${bam_dir}/../bam_list"
find -L "${bam_dir}" -type f -name "*.bam" | sort -u | head -n 10 > ${bam_list}
bam_array_length=$(wc -l < ${bam_list})
mkdir -p ${output_dir}/logs

sbatch --output "${output_dir}/logs/%A_%a.log" \
    --error "${output_dir}/logs/%A_%a.log" \
    --array "1-${bam_array_length}" \
    --time 12:00:00 \
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
