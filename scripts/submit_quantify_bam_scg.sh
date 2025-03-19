#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/quantify_test_bams.sh"
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    echo "Error: Config file $CONFIG_FILE not found!"
    exit 1
fi

bam_list="${bam_dir}/../bam_list"
find -L "${bam_dir}" -type f -name "*.bam" | sort -u > ${bam_list}
bam_array_length=$(wc -l < ${bam_list})
mkdir -p ${output_dir}/logs

sbatch --output "${output_dir}/logs/%A_%a.log" \
    --error "${output_dir}/logs/%A_%a.log" \
    --array "1-${bam_array_length}" \
    --time 6:00:00 \
    --account smontgom \
    --partition batch \
    --cpus-per-task 1 \
    --mem 64G \
    --constraint="nvme" \
    --job-name quantiy_bam \
    ${code_dir}/quantify_bam.sh \
        --reference_dir ${reference_dir} \
        --vcf_dir ${vcf_dir} \
        --output_dir ${output_dir} \
        --code_dir ${code_dir} \
        --reference_fasta ${reference_fasta} \
        --bam_list=${bam_list} \
        --genes_gtf ${genes_gtf} \
        --chr_sizes ${chr_sizes} \
        --intervals_bed ${intervals_bed}
