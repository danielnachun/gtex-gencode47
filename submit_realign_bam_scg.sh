#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# sample args
bam_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/test_samples
reference_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/references
reference_fasta=Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
rsem_ref_dir=rsem_reference_GRCh38_gencode47
star_index=STAR_genome_GRCh38_noALT_noHLA_noDecoy_v47_oh75
vcf_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs
output_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/output_kate/
code_dir=$(realpath $(dirname ${BASH_SOURCE[0]}))

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
