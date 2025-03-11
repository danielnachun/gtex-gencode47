#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# sample args
bam_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/output_francois/sample_bams/
reference_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/references
reference_fasta=Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
chr_sizes=GRCh38.chrsizes
genes_gtf=gencode.v47.genes.gtf
intervals_bed=gencode.v47.GRCh38.insert_size_intervals_geq1000bp.bed
vcf_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs
output_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/output_francois_bam
code_dir=$(realpath $(dirname ${BASH_SOURCE[0]}))

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
