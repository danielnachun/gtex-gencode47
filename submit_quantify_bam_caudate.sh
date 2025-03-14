#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# sample args
realign_bam_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/caudate_analysis/output/genome_bam
gtex_ids=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/caudate_analysis/data/gtex_ids.txt
reference_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/references
reference_fasta=Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
chr_sizes=GRCh38.chrsizes
genes_gtf=gencode.v47.genes.gtf
intervals_bed=gencode.v47.GRCh38.insert_size_intervals_geq1000bp.bed
vcf_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs
output_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/caudate_analysis/output/
code_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/code

# make a bam list
mkdir -p ${output_dir}
mkdir -p ${output_dir}/logs

bam_list="$output_dir/realigned_bam_list"
awk -v dir="$realign_bam_dir" '{print dir "/" $0 ".Aligned.sortedByCoord.out.patched.v11md.bam"}' "$gtex_ids" > "$bam_list"
bam_array_length=$(wc -l < ${bam_list})

sbatch --output "${output_dir}/logs/%A_%a.log" \
    --error "${output_dir}/logs/%A_%a.log" \
    --array "1-${bam_array_length}%250" \
    --time 1:00:00 \
    --partition normal,owners \
    --cpus-per-task 1 \
    --mem 64G \
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
