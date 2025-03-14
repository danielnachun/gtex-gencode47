#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# sample args
bam_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/raw/GTEx_Analysis_2022-06-06_v10_RNAseq_BAM_files
gtex_ids=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/caudate_analysis/data/gtex_ids.txt
reference_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/references
reference_fasta=Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
rsem_ref_dir=rsem_reference_GRCh38_gencode47
star_index=STAR_genome_GRCh38_noALT_noHLA_noDecoy_v47_oh75
vcf_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs
output_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/caudate_analysis/output/
code_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/code


# make a bam list
mkdir -p ${output_dir}
mkdir -p ${output_dir}/logs

bam_list="$output_dir/bam_list"
awk -v dir="$bam_dir" '{print dir "/" $0 ".Aligned.sortedByCoord.out.patched.md.bam"}' "$gtex_ids" > "$bam_list"
bam_array_length=$(wc -l < ${bam_list})

sbatch --output "${output_dir}/logs/%A_%a.log" \
    --error "${output_dir}/logs/%A_%a.log" \
    --array "1-${bam_array_length}%250" \
    --time 6:00:00 \
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
