#!/usr/bin/env bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o pipefail -o errexit


chromosome=chr21
bam_files="/home/klawren/oak/gtex/output/all_tissues/genome_bam/GTEX-1A8G6-0006-SM-7PC1X.Aligned.sortedByCoord.out.patched.v11md.bam \
    /home/klawren/oak/gtex/output/all_tissues/genome_bam/GTEX-1A8G6-0008-SM-EWRMK.Aligned.sortedByCoord.out.patched.v11md.bam \
    /home/klawren/oak/gtex/output/all_tissues/genome_bam/GTEX-1A8G6-0011-R1a-SM-7P8PC.Aligned.sortedByCoord.out.patched.v11md.bam"
participant_id="GTEX-1A8G6-0008"
genotype_vcf=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/realign_references/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.vcf.gz
gene_model_bed=/home/klawren/oak/gtex/data/quantify_references/gencode.v47.genes.phaser.bed
dir_prefix=/home/klawren/oak/gtex/output/phaser_test
local_reference_dir=/home/klawren/oak/gtex/output/phaser_test
haplotype_blacklist_bed=/home/klawren/oak/gtex/data/realign_references/hg38_haplo_count_blacklist.bed
code_dir=/home/klawren/oak/gtex/scripts/processing

# run phaser
bash ${code_dir}/run_phaser.sh \
    --bam_files ${bam_files} \
    --participant_id ${participant_id} \
    --chromosome ${chromosome} \
    --genotype_vcf ${genotype_vcf} \
    --gene_model_bed ${gene_model_bed} \
    --output_dir ${dir_prefix}/output/phaser \
    --working_dir ${local_reference_dir}/${participant_id} \
    --num_threads 8
#    --haplotype_blacklist_bed ${haplotype_blacklist_bed} \
#    --phaser_blacklist_bed ${phaser_blacklist_bed} 