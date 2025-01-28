#!/usr/bin/env bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o pipefail -o errexit


dir_prefix = /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow
sample_id = GTEX-1A3MV-0005-SM-7PC1O
reference_fasta = /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/references/reference.fasta
chr_sizes = ??
genes_gtf = /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/references/gencode.v47.genes.gtf
intervals_bed = ??
het_vcf = ??


# run rnaseq qc
run_rnaseq_qc.sh --duplicate_marked_bam ${dir_prefix}/output/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
     --genes_gtf ${genes_gtf} \
     --genome_fasta ${reference_fasta} \
     --sample_id ${sample_id} \
     --output_dir ${dir_prefix}/output/rnaseq_qc

# run coverage
run_bam_to_coverage.sh --duplicate_marked_bam ${dir_prefix}/output/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
    --chr_sizes ${chr_sizes} \
    --sample_id ${sample_id} \
    --intervals_bed ${intervals_bed} \
    --output_dir ${dir_prefix}/output/coverage

# run gatk
run_gatk.sh --sample_id ${sample_id} \
    --dir_prefix ${dir_prefix} \
    --genome_fasta ${reference_fasta} \
    --het_vcf ${het_vcf} \
    --duplicate_marked_bam ${dir_prefix}/output/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
    --output_dir ${dir_prefix}/output/gatk

# run leafcutter
run_leafcutter.sh --sample_id ${sample_id} \
    --dir_prefix ${dir_prefix} \
    --duplicate_marked_bam ${dir_prefix}/output/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam \
    --output_dir ${dir_prefix}/output/leafcutter

# run phaser
# TODO write wrapper

