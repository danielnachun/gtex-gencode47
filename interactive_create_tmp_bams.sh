#!/usr/bin/env bash

# with hardcoded paths just to generate some intermediate files
bam_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/raw/GTEx_Analysis_2022-06-06_v10_RNAseq_BAM_files/
reference_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/references
reference_fasta=Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
rsem_ref_dir=rsem_reference_GRCh38_gencode47
star_index=STAR_genome_GRCh38_noALT_noHLA_noDecoy_v47_oh75
vcf_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs
output_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/output_kate/tmp_new_star_2
code_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/code

sample_id=GTEX-YFCO-0626-SM-HM8UJ
participant_id=$(echo ${sample_id} | cut -d '-' -f1,2)
vcf_file=${participant_id}.snps.vcf.gz


source <(pixi shell-hook --environment realignbam --manifest-path /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/code/pixi.toml)

# process bams to fastqs
bash ${code_dir}/run_bam_to_fastq.sh --bam_file ${bam_dir}/${sample_id}.Aligned.sortedByCoord.out.patched.md.bam \
    --sample_id ${sample_id} \
    --reference_fasta ${reference_dir}/${reference_fasta} \
    --tmp_dir ${output_dir}/fastq

# align with star
bash ${code_dir}/run_fastq_to_star.sh \
    --star_index ${reference_dir}/${star_index} \
    --fastq_1 ${output_dir}/fastq/${sample_id}_1.fastq.gz \
    --fastq_2 ${output_dir}/fastq/${sample_id}_2.fastq.gz \
    --sample_id ${sample_id} \
    --vcf_file ${vcf_dir}/${vcf_file} \
    --tmp_dir ${output_dir}/star

# sync bams
bash ${code_dir}/run_bam_sync.sh \
    --initial_bam_file ${bam_dir}/${sample_id}.Aligned.sortedByCoord.out.patched.md.bam \
    --star_aligned_bam ${output_dir}/star/${sample_id}.Aligned.sortedByCoord.out.bam \
    --sample_id ${sample_id} \
    --tmp_dir ${output_dir}/bamsync \
    --output_dir ${output_dir}/flagstat

# mark duplicates, get the genome bam that we save
bash ${code_dir}/run_mark_duplicates.sh \
    --genome_bam_file ${output_dir}/bamsync/${sample_id}.Aligned.sortedByCoord.out.patched.bam \
    --output_prefix ${sample_id}.Aligned.sortedByCoord.out.patched.v11md \
    --output_dir ${output_dir}/genome_bam

# run rsem, save isoform quantification
bash ./run_rsem.sh \
    --rsem_ref_dir ${reference_dir}/${rsem_ref_dir} \
    --transcriptome_bam ${output_dir}/star/${sample_id}.Aligned.toTranscriptome.out.bam \
    --sample_id ${sample_id} \
    --output_dir ${output_dir}/rsem