#!/usr/bin/env bash

#set -o xtrace -o nounset -o errexit

# sample args
old_dir=/home/klawren/oak/gtex/test_workflow/output_fixed_rnaseqqc
new_dir=/home/klawren/oak/gtex/test_workflow/output_francois
#new_dir=/home/klawren/oak/gtex/test_workflow/output_francois








#sample_id=GTEX-1C4CL-2126-SM-7IGQC
sample_id=GTEX-YFCO-0626-SM-HM8UJ

#diff --speed-large-files <(zcat /home/klawren/oak/gtex/test_workflow/output_francois/fa-exchange2/GTEX-YFCO.snps.vcf.gz) <(zcat /home/klawren/oak/gtex/data/processed/vcfs/GTEX-YFCO.snps.vcf.gz)


# # check for gatk
# diff <(zcat ${old_dir}/gatk/${sample_id}.readcounts.txt.gz) <(zcat ${new_dir}/gatk/${sample_id}.readcounts.txt.gz)
# diff <(zcat ${old_dir}/gatk/${sample_id}.readcounts.chrX.txt.gz) <(zcat ${new_dir}/gatk/${sample_id}.readcounts.chrX.txt.gz)

# # check for leafcutter
# # ignore 4th column where junctions are named
# #diff <(zcat ${old_dir}/leafcutter/${sample_id}.regtools_junc.txt.gz) <(zcat ${new_dir}/leafcutter/${sample_id}.regtools_junc.txt.gz)
# diff <(zcat ${old_dir}/leafcutter/${sample_id}.regtools_junc.txt.gz | awk '{$4=""; print $0}' | sed 's/  *$//') <(zcat ${new_dir}/leafcutter/${sample_id}.regtools_junc.txt.gz | awk '{$4=""; print $0}' | sed 's/  *$//')


# # for rnaseq_qc
diff <(zcat ${old_dir}/rnaseq_qc/${sample_id}.exon_reads.gct.gz) <(zcat ${new_dir}/rnaseq_qc/${sample_id}.exon_reads.gct.gz)
diff ${old_dir}/rnaseq_qc/${sample_id}.fragmentSizes.txt ${new_dir}/rnaseq_qc/${sample_id}.fragmentSizes.txt
diff ${old_dir}/rnaseq_qc/${sample_id}.gc_content.tsv ${new_dir}/rnaseq_qc/${sample_id}.gc_contenct.tsv
diff <(zcat ${old_dir}/rnaseq_qc/${sample_id}.gene_reads.gct.gz) <(zcat ${new_dir}/rnaseq_qc/${sample_id}.gene_reads.gct.gz)
diff <(zcat ${old_dir}/rnaseq_qc/${sample_id}.gene_tpm.gct.gz) <(zcat ${new_dir}/rnaseq_qc/${sample_id}.gene_tpm.gct.gz)
diff ${old_dir}/rnaseq_qc/${sample_id}.metrics.tsv ${new_dir}/rnaseq_qc/${sample_id}.metrics.tsv


# # Generate MD5 checksums for the decompressed files
# old_md5=$(md5sum <(zcat "${old_dir}/rnaseq_qc/${sample_id}.gene_tpm.gct.gz") | awk '{print $1}')
# new_md5=$(md5sum <(zcat "${new_dir}/rnaseq_qc/${sample_id}.gene_tpm.gct.gz") | awk '{print $1}')

# # Print the checksums
# echo "Old MD5: $old_md5"
# echo "New MD5: $new_md5"

# # Compare the checksums
# if [[ "$old_md5" == "$new_md5" ]]; then
#     echo "The files are identical."
# else
#     echo "The files are different."
# fi



# # # for rsem
# diff ${old_dir}/rsem/${sample_id}.rsem.genes.results ${new_dir}/rsem/${sample_id}.rsem.genes.results
# diff ${old_dir}/rsem/${sample_id}.rsem.isoforms.results ${new_dir}/rsem/${sample_id}.rsem.isoforms.results




# #  # for star
# diff <(zcat ${old_dir}/star/${sample_id}.SJ.out.tab.gz) <(zcat ${new_dir}/star/${sample_id}.SJ.out.tab.gz)
# diff <(zcat ${old_dir}/star/${sample_id}.ReadsPerGene.out.tab.gz) <(zcat ${new_dir}/star/${sample_id}.ReadsPerGene.out.tab.gz)



# compare bams

# picard CompareSAMs ${old_dir}/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam ${new_dir}/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam LENIENT_HEADER=true
# n=1000
# diff <(samtools view ${old_dir}/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam | head -${n}) <(samtools view ${new_dir}/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam | head  -${n})
#diff <(samtools view ${old_dir}/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.sorted.bam | head -${n}) <(samtools view ${new_dir}/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.sorted.bam | head  -${n})

# # compare intermediate outputs
# n=1000
# diff <(samtools view ${old_dir}/star/${sample_id}.Aligned.sortedByCoord.out.bam | head -${n}) <(samtools view ${new_dir}/star/${sample_id}.Aligned.sortedByCoord.out.bam | head  -${n})


# mkdir -p "${old_dir}/sorted_bams/"
# mkdir -p "${new_dir}/sorted_bams/"

# # sort by N then compare with diff
# sample_ids=("GTEX-1C4CL-2126-SM-7IGQC" "GTEX-YFCO-0626-SM-HM8UJ" "GTEX-1NV5F-3226-SM-EXUSL")
# for sample_id in "${sample_ids[@]}"; do
#     echo "Comparing RNA-seq QC for sample ID: $sample_id"
#     echo "sorting with samtools"
#     # samtools sort ${old_dir}/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam -N -o ${old_dir}/sorted_bams/${sample_id}.Aligned.sortedByN.out.patched.v11md.bam
#     # samtools sort ${new_dir}/genome_bam/${sample_id}.Aligned.sortedByCoord.out.patched.v11md.bam -N -o ${new_dir}/sorted_bams/${sample_id}.Aligned.sortedByN.out.patched.v11md.bam

#     echo "comparing with diff"
#     diff --speed-large-files <(samtools view ${old_dir}/sorted_bams/${sample_id}.Aligned.sortedByN.out.patched.v11md.bam) <(samtools view ${new_dir}/sorted_bams/${sample_id}.Aligned.sortedByN.out.patched.v11md.bam)    
#     # Check the exit status of the diff command
#     if [ $? -eq 0 ]; then
#         echo "No differences found for sample ID: $sample_id"
#     else
#         echo "Differences found for sample ID: $sample_id"
#     fi
# done