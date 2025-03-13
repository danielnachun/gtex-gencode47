#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# sample args
# old_dir=/home/klawren/oak/gtex/test_workflow/output_kate/tmp_old
old_dir=/home/klawren/oak/gtex/test_workflow/output_dan
new_dir=/home/klawren/oak/gtex/test_workflow/output_kate


sample_id=GTEX-YFCO-0626-SM-HM8UJ

# # for rnaseq_qc
diff <(zcat ${old_dir}/rnaseq_qc/${sample_id}.exon_reads.gct.gz) <(zcat ${new_dir}/rnaseq_qc/${sample_id}.exon_reads.gct.gz)
diff ${old_dir}/rnaseq_qc/${sample_id}.fragmentSizes.txt ${new_dir}/rnaseq_qc/${sample_id}.fragmentSizes.txt
diff <(zcat ${old_dir}/rnaseq_qc/${sample_id}.gene_reads.gct.gz) <(zcat ${new_dir}/rnaseq_qc/${sample_id}.gene_reads.gct.gz)
diff <(zcat ${old_dir}/rnaseq_qc/${sample_id}.gene_tpm.gct.gz) <(zcat ${new_dir}/rnaseq_qc/${sample_id}.gene_tpm.gct.gz)
diff ${old_dir}/rnaseq_qc/${sample_id}.metrics.tsv ${new_dir}/rnaseq_qc/${sample_id}.metrics.tsv

for rsem
diff ${old_dir}/rsem/${sample_id}.rsem.genes.results ${new_dir}/rsem/${sample_id}.rsem.genes.results
diff ${old_dir}/rsem/${sample_id}.rsem.isoforms.results ${new_dir}/rsem/${sample_id}.rsem.isoforms.results

# for leafcutter
diff <(zcat ${old_dir}/leafcutter/${sample_id}.regtools_junc.txt.gz) <(zcat ${new_dir}/leafcutter/${sample_id}.regtools_junc.txt.gz)

# # for star
# diff <(zcat ${old_dir}/star/${sample_id}.SJ.out.tab.gz) <(zcat ${new_dir}/star/${sample_id}.SJ.out.tab.gz)
# diff <(zcat ${old_dir}/star/${sample_id}.ReadsPerGene.out.tab.gz) <(zcat ${new_dir}/star/${sample_id}.ReadsPerGene.out.tab.gz)
