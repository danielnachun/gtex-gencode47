
# sample args
dir_prefix = /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow
sample_id = GTEX-1A3MV-0005-SM-7PC1O
reference_fasta = /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
chr_sizes = /oak/stanford/groups/smontgom/dnachun/data/gtex/data/references/GRCh38.chrsizes
genes_gtf = /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/references/gencode.v47.annotation.gtf
# fix interval_bed download
# intervals_bed = ??
# het_vcf pulled from folder of participant vcfs


# het_vcf = /oak/stanford/groups/smontgom/dnachun/data/gtex/data/references/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.vcf.gz

./quantify_bam.sh \
    ${dir_prefix} \
    ${sample_id} \
    ${reference_fasta} \
    ${chr_sizes} \
    ${genes_gtf} \
    ${intervals_bed} \
    ${het_vcf}