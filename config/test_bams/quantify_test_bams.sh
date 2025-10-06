# sample args for quantifying 3 test bams
realign_bam_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/test_bams/output_kate/genome_bam
gtex_ids=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/other_references/test_bams/test_bams_samples.txt
reference_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/quantify_references
reference_fasta=Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
chr_sizes=GRCh38.chrsizes
genes_gtf=gencode.v47.genes.gtf
intervals_bed=gencode.v47.GRCh38.insert_size_intervals_geq1000bp.bed
ipa_annotation=IPAFinder_anno_hg38.txt
editing_bed=test_edsites.bed
vcf_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs
output_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/test_bams/output_kate
code_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing
regenerate_all=true
batch_size=1
max_array_size=1000
submit_on=sherlock