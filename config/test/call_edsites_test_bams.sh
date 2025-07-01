# sample args for quantifying 3 test bams
realign_bam_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues/genome_bam
gtex_ids=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/other_references/test_bams/test_bams_samples.txt
reference_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/edsite_references
reference_fasta=Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
dbsnp=All_20180418.vcf.gz
indels_mills=Mills_and_1000G_gold_standard.indels.hg38.sites.vcf
indels_decoy=Homo_sapiens_assembly38_1000genomes_decoy.known_indels.vcf
gene_intervals_bed=gencode.v47.merged.bed
output_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/test_bams/output_kate
code_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing
regenerate_all=false
step_size=1
max_array_size=1000
submit_on=sherlock