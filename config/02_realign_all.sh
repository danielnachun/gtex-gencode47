bam_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues_quantifications/genome_bam
gtex_ids=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/other_references/all_tissues/all_samples.txt
reference_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/realign_references
reference_fasta=Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
rsem_ref_dir=rsem_reference_GRCh38_gencode47
star_index=STAR_genome_GRCh38_noALT_noHLA_noDecoy_v47_oh75
vcf_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs
output_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues_quantifications
code_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing
regenerate_all=false
# step limited by RAM availible on most sherlock nodes (~64G/bam needed)
batch_size=1
max_array_size=1000
submit_on=scg