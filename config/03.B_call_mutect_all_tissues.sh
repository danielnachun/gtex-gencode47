realign_bam_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues_quantifications/genome_bam
gtex_ids=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/other_references/all_tissues/all_samples.txt
reference_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/mutect_references
reference_fasta=Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
dbsnp=All_20180418.vcf.gz
indels_mills=Mills_and_1000G_gold_standard.indels.hg38.sites.vcf
indels_decoy=Homo_sapiens_assembly38_1000genomes_decoy.known_indels.vcf
gene_intervals_bed=gencode.v47.merged.bed
exac_reference=small_exac_common_3.hg38.combined.vcf.gz
output_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues_quantifications/
code_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing
vcf_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs
full_vcf_file=GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.vcf.gz
regenerate_all=false
num_parallel=2
max_array_size=1000
job_time=24:00:00
job_mem=128G
submit_on=sherlock

# File processing configuration
file_type=bam_files
input_dir=${realign_bam_dir}
file_pattern=Aligned.sortedByCoord.out.patched.v11md.bam$
completion_subdir=mutect