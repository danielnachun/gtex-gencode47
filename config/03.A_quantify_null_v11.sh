# sample args
realign_bam_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues/genome_bam
gtex_ids=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/other_references/all_tissues/all_samples.txt
reference_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/realign_references
reference_fasta=Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
chr_sizes=GRCh38.chrsizes
genes_gtf=gencode.v47.matched_nongenic_null_intergenic.gtf
intervals_bed=gencode.v47.GRCh38.insert_size_intervals_geq1000bp.bed
output_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues_null/v11_intergenic
code_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing
bam_file_end=Aligned.sortedByCoord.out.patched.v11md.bam
regenerate_all=false
num_parallel=16
max_array_size=1000
job_time=6:00:00
job_mem=64G