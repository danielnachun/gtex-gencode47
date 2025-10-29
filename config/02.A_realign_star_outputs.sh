bam_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues_quantifications/genome_bam
gtex_ids=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/other_references/all_tissues/all_samples.txt
reference_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/realign_references
reference_fasta=Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
star_index=STAR_genome_GRCh38_noALT_noHLA_noDecoy_v47_oh75
vcf_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs
output_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues_quantifications_star
code_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing
regenerate_all=false
num_parallel=1
max_array_size=1000
job_time=24:00:00
job_mem=128G
submit_on=sherlock

# File processing configuration
file_type=bam_files
input_dir=${bam_dir}
file_pattern=Aligned.sortedByCoord.out.patched.v11md.bam$
completion_subdir=star