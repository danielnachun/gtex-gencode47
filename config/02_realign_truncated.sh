bam_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues_quantifications/genome_bam
gtex_ids=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/other_references/truncated_samples.txt
reference_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/realign_references
reference_fasta=Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
rsem_ref_dir=rsem_reference_GRCh38_gencode47
star_index=STAR_genome_GRCh38_noALT_noHLA_noDecoy_v47_oh75
vcf_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs
output_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/truncated
code_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing
bam_file_end=Aligned.sortedByCoord.out.patched.v11md.bam
regenerate_all=true
num_parallel=1
max_array_size=1000
job_time=12:00:00
job_mem=128G
submit_on=scg

# File processing configuration
file_pattern=Aligned.sortedByCoord.out.patched.v11md.bam$