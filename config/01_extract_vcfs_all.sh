# sample args
full_vcf=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/realign_references/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.vcf.gz
output_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs
participant_id_list=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/other_references/all_tissues/all_participants.txt
code_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing
num_parallel=1
max_array_size=1000
job_time=4:00:00
job_mem=64G
submit_on=sherlock
regenerate_all=false

# File processing configuration
file_type=participants