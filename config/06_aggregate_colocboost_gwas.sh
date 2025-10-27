# Path to a file containing one tissue_id per line
tissue_id_list="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/qtl_tissue_ids.txt"

# Base directory containing per-tissue coloc outputs. The final COLOC_DIR is ${COLOC_BASE_DIR}/${tissue_id}
coloc_base_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/coloc/single_tissue_gwas_strong_only

code_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing
num_parallel=1
max_array_size=1000
job_time=12:00:00
job_mem=128G
submit_on=scg
regenerate=false
skip_robust=true
