# Path to a file containing one tissue_id per line
tissue_id_list="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/qtl_tissue_ids.txt"

# Base directory containing per-tissue coloc outputs. The final COLOC_DIR is ${COLOC_BASE_DIR}/${tissue_id}
coloc_base_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/coloc/single_tissue_qtl_only

code_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing
submit_on=sherlock
regenerate=false
