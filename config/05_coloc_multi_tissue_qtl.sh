# Multi tissue GWAS co-localization config (strong only)
ld_region_list="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/pyrho_EUR_LD_blocks.bed"
tissue_id_list="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/qtl_tissue_ids.txt"
genotype_stem="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.MAF01"
covariate_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/qtl/GTEx_Analysis_v11_eQTL_covariates"
expression_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/qtl/GTEx_Analysis_v11_eQTL_expression_matrices"

# gwas_id_list="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/gwas_ids.txt"
# gwas_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/barbeira_gtex_imputed/imputed_gwas_hg38_1.1"
# gwas_meta="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/barbeira_gtex_imputed/barberia_full_gwas_metadata.txt"
# ld_meta="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/ld_blocks_gtex_eur/ld_metadata.tsv"
# gwas_column_matching="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/barbeira_gtex_imputed/barberia_gtex_column_matching.yml"


all_v39_genes_path="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/v39_genes.txt"
region_padding=1000000
output_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/eqtl_coloc"
code_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing"
regenerate_all=FALSE
num_parallel=1
max_array_size=3
job_time=10:00
job_mem=256G
job_cpus=8
submit_on=sherlock


# Analysis mode flags
multi_tissue=TRUE
run_v39=FALSE
run_individual=FALSE
run_xqtl_only=TRUE
run_separate_gwas=FALSE
run_joint_gwas=FALSE


# File processing configuration
file_type=ld_regions_string
