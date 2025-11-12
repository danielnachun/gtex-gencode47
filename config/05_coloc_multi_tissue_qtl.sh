# Multi tissue QTL co-localization config
ld_region_list="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/ld_regions.bed"
tissue_id_list="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/qtl_tissue_ids.txt"
genotype_stem="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/genotypes/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.MAF01"
covariate_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/covariates"
expression_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/expression"
all_v39_genes_path="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/v39_genes.txt"
region_padding=1000000
association_padding=0
output_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/coloc"
code_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing"
regenerate_all=FALSE
num_parallel=1
max_array_size=1000
job_time=24:00:00
job_mem=256G
job_cpus=4
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
