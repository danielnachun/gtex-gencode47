# Multi tissue GWAS co-localization config
ld_region_list="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/ld_regions.bed"
tissue_id_list="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/qtl_tissue_ids.txt"
gwas_id_list="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/gwas_id_list.txt"
genotype_stem="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/genotypes/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.MAF01"
covariate_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/covariates"
expression_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/expression"
gwas_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/gwas"
gwas_meta="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/gwas_meta.txt"
ld_meta="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/ld_meta.txt"
gwas_column_matching="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/gwas_column_matching.txt"
all_v39_genes_path="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/v39_genes.txt"
region_padding=1000000
association_padding=0
output_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/coloc"
code_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing"
regenerate_all=FALSE

# submission params
num_parallel=1
max_array_size=1000
job_time=24:00:00
job_mem=512G
job_cpus=1
submit_on=scg

strong_only=TRUE
p_threshold=1e-6

# coloc mode flags
multi_tissue=TRUE
run_v39=TRUE
run_individual=TRUE
run_xqtl_only=FALSE
run_separate_gwas=FALSE
run_joint_gwas=TRUE

# GWAS-specific parameters
gwas_params="--gwas_id_list ${gwas_id_list} --gwas_dir ${gwas_dir} --gwas_meta ${gwas_meta} --ld_meta ${ld_meta} --gwas_column_matching ${gwas_column_matching}"
