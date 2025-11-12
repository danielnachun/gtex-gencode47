# Single tissue GWAS co-localization config (strong only)
ld_region_list="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/ld_regions.bed"
tissue_id_list="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/qtl_tissue_ids.txt"
gwas_id_list="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/gwas_id_list.txt"
genotype_stem="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.MAF01"
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
num_parallel=1
max_array_size=1000
job_time=24:00:00
job_mem=256G
job_cpus=8
submit_on=sherlock
strong_only=TRUE
p_threshold=5e-8

# Analysis mode flags
multi_tissue=FALSE
run_v39=FALSE
run_individual=FALSE
run_xqtl_only=FALSE
run_separate_gwas=TRUE
run_joint_gwas=FALSE

# GWAS-specific parameters
gwas_params="--gwas_id_list ${gwas_id_list} --gwas_dir ${gwas_dir} --gwas_meta ${gwas_meta} --ld_meta ${ld_meta} --gwas_column_matching ${gwas_column_matching} --strong_only ${strong_only} --p_threshold ${p_threshold}"

# File processing configuration
file_type=ld_regions_string
