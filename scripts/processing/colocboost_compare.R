#!/usr/bin/env Rscript

#' Unified colocboost script for QTL and GWAS analysis modes
#' Handles QTL-only, GWAS separate, and GWAS joint analyses based on boolean flags

library(argparser)

# Parse arguments
parser <- argparser::arg_parser("Script to run colocboost with flexible analysis modes") |>
    argparser::add_argument("--tissue_id", help = "Tissue ID") |>
    argparser::add_argument("--gene_region", help = "Region of genes to analyze") |>
    argparser::add_argument("--variant_region", help = "Region of variants to include (association window region)") |>
    argparser::add_argument("--ld_region", help = "LD block (for naming)") |>
    argparser::add_argument("--genotype_stem", help = "Path to genotype BED file") |>
    argparser::add_argument("--phenotype_list", help = "Path to list of phenotypes to analyze, one path per line") |>
    argparser::add_argument("--covariate_path", help = "Path to covariates file", default = NULL) |>
    argparser::add_argument("--covariate_list", help = "Optional path to list of covariate files, one per phenotype", default = NULL) |>
    # GWAS-specific parameters (optional)
    argparser::add_argument("--gwas_id_list", help = "Path to list of GWAS IDs to analyze, one per line", default = NULL) |>
    argparser::add_argument("--gwas_phenotype_list", help = "Path to list of GWAS phenotypes to analyze, one path per line", default = NULL) |>
    argparser::add_argument("--gwas_meta", help = "Path to GWAS metadata file", default = NULL) |>
    argparser::add_argument("--ld_meta", help = "Path to LD metadata file", default = NULL) |>
    argparser::add_argument("--gwas_column_matching", help = "Path to GWAS column matching file", default = NULL) |>
    # Analysis mode flags
    argparser::add_argument("--multi_tissue", help = "Multi-tissue analysis flag", default = FALSE) |>
    argparser::add_argument("--run_v39", help = "Run v39 genes analysis", default = FALSE) |>
    argparser::add_argument("--run_individual", help = "Run individual gene analysis", default = FALSE) |>
    argparser::add_argument("--run_xqtl_only", help = "Run QTL-only analysis", default = FALSE) |>
    argparser::add_argument("--run_separate_gwas", help = "Run separate GWAS analysis", default = FALSE) |>
    argparser::add_argument("--run_joint_gwas", help = "Run joint GWAS analysis", default = FALSE) |>
    # Legacy parameters (for backward compatibility)
    argparser::add_argument("--run_single_gene", help = "Run single gene analysis", default = FALSE) |>
    argparser::add_argument("--run_v39_genes", help = "Run v39 genes analysis", default = FALSE) |>
    # Output and quality control
    argparser::add_argument("--output_dir", help = "Output directory for results") |>
    argparser::add_argument("--v39_gene_id_path", help = "Path to list of v39 gene IDs, one per line", default = NULL) |>
    argparser::add_argument("--maf_cutoff", help = "Minor allele frequency cutoff", default = 0.01) |>
    argparser::add_argument("--mac_cutoff", help = "Minor allele count cutoff", default = 0) |>
    argparser::add_argument("--xvar_cutoff", help = "Genotype variance cutoff", default = 0) |>
    argparser::add_argument("--imiss_cutoff", help = "Missingness cutoff", default = 0)

parsed_args <- argparser::parse_args(parser)

# Set DEBUG based on analysis mode
DEBUG <- if (parsed_args$run_xqtl_only) TRUE else FALSE

# Extract core parameters
tissue_id         <- parsed_args$tissue_id
gene_region       <- parsed_args$gene_region
genotype_stem     <- parsed_args$genotype_stem
phenotype_list    <- readLines(parsed_args$phenotype_list)

# Check if phenotype_list is empty and handle gracefully
if (length(phenotype_list) == 0 || all(phenotype_list == "")) {
    cat("No phenotypes found for tissue", tissue_id, "- skipping analysis\n")
    quit(status = 0)
}

covariate_path    <- parsed_args$covariate_path
covariate_list_path <- parsed_args$covariate_list

# Normalize optional args: treat "none" or empty string as NULL
normalize_opt <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.na(x)) return(NULL)
    if (!nzchar(x)) return(NULL)
    if (tolower(x) == "none") return(NULL)
    x
}

covariate_path <- normalize_opt(covariate_path)
covariate_list_path <- normalize_opt(covariate_list_path)

# Extract analysis mode flags
multi_tissue <- as.logical(parsed_args$multi_tissue)
run_v39 <- as.logical(parsed_args$run_v39)
run_individual <- as.logical(parsed_args$run_individual)
run_xqtl_only <- as.logical(parsed_args$run_xqtl_only)
run_separate_gwas <- as.logical(parsed_args$run_separate_gwas)
run_joint_gwas <- as.logical(parsed_args$run_joint_gwas)

# Legacy parameter mapping
run_single_gene <- as.logical(parsed_args$run_single_gene) || run_individual
run_v39_genes <- as.logical(parsed_args$run_v39_genes) || run_v39

# Determine analysis mode
analysis_mode <- "basic"
if (run_xqtl_only) {
    analysis_mode <- "qtl_only"
} else if (run_separate_gwas) {
    analysis_mode <- "gwas_separate"
} else if (run_joint_gwas) {
    analysis_mode <- "gwas_joint"
}

cat("Analysis mode:", analysis_mode, "\n")
cat("Multi-tissue:", multi_tissue, "\n")
cat("Run v39:", run_v39, "\n")
cat("Run individual:", run_individual, "\n")

# Extract GWAS parameters if running GWAS analysis
if (run_separate_gwas || run_joint_gwas) {
    gwas_id_list <- readLines(parsed_args$gwas_id_list)
    gwas_phenotype_list <- readLines(parsed_args$gwas_phenotype_list)
    gwas_meta           <- parsed_args$gwas_meta
    ld_meta             <- parsed_args$ld_meta
    gwas_column_matching <- parsed_args$gwas_column_matching
    
    cat("GWAS analysis enabled\n")
    cat("  GWAS IDs:", length(gwas_id_list), "\n")
    cat("  GWAS phenotypes:", length(gwas_phenotype_list), "\n")
}

# Extract other parameters
output_dir <- parsed_args$output_dir
v39_gene_id_path <- normalize_opt(parsed_args$v39_gene_id_path)
maf_cutoff <- as.numeric(parsed_args$maf_cutoff)
mac_cutoff <- as.numeric(parsed_args$mac_cutoff)
xvar_cutoff <- as.numeric(parsed_args$xvar_cutoff)
imiss_cutoff <- as.numeric(parsed_args$imiss_cutoff)

cat("Arguments:\n")
cat("  tissue_id:", tissue_id, "\n")
cat("  gene_region:", gene_region, "\n")
cat("  variant_region:", parsed_args$variant_region, "\n")
cat("  ld_region:", parsed_args$ld_region, "\n")
cat("  genotype_stem:", genotype_stem, "\n")
cat("  phenotype_list:", length(phenotype_list), "files\n")
cat("  covariate_path:", if (is.null(covariate_path)) "NULL" else covariate_path, "\n")
cat("  covariate_list:", if (is.null(covariate_list_path)) "NULL" else covariate_list_path, "\n")
cat("  output_dir:", output_dir, "\n")
cat("  maf_cutoff:", maf_cutoff, "\n")
cat("  mac_cutoff:", mac_cutoff, "\n")
cat("  xvar_cutoff:", xvar_cutoff, "\n")
cat("  imiss_cutoff:", imiss_cutoff, "\n")
cat("  run_single_gene:", run_single_gene, "\n")
cat("  run_v39_genes:", run_v39_genes, "\n")

cat("\nLoading required libraries...\n")
library(pecotmr)
library(colocboost)
cat("Libraries loaded successfully\n")

# Load v39 gene IDs if specified
v39_gene_ids <- NULL
if (!is.null(v39_gene_id_path) && run_v39_genes) {
    v39_gene_ids <- readLines(v39_gene_id_path)
    cat("Loaded", length(v39_gene_ids), "v39 gene IDs\n")
}

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load genotype data
cat("\nLoading genotype data...\n")
genotype_data <- load_genotype_data(
    genotype_stem = genotype_stem,
    variant_region = parsed_args$variant_region,
    maf_cutoff = maf_cutoff,
    mac_cutoff = mac_cutoff,
    xvar_cutoff = xvar_cutoff,
    imiss_cutoff = imiss_cutoff
)
cat("Loaded genotype data:", nrow(genotype_data$genotypes), "variants,", ncol(genotype_data$genotypes), "samples\n")

# Load phenotype data
cat("\nLoading phenotype data...\n")
phenotype_data <- load_phenotype_data(
    phenotype_list = phenotype_list,
    covariate_path = covariate_path,
    covariate_list_path = covariate_list_path
)
cat("Loaded phenotype data:", length(phenotype_data$phenotypes), "phenotypes\n")

# Load GWAS data if running GWAS analysis
gwas_data <- NULL
if (run_separate_gwas || run_joint_gwas) {
    cat("\nLoading GWAS data...\n")
    gwas_data <- load_gwas_data(
        gwas_id_list = gwas_id_list,
        gwas_phenotype_list = gwas_phenotype_list,
        gwas_meta = gwas_meta,
        ld_meta = ld_meta,
        gwas_column_matching = gwas_column_matching
    )
    cat("Loaded GWAS data:", length(gwas_data$gwas_ids), "GWAS studies\n")
}

# Run analysis based on mode
if (analysis_mode == "qtl_only") {
    cat("\nRunning QTL-only analysis...\n")
    results <- run_qtl_analysis(
        genotype_data = genotype_data,
        phenotype_data = phenotype_data,
        gene_region = gene_region,
        run_single_gene = run_single_gene,
        run_v39_genes = run_v39_genes,
        v39_gene_ids = v39_gene_ids,
        output_dir = output_dir,
        tissue_id = tissue_id,
        ld_region = parsed_args$ld_region
    )
} else if (analysis_mode == "gwas_separate") {
    cat("\nRunning separate GWAS analysis...\n")
    results <- run_gwas_separate_analysis(
        genotype_data = genotype_data,
        phenotype_data = phenotype_data,
        gwas_data = gwas_data,
        gene_region = gene_region,
        run_single_gene = run_single_gene,
        run_v39_genes = run_v39_genes,
        v39_gene_ids = v39_gene_ids,
        output_dir = output_dir,
        tissue_id = tissue_id,
        ld_region = parsed_args$ld_region
    )
} else if (analysis_mode == "gwas_joint") {
    cat("\nRunning joint GWAS analysis...\n")
    results <- run_gwas_joint_analysis(
        genotype_data = genotype_data,
        phenotype_data = phenotype_data,
        gwas_data = gwas_data,
        gene_region = gene_region,
        run_single_gene = run_single_gene,
        run_v39_genes = run_v39_genes,
        v39_gene_ids = v39_gene_ids,
        output_dir = output_dir,
        tissue_id = tissue_id,
        ld_region = parsed_args$ld_region
    )
} else {
    cat("\nUnknown analysis mode:", analysis_mode, "\n")
    quit(status = 1)
}

cat("\nAnalysis completed successfully!\n")
cat("Results saved to:", output_dir, "\n")
