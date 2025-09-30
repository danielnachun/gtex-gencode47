#!/usr/bin/env Rscript

library(argparser)

# parse arguments
parser <- argparser::arg_parser("Script to run colocboost") |>
    argparser::add_argument("--tissue_id", help = "Tissue ID") |>
    argparser::add_argument("--gene_region", help = "Region of genes to analyze") |>
    argparser::add_argument("--genotype_stem", help = "Path to genotype BED file") |>
    argparser::add_argument("--phenotype_list", help = "Path to list of phenotypes to analyze, one per line") |>
    argparser::add_argument("--covariate_path", help = "Path to covariates file") |>
    argparser::add_argument("--variant_region", help = "Region of variants to include (association window region)") |>
    argparser::add_argument("--output_path", help = "Path to save output") |>
    argparser::add_argument("--maf_cutoff", help = "Minor allele frequency cutoff", default = 0.01) |>
    argparser::add_argument("--mac_cutoff", help = "Minor allele count cutoff", default = 0) |>
    argparser::add_argument("--xvar_cutoff", help = "Genotype variance cutoff", default = 0) |>
    argparser::add_argument("--imiss_cutoff", help = "Missingness cutoff", default = 0)

parsed_args <- argparser::parse_args(parser)
tissue_id <- parsed_args$tissue_id
gene_region <- parsed_args$gene_region
genotype_stem <- parsed_args$genotype_stem
phenotype_list <- readLines(parsed_args$phenotype_list)
covariate_path <- parsed_args$covariate_path
variant_region <- parsed_args$variant_region
output_path <- parsed_args$output_path
maf_cutoff <- as.numeric(parsed_args$maf_cutoff)
mac_cutoff <- as.numeric(parsed_args$mac_cutoff)
xvar_cutoff <- as.numeric(parsed_args$xvar_cutoff)
imiss_cutoff <- as.numeric(parsed_args$imiss_cutoff)

cat("Arguments:\n")
cat("  tissue_id:", tissue_id, "\n")
cat("  gene_region:", gene_region, "\n")
cat("  variant_region:", variant_region, "\n")
cat("  genotype_stem:", genotype_stem, "\n")
cat("  phenotype_list:", paste(phenotype_list, collapse = ", "), "\n")
cat("  covariate_path:", covariate_path, "\n")
cat("  output_path:", output_path, "\n")
cat("  maf_cutoff:", maf_cutoff, "\n")
cat("  mac_cutoff:", mac_cutoff, "\n")
cat("  xvar_cutoff:", xvar_cutoff, "\n")
cat("  imiss_cutoff:", imiss_cutoff, "\n")

cat("\nLoading required libraries...\n")
library(pecotmr)
library(colocboost)
cat("Libraries loaded successfully\n")

# load in the data
cat("\nLoading regional data...\n")
data_start_time <- Sys.time()
region_data <- load_multitask_regional_data(
    region = gene_region,
    genotype_list = c(genotype_stem),
    phenotype_list = phenotype_list,
    covariate_list = rep(covariate_path, length(phenotype_list)),
    conditions_list_individual = sub("\\.bed\\.gz$", "", basename(phenotype_list)),
    mac_cutoff = mac_cutoff,
    xvar_cutoff = xvar_cutoff,
    imiss_cutoff = imiss_cutoff,
    association_window = variant_region
)
data_end_time <- Sys.time()
cat("Data loading completed in:", round(difftime(data_end_time, data_start_time, units = "mins"), 2), "minutes\n")

# run colocboost
cat("\nStarting colocboost analysis...\n")
analysis_start_time <- Sys.time()
res <- colocboost_analysis_pipeline(
    region_data,
    pip_cutoff_to_skip_ind = rep(0, length(phenotype_list))
)
analysis_end_time <- Sys.time()
cat("Colocboost analysis completed in:", round(difftime(analysis_end_time, analysis_start_time, units = "mins"), 2), "minutes\n")

# Create output directory if it doesn't exist
output_dir <- dirname(output_path)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
}

# Write the results of the colocboost analysis to a file
cat("Saving results to:", output_path, "\n")
saveRDS(res, file = output_path)
cat("Done\n")