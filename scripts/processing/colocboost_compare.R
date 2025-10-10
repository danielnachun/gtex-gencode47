#!/usr/bin/env Rscript

#' Function to run coloboost on all genes, each gene one at a time, and then on the subset of v39 genes
#' Avoids having to load data multiple times

library(argparser)

parser <- argparser::arg_parser("Script to run colocboost") |>
    argparser::add_argument("--tissue_id", help = "Tissue ID") |>
    argparser::add_argument("--gene_region", help = "Region of genes to analyze") |>
    argparser::add_argument("--variant_region", help = "Region of variants to include (association window region)") |>
    argparser::add_argument("--ld_region", help = "LD block (for naming)") |>
    argparser::add_argument("--genotype_stem", help = "Path to genotype BED file") |>
    argparser::add_argument("--phenotype_list", help = "Path to list of phenotypes to analyze, one path per line") |>
    argparser::add_argument("--covariate_path", help = "Path to covariates file", default = NULL) |>
    argparser::add_argument("--covariate_list", help = "Optional path to list of covariate files, one per phenotype", default = NULL) |>
    argparser::add_argument("--gwas_id_list", help = "Path to list of GWAS IDs to analyze, one per line") |>
    argparser::add_argument("--gwas_phenotype_list", help = "Path to list of GWAS phenotypes to analyze, one path per line") |>
    argparser::add_argument("--gwas_meta", help = "Path to GWAS metadata file") |>
    argparser::add_argument("--ld_meta", help = "Path to LD metadata file") |>
    argparser::add_argument("--gwas_column_matching", help = "Path to GWAS column matching file") |>
    argparser::add_argument("--output_dir", help = "Output directory for results") |>
    argparser::add_argument("--maf_cutoff", help = "Minor allele frequency cutoff", default = 0.01) |>
    argparser::add_argument("--mac_cutoff", help = "Minor allele count cutoff", default = 0) |>
    argparser::add_argument("--xvar_cutoff", help = "Genotype variance cutoff", default = 0) |>
    argparser::add_argument("--imiss_cutoff", help = "Missingness cutoff", default = 0) |>
    argparser::add_argument("--run_single_gene", help = "Run single gene analysis", default = FALSE) |>
    argparser::add_argument("--run_v39_genes", help = "Run v39 genes analysis", default = FALSE) |>
    argparser::add_argument("--v39_gene_id_path", help = "Path to list of v39 gene IDs, one per line", default = NULL)

parsed_args <- argparser::parse_args(parser)

DEBUG <- FALSE

tissue_id         <- parsed_args$tissue_id
gene_region       <- parsed_args$gene_region
genotype_stem     <- parsed_args$genotype_stem
phenotype_list    <- readLines(parsed_args$phenotype_list)

# Check if phenotype_list is empty and handle gracefully
if (length(phenotype_list) == 0 || all(phenotype_list == "")) {
    cat("No phenotypes found for tissue", tissue_id, "- skipping analysis\n")
    quit(status = 0)
}

gwas_id_list <- readLines(parsed_args$gwas_id_list)
gwas_phenotype_list <- readLines(parsed_args$gwas_phenotype_list)
gwas_meta           <- parsed_args$gwas_meta
ld_meta             <- parsed_args$ld_meta
gwas_column_matching <- parsed_args$gwas_column_matching

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
variant_region    <- parsed_args$variant_region
ld_region         <- parsed_args$ld_region
output_dir       <- parsed_args$output_dir
maf_cutoff        <- as.numeric(parsed_args$maf_cutoff)
mac_cutoff        <- as.numeric(parsed_args$mac_cutoff)
xvar_cutoff       <- as.numeric(parsed_args$xvar_cutoff)
imiss_cutoff      <- as.numeric(parsed_args$imiss_cutoff)
run_single_gene   <- parsed_args$run_single_gene
run_v39_genes     <- parsed_args$run_v39_genes
v39_gene_id_path  <- parsed_args$v39_gene_id_path

# Read v39 gene IDs from file if path is provided
if (!is.null(v39_gene_id_path)) {
    v39_gene_id_list <- readLines(v39_gene_id_path)
} else {
    v39_gene_id_list <- NULL
}

# # DEBUG
# gene_region <-  "chr1:0-1000000"
# variant_region <- "chr1:0-1000000"
# ld_region <- "chr1:0-1000000"

cat("Arguments:\n")
cat("  tissue_id:", tissue_id, "\n")
cat("  gene_region:", gene_region, "\n")
cat("  variant_region:", variant_region, "\n")
cat("  ld_region:", ld_region, "\n")
cat("  genotype_stem:", genotype_stem, "\n")
cat("  phenotype_list:", paste(phenotype_list, collapse = ", "), "\n")
cat("  covariate_path:", if (is.null(covariate_path)) "<NULL>" else covariate_path, "\n")
cat("  covariate_list:", if (is.null(covariate_list_path)) "<NULL>" else covariate_list_path, "\n")
cat("  gwas_id_list:", paste(gwas_id_list, collapse = ", "), "\n")
cat("  gwas_phenotype_list:", paste(gwas_phenotype_list, collapse = ", "), "\n")
cat("  gwas_meta:", gwas_meta, "\n")
cat("  ld_meta:", ld_meta, "\n")
cat("  gwas_column_matching:", gwas_column_matching, "\n")
cat("  output_path:", output_dir, "\n")
cat("  maf_cutoff:", maf_cutoff, "\n")
cat("  mac_cutoff:", mac_cutoff, "\n")
cat("  xvar_cutoff:", xvar_cutoff, "\n")
cat("  imiss_cutoff:", imiss_cutoff, "\n")
cat("  run_single_gene:", run_single_gene, "\n")
cat("  run_v39_genes:", run_v39_genes, "\n")
cat("  v39_gene_id_path:", v39_gene_id_path, "\n")


cat("\nLoading required libraries...\n")
library(devtools)
library(magrittr)
library(dplyr)
#library(pecotmr)
library(colocboost)
# TODO: make this not hardcoded, once pecotmr is updated on CRAN
devtools::load_all("/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/pecotmr")

# TODO: make this not hardcoded
# Source formatter to enable function call mode
source("/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing/format_colocboost.R")

cat("Libraries loaded successfully\n")

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
}


# ===========================
# Format GWAS Data
# ===========================

# Read the GWAS metadata file
gwas_meta_df <- read.delim(gwas_meta, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ensure the Tag column is character
gwas_meta_df$Tag <- as.character(gwas_meta_df$Tag)

# Match rows for each GWAS in gwas_id_list, preserving order
gwas_meta_matched <- gwas_meta_df[match(gwas_id_list, gwas_meta_df$Tag), ]

# Check for missing matches and warn if any
if (any(is.na(gwas_meta_matched$Tag))) {
    missing_tags <- gwas_id_list[is.na(gwas_meta_matched$Tag)]
    warning("The following GWAS IDs were not found in the metadata: ", paste(missing_tags, collapse = ", "))
}

# Extract sample size, cases, and controls as numeric vectors
n_samples <- as.numeric(gwas_meta_matched$Sample_Size)
n_cases <- as.numeric(gwas_meta_matched$Cases)
n_controls <- as.numeric(gwas_meta_matched$Controls)



# ===========================
# Load Data for All Genes
# ===========================
gene_names <- sub("\\.bed\\.gz$", "", basename(phenotype_list))

cat("\nLoading regional data...\n")
data_start_time <- Sys.time()
# Validate covariate inputs: require at least one of path or list
if (is.null(covariate_list_path) && is.null(covariate_path)) {
    stop("You must provide either --covariate_list or --covariate_path")
}

# Define the path for saving/loading sumstat_data for this LD block
sumstat_data_path <- file.path(output_dir, paste0(ld_region, ".sumstat_data.rds"))
region_data_path <- file.path(output_dir, paste0(tissue_id, ".", ld_region, ".region_data.rds"))

if (file.exists(sumstat_data_path)) {
    # If sumstat_data already exists, load region_data and overwrite sumstat_data
    cat("sumstat_data file exists for this LD block. Loading previous sumstat_data, and overwriting sumstat_data in region_data.\n")
    region_data <- load_multitask_regional_data(
        region = gene_region,
        genotype_list = c(genotype_stem),
        phenotype_list = phenotype_list,
        covariate_list = if (!is.null(covariate_list_path)) {
            cl <- readLines(covariate_list_path)
            if (length(cl) != length(phenotype_list)) {
                stop("Length of covariate_list (", length(cl), ") does not match phenotype_list (", length(phenotype_list), ")")
            }
            cl
        } else {
            rep(covariate_path, length(phenotype_list))
        },
        conditions_list_individual = gene_names,
        maf_cutoff = maf_cutoff,
        mac_cutoff = mac_cutoff,
        xvar_cutoff = xvar_cutoff,
        imiss_cutoff = imiss_cutoff,
        association_window = variant_region
    )
    sumstat_data <- readRDS(sumstat_data_path)
    region_data$sumstat_data <- sumstat_data
} else {
    # Load region data as usual
    region_data <- load_multitask_regional_data(
        region = gene_region,
        genotype_list = c(genotype_stem),
        phenotype_list = phenotype_list,
        covariate_list = if (!is.null(covariate_list_path)) {
            cl <- readLines(covariate_list_path)
            if (length(cl) != length(phenotype_list)) {
                stop("Length of covariate_list (", length(cl), ") does not match phenotype_list (", length(phenotype_list), ")")
            }
            cl
        } else {
            rep(covariate_path, length(phenotype_list))
        },
        conditions_list_individual = gene_names,
        sumstat_path_list = gwas_phenotype_list,
        column_file_path_list = rep(gwas_column_matching, length(gwas_phenotype_list)),
        LD_meta_file_path_list = rep(ld_meta, length(gwas_phenotype_list)),
        conditions_list_sumstat = gwas_id_list,
        match_LD_sumstat = lapply(gwas_id_list, function(x) c(x)),
        n_samples = n_samples,
        n_cases = n_cases,
        n_controls = n_controls,
        maf_cutoff = maf_cutoff,
        mac_cutoff = mac_cutoff,
        xvar_cutoff = xvar_cutoff,
        imiss_cutoff = imiss_cutoff,
        association_window = variant_region
    )
    # Write out sumstat_data for this LD block
    saveRDS(region_data$sumstat_data, file = sumstat_data_path)
    # Optionally, also save the full region_data for future use
    if (DEBUG) {
        saveRDS(region_data, file = region_data_path)
    }
}

data_end_time <- Sys.time()
cat("Data loading completed in:", round(difftime(data_end_time, data_start_time, units = "mins"), 2), "minutes\n")



# ===========================
# Run Colocboost on All Genes
# ===========================
cat("\nStarting colocboost analysis on all genes...\n")
analysis_start_time <- Sys.time()
res <- colocboost_analysis_pipeline(
    region_data,
    pip_cutoff_to_skip_ind = rep(-1, length(phenotype_list)),
    pip_cutoff_to_skip_sumstat = rep(-1, length(gwas_phenotype_list)),
    qc_method = c("rss_qc"), 
    maf_cutoff=maf_cutoff, 
    xqtl_coloc = FALSE,
    joint_gwas = TRUE,
    separate_gwas = FALSE
)
analysis_end_time <- Sys.time()
cat("Colocboost on all genes completed in:", round(difftime(analysis_end_time, analysis_start_time, units = "mins"), 2), "minutes\n")

output_path <- file.path(output_dir, paste0(tissue_id, ".", ld_region, ".all_genes.colocboost.rds"))
cat("Saving results to:", output_path, "\n")
saveRDS(res, file = output_path) 
# Also write formatted credible set outputs alongside the RDS
format_colocboost(output_path, sub("\\.rds$", "", output_path))




# ===========================
# Run Colocboost Analysis on individual genes
# ===========================
if (run_single_gene) {
    cat("\nRunning colocboost analysis on individual genes...\n")
    
    # Create output directory for single gene results
    single_gene_output_dir <- file.path(output_dir, "individual_genes")

    if (!dir.exists(single_gene_output_dir)) {
        dir.create(single_gene_output_dir, recursive = TRUE)
        cat("Created single gene output directory:", single_gene_output_dir, "\n")
    }
    
    # Loop through each gene
    for (i in seq_along(gene_names)) {
        gene_name <- gene_names[i]
        cat("\nProcessing gene:", gene_name, "(", i, "of", length(gene_names), ")\n")
        
        # Create single-gene region_data by extracting the i-th gene
        region_data_single <- list(
            individual_data = list(
                residual_Y = list(region_data$individual_data$residual_Y[[i]]),
                residual_X = list(region_data$individual_data$residual_X[[i]]),
                residual_Y_scalar = region_data$individual_data$residual_Y_scalar[i],
                residual_X_scalar = region_data$individual_data$residual_X_scalar[i],
                dropped_sample = list(
                    X = list(region_data$individual_data$dropped_sample$X[[i]]),
                    Y = list(region_data$individual_data$dropped_sample$Y[[i]]),
                    covar = list(region_data$individual_data$dropped_sample$covar[[i]])
                ),
                maf = list(region_data$individual_data$maf[[i]]),
                X = region_data$individual_data$X,  # Same for all genes
                chrom = region_data$individual_data$chrom,  # Same for all genes
                grange = region_data$individual_data$grange,  # Same for all genes
                X_variance = list(region_data$individual_data$X_variance[[i]])
            ),
            sumstat_data = region_data$sumstat_data  # Keep the sumstat data structure
        )
        
        # Set names for the single gene so downstream matching finds this gene
        names(region_data_single$individual_data$residual_Y) <- gene_name
        names(region_data_single$individual_data$residual_X) <- gene_name
        names(region_data_single$individual_data$dropped_sample$X) <- gene_name
        names(region_data_single$individual_data$dropped_sample$Y) <- gene_name
        names(region_data_single$individual_data$dropped_sample$covar) <- gene_name
        names(region_data_single$individual_data$maf) <- gene_name
        names(region_data_single$individual_data$X_variance) <- gene_name
        
        # Run colocboost analysis on single gene
        analysis_start_time <- Sys.time()
        res_single <- colocboost_analysis_pipeline(
            region_data_single,
            pip_cutoff_to_skip_ind = c(-1),
            pip_cutoff_to_skip_sumstat = rep(-1, length(gwas_phenotype_list)),
            qc_method = c("rss_qc"), 
            maf_cutoff=maf_cutoff, 
            xqtl_coloc = FALSE,
            joint_gwas = FALSE,
            separate_gwas = TRUE
        )
        analysis_end_time <- Sys.time()
        
        cat("Colocboost on", gene_name, "completed in:", 
            round(difftime(analysis_end_time, analysis_start_time, units = "mins"), 2), "minutes\n")
        
        # Save individual gene result
        single_gene_output_path <- file.path(single_gene_output_dir, 
                                             paste0(tissue_id, ".", ld_region, ".", sub("^.*\\.(ENSG[0-9]+).*$", "\\1", gene_name), ".colocboost.rds"))
        saveRDS(res_single, file = single_gene_output_path)
        cat("Saved result for", gene_name, "to:", single_gene_output_path, "\n")
        # Also write formatted credible set outputs alongside the RDS
        format_colocboost(single_gene_output_path, sub("\\.rds$", "", single_gene_output_path))
    }
    
    cat("\nSingle gene analysis completed for all", length(gene_names), "genes\n")
}

# ===========================
# Run Colocboost Analysis on v39 Gene List
# ===========================
if (run_v39_genes) {
    cat("\nRunning colocboost analysis on v39 gene list...\n")

    all_gene_names <- sub("\\.bed\\.gz$", "", basename(phenotype_list))

    # Match by core Ensembl IDs (strip tissue prefix and any version suffix)
    all_gene_names_core <- sub("^.*\\.(ENSG[0-9]+).*$", "\\1", all_gene_names)
    v39_core <- sub("^(ENSG[0-9]+).*$", "\\1", v39_gene_id_list)

    keep_idx <- all_gene_names_core %in% v39_core

    if (!any(keep_idx)) {
        cat("  No genes from v39_gene_id_list found in phenotype_list (after core ID matching)\n")
        cat("  Skipping v39 gene analysis\n")
    } else {
        gene_names_v39 <- all_gene_names[keep_idx]
        cat("  v39 genes matched:", paste(gene_names_v39, collapse = ", "), "\n")

        # Create v39 region_data by extracting the matching genes from existing data
        region_data_v39 <- list(
            individual_data = list(
                residual_Y = region_data$individual_data$residual_Y[keep_idx],
                residual_X = region_data$individual_data$residual_X[keep_idx],
                residual_Y_scalar = region_data$individual_data$residual_Y_scalar[keep_idx],
                residual_X_scalar = region_data$individual_data$residual_X_scalar[keep_idx],
                dropped_sample = list(
                    X = region_data$individual_data$dropped_sample$X[keep_idx],
                    Y = region_data$individual_data$dropped_sample$Y[keep_idx],
                    covar = region_data$individual_data$dropped_sample$covar[keep_idx]
                ),
                maf = region_data$individual_data$maf[keep_idx],
                X = region_data$individual_data$X,  # Same for all genes
                chrom = region_data$individual_data$chrom,  # Same for all genes
                grange = region_data$individual_data$grange,  # Same for all genes
                X_variance = region_data$individual_data$X_variance[keep_idx]
            ),
            sumstat_data = region_data$sumstat_data  # Keep the sumstat data structure
        )

        cat("Starting colocboost analysis for v39 genes...\n")
        analysis_start_time_v39 <- Sys.time()
        res_v39 <- colocboost_analysis_pipeline(
            region_data_v39,
            pip_cutoff_to_skip_ind = rep(-1, sum(keep_idx)),
            pip_cutoff_to_skip_sumstat = rep(-1, length(gwas_phenotype_list)),
            qc_method = c("rss_qc"), 
            maf_cutoff=maf_cutoff, 
            xqtl_coloc = FALSE,
            joint_gwas = TRUE,
            separate_gwas = FALSE
        )
        analysis_end_time_v39 <- Sys.time()
        cat("Colocboost on v39 genes completed in:", round(difftime(analysis_end_time_v39, analysis_start_time_v39, units = "mins"), 2), "minutes\n")

        output_path_v39 <- file.path(output_dir, paste0(tissue_id, ".", ld_region, ".v39_genes.colocboost.rds"))
        cat("Saving v39 results to:", output_path_v39, "\n")
        saveRDS(res_v39, file = output_path_v39)
        # Also write formatted credible set outputs alongside the RDS
        format_colocboost(output_path_v39, sub("\\.rds$", "", output_path_v39))
    }
}












