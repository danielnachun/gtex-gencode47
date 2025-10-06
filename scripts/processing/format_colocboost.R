#!/usr/bin/env Rscript

#' format_colocboost.R
#'
#' Read a colocboost result .rds file and output credible sets tables.
#'
#' - Accepts an input RDS path containing colocboost results produced by colocboost_analysis_pipeline
#'   within a list structure, where elements like xqtl_coloc and joint_gwas may be present.
#' - Extracts trait-specific (uCoS) and trait-shared (CoS) credible sets and associated fields.
#' - Writes a tab-delimited table to <output_prefix>.(xqtl_coloc|joint_gwas).txt for each present element.

# NOTE: argparser is only needed for CLI mode; we import it conditionally below.

format_colocboost <- function(pecotmr_colocboost, output_prefix, verbose = TRUE) {
  if (is.null(pecotmr_colocboost) || identical(pecotmr_colocboost, "")) {
    stop("Missing required argument: pecotmr_colocboost", call. = FALSE)
  }
  if (is.null(output_prefix) || identical(output_prefix, "")) {
    stop("Missing required argument: output_prefix", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    cat("Input RDS:", pecotmr_colocboost, "\n")
    cat("Output prefix:", output_prefix, "\n")
  }

  # ----------------------------
  # Load object
  # ----------------------------
  obj <- readRDS(pecotmr_colocboost)

  # ----------------------------
  # Process available components
  # ----------------------------
  processed_any <- FALSE

  if (is.list(obj) && !is.null(obj$xqtl_coloc)) {
    cs_df <- extract_credible_sets(obj$xqtl_coloc)
    outfile <- paste0(output_prefix, ".xqtl_coloc.txt")
    write_cs_table(cs_df, outfile, verbose = verbose)
    processed_any <- TRUE
  }

  if (is.list(obj) && !is.null(obj$joint_gwas)) {
    cs_df <- extract_credible_sets(obj$joint_gwas)
    outfile <- paste0(output_prefix, ".joint_gwas.txt")
    write_cs_table(cs_df, outfile, verbose = verbose)
    processed_any <- TRUE
  }

  # separate_gwas: a named list of colocboost objects; process each entry
  if (is.list(obj) && !is.null(obj$separate_gwas) && length(obj$separate_gwas) > 0) {
    sg <- obj$separate_gwas
    # If unnamed, assign indices as names
    if (is.null(names(sg)) || any(names(sg) == "")) {
      names(sg) <- if (is.null(names(sg))) seq_along(sg) else ifelse(names(sg) == "", seq_along(sg), names(sg))
    }
    for (key in names(sg)) {
      comp <- sg[[key]]
      if (is.list(comp) && !is.null(comp$data_info)) {
        cs_df <- extract_credible_sets(comp)
        outfile <- paste0(output_prefix, ".separate_gwas.", key, ".txt")
        write_cs_table(cs_df, outfile, verbose = verbose)
        processed_any <- TRUE
      }
    }
  }

  if (!processed_any) {
    stop("No recognizable colocboost components found (xqtl_coloc, joint_gwas, or separate_gwas entries)", call. = FALSE)
  }

  invisible(TRUE)
}

# ----------------------------
# Helpers
# ----------------------------
#' Extract credible sets information from a colocboost result object
#'
#' Expects an object with fields similar to those produced by colocboost_analysis_pipeline:
#' - data_info$variables
#' - data_info$z (list of z-score vectors per outcome)
#' - data_info$outcome_info$outcome_names
#' - ucos_details and/or cos_details
#' - cos_summary (for CoS-level metrics such as cos_npc)
extract_credible_sets <- function(colocboost_result) {
  trait_specific_cs <- list()
  trait_shared_cs <- list()

  # Precompute CoS-level map for cos_npc if available
  cos_npc_map <- NULL
  if (!is.null(colocboost_result$cos_summary)) {
    cs <- colocboost_result$cos_summary
    if (!is.null(cs$cos_id) && !is.null(cs$cos_npc)) {
      cos_npc_map <- setNames(as.numeric(cs$cos_npc), cs$cos_id)
    }
  }

  # Trait-specific credible sets (uCoS)
  if (!is.null(colocboost_result$ucos_details)) {
    ucos_variables <- colocboost_result$ucos_details$ucos$ucos_variables
    ucos_outcomes  <- colocboost_result$ucos_details$ucos_outcomes$outcome_index
    ucos_weights   <- colocboost_result$ucos_details$ucos_weight

    for (i in seq_along(ucos_variables)) {
      cs_name     <- names(ucos_variables)[i]
      variants    <- ucos_variables[[i]]
      outcome_idx <- ucos_outcomes[[i]]

      phenotype_names <- colocboost_result$data_info$outcome_info$outcome_names
      phenotype_id    <- phenotype_names[outcome_idx]

      all_variants   <- colocboost_result$data_info$variables
      variant_indices <- match(variants, all_variants)

      pip_values <- ucos_weights[[i]][variant_indices]

      # -log10(p) from z-scores
      z_scores           <- colocboost_result$data_info$z[[outcome_idx]][variant_indices]
      neg_log10_p_values <- -log10(2 * pnorm(-abs(z_scores)))

      cs_df <- data.frame(
        phenotype_id = rep(phenotype_id, length(variants)),
        variant_id    = variants,
        pip           = pip_values,
        neg_log10_p_value = neg_log10_p_values,
        cs_id         = rep(cs_name, length(variants)),
        cs_type       = "trait_specific",
        cos_npc       = NA_real_,
        stringsAsFactors = FALSE
      )

      trait_specific_cs[[length(trait_specific_cs) + 1]] <- cs_df
    }
  }

  # Trait-shared credible sets (CoS)
  if (!is.null(colocboost_result$cos_details)) {
    cos_variables <- colocboost_result$cos_details$cos$cos_variables
    cos_outcomes  <- colocboost_result$cos_details$cos_outcomes$outcome_index
    cos_weights   <- colocboost_result$cos_details$cos_vcp

    for (i in seq_along(cos_variables)) {
      cs_name         <- names(cos_variables)[i]
      variants        <- cos_variables[[i]]
      outcome_indices <- cos_outcomes[[i]]

      phenotype_names <- colocboost_result$data_info$outcome_info$outcome_names
      phenotype_ids   <- phenotype_names[outcome_indices]

      all_variants     <- colocboost_result$data_info$variables
      variant_indices  <- match(variants, all_variants)
      pip_values       <- cos_weights[[i]][variant_indices]

      # For each variant, collect -log10(p) across outcomes in this CoS
      variant_neg_log10_p_values_list <- vector("list", length(variants))
      for (v in seq_along(variants)) {
        values_for_variant <- numeric(length(outcome_indices))
        for (j in seq_along(outcome_indices)) {
          outcome_idx <- outcome_indices[j]
          z_scores <- colocboost_result$data_info$z[[outcome_idx]][variant_indices]
          neg_log10_p_values <- -log10(2 * pnorm(-abs(z_scores)))
          values_for_variant[j] <- neg_log10_p_values[v]
        }
        variant_neg_log10_p_values_list[[v]] <- values_for_variant
      }

      # Map cos_npc from cos_summary by cs_name if available
      cs_cos_npc <- NA_real_
      if (!is.null(cos_npc_map) && !is.null(cos_npc_map[[cs_name]])) {
        cs_cos_npc <- as.numeric(cos_npc_map[[cs_name]])
      }

      cs_df <- data.frame(
        phenotype_id = rep(paste(phenotype_ids, collapse = ","), length(variants)),
        variant_id    = variants,
        pip           = pip_values,
        neg_log10_p_value = I(variant_neg_log10_p_values_list),
        cs_id         = rep(cs_name, length(variants)),
        cs_type       = "trait_shared",
        cos_npc       = rep(cs_cos_npc, length(variants)),
        stringsAsFactors = FALSE
      )

      trait_shared_cs[[length(trait_shared_cs) + 1]] <- cs_df
    }
  }

  all_cs <- c(trait_specific_cs, trait_shared_cs)
  if (length(all_cs) > 0) {
    result_df <- do.call(rbind, all_cs)
    return(result_df)
  } else {
    return(data.frame(
      phenotype_id = character(0),
      variant_id = character(0),
      pip = numeric(0),
      neg_log10_p_value = list(),
      cs_id = character(0),
      cs_type = character(0),
      cos_npc = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
}

## The helpers below are used by both the function and CLI

# ----------------------------
# Writer
# ----------------------------
write_cs_table <- function(cs_df, outfile, verbose = TRUE) {
  if (nrow(cs_df) > 0) {
    if (is.list(cs_df$neg_log10_p_value)) {
      cs_df$neg_log10_p_value <- vapply(
        cs_df$neg_log10_p_value,
        function(x) toString(signif(as.numeric(x), 6)),
        character(1)
      )
    } else {
      cs_df$neg_log10_p_value <- signif(as.numeric(cs_df$neg_log10_p_value), 6)
    }
    cs_df$pip <- signif(as.numeric(cs_df$pip), 8)
  }

  outdir <- dirname(outfile)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  write.table(cs_df, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)
  if (isTRUE(verbose)) {
    cat("Wrote credible sets table to:", outfile, "\n")
  }
}

# ----------------------------
# CLI entrypoint (only when run via Rscript)
# ----------------------------
args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grep("^--file=", args_all)])

# Only run CLI when this file is executed directly (not when sourced)
if (length(script_path) > 0 && identical(basename(script_path), "format_colocboost.R")) {
  suppressPackageStartupMessages({
    library(argparser)
  })

  parser <- argparser::arg_parser("Format colocboost results into credible sets tables") |>
    argparser::add_argument("--pecotmr_colocboost", help = "Path to colocboost .rds file") |>
    argparser::add_argument("--output_prefix", help = "Output file prefix (directory + stem)") |>
    argparser::add_argument("--quiet", help = "Suppress all printing (no stdout messages)", flag = TRUE)

  args <- argparser::parse_args(parser)

  pecotmr_colocboost <- args$pecotmr_colocboost
  output_prefix <- args$output_prefix
  quiet <- if (!is.null(args$quiet)) args$quiet else FALSE

  invisible(format_colocboost(pecotmr_colocboost, output_prefix, verbose = !quiet))
}