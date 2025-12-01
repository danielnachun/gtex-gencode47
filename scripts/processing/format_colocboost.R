#!/usr/bin/env Rscript

#' format_colocboost.R
#'
#' Read a colocboost result .rds file and output robust credible sets tables.
#'
#' - Accepts an input RDS path containing colocboost results produced by colocboost_analysis_pipeline
#'   within a list structure, where elements like xqtl_coloc and joint_gwas may be present.
#' - Applies get_robust_colocalization() to filter results.
#' - Extracts trait-specific (uCoS) and trait-shared (CoS) credible sets and associated fields.
#' - Writes a tab-delimited table to <output_prefix>.(xqtl_coloc|joint_gwas).robust.txt for each present element.
#'

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

  # Helper to get robust and write out
  process_and_write <- function(res, prefix, verbose = TRUE) {
    # Extract credible sets from original (non-robust) results and write out
    cs_df_original <- extract_colocboost_sets(res)
    outfile_original <- paste0(prefix, ".txt")
    write_cs_table(cs_df_original, outfile_original, verbose = verbose)
    
    # Get robust colocalization
    robust_res <- suppressWarnings(colocboost::get_robust_colocalization(
      res,
      cos_npc_cutoff = 0.5,
      npc_outcome_cutoff = 0.2,
      pvalue_cutoff = NULL,
      weight_fudge_factor = 1.5,
      coverage = 0.95
    ))
    
    # Extract credible sets from robust results and write out
    cs_df_robust <- extract_colocboost_sets(robust_res)
    outfile_robust <- paste0(prefix, ".robust.txt")
    write_cs_table(cs_df_robust, outfile_robust, verbose = verbose)
  }

  if (is.list(obj) && !is.null(obj$xqtl_coloc)) {
    process_and_write(obj$xqtl_coloc, paste0(output_prefix, ".xqtl_coloc"), verbose = verbose)
    processed_any <- TRUE
  }

  if (is.list(obj) && !is.null(obj$joint_gwas)) {
    process_and_write(obj$joint_gwas, paste0(output_prefix, ".joint_gwas"), verbose = verbose)
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
        process_and_write(comp, paste0(output_prefix, ".separate_gwas.", key), verbose = verbose)
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
#' Extract colocboost credible sets into a long-format data.frame
#'
#' @param cb_output A "colocboost" object, ideally after get_robust_colocalization()
#' @param outcome_names Optional character vector giving outcome names in the
#'   same order as Y used in the original analysis (overrides cb_output$data_info$outcome_info$outcome_names).
#'
#' @return data.frame with columns:
#'   phenotype_id, variant_id, cs_id, cs_type ("cos" or "ucos"),
#'   vcp (NA for UCoS), cos_npc (NA for UCoS), npc_outcome (NA for UCoS),
#'   ucos_weight (NA for CoS), cs_change (list for CoS with one value per outcome,
#'   numeric for UCoS with single value), neg_log10_p_value (list for CoS, numeric for UCoS).
extract_colocboost_sets <- function(cb_output, outcome_names = NULL) {
  if (!inherits(cb_output, "colocboost")) {
    stop("cb_output must be a 'colocboost' object")
  }

  # Global variant + outcome names
  variables <- cb_output$data_info$variables
  outcomes  <- cb_output$data_info$outcome_info$outcome_names
  if (!is.null(outcome_names)) {
    if (length(outcome_names) != length(outcomes)) {
      stop("outcome_names length must match number of outcomes in cb_output$data_info$outcome_info$outcome_names")
    }
    outcomes <- outcome_names
  }

  rows <- list()

  ## -----------------------------
  ## 1. Colocalized CoS (multi-trait)
  ## -----------------------------
  cd <- cb_output$cos_details
  if (!is.null(cd) && !is.null(cd$cos$cos_index) && length(cd$cos$cos_index) > 0) {
    cos_indices   <- cd$cos$cos_index
    cos_variables <- cd$cos$cos_variables
    cos_outcomes  <- cd$cos_outcomes$outcome_index
    cos_npc       <- cd$cos_npc
    cos_min_npc   <- cd$cos_min_npc_outcome
    cos_vcp_list  <- cd$cos_vcp
    cos_cs_change <- cd$cos_cs_change

    for (i in seq_along(cos_indices)) {
      cs_name <- names(cos_indices)[i]
      var_idx <- cos_indices[[i]]
      if (length(var_idx) == 0) next

      variant_id <- variables[var_idx]

      # Outcome IDs for this CoS (collapse to a single string)
      out_idx <- cos_outcomes[[i]]
      phenotype_id <- paste(outcomes[out_idx], collapse = ";")

      # Variant-level colocalization probability (VCP) for this CoS
      vcp_vec <- cos_vcp_list[[i]]
      # vcp_vec is over all variants; subset to this CoS' indices if necessary
      if (length(vcp_vec) == length(variables)) {
        vcp <- vcp_vec[var_idx]
      } else if (length(vcp_vec) == length(var_idx)) {
        vcp <- vcp_vec
      } else {
        # Fallback: length mismatch, fill with NA
        vcp <- rep(NA_real_, length(var_idx))
      }

      cos_npc_val      <- as.numeric(cos_npc[i])
      npc_outcome_val  <- as.numeric(cos_min_npc[i])

      # Extract cs_change values for all outcomes in this CoS
      cs_change_values <- rep(NA_real_, length(out_idx))
      if (!is.null(cos_cs_change) && cs_name %in% rownames(cos_cs_change)) {
        cs_change_row <- cos_cs_change[cs_name, , drop = FALSE]
        for (j in seq_along(out_idx)) {
          outcome_idx <- out_idx[j]
          outcome_name <- outcomes[outcome_idx]
          # Try to get cs_change by outcome name first, then by index
          if (outcome_name %in% colnames(cs_change_row)) {
            cs_change_values[j] <- as.numeric(cs_change_row[1, outcome_name])
          } else if (as.character(outcome_idx) %in% colnames(cs_change_row)) {
            cs_change_values[j] <- as.numeric(cs_change_row[1, as.character(outcome_idx)])
          } else if (outcome_idx <= ncol(cs_change_row)) {
            cs_change_values[j] <- as.numeric(cs_change_row[1, outcome_idx])
          }
        }
      }
      # For CoS, cs_change is a list with one value per outcome (same order as phenotype_id)
      cs_change_list <- rep(list(cs_change_values), length(var_idx))

      # Compute -log10(p) for each variant across outcomes in this CoS
      # For CoS, we have multiple outcomes, so store as a list
      variant_neg_log10_p_values_list <- vector("list", length(var_idx))
      for (v in seq_along(var_idx)) {
        values_for_variant <- numeric(length(out_idx))
        for (j in seq_along(out_idx)) {
          outcome_idx <- out_idx[j]
          z_scores <- cb_output$data_info$z[[outcome_idx]][var_idx]
          neg_log10_p_values <- -log10(2 * pnorm(-abs(z_scores)))
          values_for_variant[j] <- neg_log10_p_values[v]
        }
        variant_neg_log10_p_values_list[[v]] <- values_for_variant
      }

      rows[[length(rows) + 1L]] <- data.frame(
        phenotype_id = rep(phenotype_id, length(var_idx)),
        variant_id   = variant_id,
        cs_id        = rep(cs_name, length(var_idx)),
        cs_type      = rep("cos",  length(var_idx)),
        vcp          = vcp,
        cos_npc      = rep(cos_npc_val,     length(var_idx)),
        npc_outcome  = rep(npc_outcome_val, length(var_idx)),
        ucos_weight  = rep(NA_real_,        length(var_idx)),
        cs_change    = I(cs_change_list),
        neg_log10_p_value = I(variant_neg_log10_p_values_list),
        stringsAsFactors = FALSE
      )
    }
  }

  ## -----------------------------
  ## 2. Trait-specific UCoS
  ## -----------------------------
  ud <- cb_output$ucos_details
  if (!is.null(ud) &&
      !is.null(ud$ucos$ucos_index) &&
      length(ud$ucos$ucos_index) > 0) {

    ucos_indices   <- ud$ucos$ucos_index
    ucos_variables <- ud$ucos$ucos_variables
    ucos_outcomes  <- ud$ucos_outcomes
    ucos_weights   <- ud$ucos_weight
    ucos_cs_change <- ud$ucos_outcomes_delta

    for (j in seq_along(ucos_indices)) {
      u_name  <- names(ucos_indices)[j]
      var_idx <- ucos_indices[[j]]
      if (length(var_idx) == 0) next

      variant_id <- variables[var_idx]

      # Each UCoS should be for a single outcome; extract outcome name or index
      # outcome_index and outcome_name are both lists where each element is a list
      out_idx_list <- ucos_outcomes$outcome_index[[j]]
      out_name_list <- ucos_outcomes$outcome_name[[j]]
      
      # Try to get outcome name first (more reliable)
      if (!is.null(out_name_list) && length(out_name_list) > 0) {
        phenotype_id <- out_name_list[[1]]
        # Get corresponding index
        if (!is.null(out_idx_list) && length(out_idx_list) > 0) {
          out_idx <- out_idx_list[[1]]
        } else {
          # Try to find index from name
          out_idx <- match(phenotype_id, outcomes)
        }
      } else if (!is.null(out_idx_list) && length(out_idx_list) > 0) {
        # Fall back to index if name not available
        out_idx <- out_idx_list[[1]]
        if (out_idx >= 1 && out_idx <= length(outcomes)) {
          phenotype_id <- outcomes[out_idx]
        } else {
          phenotype_id <- NA_character_
        }
      } else {
        phenotype_id <- NA_character_
        out_idx <- NA_integer_
      }

      # UCoS weights: typically a matrix (variants x 1 outcome)
      w <- ucos_weights[[j]]
      w_vec <- as.numeric(w)
      if (length(w_vec) != length(var_idx)) {
        # Length mismatch: be explicit
        w_vec <- rep(NA_real_, length(var_idx))
      }

      # Extract cs_change for this UCoS (single value)
      cs_change_val <- NA_real_
      if (!is.null(ucos_cs_change) && is.data.frame(ucos_cs_change) && 
          "ucos_outcome" %in% colnames(ucos_cs_change) && 
          "ucos_delta" %in% colnames(ucos_cs_change)) {
        # Match by outcome name (phenotype_id)
        if (!is.na(phenotype_id) && phenotype_id %in% ucos_cs_change$ucos_outcome) {
          match_idx <- which(ucos_cs_change$ucos_outcome == phenotype_id)
          if (length(match_idx) > 0) {
            # Take the first match (or if multiple UCoS for same outcome, match by position)
            # Since we iterate in order, try to match by position j
            if (j <= nrow(ucos_cs_change) && length(match_idx) > 1) {
              # Try to find a match that hasn't been used yet by checking previous iterations
              used_indices <- c()
              for (k in seq_len(j - 1)) {
                prev_outcome <- NULL
                # Try to get outcome name first
                if (!is.null(ucos_outcomes$outcome_name[[k]]) && length(ucos_outcomes$outcome_name[[k]]) > 0) {
                  prev_outcome <- ucos_outcomes$outcome_name[[k]][[1]]
                } else if (!is.null(ucos_outcomes$outcome_index[[k]]) && length(ucos_outcomes$outcome_index[[k]]) > 0) {
                  prev_idx <- ucos_outcomes$outcome_index[[k]][[1]]
                  if (prev_idx >= 1 && prev_idx <= length(outcomes)) {
                    prev_outcome <- outcomes[prev_idx]
                  }
                }
                if (!is.null(prev_outcome) && prev_outcome == phenotype_id) {
                  prev_match <- which(ucos_cs_change$ucos_outcome == prev_outcome)
                  if (length(prev_match) > 0) {
                    used_indices <- c(used_indices, prev_match[1])
                  }
                }
              }
              available_idx <- setdiff(match_idx, used_indices)
              if (length(available_idx) > 0) {
                cs_change_val <- as.numeric(ucos_cs_change$ucos_delta[available_idx[1]])
              } else {
                cs_change_val <- as.numeric(ucos_cs_change$ucos_delta[match_idx[1]])
              }
            } else {
              cs_change_val <- as.numeric(ucos_cs_change$ucos_delta[match_idx[1]])
            }
          }
        }
      }

      # Compute -log10(p) for uCoS (single outcome)
      neg_log10_p_values <- rep(NA_real_, length(var_idx))
      if (!is.null(cb_output$data_info$z) && !is.na(out_idx) && length(cb_output$data_info$z) >= out_idx) {
        z_scores <- cb_output$data_info$z[[out_idx]][var_idx]
        neg_log10_p_values <- -log10(2 * pnorm(-abs(z_scores)))
      }

      rows[[length(rows) + 1L]] <- data.frame(
        phenotype_id = rep(phenotype_id, length(var_idx)),
        variant_id   = variant_id,
        cs_id        = rep(u_name, length(var_idx)),
        cs_type      = rep("ucos", length(var_idx)),
        vcp          = rep(NA_real_, length(var_idx)),
        cos_npc      = rep(NA_real_, length(var_idx)),
        npc_outcome  = rep(NA_real_, length(var_idx)),
        ucos_weight  = w_vec,
        cs_change    = rep(cs_change_val, length(var_idx)),
        neg_log10_p_value = neg_log10_p_values,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    return(data.frame(
      phenotype_id = character(0),
      variant_id   = character(0),
      cs_id        = character(0),
      cs_type      = character(0),
      vcp          = numeric(0),
      cos_npc      = numeric(0),
      npc_outcome  = numeric(0),
      ucos_weight  = numeric(0),
      cs_change    = numeric(0),
      neg_log10_p_value = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, rows)
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
    
    # Handle npc_outcome if it's a list (for trait-shared credible sets)
    if (is.list(cs_df$npc_outcome)) {
      cs_df$npc_outcome <- vapply(
        cs_df$npc_outcome,
        function(x) toString(signif(as.numeric(x), 6)),
        character(1)
      )
    } else {
      cs_df$npc_outcome <- signif(as.numeric(cs_df$npc_outcome), 6)
    }
    
    # Format vcp and ucos_weight if present
    if ("vcp" %in% names(cs_df)) {
      cs_df$vcp <- signif(as.numeric(cs_df$vcp), 8)
    }
    if ("ucos_weight" %in% names(cs_df)) {
      cs_df$ucos_weight <- signif(as.numeric(cs_df$ucos_weight), 8)
    }
    
    # Handle cs_change if it's a list (for CoS)
    if ("cs_change" %in% names(cs_df)) {
      if (is.list(cs_df$cs_change)) {
        cs_df$cs_change <- vapply(
          cs_df$cs_change,
          function(x) toString(signif(as.numeric(x), 6)),
          character(1)
        )
      } else {
        cs_df$cs_change <- signif(as.numeric(cs_df$cs_change), 6)
      }
    }
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