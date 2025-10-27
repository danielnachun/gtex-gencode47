#!/usr/bin/env Rscript

#' format_colocboost.R
#'
#' Read a colocboost result .rds file and output credible sets tables.
#'
#' - Accepts an input RDS path containing colocboost results produced by colocboost_analysis_pipeline
#'   within a list structure, where elements like xqtl_coloc and joint_gwas may be present.
#' - Extracts trait-specific (uCoS) and trait-shared (CoS) credible sets and associated fields.
#' - Writes a tab-delimited table to <output_prefix>.(xqtl_coloc|joint_gwas).txt for each present element.
#' - Also writes a "robust" version of each result, filtered by get_robust_colocalization(...).
#'

# NOTE: argparser is only needed for CLI mode; we import it conditionally below.

# Custom function to convert filtered trait_shared entries to trait_specific
get_robust_colocalization_with_conversion <- function(cb_output,
                                                      original_obj = NULL,  # Not used in simplified approach
                                                      cos_npc_cutoff = 0.5,
                                                      npc_outcome_cutoff = 0.2,
                                                      pvalue_cutoff = NULL,
                                                      weight_fudge_factor = 1.5,
                                                      coverage = 0.95,
                                                      verbose = TRUE) {
  # Check if cb_output is a valid colocboost object
  if (is.null(cb_output)) {
    if (verbose) cat("Warning: colocboost object is NULL, returning original object\n")
    return(cb_output)
  }
  if (!inherits(cb_output, "colocboost")) {
    if (verbose) cat("Warning: Object is not a valid colocboost object (class:", class(cb_output), "), returning original object\n")
    return(cb_output)
  }
  
  # First apply standard robust filtering
  if (verbose) cat("DEBUG: Applying robust filtering with cos_npc_cutoff =", cos_npc_cutoff, "and npc_outcome_cutoff =", npc_outcome_cutoff, "\n")
  tryCatch({
    # Suppress the "No colocalization results in this region!" warning
    robust_res <- suppressWarnings(colocboost::get_robust_colocalization(
      cb_output, 
      cos_npc_cutoff = cos_npc_cutoff, 
      npc_outcome_cutoff = npc_outcome_cutoff,
      pvalue_cutoff = pvalue_cutoff,
      weight_fudge_factor = weight_fudge_factor,
      coverage = coverage
    ))
    if (verbose) cat("DEBUG: get_robust_colocalization completed successfully\n")
  }, error = function(e) {
    if (verbose) cat("Warning: get_robust_colocalization failed:", e$message, "\n")
    return(cb_output)
  })
  
  # Check if robust_res has the expected structure
  if (is.null(robust_res) || is.null(robust_res$cos_details) || is.null(robust_res$cos_details$cos)) {
    if (verbose) cat("Warning: Robust results do not have expected structure, returning original object\n")
    if (verbose) cat("DEBUG: robust_res is null:", is.null(robust_res), "\n")
    if (!is.null(robust_res)) {
      if (verbose) cat("DEBUG: robust_res$cos_details is null:", is.null(robust_res$cos_details), "\n")
      if (!is.null(robust_res$cos_details)) {
        if (verbose) cat("DEBUG: robust_res$cos_details$cos is null:", is.null(robust_res$cos_details$cos), "\n")
      }
    }
    return(cb_output)
  }
  
  if (verbose) cat("DEBUG: Robust filtering successful, checking results...\n")
  
  # Identify filtered trait_shared entries by comparing original and robust results
  original_cos_ids <- names(cb_output$cos_details$cos$cos_index)
  robust_cos_ids <- names(robust_res$cos_details$cos$cos_index)
  filtered_cos_ids <- setdiff(original_cos_ids, robust_cos_ids)
  
  if (verbose) cat("DEBUG: Original COS sets:", length(original_cos_ids), "\n")
  if (verbose) cat("DEBUG: Robust COS sets:", length(robust_cos_ids), "\n")
  if (verbose) cat("DEBUG: Filtered COS sets:", length(filtered_cos_ids), "\n")
  if (verbose) cat("Found", length(filtered_cos_ids), "filtered colocalization sets to convert\n")
  
  if (length(filtered_cos_ids) > 0) {
    # Get original cos_npc values for filtered entries
    original_cos_npc_map <- NULL
    if (!is.null(cb_output$cos_summary)) {
      cs <- cb_output$cos_summary
      if (!is.null(cs$cos_id) && !is.null(cs$cos_npc)) {
        original_cos_npc_map <- setNames(as.numeric(cs$cos_npc), cs$cos_id)
      }
    }
    
    # Check which phenotype-variant combinations are already covered in robust results
    # We need to track which specific phenotype-variant pairs are already in robust CoS
    robust_phenotype_variant_pairs <- c()
    if (!is.null(robust_res$cos_details$cos$cos_variables) && !is.null(robust_res$cos_details$cos_outcomes$outcome_index)) {
      for (robust_cos_idx in seq_along(robust_res$cos_details$cos$cos_variables)) {
        robust_variants <- robust_res$cos_details$cos$cos_variables[[robust_cos_idx]]
        robust_outcomes <- robust_res$cos_details$cos_outcomes$outcome_index[[robust_cos_idx]]
        outcome_names <- robust_res$cos_details$cos_outcomes$outcome_name[[robust_cos_idx]]
        
        # Create phenotype-variant pairs for this robust CoS
        for (outcome_idx in robust_outcomes) {
          phenotype_name <- outcome_names[which(robust_outcomes == outcome_idx)][1]
          for (variant in robust_variants) {
            pair <- paste(phenotype_name, variant, sep = "|")
            robust_phenotype_variant_pairs <- c(robust_phenotype_variant_pairs, pair)
          }
        }
      }
    }
    
    # Create converted ucos entries from filtered cos entries
    # Only convert outcomes that aren't already in smaller joint CoS sets
    converted_ucos <- list()
    
    for (cos_id in filtered_cos_ids) {
      # Get the original colocalization details
      cos_idx <- which(names(cb_output$cos_details$cos$cos_index) == cos_id)
      if (length(cos_idx) > 0) {
        # Get the outcomes for this colocalization
        outcomes <- cb_output$cos_details$cos_outcomes$outcome_index[[cos_idx]]
        outcome_names <- cb_output$cos_details$cos_outcomes$outcome_name[[cos_idx]]
        cos_variables <- cb_output$cos_details$cos$cos_variables[[cos_idx]]
        cos_weights <- cb_output$cos_details$cos_weights[[cos_idx]]
        
        # Get original cos_npc for this colocalization
        original_cos_npc <- if (!is.null(original_cos_npc_map) && !is.null(original_cos_npc_map[[cos_id]])) {
          original_cos_npc_map[[cos_id]]
        } else {
          NA_real_
        }
        
        # Get npc_outcome values for this CoS if available
        cos_npc_outcome_values <- NULL
        if (!is.null(cb_output$cos_details$cos_outcomes_npc) && 
            !is.null(cb_output$cos_details$cos_outcomes_npc[[cos_id]])) {
          cos_npc_data <- cb_output$cos_details$cos_outcomes_npc[[cos_id]]
          if (!is.null(cos_npc_data$npc_outcome)) {
            # Reorder npc_outcome values to match the phenotype order in outcomes
            cos_npc_outcome_values <- cos_npc_data$npc_outcome[match(outcomes, cos_npc_data$outcomes_index)]
          }
        }
        
        if (verbose) cat("Processing", length(outcomes), "traits from", cos_id, "(original cos_npc:", original_cos_npc, ")\n")
        
        # For each outcome, create a trait_specific entry only if not already in smaller joint CoS
        for (i in seq_along(outcomes)) {
          outcome_idx <- outcomes[i]
          outcome_name <- outcome_names[i]
          
          # Check if this phenotype-variant combination is already covered in robust results
          # Create phenotype-variant pairs for this outcome
          phenotype_variant_pairs <- paste(rep(outcome_name, length(cos_variables)), cos_variables, sep = "|")
          
          # Check if any of these pairs are already covered
          already_covered <- any(phenotype_variant_pairs %in% robust_phenotype_variant_pairs)
          
          if (already_covered) {
            if (verbose) cat("  Skipping", outcome_name, "- phenotype-variant combinations already covered in robust results\n")
            next
          }
          
          # Get the trait-specific npc_outcome value for this outcome
          trait_specific_npc <- if (!is.null(cos_npc_outcome_values) && length(cos_npc_outcome_values) >= i) {
            cos_npc_outcome_values[i]
          } else {
            original_cos_npc  # Fallback to overall cos_npc if trait-specific not available
          }
          
          # Create a new trait_specific entry using the same variants but individual trait weights
          ucos_entry <- list(
            ucos_index = cos_variables,
            ucos_variables = cos_variables,
            outcome_index = outcome_idx,
            outcome_name = outcome_name,
            ucos_weight = cos_weights[, i, drop = FALSE],
            ucos_purity = list(
              min_abs_cor = 1.0,  # Single trait, perfect correlation with itself
              median_abs_cor = 1.0,
              max_abs_cor = 1.0
            ),
            # Add flag to indicate this came from a filtered trait_shared entry
            converted_from_shared = TRUE,
            original_cos_id = cos_id,
            original_cos_npc = trait_specific_npc,  # Store the trait-specific npc_outcome value
            original_npc_outcome = trait_specific_npc  # Store the original npc_outcome value for this trait
          )
          
          converted_ucos[[length(converted_ucos) + 1]] <- ucos_entry
        }
      }
    }
    
    # Add converted entries to ucos_details
    if (length(converted_ucos) > 0) {
      if (verbose) cat("Created", length(converted_ucos), "converted trait_specific entries\n")
      
      # Initialize ucos_details if it doesn't exist
      if (is.null(robust_res$ucos_details)) {
        robust_res$ucos_details <- list(
          ucos = list(ucos_index = list(), ucos_variables = list()),
          ucos_outcomes = list(outcome_index = list(), outcome_name = list()),
          ucos_weight = list(),
          ucos_purity = list(min_abs_cor = matrix(1), median_abs_cor = matrix(1), max_abs_cor = matrix(1)),
          ucos_top_variables = data.frame(),
          converted_from_shared = list(),
          original_cos_id = list(),
          original_cos_npc = list(),
          original_npc_outcome = list()
        )
      } else {
        # If ucos_details already exists, initialize conversion metadata for existing entries
        existing_count <- length(robust_res$ucos_details$ucos$ucos_variables)
        if (length(robust_res$ucos_details$converted_from_shared) == 0) {
          # Initialize conversion metadata for existing entries
          for (i in seq_len(existing_count)) {
            robust_res$ucos_details$converted_from_shared[[i]] <- FALSE
            robust_res$ucos_details$original_cos_id[[i]] <- NA_character_
            robust_res$ucos_details$original_cos_npc[[i]] <- NA_real_
            robust_res$ucos_details$original_npc_outcome[[i]] <- NA_real_
          }
        }
      }
      
      # Add converted entries
      for (i in seq_along(converted_ucos)) {
        ucos_entry <- converted_ucos[[i]]
        entry_name <- paste0("ucos", i, ":y", ucos_entry$outcome_index)
        
        robust_res$ucos_details$ucos$ucos_index[[length(robust_res$ucos_details$ucos$ucos_index) + 1]] <- ucos_entry$ucos_index
        robust_res$ucos_details$ucos$ucos_variables[[length(robust_res$ucos_details$ucos$ucos_variables) + 1]] <- ucos_entry$ucos_variables
        robust_res$ucos_details$ucos_outcomes$outcome_index[[length(robust_res$ucos_details$ucos_outcomes$outcome_index) + 1]] <- ucos_entry$outcome_index
        robust_res$ucos_details$ucos_outcomes$outcome_name[[length(robust_res$ucos_details$ucos_outcomes$outcome_name) + 1]] <- ucos_entry$outcome_name
        robust_res$ucos_details$ucos_weight[[length(robust_res$ucos_details$ucos_weight) + 1]] <- ucos_entry$ucos_weight
        
        # Set names for the lists
        names(robust_res$ucos_details$ucos$ucos_variables)[length(robust_res$ucos_details$ucos$ucos_variables)] <- entry_name
        names(robust_res$ucos_details$ucos_outcomes$outcome_index)[length(robust_res$ucos_details$ucos_outcomes$outcome_index)] <- entry_name
        
        # Add conversion flags and original cos_npc
        robust_res$ucos_details$converted_from_shared[[length(robust_res$ucos_details$converted_from_shared) + 1]] <- ucos_entry$converted_from_shared
        robust_res$ucos_details$original_cos_id[[length(robust_res$ucos_details$original_cos_id) + 1]] <- ucos_entry$original_cos_id
        robust_res$ucos_details$original_cos_npc[[length(robust_res$ucos_details$original_cos_npc) + 1]] <- ucos_entry$original_cos_npc
        robust_res$ucos_details$original_npc_outcome[[length(robust_res$ucos_details$original_npc_outcome) + 1]] <- ucos_entry$original_npc_outcome
      }
    }
  }
  
  return(robust_res)
}



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

  # Helper to process and write both normal and robust versions
  process_and_write <- function(res, prefix, verbose = TRUE) {
    # Standard
    cs_df <- extract_credible_sets(res, verbose = verbose, include_conversion_metadata = FALSE)
    outfile <- paste0(prefix, ".txt")
    write_cs_table(cs_df, outfile, verbose = verbose)

    # Robust with trait_shared to trait_specific conversion
    # Pass the original colocboost object for proper ucos computation
    robust_res <- get_robust_colocalization_with_conversion(res, obj, cos_npc_cutoff = 0.5, npc_outcome_cutoff = 0.2, verbose = verbose)
    cs_df_robust <- extract_credible_sets(robust_res, verbose = verbose, include_conversion_metadata = TRUE)
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
#' Extract credible sets information from a colocboost result object
#'
#' Expects an object with fields similar to those produced by colocboost_analysis_pipeline:
#' - data_info$variables
#' - data_info$z (list of z-score vectors per outcome)
#' - data_info$outcome_info$outcome_names
#' - ucos_details and/or cos_details
#' - cos_summary (for CoS-level metrics such as cos_npc)
extract_credible_sets <- function(colocboost_result, verbose = TRUE, include_conversion_metadata = FALSE) {
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
    if (verbose) cat("Available ucos_details fields:", names(colocboost_result$ucos_details), "\n")
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
      
      # npc_outcome is not available for trait-specific credible sets
      npc_outcome_values <- NA_real_

      # Check if this is a converted entry from a filtered trait_shared
      cs_type <- "trait_specific"
      original_cos_npc <- NA_real_
      if (!is.null(colocboost_result$ucos_details$converted_from_shared) && 
          length(colocboost_result$ucos_details$converted_from_shared) >= i &&
          colocboost_result$ucos_details$converted_from_shared[[i]]) {
        cs_type <- "trait_specific_from_shared"
        # Get the original cos_npc value for converted entries
        if (!is.null(colocboost_result$ucos_details$original_cos_npc) && 
            length(colocboost_result$ucos_details$original_cos_npc) >= i) {
          original_cos_npc <- colocboost_result$ucos_details$original_cos_npc[[i]]
        }
      }
      
      # Check for mismatched vector lengths before creating data.frame
      if (length(variants) == 0) {
        warning(paste("Skipping credible set", cs_name, "- empty variants vector"))
        next
      }
      if (is.null(cs_name) || length(cs_name) == 0 || cs_name == "") {
        warning(paste("Skipping credible set - empty or null cs_name"))
        next
      }
      if (length(pip_values) != length(variants)) {
        warning(paste("Skipping credible set", cs_name, "- pip_values length (", length(pip_values), ") != variants length (", length(variants), ")"))
        next
      }
      if (length(neg_log10_p_values) != length(variants)) {
        warning(paste("Skipping credible set", cs_name, "- neg_log10_p_values length (", length(neg_log10_p_values), ") != variants length (", length(variants), ")"))
        next
      }
      
      # Get additional metadata for converted entries
      original_cos_id_value <- NA_character_
      original_npc_outcome_value <- NA_real_
      if (cs_type == "trait_specific_from_shared") {
        if (!is.null(colocboost_result$ucos_details$original_cos_id) && 
            length(colocboost_result$ucos_details$original_cos_id) >= i) {
          original_cos_id_value <- colocboost_result$ucos_details$original_cos_id[[i]]
        }
        # Get the original npc_outcome value for converted entries
        if (!is.null(colocboost_result$ucos_details$original_npc_outcome) && 
            length(colocboost_result$ucos_details$original_npc_outcome) >= i) {
          original_npc_outcome_value <- colocboost_result$ucos_details$original_npc_outcome[[i]]
        }
      }
      
      if (include_conversion_metadata) {
        cs_df <- data.frame(
          phenotype_id = rep(phenotype_id, length(variants)),
          variant_id    = variants,
          pip           = pip_values,
          neg_log10_p_value = neg_log10_p_values,
          cs_id         = rep(cs_name, length(variants)),
          cs_type       = cs_type,
          cos_npc       = if (cs_type == "trait_specific_from_shared") rep(original_cos_npc, length(variants)) else rep(NA_real_, length(variants)),
          npc_outcome   = rep(npc_outcome_values, length(variants)),
          original_cos_id = if (cs_type == "trait_specific_from_shared") rep(original_cos_id_value, length(variants)) else rep(NA_character_, length(variants)),
          original_cos_npc = if (cs_type == "trait_specific_from_shared") rep(original_cos_npc, length(variants)) else rep(NA_real_, length(variants)),
          original_npc_outcome = if (cs_type == "trait_specific_from_shared") rep(original_npc_outcome_value, length(variants)) else rep(NA_real_, length(variants)),
          converted_from_shared = rep(cs_type == "trait_specific_from_shared", length(variants)),
          stringsAsFactors = FALSE
        )
      } else {
        cs_df <- data.frame(
          phenotype_id = rep(phenotype_id, length(variants)),
          variant_id    = variants,
          pip           = pip_values,
          neg_log10_p_value = neg_log10_p_values,
          cs_id         = rep(cs_name, length(variants)),
          cs_type       = "trait_specific",  # Always use standard cs_type for non-robust
          cos_npc       = NA_real_,
          npc_outcome   = rep(npc_outcome_values, length(variants)),
          original_npc_outcome = NA_real_,  # Not available for non-converted entries
          stringsAsFactors = FALSE
        )
      }

      trait_specific_cs[[length(trait_specific_cs) + 1]] <- cs_df
    }
  }

  # Trait-shared credible sets (CoS)
  if (!is.null(colocboost_result$cos_details)) {
    if (verbose) cat("Available cos_details fields:", names(colocboost_result$cos_details), "\n")
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
      
      # Extract npc_outcome for trait-shared credible sets from cos_outcomes_npc
      # For trait-shared credible sets, we need to match each variant to its corresponding npc_outcome values
      # The npc_outcome values must be ordered to match the phenotype order in phenotype_ids
      variant_npc_outcome_list <- vector("list", length(variants))
      
      # Try to find npc_outcome data - first try the current cs_name, then try original names
      cos_npc_data <- NULL
      npc_source_name <- cs_name
      
      if (!is.null(colocboost_result$cos_details$cos_outcomes_npc) && 
          !is.null(colocboost_result$cos_details$cos_outcomes_npc[[cs_name]])) {
        cos_npc_data <- colocboost_result$cos_details$cos_outcomes_npc[[cs_name]]
      } else {
        # Try to find npc_outcome data under original colocalization names
        # Look for names that contain the same cos number and outcomes
        if (!is.null(colocboost_result$cos_details$cos_outcomes_npc)) {
          available_npc_names <- names(colocboost_result$cos_details$cos_outcomes_npc)
          # Extract cos number from current name (e.g., "cos2:y85_y86_y87" -> "cos2")
          cos_number <- sub(":.*", "", cs_name)
          # Find original names that start with the same cos number
          matching_names <- available_npc_names[grepl(paste0("^", cos_number, ":"), available_npc_names)]
          
          if (length(matching_names) > 0) {
            # Use the first matching name (there should typically be only one)
            npc_source_name <- matching_names[1]
            cos_npc_data <- colocboost_result$cos_details$cos_outcomes_npc[[npc_source_name]]
            if (verbose) cat("Found npc_outcome data under original name:", npc_source_name, "\n")
          }
        }
      }
      
      if (!is.null(cos_npc_data) && !is.null(cos_npc_data$npc_outcome)) {
        # Reorder npc_outcome values to match the phenotype order in outcome_indices
        # cos_npc_data$outcomes_index gives the original order, we need to reorder to match outcome_indices
        npc_reordered <- cos_npc_data$npc_outcome[match(outcome_indices, cos_npc_data$outcomes_index)]
        
        # For trait-shared credible sets, each variant gets the reordered npc_outcome values
        for (v in seq_along(variants)) {
          variant_npc_outcome_list[[v]] <- npc_reordered
        }
        if (verbose) cat("Found npc_outcome for cos", cs_name, ":", paste(npc_reordered, collapse = ", "), "\n")
      } else {
        if (verbose) cat("npc_outcome not found for cos", cs_name, "\n")
        # Fill with NA values if not found
        for (v in seq_along(variants)) {
          variant_npc_outcome_list[[v]] <- rep(NA_real_, length(outcome_indices))
        }
      }

      # Check for mismatched vector lengths before creating data.frame
      if (length(variants) == 0) {
        warning(paste("Skipping trait-shared credible set", cs_name, "- empty variants vector"))
        next
      }
      if (is.null(cs_name) || length(cs_name) == 0 || cs_name == "") {
        warning(paste("Skipping trait-shared credible set - empty or null cs_name"))
        next
      }
      if (length(pip_values) != length(variants)) {
        warning(paste("Skipping trait-shared credible set", cs_name, "- pip_values length (", length(pip_values), ") != variants length (", length(variants), ")"))
        next
      }
      if (length(variant_neg_log10_p_values_list) != length(variants)) {
        warning(paste("Skipping trait-shared credible set", cs_name, "- variant_neg_log10_p_values_list length (", length(variant_neg_log10_p_values_list), ") != variants length (", length(variants), ")"))
        next
      }
      
      if (include_conversion_metadata) {
        cs_df <- data.frame(
          phenotype_id = rep(paste(phenotype_ids, collapse = ","), length(variants)),
          variant_id    = variants,
          pip           = pip_values,
          neg_log10_p_value = I(variant_neg_log10_p_values_list),
          cs_id         = rep(cs_name, length(variants)),
          cs_type       = "trait_shared",
          cos_npc       = rep(cs_cos_npc, length(variants)),
          npc_outcome   = I(variant_npc_outcome_list),
          original_cos_id = rep(NA_character_, length(variants)),
          original_cos_npc = rep(NA_real_, length(variants)),
          converted_from_shared = rep(FALSE, length(variants)),
          stringsAsFactors = FALSE
        )
      } else {
        cs_df <- data.frame(
          phenotype_id = rep(paste(phenotype_ids, collapse = ","), length(variants)),
          variant_id    = variants,
          pip           = pip_values,
          neg_log10_p_value = I(variant_neg_log10_p_values_list),
          cs_id         = rep(cs_name, length(variants)),
          cs_type       = "trait_shared",
          cos_npc       = rep(cs_cos_npc, length(variants)),
          npc_outcome   = I(variant_npc_outcome_list),
          stringsAsFactors = FALSE
        )
      }

      trait_shared_cs[[length(trait_shared_cs) + 1]] <- cs_df
    }
  }

  all_cs <- c(trait_specific_cs, trait_shared_cs)
  if (length(all_cs) > 0) {
    # Ensure all data frames have the same column structure before rbinding
    if (include_conversion_metadata) {
      # Standardize columns for conversion metadata case
      all_cs <- lapply(all_cs, function(df) {
        # Ensure all required columns exist
        required_cols <- c("phenotype_id", "variant_id", "pip", "neg_log10_p_value", 
                          "cs_id", "cs_type", "cos_npc", "npc_outcome", 
                          "original_cos_id", "original_cos_npc", "original_npc_outcome", 
                          "converted_from_shared")
        
        for (col in required_cols) {
          if (!col %in% names(df)) {
            if (col == "original_cos_id") {
              df[[col]] <- NA_character_
            } else if (col == "original_cos_npc" || col == "original_npc_outcome") {
              df[[col]] <- NA_real_
            } else if (col == "converted_from_shared") {
              df[[col]] <- FALSE
            } else {
              df[[col]] <- NA
            }
          }
        }
        
        # Reorder columns to match expected order
        df <- df[, required_cols, drop = FALSE]
        return(df)
      })
    } else {
      # Standardize columns for non-conversion metadata case
      all_cs <- lapply(all_cs, function(df) {
        # Ensure all required columns exist
        required_cols <- c("phenotype_id", "variant_id", "pip", "neg_log10_p_value", 
                          "cs_id", "cs_type", "cos_npc", "npc_outcome")
        
        for (col in required_cols) {
          if (!col %in% names(df)) {
            if (col == "cos_npc") {
              df[[col]] <- NA_real_
            } else {
              df[[col]] <- NA
            }
          }
        }
        
        # Reorder columns to match expected order
        df <- df[, required_cols, drop = FALSE]
        return(df)
      })
    }
    
    result_df <- do.call(rbind, all_cs)
    return(result_df)
  } else {
    if (include_conversion_metadata) {
      return(data.frame(
        phenotype_id = character(0),
        variant_id = character(0),
        pip = numeric(0),
        neg_log10_p_value = list(),
        cs_id = character(0),
        cs_type = character(0),
        cos_npc = numeric(0),
        npc_outcome = numeric(0),
        original_cos_id = character(0),
        original_cos_npc = numeric(0),
        converted_from_shared = logical(0),
        stringsAsFactors = FALSE
      ))
    } else {
      return(data.frame(
        phenotype_id = character(0),
        variant_id = character(0),
        pip = numeric(0),
        neg_log10_p_value = list(),
        cs_id = character(0),
        cs_type = character(0),
        cos_npc = numeric(0),
        npc_outcome = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
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