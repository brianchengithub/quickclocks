#' ============================================================================
#' Clock Calculator - Main Orchestration Module
#' ============================================================================
#' 
#' This module orchestrates the calculation of all epigenetic clocks
#' 

#' Calculate all requested epigenetic clocks
#' 
#' @param betas Preprocessed beta matrix (probes x samples)
#' @param config Configuration list
#' @param reference_betas Named vector of reference beta values for imputation
#' @param verbose Print progress messages
#' @return List containing all clock results
#' @export
calculate_all_clocks <- function(betas, config, reference_betas = NULL, verbose = TRUE) {
  
  log_message("========================================", verbose)
  log_message("EPIGENETIC CLOCK CALCULATION", verbose)
  log_message("========================================", verbose)
  
  results <- list()
  
  # Get probe information
  available_probes <- rownames(betas)
  n_samples <- ncol(betas)
  
  log_message(sprintf("Input: %d probes x %d samples", 
                      length(available_probes), n_samples), verbose)
  
  # =========================================================================
  # Check clock availability
  # =========================================================================
  log_message("\n--- Checking Clock Availability ---", verbose)
  
  availability <- check_all_clock_availability(available_probes, verbose = verbose)
  results$availability <- availability
  
  # =========================================================================
  # SeSAMe-based calculations
  # =========================================================================
  if (isTRUE(config$calculate_sesame)) {
    log_message("\n--- SeSAMe Calculations ---", verbose)
    
    # Cell composition
    tryCatch({
      log_message("Estimating cell composition...", verbose)
      results$cell_composition <- estimate_cell_composition_wrapper(betas, verbose)
    }, error = function(e) {
      log_message(sprintf("Cell composition failed: %s", e$message), verbose, "warning")
      results$cell_composition <- NULL
    })
    
    # Sex inference
    tryCatch({
      log_message("Inferring sex...", verbose)
      results$sex <- infer_sex_wrapper(betas, verbose)
    }, error = function(e) {
      log_message(sprintf("Sex inference failed: %s", e$message), verbose, "warning")
      results$sex <- NULL
    })
    
    # SeSAMe age clocks (if available)
    tryCatch({
      log_message("Calculating SeSAMe age clocks...", verbose)
      results$sesame_age <- calculate_sesame_age(betas, verbose)
    }, error = function(e) {
      log_message(sprintf("SeSAMe age clocks failed: %s", e$message), verbose, "warning")
      results$sesame_age <- NULL
    })
  }
  
  # =========================================================================
  # DunedinPACE
  # =========================================================================
  if (isTRUE(config$calculate_dunedin_pace)) {
    log_message("\n--- DunedinPACE ---", verbose)
    
    tryCatch({
      results$dunedin_pace <- calculate_dunedin_pace(betas, verbose = verbose)
    }, error = function(e) {
      log_message(sprintf("DunedinPACE failed: %s", e$message), verbose, "warning")
      results$dunedin_pace <- NULL
    })
  }
  
  # =========================================================================
  # PC Clocks (Morgan Levine Lab)
  # =========================================================================
  if (isTRUE(config$calculate_pc_clocks)) {
    log_message("\n--- PC Clocks (Levine Lab) ---", verbose)
    
    tryCatch({
      results$pc_clocks <- calculate_pc_clocks(betas, verbose = verbose)
    }, error = function(e) {
      log_message(sprintf("PC-Clocks failed: %s", e$message), verbose, "warning")
      results$pc_clocks <- NULL
    })
  }
  
  # =========================================================================
  # EpiTOC2 / Mitotic Clocks
  # =========================================================================
  if (isTRUE(config$calculate_epitoc2)) {
    log_message("\n--- Mitotic Clocks (epiTOC2) ---", verbose)
    
    tryCatch({
      results$mitotic_clocks <- calculate_all_mitotic_clocks(betas, verbose = verbose)
    }, error = function(e) {
      log_message(sprintf("Mitotic clocks failed: %s", e$message), verbose, "warning")
      results$mitotic_clocks <- NULL
    })
  }
  
  # =========================================================================
  # methylCIPHER Clocks
  # =========================================================================
  if (isTRUE(config$calculate_methylcipher)) {
    log_message("\n--- methylCIPHER Clocks ---", verbose)
    
    # Get list of clocks to calculate
    methylcipher_clocks <- config$methylcipher_clocks
    if (is.null(methylcipher_clocks) || length(methylcipher_clocks) == 0) {
      methylcipher_clocks <- get_methylcipher_clock_list()
    }
    
    tryCatch({
      results$methylcipher <- calculate_methylcipher_clocks(
        betas, 
        clocks = methylcipher_clocks,
        verbose = verbose
      )
    }, error = function(e) {
      log_message(sprintf("methylCIPHER failed: %s", e$message), verbose, "warning")
      results$methylcipher <- NULL
    })
  }
  
  # =========================================================================
  # Compile Summary
  # =========================================================================
  log_message("\n========================================", verbose)
  log_message("CALCULATION SUMMARY", verbose)
  log_message("========================================", verbose)
  
  results$summary <- compile_clock_summary(results, verbose)
  
  return(results)
}


#' Compile summary of all clock calculations
#' 
#' @param results List of all results
#' @param verbose Print messages
#' @return Summary data frame
#' @keywords internal
compile_clock_summary <- function(results, verbose = TRUE) {
  
  summary_rows <- list()
  
  # SeSAMe results
  if (!is.null(results$cell_composition)) {
    summary_rows$cell_composition <- list(
      category = "QC",
      clock = "Cell Composition",
      status = "Success",
      n_samples = nrow(results$cell_composition)
    )
    log_message("  [OK] Cell Composition", verbose)
  }
  
  if (!is.null(results$sex)) {
    summary_rows$sex <- list(
      category = "QC",
      clock = "Sex Inference",
      status = "Success",
      n_samples = nrow(results$sex)
    )
    log_message("  [OK] Sex Inference", verbose)
  }
  
  # DunedinPACE
  if (!is.null(results$dunedin_pace)) {
    summary_rows$dunedin_pace <- list(
      category = "Pace of Aging",
      clock = "DunedinPACE",
      status = "Success",
      n_samples = nrow(results$dunedin_pace)
    )
    log_message("  [OK] DunedinPACE", verbose)
  }
  
  # PC Clocks
  if (!is.null(results$pc_clocks)) {
    pc_clock_names <- setdiff(colnames(results$pc_clocks), "sample_id")
    for (clock in pc_clock_names) {
      summary_rows[[paste0("pc_", clock)]] <- list(
        category = "PC Clocks",
        clock = clock,
        status = "Success",
        n_samples = nrow(results$pc_clocks)
      )
    }
    log_message(sprintf("  [OK] PC Clocks (%d clocks)", length(pc_clock_names)), verbose)
  }
  
  # Mitotic clocks
  if (!is.null(results$mitotic_clocks)) {
    mitotic_names <- setdiff(colnames(results$mitotic_clocks), "sample_id")
    for (clock in mitotic_names) {
      summary_rows[[paste0("mitotic_", clock)]] <- list(
        category = "Mitotic",
        clock = clock,
        status = "Success",
        n_samples = nrow(results$mitotic_clocks)
      )
    }
    log_message(sprintf("  [OK] Mitotic Clocks (%d clocks)", length(mitotic_names)), verbose)
  }
  
  # methylCIPHER
  if (!is.null(results$methylcipher)) {
    mc_names <- setdiff(colnames(results$methylcipher), "sample_id")
    for (clock in mc_names) {
      category <- get_methylcipher_category(clock)
      summary_rows[[paste0("mc_", clock)]] <- list(
        category = category,
        clock = clock,
        status = "Success",
        n_samples = nrow(results$methylcipher)
      )
    }
    log_message(sprintf("  [OK] methylCIPHER (%d clocks)", length(mc_names)), verbose)
  }
  
  # Convert to data frame
  if (length(summary_rows) > 0) {
    summary_df <- do.call(rbind, lapply(summary_rows, as.data.frame, 
                                        stringsAsFactors = FALSE))
    rownames(summary_df) <- NULL
    return(summary_df)
  }
  
  return(data.frame(
    category = character(),
    clock = character(),
    status = character(),
    n_samples = integer(),
    stringsAsFactors = FALSE
  ))
}


#' Get methylCIPHER clock category
#' @keywords internal
get_methylcipher_category <- function(clock_name) {
  
  age_clocks <- c("Horvath1", "Horvath2", "Hannum", "PhenoAge", 
                  "GrimAge1", "GrimAge2", "Zhang", "Zhang2019", 
                  "Lin", "DNAmTL")
  
  specialized <- c("AdaptAge", "CausAge", "DamAge", "HypoClock",
                   "MiAge", "SystemsAge", "RetroelementAge450K")
  
  lifestyle <- c("Alcohol", "BMI", "Smoking")
  
  organ <- c("Blood", "Brain", "Heart", "Kidney", 
             "Liver", "Lung", "MusculoSkeletal")
  
  biological <- c("Inflammation", "Hormone", "Immune", "Metabolic")
  
  if (clock_name %in% age_clocks) return("Age")
  if (clock_name %in% specialized) return("Specialized")
  if (clock_name %in% lifestyle) return("Lifestyle")
  if (clock_name %in% organ) return("Organ/System")
  if (clock_name %in% biological) return("Biological")
  
  return("Other")
}


#' Combine all clock results into a single data frame
#' 
#' @param results List of clock results
#' @return Combined data frame with all clock values per sample
#' @export
combine_clock_results <- function(results) {
  
  # Start with sample IDs
  sample_ids <- NULL
  
  # Find sample IDs from first available result
  for (name in names(results)) {
    if (!is.null(results[[name]]) && is.data.frame(results[[name]])) {
      if ("sample_id" %in% colnames(results[[name]])) {
        sample_ids <- results[[name]]$sample_id
        break
      }
    }
  }
  
  if (is.null(sample_ids)) {
    warning("No sample IDs found in results")
    return(NULL)
  }
  
  combined <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
  
  # Add each result set
  result_names <- c("dunedin_pace", "pc_clocks", "mitotic_clocks", "methylcipher",
                    "sesame_age", "cell_composition", "sex")
  
  for (name in result_names) {
    if (!is.null(results[[name]]) && is.data.frame(results[[name]])) {
      df <- results[[name]]
      
      # Get value columns (not sample_id)
      value_cols <- setdiff(colnames(df), "sample_id")
      
      if (length(value_cols) > 0) {
        # Merge by sample_id
        combined <- merge(combined, df[, c("sample_id", value_cols), drop = FALSE],
                         by = "sample_id", all.x = TRUE)
      }
    }
  }
  
  return(combined)
}


#' Run clock calculation pipeline
#' 
#' @param betas Beta matrix
#' @param config Configuration
#' @param reference_betas Reference betas for imputation
#' @param output_dir Output directory
#' @param verbose Verbose output
#' @return List of results
#' @export
run_clock_pipeline <- function(betas, config, reference_betas = NULL, 
                               output_dir = ".", verbose = TRUE) {
  
  start_time <- Sys.time()
  
  # Validate inputs
  if (!is.matrix(betas) && !is.data.frame(betas)) {
    stop("betas must be a matrix or data frame")
  }
  
  if (is.data.frame(betas)) {
    betas <- as.matrix(betas)
  }
  
  # Ensure numeric
  if (!is.numeric(betas)) {
    stop("betas must contain numeric values")
  }
  
  log_message("Starting clock calculation pipeline...", verbose)
  log_message(sprintf("  Samples: %d", ncol(betas)), verbose)
  log_message(sprintf("  Probes: %d", nrow(betas)), verbose)
  
  # Calculate all clocks
  results <- calculate_all_clocks(betas, config, reference_betas, verbose)
  
  # Combine results
  results$combined <- combine_clock_results(results)
  
  # Timing
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  
  log_message(sprintf("\nTotal time: %.2f minutes", as.numeric(elapsed)), verbose)
  
  return(results)
}


#' Quick calculation of specific clocks only
#' 
#' @param betas Beta matrix
#' @param clocks Character vector of clock names
#' @param verbose Verbose output
#' @return Data frame with requested clock values
#' @export
calculate_specific_clocks <- function(betas, clocks, verbose = TRUE) {
  
  # Categorize requested clocks
  pc_clocks <- c("PCHorvath1", "PCHorvath2", "PCHannum", "PCPhenoAge", "PCGrimAge")
  mitotic_clocks <- c("epiTOC2", "HypoClock", "MiAge")
  
  results <- data.frame(sample_id = colnames(betas), stringsAsFactors = FALSE)
  
  # DunedinPACE
  if ("DunedinPACE" %in% clocks) {
    tryCatch({
      pace <- calculate_dunedin_pace(betas, verbose = verbose)
      if (!is.null(pace) && "DunedinPACE" %in% colnames(pace)) {
        results$DunedinPACE <- pace$DunedinPACE
      }
    }, error = function(e) {
      log_message(sprintf("DunedinPACE failed: %s", e$message), verbose, "warning")
    })
  }
  
  # PC clocks
  requested_pc <- intersect(clocks, pc_clocks)
  if (length(requested_pc) > 0) {
    tryCatch({
      pc_results <- calculate_pc_clocks(betas, verbose = verbose)
      if (!is.null(pc_results)) {
        for (clock in requested_pc) {
          if (clock %in% colnames(pc_results)) {
            results[[clock]] <- pc_results[[clock]]
          }
        }
      }
    }, error = function(e) {
      log_message(sprintf("PC clocks failed: %s", e$message), verbose, "warning")
    })
  }
  
  # Mitotic clocks
  requested_mitotic <- intersect(clocks, mitotic_clocks)
  if (length(requested_mitotic) > 0) {
    tryCatch({
      mitotic_results <- calculate_all_mitotic_clocks(betas, verbose = verbose)
      if (!is.null(mitotic_results)) {
        for (clock in requested_mitotic) {
          if (clock %in% colnames(mitotic_results)) {
            results[[clock]] <- mitotic_results[[clock]]
          }
        }
      }
    }, error = function(e) {
      log_message(sprintf("Mitotic clocks failed: %s", e$message), verbose, "warning")
    })
  }
  
  # methylCIPHER clocks
  mc_clocks <- get_methylcipher_clock_list()
  requested_mc <- intersect(clocks, mc_clocks)
  if (length(requested_mc) > 0) {
    tryCatch({
      mc_results <- calculate_methylcipher_clocks(betas, clocks = requested_mc, 
                                                  verbose = verbose)
      if (!is.null(mc_results)) {
        for (clock in requested_mc) {
          if (clock %in% colnames(mc_results)) {
            results[[clock]] <- mc_results[[clock]]
          }
        }
      }
    }, error = function(e) {
      log_message(sprintf("methylCIPHER clocks failed: %s", e$message), verbose, "warning")
    })
  }
  
  return(results)
}
