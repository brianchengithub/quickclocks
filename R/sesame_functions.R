#' ============================================================================
#' SeSAMe Functions Wrapper for Epigenetic Clock Calculator
#' ============================================================================
#' 
#' Wrapper functions for SeSAMe's cell composition and sex inference.

#' Estimate cell composition using SeSAMe
#' 
#' @param betas Matrix of beta values (probes x samples)
#' @param reference Reference tissue type ("blood", "saliva", etc.)
#' @param verbose Print progress
#' @return Data frame with cell type proportions per sample
#' @export
estimate_cell_composition_wrapper <- function(betas, reference = "blood",
                                               verbose = TRUE) {
  
  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("Package 'sesame' required for cell composition estimation")
  }
  
  log_message("Estimating cell type composition...", verbose)
  
  # Get sample names
  sample_ids <- colnames(betas)
  n_samples <- ncol(betas)
  
  # Initialize results
  results <- NULL
  
  # Process each sample
  for (i in seq_len(n_samples)) {
    sample_betas <- betas[, i, drop = TRUE]
    names(sample_betas) <- rownames(betas)
    
    # Remove NAs
    sample_betas <- sample_betas[!is.na(sample_betas)]
    
    tryCatch({
      # Use SeSAMe's cell composition estimation
      cell_comp <- sesame::estimateCellComposition(
        sample_betas,
        reference = reference
      )
      
      if (!is.null(cell_comp)) {
        cell_df <- as.data.frame(t(cell_comp))
        cell_df$sample_id <- sample_ids[i]
        
        if (is.null(results)) {
          results <- cell_df
        } else {
          results <- rbind(results, cell_df)
        }
      }
      
    }, error = function(e) {
      if (verbose) {
        log_message(sprintf("  Warning: Cell composition failed for %s: %s",
                            sample_ids[i], e$message), verbose, "warning")
      }
    })
    
    if (verbose && i %% 50 == 0) {
      log_message(sprintf("  Processed %d / %d samples", i, n_samples), verbose)
    }
  }
  
  if (is.null(results)) {
    warning("Cell composition estimation failed for all samples")
    return(NULL)
  }
  
  # Reorder columns
  cell_types <- setdiff(colnames(results), "sample_id")
  results <- results[, c("sample_id", cell_types), drop = FALSE]
  
  rownames(results) <- NULL
  
  return(results)
}


#' Infer sex using SeSAMe
#' 
#' @param betas Matrix of beta values (probes x samples)
#' @param verbose Print progress
#' @return Data frame with inferred sex per sample
#' @export
infer_sex_wrapper <- function(betas, verbose = TRUE) {
  
  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("Package 'sesame' required for sex inference")
  }
  
  log_message("Inferring sex from methylation data...", verbose)
  
  # Get sample names
  sample_ids <- colnames(betas)
  n_samples <- ncol(betas)
  
  # Initialize results
  results <- data.frame(
    sample_id = sample_ids,
    inferred_sex = character(n_samples),
    sex_score = numeric(n_samples),
    stringsAsFactors = FALSE
  )
  
  # Process each sample
  for (i in seq_len(n_samples)) {
    sample_betas <- betas[, i, drop = TRUE]
    names(sample_betas) <- rownames(betas)
    
    # Remove NAs
    sample_betas <- sample_betas[!is.na(sample_betas)]
    
    tryCatch({
      # Use SeSAMe's sex inference
      sex_result <- sesame::inferSex(sample_betas)
      
      # SeSAMe returns a list or character
      if (is.list(sex_result)) {
        results$inferred_sex[i] <- sex_result$sex
        results$sex_score[i] <- sex_result$score %||% NA_real_
      } else {
        results$inferred_sex[i] <- as.character(sex_result)
        results$sex_score[i] <- NA_real_
      }
      
    }, error = function(e) {
      if (verbose) {
        log_message(sprintf("  Warning: Sex inference failed for %s: %s",
                            sample_ids[i], e$message), verbose, "warning")
      }
      results$inferred_sex[i] <- NA_character_
      results$sex_score[i] <- NA_real_
    })
  }
  
  # Standardize sex labels
  results$inferred_sex <- toupper(results$inferred_sex)
  results$inferred_sex[results$inferred_sex == "MALE"] <- "M"
  results$inferred_sex[results$inferred_sex == "FEMALE"] <- "F"
  
  log_message(sprintf("  Inferred: %d Male, %d Female, %d Unknown",
                      sum(results$inferred_sex == "M", na.rm = TRUE),
                      sum(results$inferred_sex == "F", na.rm = TRUE),
                      sum(is.na(results$inferred_sex))), verbose)
  
  return(results)
}


#' Calculate age using SeSAMe's built-in clocks
#' 
#' @param betas Matrix of beta values
#' @param clock_name Name of clock ("Horvath", "Hannum", etc.)
#' @param verbose Print progress
#' @return Named numeric vector of age estimates
#' @export
calculate_sesame_age <- function(betas, clock_name = "Horvath", verbose = TRUE) {
  
  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("Package 'sesame' required")
  }
  
  log_message(sprintf("Calculating %s age using SeSAMe...", clock_name), verbose)
  
  sample_ids <- colnames(betas)
  n_samples <- ncol(betas)
  
  ages <- numeric(n_samples)
  names(ages) <- sample_ids
  
  for (i in seq_len(n_samples)) {
    sample_betas <- betas[, i, drop = TRUE]
    names(sample_betas) <- rownames(betas)
    sample_betas <- sample_betas[!is.na(sample_betas)]
    
    tryCatch({
      # Try to use predictAge function (if available in sesame version)
      age <- sesame::predictAge(sample_betas, clock = clock_name)
      ages[i] <- age
    }, error = function(e) {
      ages[i] <- NA_real_
    })
  }
  
  return(ages)
}


#' Get SeSAMe reference cell types
#' 
#' @param reference Reference type
#' @return Character vector of cell types
#' @export
get_sesame_cell_types <- function(reference = "blood") {
  
  if (!requireNamespace("sesame", quietly = TRUE)) {
    return(NULL)
  }
  
  # Blood cell types
  blood_types <- c(
    "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"
  )
  
  # Extended blood types
  blood_extended <- c(
    "CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"
  )
  
  switch(reference,
    "blood" = blood_types,
    "blood_extended" = blood_extended,
    blood_types  # default
  )
}


#' Null coalescing operator
#' @keywords internal
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}
