#' ============================================================================
#' Imputation Module for Epigenetic Clock Calculator
#' ============================================================================
#' 
#' Implements mean/median imputation for missing values and zero-shot
#' imputation for entirely missing probes using reference data.

#' Perform imputation on beta value matrix
#' 
#' @param betas Matrix of beta values (probes x samples)
#' @param reference_betas Named numeric vector of reference beta values
#' @param method Imputation method: "median" or "mean"
#' @param zero_shot Enable zero-shot imputation for missing probes
#' @param verbose Print progress messages
#' @return List with imputed betas and imputation statistics
#' @export
perform_imputation <- function(betas, reference_betas, method = "median",
                                zero_shot = TRUE, verbose = TRUE) {
  
  log_message("Starting imputation...", verbose)
  
  stats <- list(
    original_probes = nrow(betas),
    original_samples = ncol(betas),
    original_missing = sum(is.na(betas)),
    values_imputed = 0,
    probes_added = 0
  )
  
  # Step 1: Impute missing values within existing probes
  impute_result <- impute_missing_values(
    betas = betas,
    reference_betas = reference_betas,
    method = method,
    verbose = verbose
  )
  
  betas <- impute_result$betas
  stats$values_imputed <- impute_result$n_imputed
  
  # Step 2: Zero-shot imputation for entirely missing probes
  if (zero_shot && !is.null(reference_betas)) {
    zeroshot_result <- zero_shot_imputation(
      betas = betas,
      reference_betas = reference_betas,
      verbose = verbose
    )
    
    betas <- zeroshot_result$betas
    stats$probes_added <- zeroshot_result$n_added
  }
  
  stats$final_probes <- nrow(betas)
  stats$final_missing <- sum(is.na(betas))
  
  return(list(betas = betas, stats = stats))
}


#' Impute missing values within existing probes
#' 
#' For each probe with missing values, imputes using:
#' 1. Reference beta value (if available)
#' 2. Row median/mean across samples (fallback)
#' 
#' @param betas Matrix of beta values
#' @param reference_betas Named vector of reference values
#' @param method "median" or "mean"
#' @param verbose Print messages
#' @return List with imputed betas and count
#' @export
impute_missing_values <- function(betas, reference_betas = NULL, 
                                   method = "median", verbose = TRUE) {
  
  n_imputed <- 0
  
  # Find probes with missing values
  probe_has_missing <- rowSums(is.na(betas)) > 0
  probes_to_impute <- names(which(probe_has_missing))
  
  if (length(probes_to_impute) == 0) {
    log_message("  No missing values to impute", verbose)
    return(list(betas = betas, n_imputed = 0))
  }
  
  log_message(sprintf("  Imputing missing values in %d probes...", 
                      length(probes_to_impute)), verbose)
  
  for (probe in probes_to_impute) {
    missing_idx <- which(is.na(betas[probe, ]))
    
    if (length(missing_idx) == 0) next
    
    # Determine imputation value
    impute_value <- NULL
    
    # First choice: reference value
    if (!is.null(reference_betas) && probe %in% names(reference_betas)) {
      impute_value <- reference_betas[probe]
    }
    
    # Second choice: calculate from available data
    if (is.null(impute_value) || is.na(impute_value)) {
      available_values <- betas[probe, -missing_idx]
      
      if (length(available_values) > 0 && !all(is.na(available_values))) {
        impute_value <- switch(method,
          "median" = median(available_values, na.rm = TRUE),
          "mean" = mean(available_values, na.rm = TRUE),
          median(available_values, na.rm = TRUE)
        )
      }
    }
    
    # Third choice: use 0.5 (neutral beta value) as last resort
    if (is.null(impute_value) || is.na(impute_value)) {
      impute_value <- 0.5
    }
    
    # Apply imputation
    betas[probe, missing_idx] <- impute_value
    n_imputed <- n_imputed + length(missing_idx)
  }
  
  log_message(sprintf("  Imputed %d values", n_imputed), verbose)
  
  return(list(betas = betas, n_imputed = n_imputed))
}


#' Zero-shot imputation for entirely missing probes
#' 
#' Adds probes that are needed for clock calculations but are entirely
#' missing from the data, using reference values.
#' 
#' @param betas Matrix of beta values
#' @param reference_betas Named vector of reference values
#' @param required_probes Optional vector of probes to add
#' @param verbose Print messages
#' @return List with expanded betas and count of added probes
#' @export
zero_shot_imputation <- function(betas, reference_betas, 
                                  required_probes = NULL, verbose = TRUE) {
  
  if (is.null(reference_betas) || length(reference_betas) == 0) {
    return(list(betas = betas, n_added = 0))
  }
  
  current_probes <- rownames(betas)
  
  # Determine which probes to add
  if (!is.null(required_probes)) {
    # Add only specified required probes
    probes_to_add <- setdiff(required_probes, current_probes)
    probes_to_add <- intersect(probes_to_add, names(reference_betas))
  } else {
    # Add all reference probes not in current data
    # (This is typically not done as it would add too many probes)
    probes_to_add <- character(0)
  }
  
  if (length(probes_to_add) == 0) {
    return(list(betas = betas, n_added = 0))
  }
  
  log_message(sprintf("  Adding %d probes via zero-shot imputation...", 
                      length(probes_to_add)), verbose)
  
  # Create new rows for missing probes
  n_samples <- ncol(betas)
  
  new_rows <- matrix(
    NA_real_,
    nrow = length(probes_to_add),
    ncol = n_samples,
    dimnames = list(probes_to_add, colnames(betas))
  )
  
  # Fill with reference values (same value for all samples)
  for (probe in probes_to_add) {
    new_rows[probe, ] <- reference_betas[probe]
  }
  
  # Combine matrices
  betas <- rbind(betas, new_rows)
  
  return(list(betas = betas, n_added = length(probes_to_add)))
}


#' Get all probes required by specified clocks
#' 
#' @param clock_names Vector of clock names
#' @return Character vector of required probe IDs
#' @export
get_required_probes_for_clocks <- function(clock_names = NULL) {
  # This would typically load probe requirements from each clock package
  # For now, return an empty vector - actual requirements are loaded
  # in clock-specific modules
  
  all_probes <- character(0)
  
  # The actual implementation would query each clock package
  # for its required probes and combine them
  
  return(unique(all_probes))
}


#' KNN imputation (alternative method)
#' 
#' K-nearest neighbors imputation using sample similarity.
#' More sophisticated but computationally expensive.
#' 
#' @param betas Matrix of beta values
#' @param k Number of neighbors
#' @param verbose Print messages
#' @return Imputed beta matrix
#' @export
knn_imputation <- function(betas, k = 10, verbose = TRUE) {
  if (!requireNamespace("impute", quietly = TRUE)) {
    stop("Package 'impute' required for KNN imputation. ",
         "Install with: BiocManager::install('impute')")
  }
  
  log_message(sprintf("  Performing KNN imputation (k=%d)...", k), verbose)
  
  # impute.knn works on probes as rows
  result <- impute::impute.knn(betas, k = k, rowmax = 0.5, colmax = 0.8)
  
  return(result$data)
}


#' Random forest imputation (alternative method)
#' 
#' Uses random forest to predict missing values.
#' Most sophisticated but slowest.
#' 
#' @param betas Matrix of beta values  
#' @param verbose Print messages
#' @return Imputed beta matrix
#' @export
rf_imputation <- function(betas, verbose = TRUE) {
  if (!requireNamespace("missForest", quietly = TRUE)) {
    stop("Package 'missForest' required for RF imputation. ",
         "Install with: install.packages('missForest')")
  }
  
  log_message("  Performing random forest imputation...", verbose)
  
  # missForest works on transposed data (samples as rows)
  result <- missForest::missForest(t(betas), verbose = verbose)
  
  return(t(result$ximp))
}


#' Check imputation quality
#' 
#' @param original Original beta matrix
#' @param imputed Imputed beta matrix
#' @return Data frame with quality metrics
#' @export
check_imputation_quality <- function(original, imputed) {
  # Compare distributions
  orig_mean <- mean(original, na.rm = TRUE)
  orig_sd <- sd(original, na.rm = TRUE)
  
  imp_mean <- mean(imputed, na.rm = TRUE)
  imp_sd <- sd(imputed, na.rm = TRUE)
  
  # Find which values were imputed
  was_missing <- is.na(original)
  imputed_values <- imputed[was_missing]
  
  data.frame(
    original_mean = orig_mean,
    original_sd = orig_sd,
    imputed_mean = imp_mean,
    imputed_sd = imp_sd,
    mean_of_imputed_values = mean(imputed_values),
    sd_of_imputed_values = sd(imputed_values),
    n_imputed = sum(was_missing),
    pct_imputed = 100 * mean(was_missing)
  )
}


#' Null coalescing operator (if not already defined)
#' @keywords internal
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}
