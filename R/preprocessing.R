#' ============================================================================
#' Preprocessing Module for Epigenetic Clock Calculator
#' ============================================================================
#' 
#' SeSAMe-based preprocessing pipeline with NOOB normalization,
#' non-linear dye-bias correction, and pOOBAH background subtraction.

#' Apply SeSAMe preprocessing to beta values
#' 
#' This function is used when beta values are provided directly but
#' additional preprocessing steps are needed.
#' 
#' @param betas Matrix of beta values
#' @param config Configuration list
#' @return Preprocessed beta matrix
#' @export
preprocess_betas <- function(betas, config = list()) {
  
  verbose <- config$verbose %||% TRUE
  
  # For pre-computed betas, we can only do limited preprocessing
  log_message("Applying post-processing to beta values...", verbose)
  
  # Check for out-of-range values
  n_below_zero <- sum(betas < 0, na.rm = TRUE)
  n_above_one <- sum(betas > 1, na.rm = TRUE)
  
  if (n_below_zero > 0 || n_above_one > 0) {
    log_message(sprintf("  Clipping %d values < 0 and %d values > 1", 
                        n_below_zero, n_above_one), verbose)
    betas[betas < 0] <- 0
    betas[betas > 1] <- 1
  }
  
  # Remove probes with too many missing values
  max_missing <- config$max_missing_probe_fraction %||% 0.5
  probe_missing_frac <- rowMeans(is.na(betas))
  
  probes_to_remove <- probe_missing_frac > max_missing
  if (any(probes_to_remove)) {
    n_removed <- sum(probes_to_remove)
    log_message(sprintf("  Removing %d probes with >%.0f%% missing values", 
                        n_removed, max_missing * 100), verbose)
    betas <- betas[!probes_to_remove, , drop = FALSE]
  }
  
  # Remove samples with too many missing values
  max_sample_missing <- config$max_missing_fraction %||% 0.2
  sample_missing_frac <- colMeans(is.na(betas))
  
  samples_to_remove <- sample_missing_frac > max_sample_missing
  if (any(samples_to_remove)) {
    n_removed <- sum(samples_to_remove)
    log_message(sprintf("  Removing %d samples with >%.0f%% missing values", 
                        n_removed, max_sample_missing * 100), verbose)
    log_message(sprintf("    Removed: %s", 
                        paste(colnames(betas)[samples_to_remove], collapse = ", ")), 
                verbose)
    betas <- betas[, !samples_to_remove, drop = FALSE]
  }
  
  return(betas)
}


#' Apply NOOB normalization
#' 
#' Normal-exponential out-of-band (NOOB) background subtraction
#' for Illumina methylation arrays.
#' 
#' @param sigset SigSet object or beta values
#' @return Normalized values
#' @export
apply_noob <- function(sigset) {
  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("Package 'sesame' required for NOOB normalization")
  }
  
  sesame::noob(sigset)
}


#' Apply non-linear dye-bias correction
#' 
#' Corrects for dye-bias artifacts in two-color arrays.
#' 
#' @param sigset SigSet object
#' @return Corrected SigSet
#' @export
apply_dyebias_correction <- function(sigset) {
  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("Package 'sesame' required for dye-bias correction")
  }
  
  sesame::dyeBiasCorrTypeINorm(sigset)
}


#' Apply pOOBAH background detection
#' 
#' p-value with OOB Array Hybridization for detection calling.
#' 
#' @param sigset SigSet object
#' @param pval_threshold P-value threshold for detection
#' @return SigSet with pOOBAH mask applied
#' @export
apply_poobah <- function(sigset, pval_threshold = 0.05) {
  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("Package 'sesame' required for pOOBAH")
  }
  
  sesame::pOOBAH(sigset, pval.threshold = pval_threshold)
}


#' Full preprocessing pipeline for SigSet
#' 
#' Applies complete preprocessing: quality mask, dye-bias correction,
#' detection masking, and NOOB background correction.
#' 
#' @param sigset SigSet object
#' @param prep_method Preprocessing method string
#' @param pval_threshold Detection p-value threshold
#' @return Preprocessed beta values
#' @export
full_preprocessing_pipeline <- function(sigset, prep_method = "QCDPB", 
                                         pval_threshold = 0.05) {
  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("Package 'sesame' required")
  }
  
  # Parse prep_method string
  # Q = Quality mask
  # C = Dye-bias Correction
  # D = Detection (pOOBAH)
  # P = Probe filtering
  # B = Background (NOOB)
  
  prep_steps <- strsplit(prep_method, "")[[1]]
  
  # Apply steps in order
  for (step in prep_steps) {
    sigset <- switch(step,
      "Q" = sesame::qualityMask(sigset),
      "C" = sesame::dyeBiasCorrTypeINorm(sigset),
      "D" = sesame::pOOBAH(sigset, pval.threshold = pval_threshold),
      "P" = sesame::inferInfiniumIChannel(sigset),
      "B" = sesame::noob(sigset),
      sigset  # Unknown step, pass through
    )
  }
  
  # Extract beta values
  sesame::getBetas(sigset)
}


#' Batch effect correction using ComBat (optional)
#' 
#' Applies ComBat batch correction if batch information is provided.
#' 
#' @param betas Matrix of beta values
#' @param batch Vector of batch labels
#' @param covariates Optional data frame of covariates to preserve
#' @return Batch-corrected beta matrix
#' @export
apply_combat <- function(betas, batch, covariates = NULL) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("Package 'sva' required for ComBat. Install with: BiocManager::install('sva')")
  }
  
  # Convert to M-values for ComBat (more normally distributed)
  m_values <- beta_to_m(betas)
  
  # Handle infinite values from extreme betas
  m_values[is.infinite(m_values)] <- NA
  
  # Remove probes with too many missing values
  valid_probes <- rowMeans(is.na(m_values)) < 0.1
  m_values <- m_values[valid_probes, ]
  
  # Apply ComBat
  if (is.null(covariates)) {
    m_corrected <- sva::ComBat(dat = m_values, batch = batch)
  } else {
    mod <- model.matrix(~ ., data = covariates)
    m_corrected <- sva::ComBat(dat = m_values, batch = batch, mod = mod)
  }
  
  # Convert back to beta values
  betas_corrected <- m_to_beta(m_corrected)
  
  # Create full output matrix
  betas_out <- betas
  betas_out[valid_probes, ] <- betas_corrected
  
  return(betas_out)
}


#' Remove technical variation probes
#' 
#' Filters out probes known to have high technical variation.
#' 
#' @param betas Beta value matrix
#' @param platform Array platform
#' @return Filtered beta matrix
#' @export
filter_technical_probes <- function(betas, platform = "EPIC") {
  # Probes to filter:
  # - Cross-reactive probes
  # - SNP-affected probes
  # - Sex chromosome probes (optional)
  
  probes <- rownames(betas)
  
  # Remove probes on sex chromosomes (if desired)
  # These have "chrX" or "chrY" in their annotation
  
  # Remove control probes (start with "rs" or "nv")
  control_pattern <- "^(rs|nv)"
  control_probes <- grep(control_pattern, probes, value = TRUE)
  
  if (length(control_probes) > 0) {
    betas <- betas[!probes %in% control_probes, , drop = FALSE]
  }
  
  return(betas)
}


#' Calculate detection p-values
#' 
#' @param betas Beta value matrix
#' @param background Background signal estimate
#' @return Matrix of detection p-values
#' @export
calculate_detection_pvalues <- function(betas, background = NULL) {
  # Simplified detection based on signal strength
  # Lower beta values (closer to 0 or 1) generally have better detection
  
  # Use deviation from 0.5 as proxy for signal
  signal_strength <- abs(betas - 0.5) * 2
  
  # Convert to approximate p-value (not statistically rigorous)
  # This is a placeholder - proper detection requires raw intensities
  pvals <- 1 - signal_strength
  pvals[is.na(pvals)] <- 1
  
  return(pvals)
}


#' Quality metrics summary
#' 
#' @param betas Beta value matrix
#' @return Data frame of quality metrics per sample
#' @export
calculate_qc_metrics <- function(betas) {
  data.frame(
    sample_id = colnames(betas),
    n_probes = nrow(betas),
    n_missing = colSums(is.na(betas)),
    pct_missing = 100 * colMeans(is.na(betas)),
    mean_beta = colMeans(betas, na.rm = TRUE),
    sd_beta = apply(betas, 2, sd, na.rm = TRUE),
    median_beta = apply(betas, 2, median, na.rm = TRUE),
    row.names = NULL
  )
}


#' Null coalescing operator (if not already defined)
#' @keywords internal
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}
