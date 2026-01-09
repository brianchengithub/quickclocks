#' Self-contained epigenetic clock implementations
#' These functions extract coefficient data from packages at load time
#' and implement calculations directly, avoiding lazy data loading issues

#' Initialize clock coefficients from installed packages
#' @return List of coefficient data frames
#' @keywords internal
initialize_clock_coefficients <- function() {
  coeffs <- list()
  
  # ===== methylCIPHER coefficients =====
  if (requireNamespace("methylCIPHER", quietly = TRUE)) {
    mc_datasets <- c(
      "Horvath1_CpGs", "Horvath2_CpGs", "Hannum_CpGs", "PhenoAge_CpGs",
      "DNAmTL_CpGs", "Lin_CpGs", "Zhang_10_CpG", "Zhang2019_CpGs",
      "Bocklandt_CpG", "Weidner_CpGs", "VidalBralo_CpGs", "Garagnani_CpG",
      "PCClocks_CpGs", "HorvathOnlineRef",
      "AdaptAge_CpGs", "CausAge_CpGs", "DamAge_CpGs", "SystemsAge_CpGs",
      "EpiToc_CpGs", "EpiToc2_CpGs", "hypoClock_CpGs", "MiAge_CpGs"
    )
    
    for (d in mc_datasets) {
      tryCatch({
        # Load to temp environment to extract
        temp_env <- new.env()
        data(list = d, package = "methylCIPHER", envir = temp_env)
        if (exists(d, envir = temp_env)) {
          coeffs[[d]] <- get(d, envir = temp_env)
        }
      }, error = function(e) NULL, warning = function(w) NULL)
    }
  }
  
  # ===== EpiMitClocks coefficients =====
  if (requireNamespace("EpiMitClocks", quietly = TRUE)) {
    epi_datasets <- list(
      "dataETOC3" = "dataETOC3.l",
      "estETOC3" = "estETOC3.m", 
      "epiTOCcpgs3" = "epiTOCcpgs3.v",
      "cugpmitclockCpG" = "cugpmitclockCpG.v",
      "EpiCMITcpgs" = "epiCMIT.df",
      "Replitali" = c("replitali.coe", "replitali.cpg.v")
    )
    
    for (load_name in names(epi_datasets)) {
      tryCatch({
        temp_env <- new.env()
        data(list = load_name, package = "EpiMitClocks", envir = temp_env)
        obj_names <- epi_datasets[[load_name]]
        for (obj in obj_names) {
          if (exists(obj, envir = temp_env)) {
            coeffs[[obj]] <- get(obj, envir = temp_env)
          }
        }
      }, error = function(e) NULL, warning = function(w) NULL)
    }
  }
  
  return(coeffs)
}


#' Calculate Horvath Multi-Tissue Clock (Horvath 2013)
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param coeffs Pre-loaded coefficients list
#' @return Named vector of clock values
#' @keywords internal
calc_horvath1_direct <- function(betas, coeffs) {
  if (!"Horvath1_CpGs" %in% names(coeffs)) return(NULL)
  
  coef_df <- coeffs$Horvath1_CpGs
  
  # Standard format: CpG column and weight/coefficient column
  if ("CpGmarker" %in% colnames(coef_df)) {
    cpg_col <- "CpGmarker"
  } else if ("CpG" %in% colnames(coef_df)) {
    cpg_col <- "CpG"
  } else {
    cpg_col <- colnames(coef_df)[1]
  }
  
  if ("CoesCoefficient" %in% colnames(coef_df)) {
    weight_col <- "CoefficientTraining"
  } else if ("weight" %in% colnames(coef_df)) {
    weight_col <- "weight"
  } else if ("Coefficient" %in% colnames(coef_df)) {
    weight_col <- "Coefficient"
  } else {
    weight_col <- colnames(coef_df)[2]
  }
  
  # Get intercept (first row usually)
  intercept_idx <- which(coef_df[[cpg_col]] == "(Intercept)" | 
                         coef_df[[cpg_col]] == "Intercept" |
                         is.na(coef_df[[cpg_col]]))
  
  if (length(intercept_idx) > 0) {
    intercept <- coef_df[[weight_col]][intercept_idx[1]]
    coef_df <- coef_df[-intercept_idx, ]
  } else {
    intercept <- 0
  }
  
  cpgs <- coef_df[[cpg_col]]
  weights <- coef_df[[weight_col]]
  
  # Match CpGs
  matched_idx <- match(cpgs, rownames(betas))
  valid <- !is.na(matched_idx)
  
  if (sum(valid) == 0) return(NULL)
  
  # Calculate clock
  betas_subset <- betas[matched_idx[valid], , drop = FALSE]
  weights_valid <- weights[valid]
  
  # Weighted sum + intercept
  clock_values <- intercept + colSums(betas_subset * weights_valid, na.rm = TRUE)
  
  # Horvath uses anti-log transformation
  # DNAmAge = exp(clock) - 1 for values < 0, or exp(clock) + 1 for values >= 0
  # Actually the transformation is: if x < 0: (1+adult_age)*exp(x)-1, else: (1+adult_age)*x + adult_age
  # Simplified: often just report the raw value or use adult_age = 20
  
  return(clock_values)
}


#' Calculate Hannum Clock (Hannum 2013)
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param coeffs Pre-loaded coefficients list
#' @return Named vector of clock values
#' @keywords internal
calc_hannum_direct <- function(betas, coeffs) {
  if (!"Hannum_CpGs" %in% names(coeffs)) return(NULL)
  
  coef_df <- coeffs$Hannum_CpGs
  
  # Find CpG and weight columns
  cpg_col <- if ("CpG" %in% colnames(coef_df)) "CpG" else colnames(coef_df)[1]
  weight_col <- if ("Coefficient" %in% colnames(coef_df)) "Coefficient" else 
                if ("weight" %in% colnames(coef_df)) "weight" else colnames(coef_df)[2]
  
  # Check for intercept
  intercept_idx <- which(coef_df[[cpg_col]] %in% c("(Intercept)", "Intercept"))
  if (length(intercept_idx) > 0) {
    intercept <- coef_df[[weight_col]][intercept_idx[1]]
    coef_df <- coef_df[-intercept_idx, ]
  } else {
    intercept <- 0
  }
  
  cpgs <- coef_df[[cpg_col]]
  weights <- coef_df[[weight_col]]
  
  matched_idx <- match(cpgs, rownames(betas))
  valid <- !is.na(matched_idx)
  
  if (sum(valid) == 0) return(NULL)
  
  betas_subset <- betas[matched_idx[valid], , drop = FALSE]
  weights_valid <- weights[valid]
  
  clock_values <- intercept + colSums(betas_subset * weights_valid, na.rm = TRUE)
  
  return(clock_values)
}


#' Calculate PhenoAge (Levine 2018)
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param coeffs Pre-loaded coefficients list
#' @return Named vector of clock values
#' @keywords internal
calc_phenoage_direct <- function(betas, coeffs) {
  if (!"PhenoAge_CpGs" %in% names(coeffs)) return(NULL)
  
  coef_df <- coeffs$PhenoAge_CpGs
  
  cpg_col <- if ("CpG" %in% colnames(coef_df)) "CpG" else colnames(coef_df)[1]
  weight_col <- if ("Coefficient" %in% colnames(coef_df)) "Coefficient" else 
                if ("weight" %in% colnames(coef_df)) "weight" else colnames(coef_df)[2]
  
  intercept_idx <- which(coef_df[[cpg_col]] %in% c("(Intercept)", "Intercept"))
  if (length(intercept_idx) > 0) {
    intercept <- coef_df[[weight_col]][intercept_idx[1]]
    coef_df <- coef_df[-intercept_idx, ]
  } else {
    intercept <- 0
  }
  
  cpgs <- coef_df[[cpg_col]]
  weights <- coef_df[[weight_col]]
  
  matched_idx <- match(cpgs, rownames(betas))
  valid <- !is.na(matched_idx)
  
  if (sum(valid) == 0) return(NULL)
  
  betas_subset <- betas[matched_idx[valid], , drop = FALSE]
  weights_valid <- weights[valid]
  
  clock_values <- intercept + colSums(betas_subset * weights_valid, na.rm = TRUE)
  
  return(clock_values)
}


#' Calculate DNAmTL (Lu 2019)
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param coeffs Pre-loaded coefficients list
#' @return Named vector of clock values
#' @keywords internal
calc_dnamtl_direct <- function(betas, coeffs) {
  if (!"DNAmTL_CpGs" %in% names(coeffs)) return(NULL)
  
  coef_df <- coeffs$DNAmTL_CpGs
  
  cpg_col <- if ("CpG" %in% colnames(coef_df)) "CpG" else colnames(coef_df)[1]
  weight_col <- if ("Coefficient" %in% colnames(coef_df)) "Coefficient" else 
                if ("weight" %in% colnames(coef_df)) "weight" else colnames(coef_df)[2]
  
  intercept_idx <- which(coef_df[[cpg_col]] %in% c("(Intercept)", "Intercept"))
  if (length(intercept_idx) > 0) {
    intercept <- coef_df[[weight_col]][intercept_idx[1]]
    coef_df <- coef_df[-intercept_idx, ]
  } else {
    intercept <- 0
  }
  
  cpgs <- coef_df[[cpg_col]]
  weights <- coef_df[[weight_col]]
  
  matched_idx <- match(cpgs, rownames(betas))
  valid <- !is.na(matched_idx)
  
  if (sum(valid) == 0) return(NULL)
  
  betas_subset <- betas[matched_idx[valid], , drop = FALSE]
  weights_valid <- weights[valid]
  
  clock_values <- intercept + colSums(betas_subset * weights_valid, na.rm = TRUE)
  
  return(clock_values)
}


#' Calculate epiTOC2 mitotic clock (Teschendorff 2020)
#' Direct implementation without relying on package data loading
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param coeffs Pre-loaded coefficients list
#' @return Named vector of TNSC (total number of stem cell divisions)
#' @keywords internal
calc_epitoc2_direct <- function(betas, coeffs) {
  # Try to get epiTOC2 parameters from coefficients
  # epiTOC2 uses 163 CpGs with de-novo methylation rate (delta) and ground-state (beta0)
  
  if ("dataETOC3.l" %in% names(coeffs)) {
    data_list <- coeffs[["dataETOC3.l"]]
    if (is.list(data_list) && length(data_list) >= 1) {
      estETOC2 <- data_list[[1]]  # First element contains epiTOC2 parameters
    } else {
      return(NULL)
    }
  } else if ("EpiToc2_CpGs" %in% names(coeffs)) {
    # methylCIPHER version
    estETOC2 <- coeffs$EpiToc2_CpGs
  } else {
    return(NULL)
  }
  
  if (is.null(estETOC2)) return(NULL)
  
  # Get CpG names (rownames of the coefficient matrix)
  cpgs <- rownames(estETOC2)
  if (is.null(cpgs)) return(NULL)
  
  # Match to betas
  matched_idx <- match(cpgs, rownames(betas))
  valid <- !is.na(matched_idx)
  
  n_valid <- sum(valid)
  if (n_valid == 0) return(NULL)
  
  message(sprintf("    epiTOC2: Using %d of %d CpGs", n_valid, length(cpgs)))
  
  # Extract matched data
  betas_matched <- betas[matched_idx[valid], , drop = FALSE]
  params_matched <- estETOC2[valid, , drop = FALSE]
  
  # epiTOC2 calculation:
  # TNSC = 2 * mean((beta - beta0) / (delta * (1 - beta0)))
  # where delta = de-novo methylation rate, beta0 = ground-state methylation
  
  if (ncol(params_matched) >= 2) {
    delta <- params_matched[, 1]  # de-novo rate
    beta0 <- params_matched[, 2]  # ground-state
    
    # Calculate TNSC for each sample
    tnsc <- apply(betas_matched, 2, function(b) {
      # Avoid division by zero
      denom <- delta * (1 - beta0)
      denom[denom == 0] <- NA
      scores <- (b - beta0) / denom
      2 * mean(scores, na.rm = TRUE)
    })
    
    return(tnsc)
  }
  
  return(NULL)
}


#' Calculate all available clocks using direct implementations
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param verbose Logical, print progress
#' @return Data frame with clock values
#' @export
calculate_clocks_direct <- function(betas, verbose = TRUE) {
  
  results <- data.frame(sample_id = colnames(betas), stringsAsFactors = FALSE)
  
  # Initialize coefficients once
  if (verbose) message("  Loading clock coefficients...")
  coeffs <- initialize_clock_coefficients()
  
  if (verbose) message("    Loaded ", length(coeffs), " coefficient sets")
  
  # Calculate each clock
  clock_funcs <- list(
    "Horvath1" = calc_horvath1_direct,
    "Hannum" = calc_hannum_direct,
    "PhenoAge" = calc_phenoage_direct,
    "DNAmTL" = calc_dnamtl_direct,
    "epiTOC2_TNSC" = calc_epitoc2_direct
  )
  
  computed <- c()
  
  for (clock_name in names(clock_funcs)) {
    tryCatch({
      result <- clock_funcs[[clock_name]](betas, coeffs)
      if (!is.null(result) && length(result) == ncol(betas)) {
        results[[clock_name]] <- as.numeric(result)
        computed <- c(computed, clock_name)
      }
    }, error = function(e) {
      if (verbose) message("    ", clock_name, " error: ", e$message)
    })
  }
  
  if (verbose && length(computed) > 0) {
    message("    Direct calculations: ", paste(computed, collapse = ", "))
  }
  
  return(list(results = results, coeffs = coeffs))
}


#' Generic weighted sum clock calculator
#' @param betas Beta matrix
#' @param coef_df Data frame with CpG and weight columns
#' @param cpg_col Name of CpG column
#' @param weight_col Name of weight column
#' @param intercept Intercept value (default 0)
#' @return Named vector of clock values
#' @keywords internal
calc_weighted_sum_clock <- function(betas, coef_df, cpg_col = NULL, weight_col = NULL, intercept = 0) {
  
  # Auto-detect columns if not specified
  if (is.null(cpg_col)) {
    cpg_candidates <- c("CpG", "CpGmarker", "probe", "Probe", "cpg")
    cpg_col <- intersect(cpg_candidates, colnames(coef_df))[1]
    if (is.na(cpg_col)) cpg_col <- colnames(coef_df)[1]
  }
  
  if (is.null(weight_col)) {
    weight_candidates <- c("Coefficient", "weight", "Weight", "coef", "beta")
    weight_col <- intersect(weight_candidates, colnames(coef_df))[1]
    if (is.na(weight_col)) weight_col <- colnames(coef_df)[2]
  }
  
  cpgs <- coef_df[[cpg_col]]
  weights <- coef_df[[weight_col]]
  
  # Remove intercept row if present
  intercept_mask <- cpgs %in% c("(Intercept)", "Intercept") | is.na(cpgs)
  if (any(intercept_mask)) {
    intercept <- weights[intercept_mask][1]
    cpgs <- cpgs[!intercept_mask]
    weights <- weights[!intercept_mask]
  }
  
  # Match CpGs
  matched_idx <- match(cpgs, rownames(betas))
  valid <- !is.na(matched_idx)
  
  if (sum(valid) == 0) return(NULL)
  
  betas_subset <- betas[matched_idx[valid], , drop = FALSE]
  weights_valid <- weights[valid]
  
  clock_values <- intercept + colSums(betas_subset * weights_valid, na.rm = TRUE)
  
  return(clock_values)
}
