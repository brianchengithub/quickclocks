#' ============================================================================
#' DunedinPACE Wrapper for Epigenetic Clock Calculator
#' ============================================================================
#' 
#' Wrapper for DunedinPACE (Pace of Aging Calculated from Epigenome)
#' https://github.com/danbelsky/DunedinPACE

#' Calculate DunedinPACE
#' 
#' @param betas Matrix of beta values (probes x samples)
#' @param verbose Print progress
#' @return Data frame with DunedinPACE results
#' @export
calculate_dunedin_pace <- function(betas, verbose = TRUE) {
  
  if (!requireNamespace("DunedinPACE", quietly = TRUE)) {
    warning("Package 'DunedinPACE' not installed. ",
            "Install with: devtools::install_github('danbelsky/DunedinPACE')")
    return(NULL)
  }
  
  log_message("Calculating DunedinPACE...", verbose)
  
  tryCatch({
    # DunedinPACE expects a data frame with samples as rows
    # Need to transpose the beta matrix
    betas_df <- as.data.frame(t(betas))
    
    # Calculate DunedinPACE
    # The function signature may vary by version
    pace_result <- DunedinPACE::DunedinPACE(betas_df)
    
    # Format results
    if (is.data.frame(pace_result)) {
      results <- pace_result
      results$sample_id <- rownames(betas_df)
    } else if (is.numeric(pace_result)) {
      results <- data.frame(
        sample_id = rownames(betas_df),
        DunedinPACE = pace_result,
        stringsAsFactors = FALSE
      )
    } else {
      # Try to convert to data frame
      results <- data.frame(
        sample_id = colnames(betas),
        DunedinPACE = as.numeric(pace_result),
        stringsAsFactors = FALSE
      )
    }
    
    log_message(sprintf("  Calculated DunedinPACE for %d samples", nrow(results)), 
                verbose)
    log_message(sprintf("  Mean: %.3f, SD: %.3f", 
                        mean(results$DunedinPACE, na.rm = TRUE),
                        sd(results$DunedinPACE, na.rm = TRUE)), verbose)
    
    return(results)
    
  }, error = function(e) {
    warning("DunedinPACE calculation failed: ", e$message)
    return(NULL)
  })
}


#' Calculate DunedinPACE using PoAm45 method
#' 
#' Alternative calculation using the PoAm45 variant if available.
#' 
#' @param betas Matrix of beta values
#' @param verbose Print progress
#' @return Data frame with results
#' @export
calculate_dunedin_poam45 <- function(betas, verbose = TRUE) {
  
  if (!requireNamespace("DunedinPACE", quietly = TRUE)) {
    return(NULL)
  }
  
  log_message("Calculating DunedinPoAm45...", verbose)
  
  tryCatch({
    betas_df <- as.data.frame(t(betas))
    
    # Try PoAm45 if available
    if ("PoAm45" %in% ls(getNamespace("DunedinPACE"))) {
      poam_result <- DunedinPACE::PoAm45(betas_df)
    } else {
      log_message("  PoAm45 not available in installed version", verbose)
      return(NULL)
    }
    
    results <- data.frame(
      sample_id = rownames(betas_df),
      DunedinPoAm45 = as.numeric(poam_result),
      stringsAsFactors = FALSE
    )
    
    return(results)
    
  }, error = function(e) {
    warning("DunedinPoAm45 calculation failed: ", e$message)
    return(NULL)
  })
}


#' Get DunedinPACE required probes
#' 
#' @return Character vector of required probe IDs
#' @export
get_dunedin_pace_probes <- function() {
  
  if (!requireNamespace("DunedinPACE", quietly = TRUE)) {
    return(NULL)
  }
  
  tryCatch({
    # Try to get probes from package data
    # The exact location depends on package version
    
    # Method 1: Try package data
    probes <- tryCatch({
      data("DunedinPACE_probes", package = "DunedinPACE", envir = environment())
      get("DunedinPACE_probes", envir = environment())
    }, error = function(e) NULL)
    
    if (!is.null(probes)) return(probes)
    
    # Method 2: Try to extract from model
    probes <- tryCatch({
      model <- DunedinPACE::DunedinPACE_model
      if (is.list(model) && "probes" %in% names(model)) {
        model$probes
      } else {
        names(model)
      }
    }, error = function(e) NULL)
    
    return(probes)
    
  }, error = function(e) {
    NULL
  })
}


#' Check DunedinPACE probe availability
#' 
#' @param available_probes Character vector of available probe IDs
#' @return List with availability statistics
#' @export
check_dunedin_probes <- function(available_probes) {
  
  required <- get_dunedin_pace_probes()
  
  if (is.null(required)) {
    return(list(
      n_required = NA,
      n_available = NA,
      pct_available = NA,
      missing = NULL
    ))
  }
  
  n_required <- length(required)
  n_available <- sum(required %in% available_probes)
  missing <- setdiff(required, available_probes)
  
  list(
    n_required = n_required,
    n_available = n_available,
    pct_available = 100 * n_available / n_required,
    missing = missing
  )
}


#' Null coalescing operator
#' @keywords internal
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}
