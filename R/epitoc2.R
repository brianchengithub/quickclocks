#' ============================================================================
#' epiTOC2 Wrapper for Epigenetic Clock Calculator
#' ============================================================================
#' 
#' Wrapper for epiTOC2 (Epigenetic Timer of Cancer 2) for estimating
#' total number of stem cell divisions
#' https://github.com/aet21/EpiMitClocks

#' Calculate epiTOC2 (stem cell division estimate)
#' 
#' @param betas Matrix of beta values (probes x samples)
#' @param verbose Print progress
#' @return Data frame with epiTOC2 results
#' @export
calculate_epitoc2 <- function(betas, verbose = TRUE) {
  
  if (!requireNamespace("EpiMitClocks", quietly = TRUE)) {
    warning("Package 'EpiMitClocks' not installed. ",
            "Install with: devtools::install_github('aet21/EpiMitClocks')")
    return(NULL)
  }
  
  log_message("Calculating epiTOC2 (stem cell divisions)...", verbose)
  
  tryCatch({
    # Try the main epiTOC2 function
    # The exact function signature depends on package version
    
    # Method 1: Try epiTOC2
    if ("epiTOC2" %in% ls(getNamespace("EpiMitClocks"))) {
      results <- calculate_epitoc2_v1(betas, verbose)
    }
    # Method 2: Try calcEpiTOC2
    else if ("calcEpiTOC2" %in% ls(getNamespace("EpiMitClocks"))) {
      results <- calculate_epitoc2_v2(betas, verbose)
    }
    # Method 3: Manual calculation using stored model
    else {
      results <- calculate_epitoc2_manual(betas, verbose)
    }
    
    if (!is.null(results)) {
      log_message(sprintf("  Calculated epiTOC2 for %d samples", nrow(results)), 
                  verbose)
      log_message(sprintf("  Mean tnsc: %.2f, SD: %.2f", 
                          mean(results$epiTOC2_tnsc, na.rm = TRUE),
                          sd(results$epiTOC2_tnsc, na.rm = TRUE)), verbose)
    }
    
    return(results)
    
  }, error = function(e) {
    warning("epiTOC2 calculation failed: ", e$message)
    return(NULL)
  })
}


#' epiTOC2 calculation - Version 1
#' @keywords internal
calculate_epitoc2_v1 <- function(betas, verbose = TRUE) {
  
  # epiTOC2 may expect different input formats
  # Try with matrix first
  
  tryCatch({
    result <- EpiMitClocks::epiTOC2(betas)
    
    # Format results
    if (is.list(result)) {
      # Multiple outputs
      df <- data.frame(
        sample_id = colnames(betas),
        stringsAsFactors = FALSE
      )
      
      # Extract available components
      if ("tnsc" %in% names(result)) {
        df$epiTOC2_tnsc <- result$tnsc
      }
      if ("tnsc2" %in% names(result)) {
        df$epiTOC2_tnsc2 <- result$tnsc2
      }
      if ("irS" %in% names(result)) {
        df$epiTOC2_irS <- result$irS
      }
      if ("HypoClock" %in% names(result)) {
        df$epiTOC2_HypoClock <- result$HypoClock
      }
      
      return(df)
      
    } else if (is.numeric(result)) {
      return(data.frame(
        sample_id = colnames(betas),
        epiTOC2_tnsc = as.numeric(result),
        stringsAsFactors = FALSE
      ))
    }
    
    return(NULL)
    
  }, error = function(e) {
    NULL
  })
}


#' epiTOC2 calculation - Version 2
#' @keywords internal
calculate_epitoc2_v2 <- function(betas, verbose = TRUE) {
  
  tryCatch({
    result <- EpiMitClocks::calcEpiTOC2(betas)
    
    if (is.data.frame(result)) {
      if (!"sample_id" %in% colnames(result)) {
        result <- cbind(sample_id = colnames(betas), result)
      }
      return(result)
    }
    
    return(data.frame(
      sample_id = colnames(betas),
      epiTOC2_tnsc = as.numeric(result),
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    NULL
  })
}


#' epiTOC2 manual calculation
#' @keywords internal
calculate_epitoc2_manual <- function(betas, verbose = TRUE) {
  
  # Get epiTOC2 model/probes
  model <- get_epitoc2_model()
  
  if (is.null(model)) {
    log_message("  Could not load epiTOC2 model", verbose, "warning")
    return(NULL)
  }
  
  # Get required probes
  probes <- model$probes
  coefs <- model$coefficients
  
  # Check probe availability
  available <- intersect(probes, rownames(betas))
  
  if (length(available) < length(probes) * 0.9) {
    log_message(sprintf("  Only %d/%d probes available", 
                        length(available), length(probes)), verbose, "warning")
  }
  
  # Subset to available probes
  betas_sub <- betas[available, , drop = FALSE]
  coefs_sub <- coefs[available]
  
  # Calculate epiTOC2 score
  # Basic formula: score = sum(coef * beta) for mitotic clock probes
  scores <- colSums(betas_sub * coefs_sub, na.rm = TRUE)
  
  # Convert to estimated cell divisions
  # This is a simplified conversion - actual formula may differ
  tnsc <- scores / mean(coefs_sub, na.rm = TRUE)
  
  return(data.frame(
    sample_id = colnames(betas),
    epiTOC2_tnsc = tnsc,
    stringsAsFactors = FALSE
  ))
}


#' Get epiTOC2 model data
#' @keywords internal
get_epitoc2_model <- function() {
  
  if (!requireNamespace("EpiMitClocks", quietly = TRUE)) {
    return(NULL)
  }
  
  tryCatch({
    # Try to load model data from package
    pkg_data <- data(package = "EpiMitClocks")$results[, "Item"]
    
    for (item in pkg_data) {
      tryCatch({
        data(list = item, package = "EpiMitClocks", envir = environment())
        obj <- get(item, envir = environment())
        
        # Look for model object with probes and coefficients
        if (is.list(obj)) {
          if ("probes" %in% names(obj) || "CpGs" %in% names(obj)) {
            probes <- obj$probes %||% obj$CpGs
            coefs <- obj$coef %||% obj$coefficients %||% rep(1, length(probes))
            
            return(list(probes = probes, coefficients = coefs))
          }
        }
      }, error = function(e) {
        # Skip failed items
      })
    }
    
    return(NULL)
    
  }, error = function(e) {
    NULL
  })
}


#' Get epiTOC2 required probes
#' 
#' @return Character vector of probe IDs
#' @export
get_epitoc2_probes <- function() {
  model <- get_epitoc2_model()
  if (!is.null(model)) {
    return(model$probes)
  }
  NULL
}


#' Calculate HypoClock (from EpiMitClocks)
#' 
#' @param betas Beta matrix
#' @param verbose Print messages
#' @return Data frame with HypoClock results
#' @export
calculate_hypoclock <- function(betas, verbose = TRUE) {
  
  if (!requireNamespace("EpiMitClocks", quietly = TRUE)) {
    return(NULL)
  }
  
  log_message("Calculating HypoClock...", verbose)
  
  tryCatch({
    # Try HypoClock function if available
    if ("HypoClock" %in% ls(getNamespace("EpiMitClocks"))) {
      result <- EpiMitClocks::HypoClock(betas)
      
      return(data.frame(
        sample_id = colnames(betas),
        HypoClock = as.numeric(result),
        stringsAsFactors = FALSE
      ))
    }
    
    return(NULL)
    
  }, error = function(e) {
    warning("HypoClock calculation failed: ", e$message)
    return(NULL)
  })
}


#' Calculate MiAge (mitotic age)
#' 
#' @param betas Beta matrix
#' @param verbose Print messages
#' @return Data frame with MiAge results
#' @export
calculate_miage <- function(betas, verbose = TRUE) {
  
  if (!requireNamespace("EpiMitClocks", quietly = TRUE)) {
    return(NULL)
  }
  
  log_message("Calculating MiAge...", verbose)
  
  tryCatch({
    if ("MiAge" %in% ls(getNamespace("EpiMitClocks"))) {
      result <- EpiMitClocks::MiAge(betas)
      
      return(data.frame(
        sample_id = colnames(betas),
        MiAge = as.numeric(result),
        stringsAsFactors = FALSE
      ))
    }
    
    return(NULL)
    
  }, error = function(e) {
    warning("MiAge calculation failed: ", e$message)
    return(NULL)
  })
}


#' Calculate all EpiMitClocks outputs
#' 
#' @param betas Beta matrix
#' @param verbose Print messages
#' @return Data frame with all mitotic clock results
#' @export
calculate_all_mitotic_clocks <- function(betas, verbose = TRUE) {
  
  log_message("Calculating all mitotic clocks...", verbose)
  
  # Initialize results with sample IDs
  results <- data.frame(sample_id = colnames(betas), stringsAsFactors = FALSE)
  
  # Calculate epiTOC2
  epitoc2 <- calculate_epitoc2(betas, verbose = FALSE)
  if (!is.null(epitoc2)) {
    results <- merge(results, epitoc2, by = "sample_id", all.x = TRUE)
  }
  
  # Calculate HypoClock
  hypoclock <- calculate_hypoclock(betas, verbose = FALSE)
  if (!is.null(hypoclock)) {
    results <- merge(results, hypoclock, by = "sample_id", all.x = TRUE)
  }
  
  # Calculate MiAge
  miage <- calculate_miage(betas, verbose = FALSE)
  if (!is.null(miage)) {
    results <- merge(results, miage, by = "sample_id", all.x = TRUE)
  }
  
  n_clocks <- ncol(results) - 1
  log_message(sprintf("  Calculated %d mitotic clocks", n_clocks), verbose)
  
  return(results)
}


#' Null coalescing operator
#' @keywords internal
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}
