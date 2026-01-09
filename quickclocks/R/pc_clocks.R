#' ============================================================================
#' PC-Clocks Wrapper for Epigenetic Clock Calculator
#' ============================================================================
#' 
#' Wrapper for PC-Clocks (Principal Component Clocks)
#' https://github.com/MorganLevineLab/PC-Clocks

#' Calculate all PC-Clocks
#' 
#' @param betas Matrix of beta values (probes x samples)
#' @param verbose Print progress
#' @return Data frame with all PC-Clock results
#' @export
calculate_pc_clocks <- function(betas, verbose = TRUE) {
  
  if (!requireNamespace("PCClocks", quietly = TRUE)) {
    warning("Package 'PCClocks' not installed. ",
            "Install with: devtools::install_github('MorganLevineLab/PC-Clocks')")
    return(NULL)
  }
  
  log_message("Calculating PC-Clocks...", verbose)
  
  tryCatch({
    # PC-Clocks may expect a specific input format
    # Check the package documentation for exact requirements
    
    # Try the main calculation function
    # The exact function name may vary by package version
    
    # Method 1: Try PredictAge if available
    if ("PredictAge" %in% ls(getNamespace("PCClocks"))) {
      results <- calculate_pc_clocks_v1(betas, verbose)
    }
    # Method 2: Try calcPCClocks
    else if ("calcPCClocks" %in% ls(getNamespace("PCClocks"))) {
      results <- calculate_pc_clocks_v2(betas, verbose)
    }
    # Method 3: Generic approach
    else {
      results <- calculate_pc_clocks_generic(betas, verbose)
    }
    
    return(results)
    
  }, error = function(e) {
    warning("PC-Clocks calculation failed: ", e$message)
    return(NULL)
  })
}


#' PC-Clocks calculation - Version 1 (using PredictAge)
#' @keywords internal
calculate_pc_clocks_v1 <- function(betas, verbose = TRUE) {
  
  # Transpose to samples as rows
  betas_t <- t(betas)
  
  # Calculate each PC clock
  clocks <- c("PCHorvath1", "PCHorvath2", "PCHannum", "PCPhenoAge", "PCGrimAge")
  
  results <- data.frame(sample_id = rownames(betas_t), stringsAsFactors = FALSE)
  
  for (clock in clocks) {
    log_message(sprintf("  Calculating %s...", clock), verbose)
    
    tryCatch({
      clock_values <- PCClocks::PredictAge(betas_t, clock = clock)
      results[[clock]] <- as.numeric(clock_values)
    }, error = function(e) {
      log_message(sprintf("    Warning: %s failed: %s", clock, e$message), 
                  verbose, "warning")
      results[[clock]] <- NA_real_
    })
  }
  
  return(results)
}


#' PC-Clocks calculation - Version 2 (using calcPCClocks)
#' @keywords internal
calculate_pc_clocks_v2 <- function(betas, verbose = TRUE) {
  
  betas_t <- t(betas)
  
  tryCatch({
    results <- PCClocks::calcPCClocks(betas_t)
    
    # Ensure sample_id column
    if (!"sample_id" %in% colnames(results)) {
      results <- cbind(sample_id = rownames(betas_t), results)
    }
    
    return(as.data.frame(results))
    
  }, error = function(e) {
    warning("calcPCClocks failed: ", e$message)
    return(NULL)
  })
}


#' PC-Clocks calculation - Generic approach
#' @keywords internal
calculate_pc_clocks_generic <- function(betas, verbose = TRUE) {
  
  log_message("  Using generic PC-Clocks approach...", verbose)
  
  # Get exported functions from PCClocks
  pc_funcs <- ls(getNamespace("PCClocks"))
  
  # Look for age prediction functions
  age_funcs <- grep("(age|Age|clock|Clock)", pc_funcs, value = TRUE)
  
  if (length(age_funcs) == 0) {
    warning("No recognizable age prediction functions in PCClocks")
    return(NULL)
  }
  
  log_message(sprintf("  Found functions: %s", 
                      paste(age_funcs, collapse = ", ")), verbose)
  
  # Try each function
  betas_t <- t(betas)
  results <- data.frame(sample_id = rownames(betas_t), stringsAsFactors = FALSE)
  
  for (func_name in age_funcs) {
    tryCatch({
      func <- get(func_name, envir = asNamespace("PCClocks"))
      if (is.function(func)) {
        clock_result <- func(betas_t)
        if (is.numeric(clock_result) && length(clock_result) == nrow(betas_t)) {
          results[[func_name]] <- clock_result
        }
      }
    }, error = function(e) {
      # Skip failed functions
    })
  }
  
  if (ncol(results) == 1) {
    warning("No PC-Clocks could be calculated")
    return(NULL)
  }
  
  return(results)
}


#' Calculate individual PC-Horvath1 clock
#' 
#' @param betas Beta matrix
#' @param verbose Print messages
#' @return Numeric vector of ages
#' @export
calculate_pc_horvath1 <- function(betas, verbose = TRUE) {
  
  if (!requireNamespace("PCClocks", quietly = TRUE)) {
    return(NULL)
  }
  
  tryCatch({
    PCClocks::PredictAge(t(betas), clock = "PCHorvath1")
  }, error = function(e) {
    NULL
  })
}


#' Calculate individual PC-Horvath2 clock
#' @export
calculate_pc_horvath2 <- function(betas, verbose = TRUE) {
  
  if (!requireNamespace("PCClocks", quietly = TRUE)) {
    return(NULL)
  }
  
  tryCatch({
    PCClocks::PredictAge(t(betas), clock = "PCHorvath2")
  }, error = function(e) {
    NULL
  })
}


#' Calculate individual PC-Hannum clock
#' @export
calculate_pc_hannum <- function(betas, verbose = TRUE) {
  
  if (!requireNamespace("PCClocks", quietly = TRUE)) {
    return(NULL)
  }
  
  tryCatch({
    PCClocks::PredictAge(t(betas), clock = "PCHannum")
  }, error = function(e) {
    NULL
  })
}


#' Calculate individual PC-PhenoAge clock
#' @export
calculate_pc_phenoage <- function(betas, verbose = TRUE) {
  
  if (!requireNamespace("PCClocks", quietly = TRUE)) {
    return(NULL)
  }
  
  tryCatch({
    PCClocks::PredictAge(t(betas), clock = "PCPhenoAge")
  }, error = function(e) {
    NULL
  })
}


#' Calculate individual PC-GrimAge clock
#' @export
calculate_pc_grimage <- function(betas, verbose = TRUE) {
  
  if (!requireNamespace("PCClocks", quietly = TRUE)) {
    return(NULL)
  }
  
  tryCatch({
    PCClocks::PredictAge(t(betas), clock = "PCGrimAge")
  }, error = function(e) {
    NULL
  })
}


#' Get PC-Clocks required probes
#' 
#' @param clock_name Optional specific clock name
#' @return Character vector of required probe IDs
#' @export
get_pc_clock_probes <- function(clock_name = NULL) {
  
  if (!requireNamespace("PCClocks", quietly = TRUE)) {
    return(NULL)
  }
  
  tryCatch({
    # Try to access probe lists from package
    probes <- NULL
    
    # Check for stored coefficients or probe lists
    pkg_data <- data(package = "PCClocks")$results[, "Item"]
    
    for (item in pkg_data) {
      tryCatch({
        data(list = item, package = "PCClocks", envir = environment())
        obj <- get(item, envir = environment())
        
        if (is.data.frame(obj) && "CpG" %in% colnames(obj)) {
          probes <- c(probes, obj$CpG)
        } else if (is.list(obj) && "CpGs" %in% names(obj)) {
          probes <- c(probes, obj$CpGs)
        } else if (is.numeric(obj) && !is.null(names(obj))) {
          probes <- c(probes, names(obj))
        }
      }, error = function(e) {
        # Skip items that can't be loaded
      })
    }
    
    return(unique(probes))
    
  }, error = function(e) {
    NULL
  })
}


#' Check if PC-Clocks can be calculated
#' 
#' @param available_probes Available probe IDs
#' @return List with availability info
#' @export
check_pc_clocks_availability <- function(available_probes) {
  
  required <- get_pc_clock_probes()
  
  if (is.null(required) || length(required) == 0) {
    # Can't determine requirements - check if package works
    pkg_available <- requireNamespace("PCClocks", quietly = TRUE)
    
    return(list(
      available = pkg_available,
      n_required = NA,
      n_available = NA,
      pct_available = if (pkg_available) 100 else 0,
      reason = if (pkg_available) "Package available" else "Package not installed"
    ))
  }
  
  n_required <- length(required)
  n_available <- sum(required %in% available_probes)
  pct <- 100 * n_available / n_required
  
  list(
    available = pct >= 95,
    n_required = n_required,
    n_available = n_available,
    pct_available = pct,
    reason = sprintf("%.1f%% probes available", pct)
  )
}


#' Null coalescing operator
#' @keywords internal
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}
