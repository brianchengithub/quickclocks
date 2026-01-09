#' ============================================================================
#' methylCIPHER Wrapper for Epigenetic Clock Calculator
#' ============================================================================
#' 
#' Wrapper for methylCIPHER (Comprehensive Interrogation of Phenotypes 
#' through Human Epigenomic Research)
#' https://github.com/HigginsChenLab/methylCIPHER

#' Calculate all methylCIPHER clocks
#' 
#' @param betas Matrix of beta values (probes x samples)
#' @param clocks Character vector of clock names to calculate (NULL = all)
#' @param verbose Print progress
#' @return Data frame with all requested clock results
#' @export
calculate_methylcipher_clocks <- function(betas, clocks = NULL, verbose = TRUE) {
  
  if (!requireNamespace("methylCIPHER", quietly = TRUE)) {
    warning("Package 'methylCIPHER' not installed. ",
            "Install with: devtools::install_github('HigginsChenLab/methylCIPHER')")
    return(NULL)
  }
  
  log_message("Calculating methylCIPHER clocks...", verbose)
  
  # Define all available clocks
  all_clocks <- get_methylcipher_clock_list()
  
  # Filter to requested clocks
  if (!is.null(clocks)) {
    available_requested <- intersect(clocks, all_clocks)
    if (length(available_requested) < length(clocks)) {
      not_found <- setdiff(clocks, all_clocks)
      warning("Some requested clocks not found: ", paste(not_found, collapse = ", "))
    }
    clocks_to_calc <- available_requested
  } else {
    clocks_to_calc <- all_clocks
  }
  
  log_message(sprintf("  Calculating %d clocks...", length(clocks_to_calc)), verbose)
  
  tryCatch({
    # Try the main methylCIPHER calculation function
    # Different package versions may have different interfaces
    
    if ("calcAllClocks" %in% ls(getNamespace("methylCIPHER"))) {
      results <- methylCIPHER::calcAllClocks(betas, clocks = clocks_to_calc)
    } else if ("DNAmAge" %in% ls(getNamespace("methylCIPHER"))) {
      results <- calculate_methylcipher_individual(betas, clocks_to_calc, verbose)
    } else {
      results <- calculate_methylcipher_generic(betas, clocks_to_calc, verbose)
    }
    
    # Ensure sample_id column
    if (!is.null(results) && !"sample_id" %in% colnames(results)) {
      results <- cbind(sample_id = colnames(betas), results)
    }
    
    return(results)
    
  }, error = function(e) {
    warning("methylCIPHER calculation failed: ", e$message)
    return(NULL)
  })
}


#' Get list of all methylCIPHER clocks
#' @return Character vector of clock names
#' @export
get_methylcipher_clock_list <- function() {
  c(
    # Age clocks
    "Horvath1", "Horvath2", "Hannum", "PhenoAge", 
    "GrimAge1", "GrimAge2", "Zhang", "Zhang2019", 
    "Lin", "DNAmTL",
    
    # Specialized clocks
    "AdaptAge", "CausAge", "DamAge", "HypoClock",
    "MiAge", "SystemsAge", "RetroelementAge450K",
    
    # Lifestyle predictors
    "Alcohol", "BMI", "Smoking",
    
    # Organ/system clocks
    "Blood", "Brain", "Heart", "Kidney", 
    "Liver", "Lung", "MusculoSkeletal",
    
    # Biological clocks
    "Inflammation", "Hormone", "Immune", "Metabolic"
  )
}


#' Calculate individual methylCIPHER clocks
#' @keywords internal
calculate_methylcipher_individual <- function(betas, clocks, verbose = TRUE) {
  
  results <- data.frame(sample_id = colnames(betas), stringsAsFactors = FALSE)
  
  for (clock in clocks) {
    if (verbose) {
      log_message(sprintf("    Calculating %s...", clock), verbose)
    }
    
    clock_result <- tryCatch({
      # Try to get the calculation function
      func_name <- paste0("calc", clock)
      
      if (exists(func_name, envir = asNamespace("methylCIPHER"))) {
        func <- get(func_name, envir = asNamespace("methylCIPHER"))
        func(betas)
      } else if (exists(clock, envir = asNamespace("methylCIPHER"))) {
        func <- get(clock, envir = asNamespace("methylCIPHER"))
        func(betas)
      } else {
        # Try generic DNAmAge function with clock parameter
        methylCIPHER::DNAmAge(betas, clock = clock)
      }
    }, error = function(e) {
      if (verbose) {
        log_message(sprintf("      Warning: %s", e$message), verbose, "warning")
      }
      return(rep(NA_real_, ncol(betas)))
    })
    
    if (length(clock_result) == ncol(betas)) {
      results[[clock]] <- as.numeric(clock_result)
    }
  }
  
  return(results)
}


#' Calculate methylCIPHER clocks using generic approach
#' @keywords internal
calculate_methylcipher_generic <- function(betas, clocks, verbose = TRUE) {
  
  # Get all exported functions
  pkg_funcs <- ls(getNamespace("methylCIPHER"))
  
  results <- data.frame(sample_id = colnames(betas), stringsAsFactors = FALSE)
  
  for (clock in clocks) {
    # Try different function naming patterns
    patterns <- c(
      paste0("calc", clock),
      paste0("Calculate", clock),
      clock,
      tolower(clock),
      toupper(clock)
    )
    
    found <- FALSE
    for (pattern in patterns) {
      if (pattern %in% pkg_funcs) {
        tryCatch({
          func <- get(pattern, envir = asNamespace("methylCIPHER"))
          result <- func(betas)
          
          if (length(result) == ncol(betas)) {
            results[[clock]] <- as.numeric(result)
            found <- TRUE
            break
          }
        }, error = function(e) {
          # Try next pattern
        })
      }
    }
    
    if (!found && verbose) {
      log_message(sprintf("    Warning: Could not calculate %s", clock), 
                  verbose, "warning")
    }
  }
  
  return(results)
}


#' Calculate Horvath1 (Pan-tissue) clock
#' @export
calculate_horvath1 <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "Horvath1", verbose)
}

#' Calculate Horvath2 (Skin+Blood) clock
#' @export
calculate_horvath2 <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "Horvath2", verbose)
}

#' Calculate Hannum clock
#' @export
calculate_hannum <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "Hannum", verbose)
}

#' Calculate PhenoAge (Levine) clock
#' @export
calculate_phenoage <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "PhenoAge", verbose)
}

#' Calculate GrimAge v1
#' @export
calculate_grimage1 <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "GrimAge1", verbose)
}

#' Calculate GrimAge v2
#' @export
calculate_grimage2 <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "GrimAge2", verbose)
}

#' Calculate Zhang clock
#' @export
calculate_zhang <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "Zhang", verbose)
}

#' Calculate Zhang2019 clock
#' @export
calculate_zhang2019 <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "Zhang2019", verbose)
}

#' Calculate Lin clock
#' @export
calculate_lin <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "Lin", verbose)
}

#' Calculate DNAmTL (telomere length)
#' @export
calculate_dnamtl <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "DNAmTL", verbose)
}

#' Calculate AdaptAge
#' @export
calculate_adaptage <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "AdaptAge", verbose)
}

#' Calculate CausAge
#' @export
calculate_causage <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "CausAge", verbose)
}

#' Calculate DamAge
#' @export
calculate_damage <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "DamAge", verbose)
}

#' Calculate SystemsAge
#' @export
calculate_systemsage <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "SystemsAge", verbose)
}

#' Calculate Retroelement-Age 450K
#' @export
calculate_retroelement_age <- function(betas, verbose = TRUE) {
  calculate_single_clock(betas, "RetroelementAge450K", verbose)
}


#' Helper function to calculate a single clock
#' @keywords internal
calculate_single_clock <- function(betas, clock_name, verbose = TRUE) {
  
  if (!requireNamespace("methylCIPHER", quietly = TRUE)) {
    return(NULL)
  }
  
  log_message(sprintf("Calculating %s...", clock_name), verbose)
  
  tryCatch({
    result <- methylCIPHER::calcAllClocks(betas, clocks = clock_name)
    
    if (is.data.frame(result)) {
      return(result)
    }
    
    return(data.frame(
      sample_id = colnames(betas),
      clock_name = as.numeric(result),
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    warning(sprintf("%s calculation failed: %s", clock_name, e$message))
    return(NULL)
  })
}


#' Calculate lifestyle predictors (Alcohol, BMI, Smoking)
#' 
#' @param betas Beta matrix
#' @param verbose Print messages
#' @return Data frame with lifestyle predictions
#' @export
calculate_lifestyle_clocks <- function(betas, verbose = TRUE) {
  
  log_message("Calculating lifestyle predictors...", verbose)
  
  clocks <- c("Alcohol", "BMI", "Smoking")
  calculate_methylcipher_clocks(betas, clocks = clocks, verbose = verbose)
}


#' Calculate organ/system clocks
#' 
#' @param betas Beta matrix
#' @param verbose Print messages
#' @return Data frame with organ clock results
#' @export
calculate_organ_clocks <- function(betas, verbose = TRUE) {
  
  log_message("Calculating organ/system clocks...", verbose)
  
  clocks <- c("Blood", "Brain", "Heart", "Kidney", 
              "Liver", "Lung", "MusculoSkeletal")
  calculate_methylcipher_clocks(betas, clocks = clocks, verbose = verbose)
}


#' Calculate biological clocks
#' 
#' @param betas Beta matrix
#' @param verbose Print messages
#' @return Data frame with biological clock results
#' @export
calculate_biological_clocks <- function(betas, verbose = TRUE) {
  
  log_message("Calculating biological system clocks...", verbose)
  
  clocks <- c("Inflammation", "Hormone", "Immune", "Metabolic")
  calculate_methylcipher_clocks(betas, clocks = clocks, verbose = verbose)
}


#' Get probes required for a methylCIPHER clock
#' 
#' @param clock_name Clock name
#' @return Character vector of probe IDs
#' @export
get_methylcipher_clock_probes <- function(clock_name) {
  
  if (!requireNamespace("methylCIPHER", quietly = TRUE)) {
    return(NULL)
  }
  
  tryCatch({
    # Try to get probe requirements from package data
    # The exact structure depends on package version
    
    data_name <- paste0(clock_name, "_probes")
    
    if (exists(data_name, envir = asNamespace("methylCIPHER"))) {
      return(get(data_name, envir = asNamespace("methylCIPHER")))
    }
    
    # Try loading from package data
    pkg_data <- tryCatch({
      data(list = data_name, package = "methylCIPHER", envir = environment())
      get(data_name, envir = environment())
    }, error = function(e) {
      NULL
    })
    
    return(pkg_data)
    
  }, error = function(e) {
    NULL
  })
}


#' Check methylCIPHER availability for a specific clock
#' 
#' @param clock_name Clock name
#' @param available_probes Available probe IDs
#' @return List with availability info
#' @export
check_methylcipher_clock_availability <- function(clock_name, available_probes) {
  
  required <- get_methylcipher_clock_probes(clock_name)
  
  if (is.null(required) || length(required) == 0) {
    pkg_available <- requireNamespace("methylCIPHER", quietly = TRUE)
    
    return(list(
      available = pkg_available,
      n_required = NA,
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
