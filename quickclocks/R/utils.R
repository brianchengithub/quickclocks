#' ============================================================================
#' Utility Functions for Epigenetic Clock Calculator
#' ============================================================================

#' Log a message with timestamp
#' @param message The message to log
#' @param verbose Whether to print the message
#' @param level Log level: "info", "warning", "error", "debug"
#' @export
log_message <- function(message, verbose = TRUE, level = "info") {
  if (!verbose) return(invisible(NULL))
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  prefix <- switch(level,
    "info" = "[INFO]",
    "warning" = "[WARN]",
    "error" = "[ERROR]",
    "debug" = "[DEBUG]",
    "[INFO]"
  )
  
  cat(sprintf("%s %s %s\n", timestamp, prefix, message))
}

#' Print a section header
#' @param title Section title
#' @param verbose Whether to print
#' @export
log_section <- function(title, verbose = TRUE) {
  if (!verbose) return(invisible(NULL))
  
  cat("\n")
  cat("───────────────────────────────────────────────────────────────\n")
  cat(sprintf("  %s\n", title))
  cat("───────────────────────────────────────────────────────────────\n")
  cat("\n")
}

#' Check if a package is installed
#' @param pkg Package name
#' @return Logical indicating if package is available
#' @export
is_package_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

#' Load a package with error handling
#' @param pkg Package name
#' @param verbose Print messages
#' @return Logical indicating success
#' @export
safe_load_package <- function(pkg, verbose = TRUE) {
  if (!is_package_installed(pkg)) {
    if (verbose) {
      log_message(sprintf("Package '%s' not installed", pkg), verbose, "warning")
    }
    return(FALSE)
  }
  
  tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    return(TRUE)
  }, error = function(e) {
    if (verbose) {
      log_message(sprintf("Failed to load package '%s': %s", pkg, e$message), 
                  verbose, "error")
    }
    return(FALSE)
  })
}

#' Check and load all required packages
#' @param verbose Print messages
#' @return List of package availability
#' @export
check_required_packages <- function(verbose = TRUE) {
  # Core required packages
  core_packages <- c("sesame", "sesameData")
  
  # Clock-specific packages
  clock_packages <- c("DunedinPACE", "PCClocks", "methylclock", "methylCIPHER")
  
  # Optional utility packages
  util_packages <- c("parallel", "doParallel", "foreach", "data.table")
  
  results <- list(
    core = sapply(core_packages, is_package_installed),
    clocks = sapply(clock_packages, is_package_installed),
    utils = sapply(util_packages, is_package_installed)
  )
  
  if (verbose) {
    log_message("Package availability:", verbose)
    
    cat("  Core packages:\n")
    for (pkg in names(results$core)) {
      status <- if (results$core[pkg]) "✓" else "✗"
      cat(sprintf("    %s %s\n", status, pkg))
    }
    
    cat("  Clock packages:\n")
    for (pkg in names(results$clocks)) {
      status <- if (results$clocks[pkg]) "✓" else "✗"
      cat(sprintf("    %s %s\n", status, pkg))
    }
  }
  
  return(results)
}

#' Detect the methylation array platform from probe names
#' 
#' Supports: MSA, EPIC+ (EPICv2), EPIC (EPICv1), 450K, 27K
#' 
#' @param probe_names Character vector of probe names
#' @return Platform string
#' @export
detect_platform <- function(probe_names) {
  n_probes <- length(probe_names)
  
  # Count cg and ch probes
  n_cg <- sum(grepl("^cg", probe_names))
  n_ch <- sum(grepl("^ch", probe_names))
  
  # Check for platform-specific probes
  # EPICv2/EPIC+ has probes with "_TC" suffix and "nv_" prefix
  has_epic_v2 <- any(grepl("^cg.*_TC", probe_names)) || 
                 any(grepl("^nv_", probe_names, ignore.case = TRUE))
  
  # MSA (Methylation Screening Array) has specific probe naming
  has_msa <- any(grepl("^MSA", probe_names, ignore.case = TRUE))
  
  # Approximate probe counts per platform
  # MSA: ~285k probes
  # EPIC+ (EPICv2): ~930k probes
  # EPIC v1: ~865k probes
  # 450K: ~485k probes
  # 27K: ~27k probes
  
  if (has_msa || (n_probes > 250000 && n_probes < 350000)) {
    return("MSA")
  } else if (has_epic_v2 || n_probes > 900000) {
    return("EPICv2/EPIC+")
  } else if (n_probes > 800000) {
    return("EPIC")
  } else if (n_probes > 400000) {
    return("450K")
  } else if (n_probes > 20000) {
    return("27K")
  } else {
    return("Unknown")
  }
}

#' Validate a beta value matrix
#' @param betas Matrix of beta values
#' @return Validated matrix
#' @export
validate_beta_matrix <- function(betas) {
  # Check if matrix
  if (!is.matrix(betas)) {
    if (is.data.frame(betas)) {
      betas <- as.matrix(betas)
    } else {
      stop("Input must be a matrix or data.frame")
    }
  }
  
  # Check for row names (probe IDs)
  if (is.null(rownames(betas))) {
    stop("Beta matrix must have row names (probe IDs)")
  }
  
  # Check for column names (sample IDs)
  if (is.null(colnames(betas))) {
    warning("Beta matrix has no column names. Generating sample IDs.")
    colnames(betas) <- paste0("Sample_", seq_len(ncol(betas)))
  }
  
  # Check numeric values
  if (!is.numeric(betas)) {
    stop("Beta values must be numeric")
  }
  
  # Check value range (allowing for some measurement error)
  valid_range <- betas >= -0.1 & betas <= 1.1
  if (!all(valid_range, na.rm = TRUE)) {
    n_out_of_range <- sum(!valid_range, na.rm = TRUE)
    warning(sprintf("%d values outside expected range [0, 1]", n_out_of_range))
  }
  
  # Clip extreme values
  betas[betas < 0] <- 0
  betas[betas > 1] <- 1
  
  return(betas)
}

#' Calculate missing value statistics
#' @param betas Beta value matrix
#' @return List of statistics
#' @export
calculate_missing_stats <- function(betas) {
  n_total <- length(betas)
  n_missing <- sum(is.na(betas))
  
  # Per-sample statistics
  sample_missing <- colSums(is.na(betas))
  sample_missing_pct <- 100 * sample_missing / nrow(betas)
  
  # Per-probe statistics
  probe_missing <- rowSums(is.na(betas))
  probe_missing_pct <- 100 * probe_missing / ncol(betas)
  
  list(
    n_total = n_total,
    n_missing = n_missing,
    overall_missing_pct = 100 * n_missing / n_total,
    sample_missing = sample_missing,
    sample_missing_pct = sample_missing_pct,
    probe_missing = probe_missing,
    probe_missing_pct = probe_missing_pct,
    samples_above_threshold = names(sample_missing_pct[sample_missing_pct > 20]),
    probes_above_threshold = names(probe_missing_pct[probe_missing_pct > 50])
  )
}

#' Safe wrapper to run a function and catch errors
#' @param func Function to run
#' @param ... Arguments to pass to function
#' @param error_value Value to return on error
#' @param verbose Print error messages
#' @return Function result or error_value
#' @export
safe_run <- function(func, ..., error_value = NULL, verbose = TRUE) {
  tryCatch({
    func(...)
  }, error = function(e) {
    if (verbose) {
      log_message(sprintf("Error: %s", e$message), verbose, "error")
    }
    return(error_value)
  })
}

#' Convert M-values to beta values
#' @param m_values Matrix of M-values
#' @return Matrix of beta values
#' @export
m_to_beta <- function(m_values) {
  2^m_values / (2^m_values + 1)
}

#' Convert beta values to M-values
#' @param beta_values Matrix of beta values
#' @return Matrix of M-values
#' @export
beta_to_m <- function(beta_values) {
  log2(beta_values / (1 - beta_values))
}

#' Subset matrix to specific probes
#' @param betas Beta matrix
#' @param probes Vector of probe IDs to keep
#' @param allow_missing Allow missing probes (fill with NA)
#' @return Subsetted matrix
#' @export
subset_probes <- function(betas, probes, allow_missing = TRUE) {
  available <- intersect(probes, rownames(betas))
  missing <- setdiff(probes, rownames(betas))
  
  if (length(missing) > 0 && !allow_missing) {
    stop(sprintf("Missing %d required probes", length(missing)))
  }
  
  # Create output matrix
  result <- matrix(
    NA_real_,
    nrow = length(probes),
    ncol = ncol(betas),
    dimnames = list(probes, colnames(betas))
  )
  
  # Fill in available probes
  result[available, ] <- betas[available, , drop = FALSE]
  
  return(result)
}

#' Format a data frame for display
#' @param df Data frame
#' @param max_rows Maximum rows to show
#' @return Formatted string
#' @export
format_df_summary <- function(df, max_rows = 5) {
  if (nrow(df) <= max_rows) {
    return(capture.output(print(df)))
  }
  
  top <- head(df, max_rows)
  lines <- capture.output(print(top))
  c(lines, sprintf("... and %d more rows", nrow(df) - max_rows))
}

#' Create a progress bar for loops
#' @param total Total iterations
#' @param prefix Prefix string
#' @return Progress bar function
#' @export
create_progress_bar <- function(total, prefix = "Progress") {
  current <- 0
  width <- 50
  
  function(increment = 1) {
    current <<- current + increment
    pct <- current / total
    filled <- round(pct * width)
    bar <- paste0(
      strrep("█", filled),
      strrep("░", width - filled)
    )
    cat(sprintf("\r%s: [%s] %3.0f%% (%d/%d)", 
                prefix, bar, pct * 100, current, total))
    if (current >= total) cat("\n")
    flush.console()
  }
}

#' Clean probe names for consistency
#' @param probes Vector of probe names
#' @return Cleaned probe names
#' @export
clean_probe_names <- function(probes) {
  # Remove any whitespace
  probes <- trimws(probes)
  
  # Standardize case for cg/ch probes
  probes <- gsub("^CG", "cg", probes)
  probes <- gsub("^CH", "ch", probes)
  
  return(probes)
}

#' Get memory usage
#' @return Memory usage in GB
#' @export
get_memory_usage <- function() {
  gc()
  mem <- gc()
  sum(mem[, 2]) / 1024  # Convert to GB
}

#' Check available memory
#' @param required_gb Required memory in GB
#' @param verbose Print messages
#' @return Logical indicating if sufficient memory
#' @export
check_memory <- function(required_gb = 4, verbose = TRUE) {
  current_usage <- get_memory_usage()
  
  # Try to get system memory (platform dependent)
  total_mem <- tryCatch({
    if (.Platform$OS.type == "unix") {
      as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo", 
                        intern = TRUE)) / 1024 / 1024
    } else {
      # Windows - rough estimate
      8
    }
  }, error = function(e) 8)
  
  available <- total_mem - current_usage
  
  if (verbose && available < required_gb) {
    log_message(sprintf(
      "Warning: Low memory. Available: %.1f GB, Required: %.1f GB",
      available, required_gb
    ), verbose, "warning")
  }
  
  return(available >= required_gb)
}
