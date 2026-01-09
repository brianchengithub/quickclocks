#' ============================================================================
#' Epigenetic Clock Calculator - Package Entry Point
#' ============================================================================
#' 
#' This file provides the main entry point for the package.
#' The primary function is calculate_clocks() which is defined in 
#' R/calculate_clocks.R
#' 
#' @examples
#' # From IDAT directory
#' results <- calculate_clocks("/path/to/idat/directory")
#' 
#' # From beta matrix  
#' results <- calculate_clocks(my_beta_matrix)
#' 
#' # With specific cores
#' results <- calculate_clocks(my_beta_matrix, n_cores = 4)
#' 
#' @docType package
#' @name EpigeneticClockCalculator
NULL

# Source all modules when file is sourced directly (not as installed package)
.onLoad <- function(libname, pkgname) {
  # Package is loaded - modules are already available
}

# For sourcing directly (development mode)
if (sys.nframe() == 0 || !isNamespaceLoaded("EpigeneticClockCalculator")) {
  
  # Find package root
  pkg_root <- getwd()
  
  if (file.exists(file.path(pkg_root, "R", "calculate_clocks.R"))) {
    # Source all R files
    r_files <- list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE)
    
    # Source utils first (contains helper functions)
    utils_file <- file.path(pkg_root, "R", "utils.R")
    if (file.exists(utils_file)) {
      source(utils_file, local = FALSE)
    }
    
    # Source remaining files
    for (f in r_files) {
      if (!grepl("utils\\.R$", f)) {
        tryCatch(
          source(f, local = FALSE),
          error = function(e) message("Note: Could not source ", basename(f))
        )
      }
    }
  }
}
