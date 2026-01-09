#!/usr/bin/env Rscript
#' ============================================================================
#' Install Dependencies for Epigenetic Clock Calculator
#' ============================================================================
#' 
#' This script installs all required packages from CRAN, Bioconductor, and GitHub.
#' Run this script before using the epigenetic clock calculator:
#' 
#'   source("install_dependencies.R")
#'
#' After installation, use the package with:
#' 
#'   source("R/calculate_clocks.R")
#'   results <- calculate_clocks("/path/to/idats")
#'   results <- calculate_clocks(my_beta_matrix)
#'

cat("============================================================\n")
cat("Epigenetic Clock Calculator - Dependency Installation\n")
cat("============================================================\n\n")

#' Helper function to install if missing
install_if_missing <- function(pkg, source = "CRAN", repo = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s from %s...\n", pkg, source))
    
    result <- tryCatch({
      if (source == "CRAN") {
        install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
      } else if (source == "Bioconductor") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
        }
        BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
      } else if (source == "GitHub") {
        if (!requireNamespace("devtools", quietly = TRUE)) {
          install.packages("devtools", repos = "https://cloud.r-project.org", quiet = TRUE)
        }
        devtools::install_github(repo, quiet = TRUE, upgrade = "never")
      }
      TRUE
    }, error = function(e) {
      cat(sprintf("  WARNING: Failed to install %s: %s\n", pkg, e$message))
      FALSE
    })
    
    if (result && requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("  SUCCESS: %s installed\n", pkg))
    } else {
      cat(sprintf("  FAILED: %s could not be installed\n", pkg))
    }
    return(result)
  } else {
    cat(sprintf("  SKIP: %s already installed\n", pkg))
    return(TRUE)
  }
}

# ============================================================================
# CRAN Packages
# ============================================================================
cat("\n--- Installing CRAN packages ---\n")

cran_packages <- c(
  "devtools",       # For GitHub installation
  "remotes",        # Alternative GitHub installation
  "parallel",       # Parallel processing (base R)
  "doParallel",     # Parallel backend
  "foreach",        # Parallel loops
  "data.table",     # Fast data manipulation
  "matrixStats",    # Matrix operations
  "tidyverse",      # Data manipulation
  "R.utils",        # Utility functions
  "yaml",           # Config file parsing
  "optparse",       # Command-line parsing
  "crayon",         # Colored console output
  "progress",       # Progress bars
  "impute",         # Imputation (may also be on Bioconductor)
  "qs2"             # Fast serialization (required for PC clocks data)
)

cran_results <- sapply(cran_packages, function(pkg) {
  install_if_missing(pkg, "CRAN")
})

# ============================================================================
# Bioconductor Packages
# ============================================================================
cat("\n--- Installing Bioconductor packages ---\n")

# Install BiocManager first
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
}

bioc_packages <- c(
  # Core SeSAMe packages
  "sesame",
  "sesameData",
  
  # Illumina array annotations
  "IlluminaHumanMethylation27kmanifest",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICv2manifest",
  
  "IlluminaHumanMethylation27kanno.ilmn12.hg19",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  
  # Other useful packages
  "minfi",
  "GenomicRanges",
  "SummarizedExperiment",
  "impute",  # KNN imputation
  
  # Cell type deconvolution
  "EpiDISH"  # RPC, CP, CBS methods with 7 and 12 cell type references
)

bioc_results <- sapply(bioc_packages, function(pkg) {
  install_if_missing(pkg, "Bioconductor")
})

# ============================================================================
# GitHub Packages
# ============================================================================
cat("\n--- Installing GitHub packages ---\n")
cat("Note: Some packages may require manual installation if automatic fails.\n\n")

github_packages <- list(
  list(pkg = "DunedinPACE", repo = "danbelsky/DunedinPACE"),
  list(pkg = "EpiMitClocks", repo = "aet21/EpiMitClocks"),
  list(pkg = "methylCIPHER", repo = "HigginsChenLab/methylCIPHER")
)

github_results <- sapply(github_packages, function(item) {
  install_if_missing(item$pkg, "GitHub", item$repo)
})

# ============================================================================
# Installation Summary
# ============================================================================
cat("\n============================================================\n")
cat("Installation Summary\n")
cat("============================================================\n\n")

# CRAN summary
cat("CRAN Packages:\n")
for (i in seq_along(cran_packages)) {
  status <- if (cran_results[i]) "OK" else "FAILED"
  cat(sprintf("  [%s] %s\n", status, cran_packages[i]))
}

# Bioconductor summary
cat("\nBioconductor Packages:\n")
for (i in seq_along(bioc_packages)) {
  status <- if (bioc_results[i]) "OK" else "FAILED"
  cat(sprintf("  [%s] %s\n", status, bioc_packages[i]))
}

# GitHub summary
cat("\nGitHub Packages:\n")
for (i in seq_along(github_packages)) {
  status <- if (github_results[i]) "OK" else "FAILED"
  cat(sprintf("  [%s] %s (%s)\n", status, github_packages[[i]]$pkg, github_packages[[i]]$repo))
}

# Overall status
all_ok <- all(c(cran_results, bioc_results, github_results))

cat("\n============================================================\n")
if (all_ok) {
  cat("All dependencies installed successfully!\n")
  cat("You can now use the epigenetic clock calculator.\n")
} else {
  cat("Some packages failed to install.\n")
  cat("Please check the errors above and try installing manually.\n")
  cat("\nCommon issues:\n")
  cat("  - GitHub rate limits: Use a personal access token\n")
  cat("  - Bioconductor version: Run BiocManager::install(version = '3.18')\n")
  cat("  - System dependencies: Install system libraries for compilation\n")
}
cat("============================================================\n")

# ============================================================================
# Download SeSAMe Data
# ============================================================================
cat("\n--- Downloading SeSAMe reference data ---\n")

if (requireNamespace("sesameData", quietly = TRUE)) {
  tryCatch({
    cat("Caching EPIC annotation data...\n")
    sesameData::sesameDataCache("EPIC")
    cat("Caching 450K annotation data...\n")
    sesameData::sesameDataCache("HM450")
    cat("SeSAMe data cached successfully.\n")
  }, error = function(e) {
    cat(sprintf("Warning: Could not cache SeSAMe data: %s\n", e$message))
    cat("Run sesameData::sesameDataCache() manually if needed.\n")
  })
} else {
  cat("SeSAMe not available for data caching.\n")
}

cat("\nDependency installation complete.\n")
