#!/usr/bin/env Rscript
#' Diagnostic script to check epigenetic clock packages
#' Run this to understand what's available and working

cat("=== 'QUICKCLOCKS' PACKAGE DIAGNOSTICS ===\n\n")

# ============================================================================
# 1. Check SeSAMe
# ============================================================================
cat("--- SeSAMe ---\n")
if (requireNamespace("sesame", quietly = TRUE)) {
  cat("Installed: YES\n")
  cat("Version: ", as.character(packageVersion("sesame")), "\n")
  
  # List cell composition functions
  sesame_ns <- asNamespace("sesame")
  all_funcs <- ls(sesame_ns)
  
  cell_funcs <- all_funcs[grepl("cell|leuk|compos", all_funcs, ignore.case = TRUE)]
  cat("Cell composition functions: ", paste(cell_funcs, collapse = ", "), "\n")
  
  sex_funcs <- all_funcs[grepl("sex|gender", all_funcs, ignore.case = TRUE)]
  cat("Sex inference functions: ", paste(sex_funcs, collapse = ", "), "\n")
  
  # Check function signatures
  if ("estimateCellComposition" %in% all_funcs) {
    cat("\nestimateCellComposition signature:\n")
    print(args(sesame::estimateCellComposition))
  }
  if ("inferSex" %in% all_funcs) {
    cat("\ninferSex signature:\n")
    print(args(sesame::inferSex))
  }
} else {
  cat("Installed: NO\n")
}

# ============================================================================
# 2. Check EpiMitClocks
# ============================================================================
cat("\n--- EpiMitClocks ---\n")
if (requireNamespace("EpiMitClocks", quietly = TRUE)) {
  cat("Installed: YES\n")
  cat("Version: ", as.character(packageVersion("EpiMitClocks")), "\n")
  
  # Check available data
  cat("\nAvailable datasets:\n")
  pkg_data <- data(package = "EpiMitClocks")
  if (nrow(pkg_data$results) > 0) {
    print(pkg_data$results[, c("Item", "Title")])
  } else {
    cat("  No datasets found!\n")
  }
  
  # List functions
  epi_ns <- asNamespace("EpiMitClocks")
  all_funcs <- ls(epi_ns)
  cat("\nAvailable functions: ", paste(all_funcs[1:min(20, length(all_funcs))], collapse = ", "), "\n")
  
  # Try loading data
  cat("\nTrying to load dataETOC2:\n")
  tryCatch({
    data("dataETOC2", package = "EpiMitClocks")
    cat("  SUCCESS - dataETOC2 loaded\n")
    if (exists("dataETOC2.l")) cat("  dataETOC2.l exists\n")
  }, error = function(e) cat("  FAILED:", e$message, "\n"),
     warning = function(w) cat("  WARNING:", w$message, "\n"))
  
  cat("\nTrying to load dataETOC3:\n")
  tryCatch({
    data("dataETOC3", package = "EpiMitClocks")
    cat("  SUCCESS - dataETOC3 loaded\n")
    if (exists("dataETOC3.l")) cat("  dataETOC3.l exists\n")
  }, error = function(e) cat("  FAILED:", e$message, "\n"),
     warning = function(w) cat("  WARNING:", w$message, "\n"))
  
  # Check main function signature
  if ("EpiMitClocks" %in% all_funcs) {
    cat("\nEpiMitClocks function signature:\n")
    print(args(EpiMitClocks::EpiMitClocks))
  }
  if ("epiTOC2" %in% all_funcs) {
    cat("\nepiTOC2 function signature:\n")
    print(args(EpiMitClocks::epiTOC2))
  }
} else {
  cat("Installed: NO\n")
}

# ============================================================================
# 3. Check methylCIPHER
# ============================================================================
cat("\n--- methylCIPHER ---\n")
if (requireNamespace("methylCIPHER", quietly = TRUE)) {
  cat("Installed: YES\n")
  cat("Version: ", as.character(packageVersion("methylCIPHER")), "\n")
  
  # List all calc* functions
  mc_ns <- asNamespace("methylCIPHER")
  all_funcs <- ls(mc_ns)
  calc_funcs <- all_funcs[grepl("^calc", all_funcs)]
  cat("\nAvailable calc* functions (", length(calc_funcs), "):\n")
  cat(paste(calc_funcs, collapse = ", "), "\n")
  
  # Check specific function signatures
  test_funcs <- c("calcHorvath1", "calcPhenoAge", "calcPCHorvath1", "calcCoreClocks")
  for (fn in test_funcs) {
    if (fn %in% all_funcs) {
      cat(sprintf("\n%s signature:\n", fn))
      print(args(get(fn, envir = mc_ns)))
    }
  }
  
  # Check available data
  cat("\nAvailable datasets:\n")
  pkg_data <- data(package = "methylCIPHER")
  if (nrow(pkg_data$results) > 0) {
    print(pkg_data$results[, c("Item", "Title")])
  } else {
    cat("  No datasets found - this might be a problem!\n")
  }
} else {
  cat("Installed: NO\n")
}

# ============================================================================
# 4. Check DunedinPACE
# ============================================================================
cat("\n--- DunedinPACE ---\n")
if (requireNamespace("DunedinPACE", quietly = TRUE)) {
  cat("Installed: YES\n")
  cat("Version: ", as.character(packageVersion("DunedinPACE")), "\n")
  
  dp_ns <- asNamespace("DunedinPACE")
  all_funcs <- ls(dp_ns)
  cat("Available functions: ", paste(all_funcs, collapse = ", "), "\n")
  
  if ("PACEProjector" %in% all_funcs) {
    cat("\nPACEProjector signature:\n")
    print(args(DunedinPACE::PACEProjector))
  }
} else {
  cat("Installed: NO\n")
}

# ============================================================================
# 5. Quick test with sample data
# ============================================================================
cat("\n--- Quick Function Tests ---\n")

# Create tiny test matrix
set.seed(42)
test_betas <- matrix(runif(1000 * 5, 0, 1), nrow = 1000, ncol = 5)
rownames(test_betas) <- paste0("cg", sprintf("%08d", 1:1000))
colnames(test_betas) <- paste0("Sample", 1:5)

cat("\nTest matrix: ", nrow(test_betas), " CpGs x ", ncol(test_betas), " samples\n")

# Test methylCIPHER
if (requireNamespace("methylCIPHER", quietly = TRUE)) {
  cat("\nTesting methylCIPHER::calcHorvath1...\n")
  
  # First load the required data
  tryCatch({
    data("Horvath1_CpGs", package = "methylCIPHER", envir = environment())
    cat("  Loaded Horvath1_CpGs\n")
  }, error = function(e) cat("  Failed to load Horvath1_CpGs:", e$message, "\n"),
     warning = function(w) cat("  Warning loading Horvath1_CpGs:", w$message, "\n"))
  
  tryCatch({
    # Transpose: methylCIPHER wants samples as rows
    test_t <- t(test_betas)
    result <- methylCIPHER::calcHorvath1(test_t, imputation = FALSE)
    cat("  Result type: ", class(result)[1], "\n")
    cat("  Result length/dim: ", 
        if(is.null(dim(result))) length(result) else paste(dim(result), collapse="x"), "\n")
    if (is.numeric(result)) cat("  Values: ", paste(head(result), collapse=", "), "\n")
  }, error = function(e) cat("  ERROR:", e$message, "\n"))
}

# Test EpiMitClocks
if (requireNamespace("EpiMitClocks", quietly = TRUE)) {
  cat("\nTesting EpiMitClocks main function...\n")
  
  # Load required data first
  tryCatch({
    data("dataETOC3", package = "EpiMitClocks", envir = environment())
    cat("  Loaded dataETOC3 successfully\n")
  }, error = function(e) cat("  Failed to load dataETOC3:", e$message, "\n"),
     warning = function(w) cat("  Warning:", w$message, "\n"))
  
  # Test main function with CpGs as rows (correct format)
  tryCatch({
    result <- EpiMitClocks::EpiMitClocks(data.m = test_betas, ages.v = NULL)
    cat("  Result type: ", class(result)[1], "\n")
    if (is.data.frame(result)) {
      cat("  Result dims: ", paste(dim(result), collapse="x"), "\n")
      cat("  Columns: ", paste(colnames(result), collapse=", "), "\n")
    }
  }, error = function(e) cat("  ERROR:", e$message, "\n"))
  
  cat("\nTesting EpiMitClocks::epiTOC2 with data.m parameter...\n")
  tryCatch({
    result <- EpiMitClocks::epiTOC2(data.m = test_betas)
    cat("  Result type: ", class(result)[1], "\n")
    if (is.list(result)) cat("  List names: ", paste(names(result), collapse=", "), "\n")
  }, error = function(e) cat("  ERROR:", e$message, "\n"))
}

cat("\n=== DIAGNOSTICS COMPLETE ===\n")
cat("\nPlease share this output so we can fix the issues.\n")
