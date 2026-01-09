#' ============================================================================
#' Input Handler for Epigenetic Clock Calculator
#' ============================================================================
#' 
#' Functions for loading IDAT files and beta value matrices from various sources.

#' Load IDAT files from a directory
#' 
#' Recursively searches for IDAT files and processes them using SeSAMe's
#' openSesame() function with specified preprocessing options.
#' 
#' @param idat_dir Path to directory containing IDAT files
#' @param sample_sheet Optional path to sample sheet CSV
#' @param config Configuration list
#' @return List with beta matrix, sample IDs, and platform info
#' @export
load_idat_files <- function(idat_dir, sample_sheet = NULL, config = list()) {
  
  if (!dir.exists(idat_dir)) {
    stop("IDAT directory does not exist: ", idat_dir)
  }
  
  verbose <- config$verbose %||% TRUE
  
  # Load required packages
  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("Package 'sesame' is required. Install with: BiocManager::install('sesame')")
  }
  
  log_message("Searching for IDAT files...", verbose)
  
  # Find all IDAT files recursively
  idat_files <- find_idat_files(idat_dir)
  
  if (length(idat_files$grn) == 0) {
    stop("No IDAT files found in: ", idat_dir)
  }
  
  log_message(sprintf("  Found %d IDAT pairs", length(idat_files$grn)), verbose)
  
  # Get sample IDs from file names
  sample_ids <- extract_sample_ids(idat_files$grn)
  
  # If sample sheet provided, use it to get sample info
  if (!is.null(sample_sheet) && file.exists(sample_sheet)) {
    log_message("Loading sample sheet...", verbose)
    sheet_info <- load_sample_sheet(sample_sheet)
    
    # Match samples
    matched <- match(sample_ids, sheet_info$sample_id)
    if (any(!is.na(matched))) {
      log_message(sprintf("  Matched %d samples from sheet", 
                          sum(!is.na(matched))), verbose)
    }
  }
  
  # Process IDAT files with SeSAMe
  log_message("Processing IDAT files with SeSAMe...", verbose)
  
  betas <- process_idats_sesame(
    idat_prefixes = idat_files$prefixes,
    prep_method = config$prep_method %||% "QCDPB",
    n_cores = config$n_cores %||% 1,
    verbose = verbose
  )
  
  # Detect platform
  platform <- detect_platform(rownames(betas))
  
  return(list(
    betas = betas,
    sample_ids = sample_ids,
    platform = platform,
    idat_files = idat_files
  ))
}


#' Find IDAT files recursively in a directory
#' 
#' @param dir Directory to search
#' @return List with green files, red files, and prefixes
#' @export
find_idat_files <- function(dir) {
  # Find all IDAT files
  all_idats <- list.files(
    dir,
    pattern = "\\.(idat|IDAT)(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(all_idats) == 0) {
    return(list(grn = character(), red = character(), prefixes = character()))
  }
  
  # Separate green and red channels
  grn_pattern <- "_Grn\\.(idat|IDAT)(\\.gz)?$"
  red_pattern <- "_Red\\.(idat|IDAT)(\\.gz)?$"
  
  grn_files <- grep(grn_pattern, all_idats, value = TRUE)
  red_files <- grep(red_pattern, all_idats, value = TRUE)
  
  # Extract prefixes (everything before _Grn/_Red)
  grn_prefixes <- gsub(grn_pattern, "", grn_files)
  red_prefixes <- gsub(red_pattern, "", red_files)
  
  # Find matching pairs
  common_prefixes <- intersect(grn_prefixes, red_prefixes)
  
  if (length(common_prefixes) < length(grn_prefixes)) {
    n_unmatched <- length(grn_prefixes) - length(common_prefixes)
    warning(sprintf("%d IDAT files without matching pair", n_unmatched))
  }
  
  # Get matched files
  matched_grn <- grn_files[grn_prefixes %in% common_prefixes]
  matched_red <- red_files[red_prefixes %in% common_prefixes]
  
  list(
    grn = matched_grn,
    red = matched_red,
    prefixes = common_prefixes
  )
}


#' Extract sample IDs from IDAT file paths
#' 
#' @param idat_files Vector of IDAT file paths
#' @return Character vector of sample IDs
#' @export
extract_sample_ids <- function(idat_files) {
  # Get base names without path
  basenames <- basename(idat_files)
  
  # Remove _Grn/_Red suffix and .idat extension
  sample_ids <- gsub("_(Grn|Red)\\.(idat|IDAT)(\\.gz)?$", "", basenames)
  
  # Clean up any duplicate underscores or weird characters
  sample_ids <- gsub("_+", "_", sample_ids)
  sample_ids <- gsub("^_|_$", "", sample_ids)
  
  # Check for duplicates
  if (anyDuplicated(sample_ids)) {
    # Make unique by adding index
    dup_idx <- which(duplicated(sample_ids) | duplicated(sample_ids, fromLast = TRUE))
    sample_ids[dup_idx] <- make.unique(sample_ids[dup_idx], sep = "_")
    warning("Duplicate sample IDs detected and made unique")
  }
  
  return(sample_ids)
}


#' Load and parse a sample sheet
#' 
#' @param path Path to sample sheet CSV
#' @return Data frame with sample information
#' @export
load_sample_sheet <- function(path) {
  if (!file.exists(path)) {
    stop("Sample sheet not found: ", path)
  }
  
  # Try to read the file
  sheet <- tryCatch({
    # First try standard CSV
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }, error = function(e) {
    # Try with semicolon separator (European format)
    read.csv2(path, stringsAsFactors = FALSE, check.names = FALSE)
  })
  
  # Standardize column names
  colnames(sheet) <- tolower(colnames(sheet))
  
  # Look for sample ID column
  id_cols <- c("sample_id", "sampleid", "sample", "barcode", "sentrix_id", 
               "basename", "sample_name")
  
  found_id_col <- intersect(colnames(sheet), id_cols)
  
  if (length(found_id_col) == 0) {
    # Use first column as ID
    warning("No standard sample ID column found. Using first column.")
    colnames(sheet)[1] <- "sample_id"
  } else {
    # Rename found column to sample_id
    colnames(sheet)[colnames(sheet) == found_id_col[1]] <- "sample_id"
  }
  
  return(sheet)
}


#' Process IDAT files using SeSAMe
#' 
#' @param idat_prefixes Vector of IDAT file prefixes
#' @param prep_method SeSAMe preprocessing method string
#' @param n_cores Number of parallel cores
#' @param verbose Print progress messages
#' @return Matrix of beta values
#' @export
process_idats_sesame <- function(idat_prefixes, prep_method = "QCDPB", 
                                  n_cores = 1, verbose = TRUE) {
  
  library(sesame)
  
  # Cache SeSAMe data if not already done
  tryCatch({
    sesameDataCache()
  }, error = function(e) {
    log_message("Note: Could not cache sesameData. Continuing...", verbose, "warning")
  })
  
  n_samples <- length(idat_prefixes)
  
  if (n_cores > 1 && n_samples > 1) {
    # Parallel processing
    log_message(sprintf("Processing %d samples in parallel (%d cores)...", 
                        n_samples, n_cores), verbose)
    
    if (requireNamespace("parallel", quietly = TRUE)) {
      # Setup parallel backend
      cl <- parallel::makeCluster(min(n_cores, n_samples))
      on.exit(parallel::stopCluster(cl))
      
      # Export required packages to workers
      parallel::clusterEvalQ(cl, library(sesame))
      
      # Process in parallel
      beta_list <- parallel::parLapply(cl, idat_prefixes, function(prefix) {
        tryCatch({
          betas <- openSesame(prefix, prep = prep_method)
          return(betas)
        }, error = function(e) {
          return(NULL)
        })
      })
    } else {
      warning("Parallel processing requested but 'parallel' package not available")
      n_cores <- 1
    }
  }
  
  if (n_cores == 1) {
    # Sequential processing with progress
    log_message(sprintf("Processing %d samples sequentially...", n_samples), verbose)
    
    beta_list <- vector("list", n_samples)
    
    for (i in seq_along(idat_prefixes)) {
      if (verbose && i %% 10 == 0) {
        log_message(sprintf("  Processed %d / %d samples", i, n_samples), verbose)
      }
      
      beta_list[[i]] <- tryCatch({
        openSesame(idat_prefixes[i], prep = prep_method)
      }, error = function(e) {
        warning(sprintf("Failed to process sample %d: %s", i, e$message))
        NULL
      })
    }
  }
  
  # Remove failed samples
  failed <- sapply(beta_list, is.null)
  if (any(failed)) {
    warning(sprintf("%d samples failed processing and were removed", sum(failed)))
    beta_list <- beta_list[!failed]
    idat_prefixes <- idat_prefixes[!failed]
  }
  
  if (length(beta_list) == 0) {
    stop("All samples failed processing")
  }
  
  # Combine into matrix
  log_message("Combining samples into matrix...", verbose)
  
  # Get all unique probe names
  all_probes <- unique(unlist(lapply(beta_list, names)))
  
  # Create output matrix
  betas <- matrix(
    NA_real_,
    nrow = length(all_probes),
    ncol = length(beta_list),
    dimnames = list(
      all_probes,
      extract_sample_ids(paste0(idat_prefixes, "_Grn.idat"))
    )
  )
  
  # Fill in values
  for (i in seq_along(beta_list)) {
    b <- beta_list[[i]]
    betas[names(b), i] <- b
  }
  
  log_message(sprintf("  Final matrix: %d probes × %d samples", 
                      nrow(betas), ncol(betas)), verbose)
  
  return(betas)
}


#' Load a pre-existing beta value matrix
#' 
#' @param path Path to RDS or CSV file containing beta matrix
#' @return Matrix of beta values
#' @export
load_beta_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("Beta matrix file not found: ", path)
  }
  
  ext <- tolower(tools::file_ext(path))
  
  betas <- switch(ext,
    "rds" = readRDS(path),
    "csv" = {
      df <- read.csv(path, row.names = 1, check.names = FALSE)
      as.matrix(df)
    },
    "txt" = {
      df <- read.delim(path, row.names = 1, check.names = FALSE)
      as.matrix(df)
    },
    "tsv" = {
      df <- read.delim(path, row.names = 1, check.names = FALSE)
      as.matrix(df)
    },
    stop("Unsupported file format: ", ext)
  )
  
  return(validate_beta_matrix(betas))
}


#' Load reference beta values for imputation
#' 
#' @param path Path to RDS file containing reference betas
#' @return Named numeric vector of reference beta values
#' @export
load_reference_betas <- function(path) {
  if (!file.exists(path)) {
    stop("Reference betas file not found: ", path)
  }
  
  ref <- readRDS(path)
  
  # Validate structure
  if (!is.numeric(ref)) {
    stop("Reference betas must be numeric")
  }
  
  if (is.null(names(ref))) {
    stop("Reference betas must be named (probe IDs)")
  }
  
  # Clean probe names
  names(ref) <- clean_probe_names(names(ref))
  
  # Check for valid beta values
  if (any(ref < 0 | ref > 1, na.rm = TRUE)) {
    warning("Some reference beta values outside [0, 1] range")
    ref <- pmax(0, pmin(1, ref))
  }
  
  return(ref)
}


#' Helper: null coalescing operator
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a
