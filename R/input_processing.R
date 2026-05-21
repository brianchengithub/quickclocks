# ============================================================================
# Input Processing
#
# Handles IDAT file loading via SeSAMe and beta matrix preparation.
#
# When loading from a directory of IDATs, this module also:
#   * Recursively discovers IDAT files and sample sheets in the input tree
#   * Builds a pheno data.frame mapping IDAT basenames to sample metadata,
#     synthesizing minimal rows for IDATs without sample-sheet matches
#   * Stashes the resulting pheno in .qc_env$discovered_pheno so the
#     clocker() orchestrator can use it when the user didn't supply one
#
# When loading IDATs, also computes per-sample chrX/chrY total signal
# intensities from the SigDFs (the methylQC sex caller's preferred input)
# and caches them in .qc_env so SigDFs themselves can be released after
# beta extraction. Beta-only inputs skip this step; sex inference then
# falls back to a beta-proxy that uses the same downstream algorithm.
# ============================================================================


#' Process IDAT files using SeSAMe
#'
#' Pipeline: SeSAMe "QCDPB" (QC mask, channel correction, dye-bias, pOOBAH
#' detection p-values, BMIQ).
#'
#' When given a directory, the function recursively discovers IDAT files
#' AND sample sheets, then drives openSesame with the discovered IDAT
#' basenames so that subdirectory-organized layouts work correctly. The
#' resulting per-sample pheno frame (with synthesized rows for any IDATs
#' lacking a sample-sheet match) is stashed in
#' \code{.qc_env$discovered_pheno} for the orchestrator to consume.
#'
#' Single-sample handling: when a single .idat path is supplied, the
#' underlying sample basename is preserved in colnames(betas).
#'
#' @keywords internal
process_idat_files <- function(input, n_cores, verbose) {

  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("sesame package required for IDAT processing. ",
         "Install with: BiocManager::install('sesame')")
  }

  # Reset any cached pheno discovery from a previous run
  .qc_env$discovered_pheno <- NULL

  # ---- Directory input: recursive discovery + sheet matching ----
  if (length(input) == 1L && dir.exists(input)) {
    if (verbose) message("  Processing IDAT directory: ", input)

    discovery <- build_pheno_from_idat_dir(input, verbose = verbose)
    if (length(discovery$idat_basenames) == 0L) {
      stop("No IDAT files found under: ", input)
    }
    .qc_env$discovered_pheno <- discovery$pheno

    # Build the explicit list of Sentrix basenames (with directories) so we
    # can drive openSesame on the deepest level rather than letting it walk
    # the tree itself -- this lets us guarantee colnames(betas) == Sentrix
    # keys, which is what the discovered pheno is keyed on.
    idat_paths <- discover_idat_files(input)
    basemap <- idat_basenames(idat_paths)
    sentrix_bases <- file.path(basemap$dir, basemap$basename)

    if (verbose) {
      message(sprintf("  Loading %d IDAT pair(s) via SeSAMe...",
                      length(sentrix_bases)))
    }
    sdfs <- sesame::openSesame(sentrix_bases, prep = "QCDPB", func = NULL,
                                 BPPARAM = NULL)

    # Wrap single-SigDF returns in a list for uniform handling
    if (!is.list(sdfs) || (!is.null(sdfs$Probe_ID) && !is.list(sdfs[[1]]))) {
      sdfs <- list(sdfs)
    }
    # Force names = Sentrix basenames
    if (length(sdfs) == length(sentrix_bases)) {
      names(sdfs) <- basemap$basename
    } else if (is.null(names(sdfs))) {
      names(sdfs) <- paste0("Sample_", seq_along(sdfs))
    }
  } else {
    # ---- Explicit file-path input: no sheet discovery ----
    if (verbose) message("  Processing ", length(input), " IDAT file(s)...")
    paths <- input
    if (length(paths) == 1L) {
      base <- sub("_(Grn|Red)\\.idat(\\.gz)?$", "", paths, ignore.case = TRUE)
      sample_name <- basename(base)
      sdfs <- sesame::openSesame(base, prep = "QCDPB", func = NULL)
      if (is.list(sdfs)) names(sdfs) <- sample_name
    } else {
      sdfs <- sesame::openSesame(paths, prep = "QCDPB", func = NULL)
    }

    if (!is.list(sdfs) || (!is.null(sdfs$Probe_ID) && !is.list(sdfs[[1]]))) {
      sdfs <- list(sdfs)
      names(sdfs) <- "Sample_1"
    }
    if (is.null(names(sdfs))) {
      names(sdfs) <- paste0("Sample_", seq_along(sdfs))
    }
  }

  # ---- Compute sex intensities from SigDFs (methylQC preferred input) ----
  sex_intensities <- tryCatch(
    compute_sex_signals_from_sdfs(sdfs),
    error = function(e) {
      if (verbose) message("    Sex intensity calc failed: ", e$message)
      NULL
    })
  if (!is.null(sex_intensities)) {
    cache_sex_signals(sex_intensities, scale = "intensity")
    if (verbose) {
      n_ok <- sum(!is.na(sex_intensities$chrX) & !is.na(sex_intensities$chrY))
      log_msg("  Cached methylQC-style intensity signals (%d/%d samples)",
              n_ok, nrow(sex_intensities), verbose = verbose)
    }
  }

  # ---- Extract betas ----
  betas <- do.call(cbind, lapply(sdfs, sesame::getBetas))

  if (is.null(colnames(betas)) || any(colnames(betas) == "")) {
    colnames(betas) <- names(sdfs)
  }

  # Align discovered pheno to actual betas column order, if we have one
  if (!is.null(.qc_env$discovered_pheno)) {
    common <- intersect(colnames(betas), rownames(.qc_env$discovered_pheno))
    if (length(common) != ncol(betas) && verbose) {
      missing_pheno <- setdiff(colnames(betas), rownames(.qc_env$discovered_pheno))
      if (length(missing_pheno) > 0L) {
        message(sprintf(
          "    %d sample(s) in betas have no discovered pheno row: %s",
          length(missing_pheno),
          paste(utils::head(missing_pheno, 3), collapse = ", ")))
      }
    }
    # Reorder pheno; pad missing rows with synthesized entries
    final_pheno <- .qc_env$discovered_pheno[
      match(colnames(betas), rownames(.qc_env$discovered_pheno)), ,
      drop = FALSE]
    missing_idx <- which(is.na(rownames(final_pheno)))
    for (mi in missing_idx) {
      rownames(final_pheno)[mi] <- colnames(betas)[mi]
      final_pheno[mi, "Sample_Name"] <- colnames(betas)[mi]
    }
    .qc_env$discovered_pheno <- final_pheno
  }

  # SigDFs no longer needed -- free memory
  rm(sdfs)
  invisible(gc(verbose = FALSE))

  betas
}


#' Load beta values from various input formats
#'
#' Side effects when input is an IDAT directory: stashes a discovered pheno
#' frame in \code{.qc_env$discovered_pheno}, available for use by clocker()
#' when the user didn't supply one explicitly.
#'
#' @keywords internal
load_input_data <- function(input, n_cores, verbose) {

  # Reset cached state from any previous run
  .qc_env$sex_signals      <- NULL
  .qc_env$discovered_pheno <- NULL

  if (is.matrix(input) || is.data.frame(input)) {
    betas <- as.matrix(input)
  } else if (is.character(input) && length(input) >= 1L) {
    is_idat <- grepl("\\.idat(\\.gz)?$", input, ignore.case = TRUE) |
               (length(input) == 1L && dir.exists(input))
    if (any(is_idat)) {
      betas <- process_idat_files(input, n_cores, verbose)
    } else if (length(input) == 1L) {
      f <- input
      if (grepl("\\.qs2?$", f) && requireNamespace("qs2", quietly = TRUE)) {
        betas <- qs2::qs_read(f)
      } else if (grepl("\\.rds$", f, ignore.case = TRUE)) {
        betas <- readRDS(f)
      } else if (grepl("\\.csv$", f, ignore.case = TRUE)) {
        betas <- as.matrix(utils::read.csv(f, row.names = 1, check.names = FALSE))
      } else if (grepl("\\.t(ab|sv)$", f, ignore.case = TRUE)) {
        betas <- as.matrix(utils::read.table(f, header = TRUE, sep = "\t",
                                                row.names = 1, check.names = FALSE))
      } else {
        stop("Unsupported file format: ", f)
      }
    } else {
      stop("Multiple file paths only supported for .idat input")
    }
  } else {
    stop("Unsupported input type: ", class(input)[1])
  }

  if (!is.matrix(betas)) betas <- as.matrix(betas)
  betas
}


#' Apply kNN imputation (and stash diagnostic frame for the orchestrator)
#' @keywords internal
perform_knn_imputation <- function(betas,
                                     reference_path,
                                     k,
                                     zero_shot_threshold,
                                     verbose) {

  imp <- knn_impute(
    betas,
    reference_path      = reference_path,
    k                   = k,
    zero_shot_threshold = zero_shot_threshold,
    verbose             = verbose
  )
  .qc_env$imputation_info <- imp$sample_info
  imp$betas
}
