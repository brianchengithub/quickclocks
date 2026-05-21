# ============================================================================
# Sample Sheet Discovery and Ingestion
#
# Robust to many real-world layouts:
#
#   1. Flat: idat_dir/<sample sheet>.csv + idat_dir/*.idat
#   2. Per-plate subdirs: parent/plateA/sheet.csv + parent/plateA/*.idat,
#                         parent/plateB/sheet.csv + parent/plateB/*.idat
#   3. Concatenated parent sheet + per-plate IDATs:
#                         parent/all_samples.csv,
#                         parent/plateA/*.idat, parent/plateB/*.idat
#   4. Mixed: parent/all_samples.csv (some), parent/plateB/sheet.csv (rest)
#   5. No sheets anywhere -- synthesize a minimal sheet from IDAT basenames.
#
# Matching strategy:
#   For each IDAT file (basename like "206203800149_R01C01_Grn.idat"), strip
#   the channel/extension to get a "Sentrix key" ("206203800149_R01C01").
#   Then match against sample sheet rows in priority order:
#     (a) Basename column = exact Sentrix key
#     (b) Sentrix_ID + "_" + Sentrix_Position = Sentrix key
#     (c) Sample_Name == Sentrix key (rare but seen)
#     (d) Sample_Name matches IDAT basename via fuzzy lookup
#   If none match, the IDAT gets a synthesized row with sample_id = Sentrix key,
#   Age = NA, Female = NA. The user gets a single consolidated warning listing
#   the IDATs that needed synthesizing.
#
# Sample sheets are recognized as CSVs whose header row contains at least one
# of: Sample_Name, Sentrix_ID, Sentrix_Position, Basename. Illumina-style
# header preamble (the "[Data]" delimiter and metadata lines above it) is
# auto-detected and skipped.
# ============================================================================


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

#' Recursively discover IDAT files under a directory
#'
#' @param root  Directory to search.
#' @return Character vector of full paths to *_Grn.idat / *_Red.idat (with
#'   optional .gz). Pairs are not deduplicated -- caller can do that via
#'   `idat_basenames()`.
#' @keywords internal
discover_idat_files <- function(root) {
  if (!dir.exists(root)) return(character(0))
  list.files(root, pattern = "_(Grn|Red)\\.idat(\\.gz)?$",
              full.names = TRUE, recursive = TRUE,
              ignore.case = TRUE)
}


#' Reduce a list of IDAT paths to unique Sentrix basenames
#'
#' "206203800149_R01C01_Grn.idat.gz" -> "206203800149_R01C01"
#' Returns a data.frame mapping basename -> full Grn / Red file paths.
#'
#' @keywords internal
idat_basenames <- function(idat_paths) {
  if (length(idat_paths) == 0L) {
    return(data.frame(basename = character(),
                       grn = character(), red = character(),
                       dir = character(),
                       stringsAsFactors = FALSE))
  }
  # Strip channel suffix and extension
  bases <- sub("_(Grn|Red)\\.idat(\\.gz)?$", "", basename(idat_paths),
                 ignore.case = TRUE)
  # Build per-basename Grn/Red map
  ub <- unique(bases)
  out <- data.frame(basename = ub, grn = NA_character_, red = NA_character_,
                     dir = NA_character_, stringsAsFactors = FALSE)
  for (i in seq_along(ub)) {
    matches <- which(bases == ub[i])
    grn_idx <- matches[grepl("_Grn\\.idat", idat_paths[matches],
                               ignore.case = TRUE)][1]
    red_idx <- matches[grepl("_Red\\.idat", idat_paths[matches],
                               ignore.case = TRUE)][1]
    if (!is.na(grn_idx)) out$grn[i] <- idat_paths[grn_idx]
    if (!is.na(red_idx)) out$red[i] <- idat_paths[red_idx]
    out$dir[i] <- dirname(idat_paths[matches[1]])
  }
  out
}


#' Recursively discover sample-sheet CSV files
#'
#' Candidates are .csv files whose first 100 lines contain a header row
#' with at least one of Sample_Name, Sentrix_ID, Sentrix_Position, Basename.
#' Per-plate sheets in subdirectories and concatenated sheets in the parent
#' are both picked up.
#'
#' @keywords internal
discover_sample_sheets <- function(root) {
  if (!dir.exists(root)) return(character(0))
  csv_files <- list.files(root, pattern = "\\.csv$",
                            full.names = TRUE, recursive = TRUE,
                            ignore.case = TRUE)
  Filter(looks_like_sample_sheet, csv_files)
}


#' Test whether a CSV file looks like an Illumina-style sample sheet
#' @keywords internal
looks_like_sample_sheet <- function(path) {
  res <- tryCatch({
    first_lines <- readLines(path, n = 100, warn = FALSE)
    keys <- c("Sample_Name", "Sentrix_ID", "Sentrix_Position", "Basename",
              "Sentrix_Barcode")
    any(vapply(first_lines, function(L) {
      any(vapply(keys, function(k) grepl(k, L, fixed = TRUE), logical(1)))
    }, logical(1)))
  }, error = function(e) FALSE)
  isTRUE(res)
}


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

#' Read a sample-sheet CSV, auto-handling Illumina-style "[Data]" preambles
#'
#' Returns a data.frame with column names normalized (whitespace trimmed).
#' Adds a `_source_sheet` column with the file path.
#'
#' @keywords internal
read_sample_sheet <- function(path, verbose = TRUE) {
  raw <- readLines(path, warn = FALSE)
  if (length(raw) == 0L) return(NULL)

  # Locate header: first line containing any known key column.
  keys <- c("Sample_Name", "Sentrix_ID", "Sentrix_Position", "Basename")
  header_line <- NA_integer_
  for (i in seq_along(raw)) {
    if (any(vapply(keys, function(k) grepl(k, raw[i], fixed = TRUE),
                    logical(1)))) {
      header_line <- i
      break
    }
  }
  if (is.na(header_line)) {
    if (verbose) message("    Sample sheet has no recognizable header: ", path)
    return(NULL)
  }

  # Drop the preamble; first data row is the header line
  raw_trimmed <- raw[header_line:length(raw)]
  tc <- textConnection(raw_trimmed)
  on.exit(close(tc), add = TRUE)
  df <- tryCatch(
    utils::read.csv(tc, stringsAsFactors = FALSE, check.names = FALSE,
                     na.strings = c("", "NA", "N/A", "na", "n/a")),
    error = function(e) {
      if (verbose) message("    Failed to parse sample sheet: ", e$message)
      NULL
    })
  if (is.null(df) || nrow(df) == 0L) return(NULL)

  # Normalize column names
  colnames(df) <- trimws(colnames(df))
  for (col in colnames(df)) {
    if (is.character(df[[col]])) df[[col]] <- trimws(df[[col]])
  }
  df[["_source_sheet"]] <- path
  df
}


# ---------------------------------------------------------------------------
# Sentrix-key derivation
# ---------------------------------------------------------------------------

#' Build a Sentrix key column for matching
#'
#' Tries (in order):
#'   Basename
#'   Sentrix_ID + "_" + Sentrix_Position
#'   Sentrix_Barcode + "_" + Sentrix_Position
#' Returns the data.frame with a new `_sentrix_key` column. Rows that can't
#' be keyed get NA.
#'
#' @keywords internal
add_sentrix_key <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(df)

  key <- rep(NA_character_, nrow(df))

  if ("Basename" %in% colnames(df)) {
    bn <- df[["Basename"]]
    bn <- sub("_(Grn|Red)\\.idat(\\.gz)?$", "", bn, ignore.case = TRUE)
    bn <- basename(bn)
    key <- ifelse(is.na(key) & nzchar(bn) & !is.na(bn), bn, key)
  }

  sid_col <- if ("Sentrix_ID" %in% colnames(df)) "Sentrix_ID"
              else if ("Sentrix_Barcode" %in% colnames(df)) "Sentrix_Barcode"
              else NULL
  pos_col <- if ("Sentrix_Position" %in% colnames(df)) "Sentrix_Position"
              else NULL
  if (!is.null(sid_col) && !is.null(pos_col)) {
    composite <- paste0(df[[sid_col]], "_", df[[pos_col]])
    key <- ifelse(is.na(key) & !is.na(df[[sid_col]]) & !is.na(df[[pos_col]]),
                  composite, key)
  }

  # Last resort: Sample_Name itself (some labs name samples by Sentrix key)
  if ("Sample_Name" %in% colnames(df)) {
    looks_sentrix <- grepl("^[0-9]{10,}_R[0-9]{2}C[0-9]{2}$", df[["Sample_Name"]])
    key <- ifelse(is.na(key) & looks_sentrix, df[["Sample_Name"]], key)
  }

  df[["_sentrix_key"]] <- key
  df
}


# ---------------------------------------------------------------------------
# Merging
# ---------------------------------------------------------------------------

#' Merge multiple sample sheets into one, with conflict resolution
#'
#' Rules:
#'   * Sheets are concatenated by Sentrix key.
#'   * If a Sentrix key appears in multiple sheets, the **sub-directory
#'     sheet wins over the parent-directory sheet** (per-plate sheets are
#'     typically more current than concatenated globals). Among equally-deep
#'     sheets, the later-modified file wins.
#'   * Columns are unioned: NA in one sheet, value in another -> value wins.
#'   * Conflicting non-NA values produce a warning listing the keys.
#'
#' @keywords internal
merge_sample_sheets <- function(sheets, root_dir, verbose = TRUE) {
  if (length(sheets) == 0L) return(NULL)
  if (length(sheets) == 1L) return(sheets[[1]])

  # Score each sheet by directory depth (deeper = more specific = preferred)
  source_paths <- vapply(sheets, function(s) s[["_source_sheet"]][1],
                          character(1))
  depths <- vapply(source_paths, function(p) {
    rel <- sub(paste0("^", normalizePath(root_dir, mustWork = FALSE),
                       "/?"), "", normalizePath(p, mustWork = FALSE),
                fixed = FALSE)
    length(strsplit(rel, "/")[[1]])
  }, integer(1))

  mtimes <- vapply(source_paths, function(p)
    as.numeric(file.info(p)$mtime), numeric(1))

  # Sort: highest depth first, then most recent mtime
  ord <- order(-depths, -mtimes)
  sheets_ordered <- sheets[ord]

  # Union of all columns
  all_cols <- unique(unlist(lapply(sheets_ordered, colnames)))
  all_cols <- setdiff(all_cols, "_sentrix_key")   # we'll re-add at end

  # Build a master frame indexed by sentrix key, where the FIRST sheet to
  # provide a non-NA value wins (sheets_ordered already sorted preferred-first)
  out_rows <- list()
  conflicts <- character()

  for (sh in sheets_ordered) {
    keys <- sh[["_sentrix_key"]]
    valid <- !is.na(keys) & nzchar(keys)
    for (i in which(valid)) {
      k <- keys[i]
      row_data <- sh[i, , drop = FALSE]
      # Drop the key/source columns when materializing
      row_data[["_sentrix_key"]] <- NULL
      if (is.null(out_rows[[k]])) {
        out_rows[[k]] <- as.list(row_data)
      } else {
        # Fill NA columns from this row; flag conflicts
        for (col in colnames(row_data)) {
          existing <- out_rows[[k]][[col]]
          incoming <- row_data[[col]]
          if (is.null(existing) || is.na(existing) ||
              (is.character(existing) && !nzchar(existing))) {
            out_rows[[k]][[col]] <- incoming
          } else if (!is.na(incoming) && !identical(existing, incoming) &&
                     !col %in% c("_source_sheet")) {
            conflicts <- c(conflicts, sprintf("%s:%s", k, col))
          }
        }
      }
    }
  }

  if (length(conflicts) > 0L && verbose) {
    warning(sprintf(
      "Sample sheet conflicts on %d field(s); subdirectory sheets preferred. First few: %s",
      length(conflicts),
      paste(utils::head(unique(conflicts), 5), collapse = ", ")),
      call. = FALSE)
  }

  # Assemble final data.frame
  if (length(out_rows) == 0L) return(NULL)
  result <- data.frame(stringsAsFactors = FALSE)
  for (k in names(out_rows)) {
    row <- out_rows[[k]]
    row[["_sentrix_key"]] <- k
    # Coerce all NULL to NA so rbind works
    for (n in names(row)) if (is.null(row[[n]])) row[[n]] <- NA
    result <- rbind(result,
                     as.data.frame(row, stringsAsFactors = FALSE,
                                    check.names = FALSE))
  }
  result
}


# ---------------------------------------------------------------------------
# Synthesis
# ---------------------------------------------------------------------------

#' Synthesize a minimal sheet row for an IDAT lacking sample-sheet metadata
#' @keywords internal
synthesize_row <- function(sentrix_key, source = "synthesized") {
  data.frame(
    Sample_Name     = sentrix_key,
    Sentrix_ID      = sub("_R[0-9]{2}C[0-9]{2}.*$", "", sentrix_key),
    Sentrix_Position = sub("^.*_(R[0-9]{2}C[0-9]{2}).*$", "\\1", sentrix_key),
    Age             = NA_real_,
    Female          = NA_integer_,
    `_sentrix_key`  = sentrix_key,
    `_source_sheet` = source,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


# ---------------------------------------------------------------------------
# Top-level entry: discover, parse, merge, match, synthesize
# ---------------------------------------------------------------------------

#' Build a pheno data.frame for an IDAT directory tree
#'
#' Looks under `root_dir` (and all subdirectories) for IDAT files and sample
#' sheets, matches them by Sentrix key, fills in synthesized rows for any
#' IDATs without a matching sheet, and returns a pheno data.frame keyed by
#' Sentrix basename (suitable for direct use as `clocker(... pheno = ...)`).
#'
#' Discovery rules:
#'   * IDATs found recursively from `root_dir`.
#'   * Sample sheets found recursively, identified by having a header row
#'     containing any of: Sample_Name, Sentrix_ID, Sentrix_Position, Basename.
#'   * When multiple sheets cover the same Sentrix key, the **deeper**
#'     subdirectory sheet wins (more specific). Among same-depth sheets,
#'     the most recently modified wins. Conflicts emit a single warning.
#'   * IDATs with no matching sheet row get a synthesized row with
#'     Age = NA, Female = NA, and a warning listing the affected basenames.
#'
#' @param root_dir Directory to search.
#' @param verbose Print discovery progress.
#' @return List with:
#'   * `pheno`         Combined data.frame keyed by Sentrix basename
#'                     (rownames = colnames(betas) after IDAT processing)
#'   * `idat_basenames` Vector of Sentrix keys found
#'   * `n_synthesized`  How many rows were synthesized
#'   * `synthesized_keys` Which Sentrix keys were synthesized
#'   * `sheets_used`    Vector of sample-sheet paths consumed
#' @keywords internal
build_pheno_from_idat_dir <- function(root_dir, verbose = TRUE) {

  if (verbose) {
    message("  Scanning ", root_dir, " for IDAT files and sample sheets...")
  }

  # ---- IDAT discovery ----
  idat_files <- discover_idat_files(root_dir)
  if (length(idat_files) == 0L) {
    if (verbose) message("    No IDAT files found under ", root_dir)
    return(list(pheno = NULL, idat_basenames = character(),
                n_synthesized = 0L, synthesized_keys = character(),
                sheets_used = character()))
  }
  idat_map <- idat_basenames(idat_files)
  if (verbose) {
    message(sprintf("    Found %d IDAT pair(s) across %d directory(ies)",
                    nrow(idat_map),
                    length(unique(idat_map$dir))))
  }

  # ---- Sample-sheet discovery ----
  sheet_paths <- discover_sample_sheets(root_dir)
  sheets <- list()
  for (p in sheet_paths) {
    df <- tryCatch(read_sample_sheet(p, verbose = verbose),
                    error = function(e) {
                      if (verbose) message("    Could not read ", p, ": ", e$message)
                      NULL
                    })
    if (!is.null(df) && nrow(df) > 0L) {
      df <- add_sentrix_key(df)
      sheets[[length(sheets) + 1L]] <- df
      if (verbose) {
        n_keyed <- sum(!is.na(df[["_sentrix_key"]]))
        message(sprintf("    Loaded %s : %d row(s), %d keyed",
                        p, nrow(df), n_keyed))
      }
    }
  }

  # ---- Merge sheets ----
  merged <- if (length(sheets) > 0L) {
    merge_sample_sheets(sheets, root_dir, verbose = verbose)
  } else {
    NULL
  }

  # ---- Match IDATs to merged sheet ----
  pheno_rows <- list()
  synthesized <- character()

  for (i in seq_len(nrow(idat_map))) {
    key <- idat_map$basename[i]
    matched <- if (!is.null(merged)) {
      which(merged[["_sentrix_key"]] == key)
    } else {
      integer(0)
    }
    if (length(matched) == 1L) {
      pheno_rows[[key]] <- merged[matched, , drop = FALSE]
    } else if (length(matched) > 1L) {
      # Shouldn't happen post-merge, but guard anyway
      if (verbose) message("    Multiple sheet rows for ", key,
                            "; using the first.")
      pheno_rows[[key]] <- merged[matched[1], , drop = FALSE]
    } else {
      # No match: synthesize
      pheno_rows[[key]] <- synthesize_row(key, source = "synthesized:no_sheet_match")
      synthesized <- c(synthesized, key)
    }
  }

  if (length(pheno_rows) == 0L) {
    return(list(pheno = NULL, idat_basenames = idat_map$basename,
                n_synthesized = 0L, synthesized_keys = character(),
                sheets_used = sheet_paths))
  }

  # Bind to a single data frame; pad missing columns with NA.
  all_cols <- unique(unlist(lapply(pheno_rows, colnames)))
  for (i in seq_along(pheno_rows)) {
    missing_cols <- setdiff(all_cols, colnames(pheno_rows[[i]]))
    for (mc in missing_cols) pheno_rows[[i]][[mc]] <- NA
    pheno_rows[[i]] <- pheno_rows[[i]][, all_cols, drop = FALSE]
  }
  pheno_df <- do.call(rbind, pheno_rows)
  rownames(pheno_df) <- pheno_df[["_sentrix_key"]]

  # Standardize the columns clocker downstream cares about: Age, Female / Sex
  # (don't drop the other metadata - the user may want it)
  pheno_df <- standardize_age_sex_columns(pheno_df)

  # ---- User-facing warnings ----
  if (length(synthesized) > 0L) {
    sample_show <- utils::head(synthesized, 5)
    extra <- if (length(synthesized) > 5) sprintf(" (... and %d more)",
                                                    length(synthesized) - 5)
              else ""
    warning(sprintf(
      "%d IDAT pair(s) had no matching sample-sheet row; synthesized minimal pheno rows. Affected: %s%s. Age acceleration will be NA for these samples.",
      length(synthesized),
      paste(sample_show, collapse = ", "),
      extra),
      call. = FALSE)
  }

  if (verbose) {
    n_with_age <- sum(!is.na(pheno_df$Age))
    n_with_sex <- if ("Female" %in% colnames(pheno_df))
      sum(!is.na(pheno_df$Female)) else 0L
    message(sprintf(
      "    Built pheno: %d sample(s)%s (Age provided for %d, Sex provided for %d, synthesized %d)",
      nrow(pheno_df),
      if (length(sheet_paths) > 0L)
        sprintf(" from %d sheet(s)", length(sheet_paths))
      else " (all synthesized)",
      n_with_age, n_with_sex, length(synthesized)))
  }

  list(
    pheno            = pheno_df,
    idat_basenames   = idat_map$basename,
    n_synthesized    = length(synthesized),
    synthesized_keys = synthesized,
    sheets_used      = sheet_paths
  )
}


#' Normalize an Illumina-style sheet's Age/Sex columns into clocker's
#' canonical Age (numeric) and Female (0/1 integer) columns.
#'
#' Handles the common case where a merged sheet has both a `Sex` column (with
#' values like "Female"/"M") and a `Female` column that's partly or entirely
#' NA (e.g., from synthesized rows). In that case, Female is populated from
#' Sex wherever Female is NA, preserving any pre-existing Female values.
#'
#' @keywords internal
standardize_age_sex_columns <- function(pheno) {
  if (is.null(pheno) || nrow(pheno) == 0L) return(pheno)

  # Age: accept Age / age / AGE / etc.
  if (!"Age" %in% colnames(pheno)) {
    age_col <- intersect(c("age", "AGE", "Sample_Age", "subject_age"),
                          colnames(pheno))[1]
    if (!is.na(age_col)) {
      pheno$Age <- suppressWarnings(as.numeric(pheno[[age_col]]))
    }
  } else {
    pheno$Age <- suppressWarnings(as.numeric(pheno$Age))
  }

  # Female: accept Female / Sex / sex / Gender (coerce via coerce_female_indicator)
  sex_col <- intersect(c("Sex", "sex", "SEX", "Gender", "gender"),
                        colnames(pheno))[1]

  if (!"Female" %in% colnames(pheno)) {
    if (!is.na(sex_col)) {
      pheno$Female <- coerce_female_indicator(pheno[[sex_col]],
                                                 na_on_unknown = TRUE)
    }
  } else {
    # Female exists. If any rows have NA Female but non-NA Sex, fill from Sex.
    needs_fill <- is.na(pheno$Female)
    if (any(needs_fill) && !is.na(sex_col)) {
      derived <- coerce_female_indicator(pheno[[sex_col]],
                                            na_on_unknown = TRUE)
      pheno$Female[needs_fill] <- derived[needs_fill]
    }
  }

  pheno
}
