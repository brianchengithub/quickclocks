#' ============================================================================
#' Output Handler - Results Formatting and Export
#' ============================================================================
#' 
#' Functions for formatting and saving clock calculation results
#'

#' Save all results to files
#' 
#' @param results List of clock results
#' @param output_dir Output directory
#' @param output_format Output format: "csv", "rds", or "both"
#' @param prefix File name prefix
#' @param verbose Print messages
#' @export
save_all_results <- function(results, output_dir, output_format = "both", 
                             prefix = "", verbose = TRUE) {
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    log_message(sprintf("Created output directory: %s", output_dir), verbose)
  }
  
  # File prefix
  if (nchar(prefix) > 0 && !endsWith(prefix, "_")) {
    prefix <- paste0(prefix, "_")
  }
  
  saved_files <- character()
  
  # =========================================================================
  # Combined results
  # =========================================================================
  if (!is.null(results$combined)) {
    log_message("Saving combined clock results...", verbose)
    
    combined_file <- save_dataframe(
      results$combined,
      output_dir,
      paste0(prefix, "clock_results"),
      output_format,
      verbose
    )
    saved_files <- c(saved_files, combined_file)
  }
  
  # =========================================================================
  # Cell composition
  # =========================================================================
  if (!is.null(results$cell_composition)) {
    log_message("Saving cell composition...", verbose)
    
    cc_file <- save_dataframe(
      results$cell_composition,
      output_dir,
      paste0(prefix, "cell_composition"),
      output_format,
      verbose
    )
    saved_files <- c(saved_files, cc_file)
  }
  
  # =========================================================================
  # Sex inference
  # =========================================================================
  if (!is.null(results$sex)) {
    log_message("Saving sex inference...", verbose)
    
    sex_file <- save_dataframe(
      results$sex,
      output_dir,
      paste0(prefix, "sex_inference"),
      output_format,
      verbose
    )
    saved_files <- c(saved_files, sex_file)
  }
  
  # =========================================================================
  # Individual clock results
  # =========================================================================
  clock_sets <- c("dunedin_pace", "pc_clocks", "mitotic_clocks", "methylcipher", "sesame_age")
  
  for (clock_set in clock_sets) {
    if (!is.null(results[[clock_set]]) && is.data.frame(results[[clock_set]])) {
      log_message(sprintf("Saving %s results...", clock_set), verbose)
      
      file <- save_dataframe(
        results[[clock_set]],
        output_dir,
        paste0(prefix, clock_set),
        output_format,
        verbose
      )
      saved_files <- c(saved_files, file)
    }
  }
  
  # =========================================================================
  # Clock availability report
  # =========================================================================
  if (!is.null(results$availability)) {
    log_message("Saving availability report...", verbose)
    
    avail_file <- file.path(output_dir, paste0(prefix, "clock_availability.txt"))
    write_clock_availability_report(results$availability, avail_file)
    saved_files <- c(saved_files, avail_file)
  }
  
  # =========================================================================
  # Summary
  # =========================================================================
  if (!is.null(results$summary)) {
    log_message("Saving calculation summary...", verbose)
    
    summary_file <- save_dataframe(
      results$summary,
      output_dir,
      paste0(prefix, "summary"),
      output_format,
      verbose
    )
    saved_files <- c(saved_files, summary_file)
  }
  
  # =========================================================================
  # Full results object (RDS only)
  # =========================================================================
  if (output_format %in% c("rds", "both")) {
    full_file <- file.path(output_dir, paste0(prefix, "full_results.rds"))
    saveRDS(results, full_file)
    saved_files <- c(saved_files, full_file)
    log_message(sprintf("  Saved: %s", full_file), verbose)
  }
  
  log_message(sprintf("\nSaved %d output files to: %s", 
                      length(saved_files), output_dir), verbose)
  
  return(saved_files)
}


#' Save a data frame to file(s)
#' 
#' @param df Data frame
#' @param output_dir Output directory
#' @param basename File base name (without extension)
#' @param format Format: "csv", "rds", or "both"
#' @param verbose Print messages
#' @return Vector of saved file paths
#' @keywords internal
save_dataframe <- function(df, output_dir, basename, format = "both", verbose = TRUE) {
  
  saved <- character()
  
  if (format %in% c("csv", "both")) {
    csv_file <- file.path(output_dir, paste0(basename, ".csv"))
    write.csv(df, csv_file, row.names = FALSE)
    saved <- c(saved, csv_file)
    if (verbose) log_message(sprintf("  Saved: %s", csv_file), verbose)
  }
  
  if (format %in% c("rds", "both")) {
    rds_file <- file.path(output_dir, paste0(basename, ".rds"))
    saveRDS(df, rds_file)
    saved <- c(saved, rds_file)
    if (verbose) log_message(sprintf("  Saved: %s", rds_file), verbose)
  }
  
  return(saved)
}


#' Format results as wide table (samples x clocks)
#' 
#' @param results Results list from calculate_all_clocks
#' @return Data frame in wide format
#' @export
format_results_wide <- function(results) {
  
  if (!is.null(results$combined)) {
    return(results$combined)
  }
  
  # Manually combine if combined not available
  return(combine_clock_results(results))
}


#' Format results as long table (one row per sample-clock pair)
#' 
#' @param results Results list from calculate_all_clocks
#' @return Data frame in long format
#' @export
format_results_long <- function(results) {
  
  wide <- format_results_wide(results)
  
  if (is.null(wide)) {
    return(NULL)
  }
  
  # Convert to long format
  clock_cols <- setdiff(colnames(wide), "sample_id")
  
  long_list <- lapply(clock_cols, function(clock) {
    data.frame(
      sample_id = wide$sample_id,
      clock = clock,
      value = wide[[clock]],
      stringsAsFactors = FALSE
    )
  })
  
  long_df <- do.call(rbind, long_list)
  rownames(long_df) <- NULL
  
  return(long_df)
}


#' Create HTML report of results
#' 
#' @param results Results list
#' @param output_file Output HTML file path
#' @param title Report title
#' @param verbose Print messages
#' @export
create_html_report <- function(results, output_file, 
                               title = "Epigenetic Clock Results", 
                               verbose = TRUE) {
  
  log_message("Generating HTML report...", verbose)
  
  # Start HTML
  html <- sprintf('<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>%s</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 20px; }
    h1 { color: #2c3e50; }
    h2 { color: #34495e; border-bottom: 1px solid #bdc3c7; padding-bottom: 5px; }
    table { border-collapse: collapse; margin: 15px 0; }
    th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
    th { background-color: #3498db; color: white; }
    tr:nth-child(even) { background-color: #f2f2f2; }
    .success { color: #27ae60; }
    .failed { color: #e74c3c; }
    .summary-box { background: #ecf0f1; padding: 15px; margin: 15px 0; border-radius: 5px; }
  </style>
</head>
<body>
<h1>%s</h1>
<p>Generated: %s</p>
', title, title, format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  
  # Summary section
  html <- paste0(html, '<div class="summary-box">\n')
  html <- paste0(html, '<h2>Summary</h2>\n')
  
  if (!is.null(results$combined)) {
    n_samples <- nrow(results$combined)
    n_clocks <- ncol(results$combined) - 1  # Exclude sample_id
    html <- paste0(html, sprintf('<p>Samples analyzed: %d</p>\n', n_samples))
    html <- paste0(html, sprintf('<p>Clocks calculated: %d</p>\n', n_clocks))
  }
  
  html <- paste0(html, '</div>\n')
  
  # Clock summary table
  if (!is.null(results$summary)) {
    html <- paste0(html, '<h2>Clock Summary</h2>\n')
    html <- paste0(html, '<table>\n')
    html <- paste0(html, '<tr><th>Category</th><th>Clock</th><th>Status</th></tr>\n')
    
    for (i in seq_len(nrow(results$summary))) {
      row <- results$summary[i, ]
      status_class <- if (row$status == "Success") "success" else "failed"
      html <- paste0(html, sprintf(
        '<tr><td>%s</td><td>%s</td><td class="%s">%s</td></tr>\n',
        row$category, row$clock, status_class, row$status
      ))
    }
    html <- paste0(html, '</table>\n')
  }
  
  # Results table (first 20 rows)
  if (!is.null(results$combined)) {
    html <- paste0(html, '<h2>Results Preview (First 20 Samples)</h2>\n')
    
    preview <- results$combined[1:min(20, nrow(results$combined)), 1:min(10, ncol(results$combined))]
    
    html <- paste0(html, '<table>\n')
    html <- paste0(html, '<tr>')
    for (col in colnames(preview)) {
      html <- paste0(html, sprintf('<th>%s</th>', col))
    }
    html <- paste0(html, '</tr>\n')
    
    for (i in seq_len(nrow(preview))) {
      html <- paste0(html, '<tr>')
      for (col in colnames(preview)) {
        val <- preview[i, col]
        if (is.numeric(val)) val <- round(val, 4)
        html <- paste0(html, sprintf('<td>%s</td>', val))
      }
      html <- paste0(html, '</tr>\n')
    }
    html <- paste0(html, '</table>\n')
    
    if (ncol(results$combined) > 10) {
      html <- paste0(html, sprintf('<p><em>Showing 10 of %d columns</em></p>\n', 
                                   ncol(results$combined)))
    }
  }
  
  # Close HTML
  html <- paste0(html, '</body>\n</html>')
  
  # Write file
  writeLines(html, output_file)
  log_message(sprintf("  Saved: %s", output_file), verbose)
  
  return(output_file)
}


#' Export results for use with other tools
#' 
#' @param results Results list
#' @param format Export format: "minfi", "sesame", "generic"
#' @param output_file Output file path
#' @export
export_for_tool <- function(results, format = "generic", output_file) {
  
  combined <- format_results_wide(results)
  
  if (is.null(combined)) {
    warning("No results to export")
    return(NULL)
  }
  
  if (format == "minfi") {
    # Format compatible with minfi phenotype data
    phenodata <- data.frame(
      Sample_Name = combined$sample_id,
      stringsAsFactors = FALSE
    )
    for (col in setdiff(colnames(combined), "sample_id")) {
      phenodata[[col]] <- combined[[col]]
    }
    write.csv(phenodata, output_file, row.names = FALSE)
    
  } else if (format == "sesame") {
    # Format compatible with SeSAMe
    rownames(combined) <- combined$sample_id
    combined$sample_id <- NULL
    saveRDS(combined, output_file)
    
  } else {
    # Generic format
    write.csv(combined, output_file, row.names = FALSE)
  }
  
  return(output_file)
}


#' Print results summary to console
#' 
#' @param results Results list
#' @export
print_results_summary <- function(results) {
  
  cat("\n========================================\n")
  cat("EPIGENETIC CLOCK RESULTS SUMMARY\n")
  cat("========================================\n\n")
  
  if (!is.null(results$combined)) {
    cat(sprintf("Samples analyzed: %d\n", nrow(results$combined)))
    cat(sprintf("Total clocks: %d\n\n", ncol(results$combined) - 1))
  }
  
  # Print by category
  if (!is.null(results$summary)) {
    categories <- unique(results$summary$category)
    
    for (cat in categories) {
      cat_results <- results$summary[results$summary$category == cat, ]
      cat(sprintf("%s:\n", cat))
      for (i in seq_len(nrow(cat_results))) {
        status_symbol <- if (cat_results$status[i] == "Success") "✓" else "✗"
        cat(sprintf("  [%s] %s\n", status_symbol, cat_results$clock[i]))
      }
      cat("\n")
    }
  }
  
  # Print statistics for main clocks
  if (!is.null(results$combined)) {
    cat("Clock Statistics:\n")
    cat("-----------------\n")
    
    clock_cols <- setdiff(colnames(results$combined), "sample_id")
    
    for (clock in clock_cols[1:min(10, length(clock_cols))]) {
      vals <- results$combined[[clock]]
      if (is.numeric(vals)) {
        cat(sprintf("  %s: mean=%.2f, sd=%.2f, range=[%.2f, %.2f]\n",
                    clock,
                    mean(vals, na.rm = TRUE),
                    sd(vals, na.rm = TRUE),
                    min(vals, na.rm = TRUE),
                    max(vals, na.rm = TRUE)))
      }
    }
    
    if (length(clock_cols) > 10) {
      cat(sprintf("  ... and %d more clocks\n", length(clock_cols) - 10))
    }
  }
  
  cat("\n========================================\n")
}


#' Generate QC report for clock calculations
#' 
#' @param results Results list
#' @param betas Original beta matrix
#' @param output_file Output file path
#' @export
generate_qc_report <- function(results, betas = NULL, output_file) {
  
  qc_data <- list()
  
  # Sample statistics
  if (!is.null(results$combined)) {
    qc_data$n_samples <- nrow(results$combined)
    qc_data$n_clocks <- ncol(results$combined) - 1
    
    # Missing values per clock
    clock_cols <- setdiff(colnames(results$combined), "sample_id")
    qc_data$missing_per_clock <- sapply(clock_cols, function(col) {
      sum(is.na(results$combined[[col]]))
    })
  }
  
  # Probe statistics (if betas provided)
  if (!is.null(betas)) {
    qc_data$n_probes <- nrow(betas)
    qc_data$missing_probes <- sum(is.na(betas))
    qc_data$pct_missing <- 100 * qc_data$missing_probes / (nrow(betas) * ncol(betas))
  }
  
  # Clock correlations
  if (!is.null(results$combined)) {
    clock_cols <- setdiff(colnames(results$combined), "sample_id")
    if (length(clock_cols) > 1) {
      clock_matrix <- as.matrix(results$combined[, clock_cols])
      qc_data$clock_correlations <- cor(clock_matrix, use = "pairwise.complete.obs")
    }
  }
  
  # Save QC report
  saveRDS(qc_data, output_file)
  
  return(qc_data)
}
