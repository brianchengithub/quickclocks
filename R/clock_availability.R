#' ============================================================================
#' Clock Availability Module for Epigenetic Clock Calculator
#' ============================================================================
#' 
#' Functions to check which clocks can be computed based on available probes
#' and report availability status.

#' Define all supported clocks and their categories
#' @return List of clock definitions
#' @export
get_clock_definitions <- function() {
  list(
    # EpiDISH cell deconvolution
    epidish = list(
      cell_composition_rpc = list(
        name = "Cell Composition (EpiDISH RPC)",
        category = "QC",
        package = "EpiDISH",
        probes = NULL  # Uses built-in reference
      ),
      cell_composition_cp = list(
        name = "Cell Composition (EpiDISH CP)",
        category = "QC",
        package = "EpiDISH",
        probes = NULL
      )
    ),
    
    # Sex inference (custom implementation)
    sex = list(
      sex_inference = list(
        name = "Sex Inference",
        category = "QC", 
        package = NULL,  # Custom implementation using manifest
        probes = NULL
      )
    ),
    
    # DunedinPACE
    dunedin = list(
      DunedinPACE = list(
        name = "DunedinPACE",
        category = "Pace of Aging",
        package = "DunedinPACE",
        probes = NULL  # Will be loaded from package
      )
    ),
    
    # PC-Clocks (Levine Lab)
    pc_clocks = list(
      PCHorvath1 = list(name = "PC Horvath1", category = "Age", package = "PCClocks"),
      PCHorvath2 = list(name = "PC Horvath2", category = "Age", package = "PCClocks"),
      PCHannum = list(name = "PC Hannum", category = "Age", package = "PCClocks"),
      PCPhenoAge = list(name = "PC PhenoAge", category = "Age", package = "PCClocks"),
      PCGrimAge = list(name = "PC GrimAge", category = "Age", package = "PCClocks")
    ),
    
    # epiTOC2
    epitoc = list(
      epiTOC2 = list(
        name = "epiTOC2 (Stem Cell Divisions)",
        category = "Mitotic",
        package = "EpiMitClocks",
        probes = NULL
      )
    ),
    
    # methylCIPHER clocks - Age
    methylcipher_age = list(
      Horvath1 = list(name = "Horvath1 (Pan-tissue)", category = "Age", package = "methylCIPHER"),
      Horvath2 = list(name = "Horvath2 (Skin+Blood)", category = "Age", package = "methylCIPHER"),
      Hannum = list(name = "Hannum", category = "Age", package = "methylCIPHER"),
      PhenoAge = list(name = "PhenoAge (Levine)", category = "Age", package = "methylCIPHER"),
      GrimAge1 = list(name = "GrimAge v1", category = "Age", package = "methylCIPHER"),
      GrimAge2 = list(name = "GrimAge v2", category = "Age", package = "methylCIPHER"),
      Zhang = list(name = "Zhang", category = "Age", package = "methylCIPHER"),
      Zhang2019 = list(name = "Zhang (2019)", category = "Age", package = "methylCIPHER"),
      Lin = list(name = "Lin", category = "Age", package = "methylCIPHER"),
      DNAmTL = list(name = "DNAmTL (Telomere Length)", category = "Age", package = "methylCIPHER")
    ),
    
    # methylCIPHER clocks - Specialized
    methylcipher_specialized = list(
      AdaptAge = list(name = "AdaptAge", category = "Specialized", package = "methylCIPHER"),
      CausAge = list(name = "CausAge", category = "Specialized", package = "methylCIPHER"),
      DamAge = list(name = "DamAge", category = "Specialized", package = "methylCIPHER"),
      HypoClock = list(name = "HypoClock", category = "Specialized", package = "methylCIPHER"),
      MiAge = list(name = "MiAge", category = "Specialized", package = "methylCIPHER"),
      SystemsAge = list(name = "SystemsAge", category = "Specialized", package = "methylCIPHER"),
      RetroelementAge450K = list(name = "Retroelement-Age 450K", category = "Specialized", package = "methylCIPHER")
    ),
    
    # methylCIPHER clocks - Lifestyle
    methylcipher_lifestyle = list(
      Alcohol = list(name = "Alcohol", category = "Lifestyle", package = "methylCIPHER"),
      BMI = list(name = "BMI", category = "Lifestyle", package = "methylCIPHER"),
      Smoking = list(name = "Smoking", category = "Lifestyle", package = "methylCIPHER")
    ),
    
    # methylCIPHER clocks - Organ/System
    methylcipher_organ = list(
      Blood = list(name = "Blood", category = "Organ", package = "methylCIPHER"),
      Brain = list(name = "Brain", category = "Organ", package = "methylCIPHER"),
      Heart = list(name = "Heart", category = "Organ", package = "methylCIPHER"),
      Kidney = list(name = "Kidney", category = "Organ", package = "methylCIPHER"),
      Liver = list(name = "Liver", category = "Organ", package = "methylCIPHER"),
      Lung = list(name = "Lung", category = "Organ", package = "methylCIPHER"),
      MusculoSkeletal = list(name = "MusculoSkeletal", category = "Organ", package = "methylCIPHER")
    ),
    
    # methylCIPHER clocks - Biological
    methylcipher_biological = list(
      Inflammation = list(name = "Inflammation", category = "Biological", package = "methylCIPHER"),
      Hormone = list(name = "Hormone", category = "Biological", package = "methylCIPHER"),
      Immune = list(name = "Immune", category = "Biological", package = "methylCIPHER"),
      Metabolic = list(name = "Metabolic", category = "Biological", package = "methylCIPHER")
    )
  )
}


#' Get flat list of all clock names
#' @return Character vector of clock names
#' @export
get_all_clock_names <- function() {
  defs <- get_clock_definitions()
  
  all_names <- c()
  for (category in names(defs)) {
    all_names <- c(all_names, names(defs[[category]]))
  }
  
  return(all_names)
}


#' Check availability of all clocks based on available probes
#' 
#' @param available_probes Character vector of probe IDs in the data
#' @param reference_betas Reference beta values (for zero-shot capability)
#' @param config Configuration list
#' @return Data frame with clock availability information
#' @export
check_all_clock_availability <- function(available_probes, reference_betas = NULL,
                                          config = list()) {
  
  verbose <- config$verbose %||% TRUE
  
  # Get clock definitions
  clock_defs <- get_clock_definitions()
  
  # Initialize results
  results <- data.frame(
    clock_name = character(),
    display_name = character(),
    category = character(),
    package = character(),
    probes_required = integer(),
    probes_available = integer(),
    probes_missing = integer(),
    probes_can_impute = integer(),
    pct_available = numeric(),
    can_compute = logical(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  
  # Check each clock family
  for (family in names(clock_defs)) {
    for (clock_id in names(clock_defs[[family]])) {
      clock_info <- clock_defs[[family]][[clock_id]]
      
      # Get required probes for this clock
      required_probes <- get_clock_probes(clock_id, clock_info$package)
      
      if (is.null(required_probes) || length(required_probes) == 0) {
        # Clock doesn't specify probes (uses built-in method)
        # Check if package is available
        pkg_available <- is_package_installed(clock_info$package)
        
        results <- rbind(results, data.frame(
          clock_name = clock_id,
          display_name = clock_info$name,
          category = clock_info$category,
          package = clock_info$package,
          probes_required = NA_integer_,
          probes_available = NA_integer_,
          probes_missing = NA_integer_,
          probes_can_impute = NA_integer_,
          pct_available = NA_real_,
          can_compute = pkg_available,
          reason = if (pkg_available) "Package available" else "Package not installed",
          stringsAsFactors = FALSE
        ))
        
      } else {
        # Check probe availability
        n_required <- length(required_probes)
        n_available <- sum(required_probes %in% available_probes)
        n_missing <- n_required - n_available
        
        # Check how many missing probes can be imputed from reference
        missing_probes <- setdiff(required_probes, available_probes)
        n_can_impute <- if (!is.null(reference_betas)) {
          sum(missing_probes %in% names(reference_betas))
        } else {
          0L
        }
        
        pct_available <- 100 * (n_available + n_can_impute) / n_required
        
        # Determine if clock can be computed
        # Allow computation if we have at least 95% of probes available or imputable
        min_pct <- config$min_probe_pct %||% 95
        can_compute <- pct_available >= min_pct
        
        reason <- if (can_compute) {
          "Sufficient probes"
        } else {
          sprintf("Only %.1f%% probes available", pct_available)
        }
        
        results <- rbind(results, data.frame(
          clock_name = clock_id,
          display_name = clock_info$name,
          category = clock_info$category,
          package = clock_info$package,
          probes_required = n_required,
          probes_available = n_available,
          probes_missing = n_missing,
          probes_can_impute = n_can_impute,
          pct_available = pct_available,
          can_compute = can_compute,
          reason = reason,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(results)
}


#' Get required probes for a specific clock
#' 
#' @param clock_name Name of the clock
#' @param package_name Package containing the clock
#' @return Character vector of probe IDs
#' @export
get_clock_probes <- function(clock_name, package_name) {
  
  # Try to get probes from package
  probes <- tryCatch({
    switch(package_name,
      "DunedinPACE" = get_dunedin_probes(),
      "PCClocks" = get_pc_clock_probes(clock_name),
      "EpiMitClocks" = get_epitoc_probes(),
      "methylCIPHER" = get_methylcipher_probes(clock_name),
      "sesame" = NULL,  # SeSAMe functions don't require specific probes
      NULL
    )
  }, error = function(e) {
    NULL
  })
  
  return(probes)
}


#' Get DunedinPACE probes
#' @return Character vector of probe IDs
#' @export
get_dunedin_probes <- function() {
  if (!is_package_installed("DunedinPACE")) {
    return(NULL)
  }
  
  tryCatch({
    # DunedinPACE stores probe info in package data
    data("DunedinPACE_probes", package = "DunedinPACE", envir = environment())
    if (exists("DunedinPACE_probes", envir = environment())) {
      return(get("DunedinPACE_probes", envir = environment()))
    }
    NULL
  }, error = function(e) {
    NULL
  })
}


#' Get PC-Clock probes
#' @param clock_name Specific PC clock name
#' @return Character vector of probe IDs
#' @export
get_pc_clock_probes <- function(clock_name) {
  if (!is_package_installed("PCClocks")) {
    return(NULL)
  }
  
  tryCatch({
    # PCClocks stores coefficients that contain probe names
    NULL  # Actual implementation depends on package structure
  }, error = function(e) {
    NULL
  })
}


#' Get epiTOC2 probes
#' @return Character vector of probe IDs
#' @export
get_epitoc_probes <- function() {
  if (!is_package_installed("EpiMitClocks")) {
    return(NULL)
  }
  
  tryCatch({
    NULL  # Actual implementation depends on package structure
  }, error = function(e) {
    NULL
  })
}


#' Get methylCIPHER probes for a specific clock
#' @param clock_name Clock name
#' @return Character vector of probe IDs
#' @export
get_methylcipher_probes <- function(clock_name) {
  if (!is_package_installed("methylCIPHER")) {
    return(NULL)
  }
  
  tryCatch({
    NULL  # Actual implementation depends on package structure
  }, error = function(e) {
    NULL
  })
}


#' Print clock availability report
#' 
#' @param availability Data frame from check_all_clock_availability
#' @param verbose Print to console
#' @export
print_clock_availability <- function(availability, verbose = TRUE) {
  if (!verbose) return(invisible(NULL))
  
  # Summary statistics
  n_total <- nrow(availability)
  n_computable <- sum(availability$can_compute)
  n_not_computable <- n_total - n_computable
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("                  CLOCK AVAILABILITY REPORT\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("\n")
  
  cat(sprintf("Total clocks defined: %d\n", n_total))
  cat(sprintf("Clocks CAN be computed: %d\n", n_computable))
  cat(sprintf("Clocks CANNOT be computed: %d\n", n_not_computable))
  cat("\n")
  
  # Group by category
  categories <- unique(availability$category)
  
  for (cat in categories) {
    cat_clocks <- availability[availability$category == cat, ]
    n_cat_computable <- sum(cat_clocks$can_compute)
    
    cat(sprintf("─── %s (%d/%d available) ───\n", 
                cat, n_cat_computable, nrow(cat_clocks)))
    
    for (i in seq_len(nrow(cat_clocks))) {
      clock <- cat_clocks[i, ]
      status <- if (clock$can_compute) "✓" else "✗"
      
      if (is.na(clock$probes_required)) {
        probe_info <- ""
      } else {
        probe_info <- sprintf(" [%d/%d probes]", 
                              clock$probes_available + clock$probes_can_impute,
                              clock$probes_required)
      }
      
      cat(sprintf("  %s %s%s\n", status, clock$display_name, probe_info))
    }
    cat("\n")
  }
  
  # List reasons for non-computable clocks
  not_computable <- availability[!availability$can_compute, ]
  if (nrow(not_computable) > 0) {
    cat("─── Reasons for unavailability ───\n")
    for (i in seq_len(nrow(not_computable))) {
      cat(sprintf("  • %s: %s\n", 
                  not_computable$display_name[i], 
                  not_computable$reason[i]))
    }
    cat("\n")
  }
}


#' Write clock availability report to file
#' 
#' @param availability Data frame from check_all_clock_availability
#' @param file_path Output file path
#' @export
write_clock_availability_report <- function(availability, file_path) {
  
  # Open file connection
  con <- file(file_path, "w")
  on.exit(close(con))
  
  # Write header
  writeLines("EPIGENETIC CLOCK AVAILABILITY REPORT", con)
  writeLines(paste("Generated:", Sys.time()), con)
  writeLines(paste(rep("=", 60), collapse = ""), con)
  writeLines("", con)
  
  # Summary
  n_total <- nrow(availability)
  n_computable <- sum(availability$can_compute)
  
  writeLines(sprintf("Total clocks: %d", n_total), con)
  writeLines(sprintf("Computable: %d", n_computable), con)
  writeLines(sprintf("Not computable: %d", n_total - n_computable), con)
  writeLines("", con)
  
  # Computable clocks
  writeLines("CLOCKS THAT CAN BE COMPUTED:", con)
  writeLines(paste(rep("-", 40), collapse = ""), con)
  
  computable <- availability[availability$can_compute, ]
  for (i in seq_len(nrow(computable))) {
    writeLines(sprintf("  • %s (%s)", 
                       computable$display_name[i],
                       computable$category[i]), con)
  }
  writeLines("", con)
  
  # Non-computable clocks
  writeLines("CLOCKS THAT CANNOT BE COMPUTED:", con)
  writeLines(paste(rep("-", 40), collapse = ""), con)
  
  not_computable <- availability[!availability$can_compute, ]
  for (i in seq_len(nrow(not_computable))) {
    writeLines(sprintf("  • %s: %s", 
                       not_computable$display_name[i],
                       not_computable$reason[i]), con)
  }
  
  # Write detailed table
  writeLines("", con)
  writeLines("DETAILED AVAILABILITY:", con)
  writeLines(paste(rep("-", 40), collapse = ""), con)
  
  # Write as formatted table
  write.csv(availability, file = gsub("\\.txt$", "_details.csv", file_path), 
            row.names = FALSE)
}


#' Null coalescing operator
#' @keywords internal
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}
