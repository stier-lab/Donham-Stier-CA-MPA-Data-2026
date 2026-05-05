#' ---
#' title: "Data Processing Pipeline"
#' description: "Process raw monitoring data into harmonized CSVs for the analysis repo"
#' author: "Emily Donham & Adrian Stier"
#' project: "Donham-Stier-CA-MPA-Data-2026"
#' ---
#'
#' This script runs the data processing pipeline from raw monitoring data
#' through harmonized CSV export.
#'
#' Usage:
#'   source(here::here("code", "R", "run_all.R"))
#'
#' Pipeline:
#'   00  - Load packages
#'   00c - Analysis constants
#'   01  - Utility functions
#'   02  - pBACIPS function
#'   03  - Import raw data
#'   04  - PISCO processing
#'   05  - KFM processing
#'   06  - LTER processing
#'   06b - Landsat processing
#'   07  - Combine all sources & export harmonized CSVs
#'
#' Outputs:
#'   output/harmonized/harmonized_response_ratios.csv
#'   output/harmonized/harmonized_raw_responses.csv
#'   output/harmonized/harmonized_landsat_rr.csv
#'   output/harmonized/harmonized_site_metadata.csv
#'   output/harmonized/checksums.md5

####################################################################################################
## Setup ###########################################################################################
####################################################################################################

rm(list = ls())
pipeline_start <- Sys.time()

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install it with install.packages('here').")
}

####################################################################################################
## File-based logging ##############################################################################
####################################################################################################

log_dir <- here::here("logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_filename <- paste0("pipeline_", format(pipeline_start, "%Y%m%d_%H%M%S"), ".log")
log_filepath <- file.path(log_dir, log_filename)

log_message <- function(..., sep = "", append = TRUE) {
  msg <- paste(..., sep = sep)
  cat(msg)
  cat(msg, file = log_filepath, append = append)
}

log_message("", append = FALSE)
log_message("========================================================================\n")
log_message("  CA MPA Kelp Forest: Data Processing Pipeline\n")
log_message("  Started: ", format(pipeline_start, "%Y-%m-%d %H:%M:%S"), "\n")
log_message("========================================================================\n\n")

cat("\n")
cat("========================================================================\n")
cat("  CA MPA Kelp Forest: Data Processing Pipeline\n")
cat("  Started: ", format(pipeline_start, "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n\n")

####################################################################################################
## Helpers #########################################################################################
####################################################################################################

# Pipeline results tracker. Accumulates per-script timing and status
pipeline_results <- data.frame(
  label     = character(0),
  script    = character(0),
  status    = character(0),
  elapsed_s = numeric(0),
  error_msg = character(0),
  stringsAsFactors = FALSE
)

source_module <- function(filename, label) {
  start_msg <- paste0("---- [", label, "] Sourcing ", filename, " ...\n")
  cat(start_msg)
  log_message(start_msg)

  t0 <- Sys.time()
  tryCatch(
    {
      source(here::here("code", "R", filename), local = FALSE)
      elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
      end_msg <- paste0("---- [", label, "] ", filename, " completed (",
                        elapsed, " s)\n\n")
      cat(end_msg)
      log_message(end_msg)

      # Record success
      pipeline_results <<- rbind(pipeline_results, data.frame(
        label = label, script = filename, status = "OK",
        elapsed_s = elapsed, error_msg = "",
        stringsAsFactors = FALSE
      ))
    },
    error = function(e) {
      elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
      err_msg <- paste0("\n****** ERROR in [", label, "] ", filename, " ******\n",
                        "Message: ", conditionMessage(e), "\n",
                        "Elapsed before failure: ", elapsed, " s\n",
                        "Pipeline halted.\n")
      cat(err_msg)
      log_message(err_msg)

      # Record failure before stopping
      pipeline_results <<- rbind(pipeline_results, data.frame(
        label = label, script = filename, status = "FAILED",
        elapsed_s = elapsed, error_msg = conditionMessage(e),
        stringsAsFactors = FALSE
      ))

      stop("Failed at: ", filename, " -- ", conditionMessage(e), call. = FALSE)
    }
  )
}

check_objects <- function(obj_names, module_label) {
  missing <- obj_names[!vapply(obj_names, exists, logical(1), envir = globalenv())]
  if (length(missing) > 0) {
    msg <- paste0("Checkpoint failed after [", module_label, "]. Missing: ",
                  paste(missing, collapse = ", "))
    log_message(msg, "\n")
    stop(msg, call. = FALSE)
  }

  # Warn if any expected data frames are empty (possible partial execution)
  for (obj_name in obj_names) {
    obj <- get(obj_name, envir = globalenv())
    if (is.data.frame(obj) && nrow(obj) == 0) {
      warn_msg <- paste0("  WARNING: ", obj_name, " is an empty data frame after [",
                          module_label, "]\n")
      cat(warn_msg)
      log_message(warn_msg)
      warning("Checkpoint warning: ", obj_name, " is empty after module ", module_label,
              call. = FALSE)
    }
  }

  check_msg <- paste0("  Checkpoint OK: ", paste(obj_names, collapse = ", "), "\n\n")
  cat(check_msg)
  log_message(check_msg)
}

####################################################################################################
## Run pipeline ####################################################################################
####################################################################################################

# --- 00: Libraries ---
source_module("00_libraries.R", "00")

# --- 00c: Constants ---
source_module("00c_analysis_constants.R", "00c")

# --- 01: Utility functions ---
source_module("01_utils.R", "01")

# --- 02: pBACIPS function ---
source_module("02_pBACIPS_function.R", "02")

# --- Data validation: check raw data files exist ---
cat("---- Checking raw data availability...\n")
required_data <- c(
  here::here("data", "ALL_sizefreq_2024.csv"),
  here::here("data", "Site_List_All.csv"),
  here::here("data", "MPAfeatures_subset.csv"),
  here::here("data", "PISCO", "MLPA_kelpforest_swath_2024.csv"),
  here::here("data", "MBON", "SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv"),
  here::here("data", "LTER", "Annual_fish_comb_20240307.csv")
)
missing_data <- required_data[!file.exists(required_data)]
if (length(missing_data) > 0) {
  msg <- paste0(
    "\n========================================================================\n",
    "  MISSING RAW DATA FILES\n",
    "========================================================================\n\n",
    paste("  -", basename(missing_data), collapse = "\n"), "\n\n",
    "Set up data symlinks or download from Dryad. See README.md.\n"
  )
  stop(msg, call. = FALSE)
}
cat("  All required data files present.\n\n")

# --- 03: Data import ---
source_module("03_data_import.R", "03")
check_objects(c("Site"), "03")

# --- 04: PISCO processing ---
source_module("04_pisco_processing.R", "04")
check_objects(c("Swath.join.sub", "SizeFreq.Urch.OG"), "04")

# --- 05: KFM processing ---
source_module("05_kfm_processing.R", "05")
check_objects(c("KFM.join.ave"), "05")

# --- 06: LTER processing ---
source_module("06_lter_processing.R", "06")
check_objects(c("LTER.join.ave"), "06")

# --- 06b: Landsat ---
source_module("06b_landsat_processing.R", "06b")
# Landsat is optional. Warn but don't halt if missing
if (!exists("Landsat.RR", envir = globalenv())) {
  warning("Landsat.RR not found after 06b. Landsat data was unavailable. Continuing without it.")
} else {
  cat("  Checkpoint OK: Landsat.RR found.\n")
}

# --- 07: Combine & export ---
source_module("07_combine_and_export.R", "07")
check_objects(c("All.RR.sub.trans", "All.Resp.sub"), "07")

####################################################################################################
## Validate exports ################################################################################
####################################################################################################

harmonized_dir <- here::here("output", "harmonized")
expected_files <- c("harmonized_response_ratios.csv", "harmonized_raw_responses.csv",
                    "harmonized_landsat_rr.csv", "harmonized_site_metadata.csv",
                    "checksums.md5")
missing_exports <- expected_files[!file.exists(file.path(harmonized_dir, expected_files))]

if (length(missing_exports) > 0) {
  warning("Missing expected exports: ", paste(missing_exports, collapse = ", "))
} else {
  cat("  All harmonized CSVs exported successfully.\n\n")
}

####################################################################################################
## Pipeline summary ################################################################################
####################################################################################################

pipeline_end <- Sys.time()
elapsed_total <- difftime(pipeline_end, pipeline_start, units = "mins")

# --- Per-script timing table ---
summary_header <- paste0(
  "\n========================================================================\n",
  "  PIPELINE SUMMARY\n",
  "========================================================================\n\n"
)
cat(summary_header)
log_message(summary_header)

timing_header <- sprintf("  %-6s %-30s %-8s %10s\n", "Step", "Script", "Status", "Time (s)")
timing_sep    <- paste0("  ", strrep("-", 58), "\n")
cat(timing_header)
cat(timing_sep)
log_message(timing_header)
log_message(timing_sep)

for (i in seq_len(nrow(pipeline_results))) {
  row <- pipeline_results[i, ]
  line <- sprintf("  %-6s %-30s %-8s %10.1f\n",
                  row$label, row$script, row$status, row$elapsed_s)
  cat(line)
  log_message(line)
}
cat(timing_sep)
log_message(timing_sep)

total_line <- sprintf("  %-6s %-30s %-8s %10.1f\n",
                      "", "TOTAL", "",
                      sum(pipeline_results$elapsed_s))
cat(total_line)
log_message(total_line)

# --- Output file inventory ---
output_header <- "\n  Harmonized output files:\n"
cat(output_header)
log_message(output_header)

csv_files_summary <- list.files(harmonized_dir, pattern = "\\.csv$", full.names = TRUE)
for (f in csv_files_summary) {
  size_kb <- round(file.size(f) / 1024, 1)
  # Read row count from CSV header (nrow = lines - 1 for header)
  n_lines <- length(readLines(f, warn = FALSE)) - 1L
  file_line <- sprintf("    %-42s %7.1f KB  (%d rows)\n", basename(f), size_kb, n_lines)
  cat(file_line)
  log_message(file_line)
}

# --- Final status ---
n_failed <- sum(pipeline_results$status == "FAILED")
status_icon <- if (n_failed == 0) "COMPLETE" else paste0("FAILED (", n_failed, " errors)")

completion_msg <- paste0(
  "\n========================================================================\n",
  "  Status:   ", status_icon, "\n",
  "  Finished: ", format(pipeline_end, "%Y-%m-%d %H:%M:%S"), "\n",
  "  Total:    ", round(as.numeric(elapsed_total), 2), " minutes\n",
  "  Output:   ", harmonized_dir, "\n",
  "========================================================================\n\n"
)
cat(completion_msg)
log_message(completion_msg)

summary_footer <- paste0("Log file: ", log_filepath, "\nDone.\n")
cat(summary_footer)
log_message(summary_footer)
