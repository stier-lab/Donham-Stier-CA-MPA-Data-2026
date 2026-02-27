# =============================================================================
# 07_combine_and_export.R
# =============================================================================
#
# PURPOSE:
#   Combine response ratio and raw response data from all monitoring programs
#   (LTER, PISCO, KFM) into unified datasets for the meta-analysis.
#
# WHAT THIS SCRIPT DOES:
#   1. Combines response ratio data from all programs into All.RR
#   2. Applies MPA exclusion rules (removes sites with insufficient data)
#   3. Handles special cases for sheephead-only MPAs
#   4. Standardizes species names to full scientific names
#   5. Standardizes source names (e.g., "LTER macro surveys" -> "LTER")
#   6. Combines raw response data (density/biomass) into All.Resp
#   7. Assigns Before/After labels and time since MPA to all data
#   8. Exports summary statistics
#
# EXCLUSION LOGIC:
#   Some MPAs only have data for certain taxa. We:
#   - Exclude MPAs from SHEEPHEAD_ONLY_MPAS list from main analysis
#   - Re-include those MPAs specifically for sheephead analysis
#   - This prevents sheephead from being over-represented in cross-taxa summaries
#
# INPUTS:
#   Response ratio objects:
#   - LTER.join.ave, Swath.join.sub (PISCO), KFM.join.ave
#   - LTER.lob, Short.lter.macro, LTER.fish, kfm.fish
#
#   Raw response objects:
#   - LTER.resp, KFM.resp, PISCO.resp
#   - KFM.fish.den.long, LTER.fish.resp, LTER.macro.resp, LTER.lob.resp
#
#   Site metadata: Site, Sites2
#
# OUTPUTS:
#   - All.RR.sub.trans: Combined response ratios with standardized names
#   - All.Resp.sub: Combined raw data with time and BA columns
#   - average_responses.csv: Summary statistics (written to project root)
#
# DEPENDENCIES:
#   Requires 00-06b scripts to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

# =============================================================================
# SECTION 1: SETUP
# =============================================================================

cat("\n=== 07_combine_and_export.R ===\n")

# =============================================================================
# SECTION 2: COMBINE RESPONSE RATIO DATA
# =============================================================================
# At this point we have response ratio dataframes from each monitoring program:
# - LTER.join.ave (LTER urchins), Swath.join.sub (PISCO), KFM.join.ave (KFM urchins/kelp)
# - LTER.lob (LTER lobsters), Short.lter.macro (LTER kelp), LTER.fish (LTER sheephead)
# - kfm.fish (KFM sheephead)
# Now we combine them into one master dataset for the meta-analysis.
#
# METHODOLOGICAL NOTE — PROPORTION-BASED lnRR:
# All lnRR values in these dataframes were computed on zero-corrected
# PROPORTIONS (PropCorr), not on raw density or biomass values.
# Each processing script (04-06) called calculate_proportions() to normalize
# raw values to the proportion of the time-series maximum within each
# MPA x taxon x status group, then spread() to wide format and
# calculate_log_response_ratio() to compute lnRR = ln(MPA / Reference).
#
# This proportion-based approach is non-standard compared to typical lnRR
# meta-analysis on raw means (e.g., Hedges et al. 1999). It was chosen to
# standardize across monitoring programs (PISCO, KFM, LTER) that use
# different transect areas, survey methods, and spatial extents. See the
# detailed rationale in calculate_proportions() in 01_utils.R and the
# manuscript methods section.
#
# NOTE on cross-program overlap: Some MPAs (e.g., Campus Point SMCA, Naples
# SMCA) have data from both PISCO and LTER. These are intentionally retained
# as separate observations because they represent independent estimates from
# different monitoring methods. The `source` column distinguishes them, and the
# meta-analysis accounts for this via a Source random effect.

# Ensure consistent column structure before binding
# Different sources may have extra columns that need to be removed:
#   - Swath.join.sub has site_designation
#   - kfm.fish has sample_method
#   - KFM.join.ave has area
# Use explicit column selection for robustness.

# Standard columns for all response ratio dataframes
RR_STANDARD_COLS <- c("CA_MPA_Name_Short", "year", "y", "lnDiff", "mpa", "reference",
                      "Diff", "resp", "time", "type", "Location", "Hectares", "source", "BA")

# Validate required columns exist, then subset to standard set
Swath.join.sub   <- validate_and_subset_columns(Swath.join.sub, RR_STANDARD_COLS, "Swath.join.sub")
kfm.fish         <- validate_and_subset_columns(kfm.fish, RR_STANDARD_COLS, "kfm.fish")
KFM.join.ave     <- validate_and_subset_columns(KFM.join.ave, RR_STANDARD_COLS, "KFM.join.ave")
LTER.join.ave    <- validate_and_subset_columns(LTER.join.ave, RR_STANDARD_COLS, "LTER.join.ave")
LTER.lob         <- validate_and_subset_columns(LTER.lob, RR_STANDARD_COLS, "LTER.lob")
Short.lter.macro <- validate_and_subset_columns(Short.lter.macro, RR_STANDARD_COLS, "Short.lter.macro")
LTER.fish        <- validate_and_subset_columns(LTER.fish, RR_STANDARD_COLS, "LTER.fish")

# Combine all response ratio datasets using rbind() (row bind)
# This stacks all dataframes vertically into one master dataframe
All.RR <- rbind(LTER.join.ave, Swath.join.sub, KFM.join.ave, LTER.lob,
                Short.lter.macro, LTER.fish, kfm.fish)

# Report per-source row counts
cat("  Response ratio rows by source:\n")
rr_source_counts <- table(All.RR$source)
for (src in sort(names(rr_source_counts))) {
  cat(sprintf("    %-10s %d rows\n", src, rr_source_counts[[src]]))
}
cat(sprintf("    %-10s %d rows\n", "TOTAL", nrow(All.RR)))

# =============================================================================
# SECTION 3: FILTER EXCLUDED MPAs / SHEEPHEAD RE-INCLUSION
# =============================================================================
# IMPORTANT: Some MPAs only have data for certain taxa (usually just sheephead from PISCO fish surveys).
# If we include these MPAs in cross-taxa summaries, sheephead would be over-represented.
# SOLUTION: Remove these "sheephead-only" MPAs from the main analysis, then selectively
# add them back ONLY for sheephead-specific analyses.
#
# SHEEPHEAD_ONLY_MPAS is a constant defined in 01_utils.R containing these MPA names.

# Remove sheephead-only MPAs from the main dataset
# The ! negates the %in% check, so we KEEP rows that are NOT in the exclusion list
n_before_exclusion <- nrow(All.RR)
All.RR.sub <- subset(All.RR, !(CA_MPA_Name_Short %in% SHEEPHEAD_ONLY_MPAS))
n_excluded <- n_before_exclusion - nrow(All.RR.sub)
cat(sprintf("  Sheephead-only MPA exclusion: removed %d rows (%d -> %d)\n",
            n_excluded, n_before_exclusion, nrow(All.RR.sub)))

# Re-include sheephead-only MPAs for SPUL (sheephead) taxa only
# Uses reinclude_sheephead_mpas() from 01_utils.R to avoid duplicating logic
All.spul <- reinclude_sheephead_mpas(All.RR, taxon_col = "y", sheephead_value = "SPUL")

# Assertion: all standard reinclusion rows must be SPUL (Santa Barbara Island
# SMR can have non-SPUL KFM rows, which is the expected special case)
spul_non_sbi <- All.spul[All.spul$CA_MPA_Name_Short != "Santa Barbara Island SMR", ]
if (nrow(spul_non_sbi) > 0 && !all(spul_non_sbi$y == "SPUL")) {
  non_spul <- unique(spul_non_sbi$y[spul_non_sbi$y != "SPUL"])
  warning("Sheephead reinclusion contains non-SPUL taxa outside SBI SMR: ",
          paste(non_spul, collapse = ", "), call. = FALSE)
}
cat(sprintf("  Sheephead re-inclusion: added back %d rows from sheephead-only MPAs (%d SPUL, %d SBI-KFM other)\n",
            nrow(All.spul),
            sum(All.spul$y == "SPUL"),
            sum(All.spul$y != "SPUL")))

# Add the sheephead-specific rows back to the main dataset
All.RR.sub <- rbind(All.RR.sub, All.spul)

# =============================================================================
# SECTION 4: DEDUPLICATION (RESPONSE RATIOS)
# =============================================================================
# Check for duplicates that might arise from overlapping data sources or processing errors.
# A duplicate is defined as same MPA + year + taxa + response type + source combination.

dup_key_cols <- c("CA_MPA_Name_Short", "year", "y", "resp", "source")
dup_check <- duplicated(All.RR.sub[, dup_key_cols])
n_dups <- sum(dup_check)

if (n_dups > 0) {
  warning("Found ", n_dups, " duplicate rows in All.RR.sub (by MPA/year/taxa/resp/source). ",
          "Removing duplicates, keeping first occurrence.", call. = FALSE, immediate. = TRUE)
  cat("\nDuplicate rows detected in All.RR.sub:", n_dups, "\n")
  cat("First 5 duplicates (before removal):\n")
  print(head(All.RR.sub[dup_check, dup_key_cols], 5))

  # FIXED (2026-02-06): Remove duplicates to prevent bias in meta-analysis
  # Keep first occurrence of each unique MPA/year/taxa/resp/source combination
  n_before <- nrow(All.RR.sub)
  All.RR.sub <- All.RR.sub[!duplicated(All.RR.sub[, dup_key_cols]), ]
  n_after <- nrow(All.RR.sub)
  cat("Removed", n_before - n_after, "duplicate rows. Rows remaining:", n_after, "\n")
} else {
  cat("Duplicate check passed: No duplicate MPA/year/taxa/resp/source combinations in All.RR.sub.\n")
}

# =============================================================================
# SECTION 5: SPECIES AND SOURCE NAME STANDARDIZATION
# =============================================================================
# Different monitoring programs use different species codes:
# - PISCO: STRPURAD, MESFRAAD, MACPYRAD, PANINT, SPUL
# - KFM/LTER: Full names like "Strongylocentrotus purpuratus"
# For consistent analysis and plotting, we standardize to full scientific names.

All.RR.sub.trans <- All.RR.sub
All.RR.sub.trans$y <- as.character(All.RR.sub.trans$y)

# standardize_species_names() from 01_utils.R converts codes to full names:
# STRPURAD -> Strongylocentrotus purpuratus (purple urchin)
# MESFRAAD -> Mesocentrotus franciscanus (red urchin)
# MACPYRAD -> Macrocystis pyrifera (giant kelp)
# PANINT -> Panulirus interruptus (California spiny lobster)
# SPUL -> Semicossyphus pulcher (California sheephead)
All.RR.sub.trans$y <- standardize_species_names(All.RR.sub.trans$y)

# Standardize source names to simple program names
# Some datasets have more descriptive source names that we simplify
All.RR.sub.trans$source[All.RR.sub.trans$source == "LTER macro surveys"] <- "LTER"
All.RR.sub.trans$source[All.RR.sub.trans$source == "LTER lob surveys"] <- "LTER"

# =============================================================================
# SECTION 6: COMBINE RAW RESPONSE DATA
# =============================================================================

# Standardize species names and column names before combining
# Use explicit column selection instead of fragile indices

# Helper function to prepare response dataframe for binding
prepare_resp_df <- function(df, status_col = "status") {
  # Rename y/taxon_name column to taxon_name if needed
  if ("y" %in% names(df) && !"taxon_name" %in% names(df)) {
    names(df)[names(df) == "y"] <- "taxon_name"
  }
  # Rename site_status to status if needed
  if ("site_status" %in% names(df) && !"status" %in% names(df)) {
    names(df)[names(df) == "site_status"] <- "status"
  }
  # Validate and select standard columns
  standard_cols <- c("CA_MPA_Name_Short", "year", "taxon_name", "source", "status", "value", "resp")
  validate_and_subset_columns(df, standard_cols, context = deparse(substitute(df)))
}

# Standardize species names
PISCO.resp$y <- standardize_species_names(as.character(PISCO.resp$y))
KFM.resp$y <- standardize_species_names(as.character(KFM.resp$y))
LTER.resp$y <- standardize_species_names(as.character(LTER.resp$y))
LTER.fish.resp$y <- standardize_species_names(as.character(LTER.fish.resp$y))
LTER.macro.resp$y <- standardize_species_names(as.character(LTER.macro.resp$y))
LTER.lob.resp$y <- standardize_species_names(as.character(LTER.lob.resp$y))

# Convert SPUL codes to full names
KFM.fish.den.long$taxon_name[KFM.fish.den.long$taxon_name == "SPUL"] <- "Semicossyphus pulcher"
LTER.fish.resp$y[LTER.fish.resp$y == "SPUL"] <- "Semicossyphus pulcher"

# Prepare each dataframe with standardized columns
PISCO.resp.std <- prepare_resp_df(PISCO.resp)
KFM.resp.std <- prepare_resp_df(KFM.resp)
KFM.fish.den.long.std <- prepare_resp_df(KFM.fish.den.long)
LTER.resp.std <- prepare_resp_df(LTER.resp, status_col = "site_status")
LTER.fish.resp.std <- prepare_resp_df(LTER.fish.resp)
LTER.macro.resp.std <- prepare_resp_df(LTER.macro.resp)
LTER.lob.resp.std <- prepare_resp_df(LTER.lob.resp)

# Combine all raw response datasets
All.Resp <- rbind(KFM.fish.den.long.std, LTER.fish.resp.std, LTER.macro.resp.std,
                  LTER.lob.resp.std, LTER.resp.std, KFM.resp.std, PISCO.resp.std)

# Report per-source row counts for raw responses
cat("  Raw response rows by source:\n")
resp_source_counts <- table(All.Resp$source)
for (src in sort(names(resp_source_counts))) {
  cat(sprintf("    %-10s %d rows\n", src, resp_source_counts[[src]]))
}
cat(sprintf("    %-10s %d rows\n", "TOTAL", nrow(All.Resp)))

# =============================================================================
# SECTION 7: FILTER RAW RESPONSE DATA (SHEEPHEAD RE-INCLUSION)
# =============================================================================

All.Resp.sub <- subset(All.Resp, !(CA_MPA_Name_Short %in% SHEEPHEAD_ONLY_MPAS))

# Re-include sheephead-only MPAs for Semicossyphus pulcher
# Uses reinclude_sheephead_mpas() from 01_utils.R to avoid duplicating logic
All.Resp.spul <- reinclude_sheephead_mpas(All.Resp, taxon_col = "taxon_name",
                                           sheephead_value = "Semicossyphus pulcher")

All.Resp.sub <- rbind(All.Resp.sub, All.Resp.spul)

# =============================================================================
# SECTION 8: DEDUPLICATION (RAW RESPONSES)
# =============================================================================
# Check for duplicates in raw response data.

dup_key_resp <- c("CA_MPA_Name_Short", "year", "taxon_name", "resp", "source", "status")
dup_check_resp <- duplicated(All.Resp.sub[, dup_key_resp])
n_dups_resp <- sum(dup_check_resp)

if (n_dups_resp > 0) {
  warning("Found ", n_dups_resp, " duplicate rows in All.Resp.sub. ",
          "Review data processing pipeline for overlapping records.")
  cat("\nDuplicate rows detected in All.Resp.sub:", n_dups_resp, "\n")
} else {
  cat("Duplicate check passed: No duplicate rows in All.Resp.sub.\n")
}

# =============================================================================
# SECTION 9: ASSIGN BEFORE/AFTER AND TIME TO RAW RESPONSE DATA
# =============================================================================

# Use vectorized utility functions instead of year-by-year if-else blocks
All.Resp.sub <- assign_ba_from_site_table(All.Resp.sub, Site)
All.Resp.sub <- assign_time_from_site_table(All.Resp.sub, Site)

# Keep only essential columns (explicit selection for robustness)
RESP_ESSENTIAL_COLS <- c("CA_MPA_Name_Short", "year", "taxon_name", "source", "status", "value", "resp", "BA", "time")
All.Resp.sub <- All.Resp.sub[, RESP_ESSENTIAL_COLS]

# =============================================================================
# SECTION 10: CALCULATE AVERAGE RESPONSES
# =============================================================================

AveResponse <- summarySE(data = All.Resp.sub, measurevar = "value",
                         groupvars = c("taxon_name", "source", "resp"))
write_csv(AveResponse, here::here("tables", "average_responses.csv"))

# =============================================================================
# SECTION 11: DATA PROVENANCE SUMMARY
# =============================================================================
# Export a summary table documenting what data is available per source, taxa, and response type.
# This is essential for the methods section and supplemental materials.

# --- Response ratio data provenance ---
rr_provenance <- All.RR.sub.trans %>%
  dplyr::group_by(source, y, resp) %>%
  dplyr::summarise(
    n_observations = dplyr::n(),
    n_MPAs = dplyr::n_distinct(CA_MPA_Name_Short),
    n_years = dplyr::n_distinct(year),
    year_range = paste(min(year), max(year), sep = "-"),
    MPAs = paste(sort(unique(CA_MPA_Name_Short)), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(source, y, resp)

write_csv(rr_provenance, here::here("tables", "table_data_provenance_rr.csv"))
cat("Data provenance (response ratios) exported to: tables/table_data_provenance_rr.csv\n")
cat("  Sources:", paste(sort(unique(rr_provenance$source)), collapse = ", "), "\n")
cat("  Total RR observations:", sum(rr_provenance$n_observations), "\n")

# --- Raw response data provenance ---
raw_provenance <- All.Resp.sub %>%
  dplyr::group_by(source, taxon_name, resp) %>%
  dplyr::summarise(
    n_observations = dplyr::n(),
    n_MPAs = dplyr::n_distinct(CA_MPA_Name_Short),
    n_years = dplyr::n_distinct(year),
    year_range = paste(min(year), max(year), sep = "-"),
    n_inside = sum(tolower(status) %in% c("inside", "mpa", "impact", "i"), na.rm = TRUE),
    n_outside = sum(tolower(status) %in% c("outside", "reference", "control", "ref", "o", "r"), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(source, taxon_name, resp)

write_csv(raw_provenance, here::here("tables", "table_data_provenance_raw.csv"))
cat("Data provenance (raw responses) exported to: tables/table_data_provenance_raw.csv\n")

# --- Filtering summary (sheephead-only and duplicate exclusions) ---
# Uses actual counts from each pipeline step (n_before_exclusion, n_excluded,
# nrow(All.spul), n_dups) rather than back-calculating from the final dataset.
filtering_summary <- data.frame(
  Step = c(
    "01_Raw_combined",
    "02_Sheephead_only_MPAs_removed",
    "03_Sheephead_MPAs_readded_for_SPUL",
    "04_Duplicates_removed",
    "05_Final_RR_dataset"
  ),
  Description = c(
    "All response ratio data combined from LTER, PISCO, KFM",
    paste("Removed sheephead-only MPAs:", paste(SHEEPHEAD_ONLY_MPAS, collapse = ", ")),
    paste("Re-included", nrow(All.spul), "rows from SHEEPHEAD_REINCLUSION_MPAS + SBI SMR"),
    paste("Removed", n_dups, "duplicate MPA/year/taxa/resp/source rows"),
    "Final dataset for effect size calculation"
  ),
  N_rows = c(
    n_before_exclusion,
    n_before_exclusion - n_excluded,
    n_before_exclusion - n_excluded + nrow(All.spul),
    n_dups,
    nrow(All.RR.sub.trans)
  ),
  stringsAsFactors = FALSE
)
write_csv(filtering_summary, here::here("outputs", "data_filtering_steps.csv"))
cat("Data filtering steps exported to: outputs/data_filtering_steps.csv\n")

# Clean up intermediate objects
rm(All.RR, All.spul, All.Resp, All.Resp.spul)


# =============================================================================
# SECTION 12: BOOTSTRAP DATA AVAILABILITY REPORT
# =============================================================================
# Document which programs provide size-frequency data for bootstrap biomass
# estimation, and which rely on cross-program PISCO size distributions.
# This makes the cross-program bootstrap assumption explicit and auditable.

bootstrap_report <- data.frame(
  Source = c("PISCO", "PISCO", "KFM", "KFM", "LTER", "LTER"),
  Taxon = c("Urchins (S. purpuratus, M. franciscanus)",
            "Lobster (P. interruptus, VRG only)",
            "Urchins (S. purpuratus, M. franciscanus)",
            "Kelp (M. pyrifera)",
            "Urchins (S. purpuratus, M. franciscanus)",
            "Lobster (P. interruptus)"),
  Size_Data_Source = c("Own (PISCO size-frequency with 25mm filter)",
                       "Own (PISCO VRG size-frequency)",
                       "Cross-program (PISCO SizeFreq.Urch.OG, no 25mm filter)",
                       "Own (KFM stipe counts per individual)",
                       "Cross-program (PISCO SizeFreq.Urch.OG, no 25mm filter)",
                       "Own (LTER measures carapace length directly)"),
  Assumption = c("None — native size data",
                 "None — native size data",
                 "Assumes PISCO urchin size structure applies to KFM sites",
                 "None — native stipe count data",
                 "Assumes PISCO urchin size structure applies to LTER sites",
                 "None — native size data"),
  stringsAsFactors = FALSE
)

write.csv(bootstrap_report, here::here("outputs", "bootstrap_data_availability.csv"), row.names = FALSE)
cat("Bootstrap data availability report exported to: outputs/bootstrap_data_availability.csv\n")

# Summary of cross-program dependencies
n_cross_program <- sum(grepl("Cross-program", bootstrap_report$Size_Data_Source))
if (n_cross_program > 0) {
  cat(sprintf("  NOTE: %d taxa/source combinations use cross-program PISCO size data for bootstrap.\n",
              n_cross_program))
  cat("  See outputs/bootstrap_data_availability.csv for details.\n")
}

# =============================================================================
# SECTION 13: EXPORT HARMONIZED CSVs
# =============================================================================
# These 4 CSVs are the handoff between the data-processing and analysis repos.
# They are also the analysis-ready dataset published on Dryad.

harmonized_dir <- here::here("output", "harmonized")
dir.create(harmonized_dir, recursive = TRUE, showWarnings = FALSE)

# --- Validate response ratio data before export ---
rr_required_cols <- c("CA_MPA_Name_Short", "year", "y", "lnDiff", "mpa", "reference",
                      "Diff", "resp", "source", "BA")
stopifnot("Expected columns missing from All.RR.sub.trans" =
            all(rr_required_cols %in% names(All.RR.sub.trans)))
cat(sprintf("  harmonized_response_ratios.csv: %d rows x %d cols\n",
            nrow(All.RR.sub.trans), ncol(All.RR.sub.trans)))

write.csv(All.RR.sub.trans, file.path(harmonized_dir, "harmonized_response_ratios.csv"), row.names = FALSE)

# --- Validate raw response data before export ---
resp_required_cols <- c("CA_MPA_Name_Short", "year", "taxon_name", "source",
                        "status", "value", "resp", "BA", "time")
stopifnot("Expected columns missing from All.Resp.sub" =
            all(resp_required_cols %in% names(All.Resp.sub)))
cat(sprintf("  harmonized_raw_responses.csv: %d rows x %d cols\n",
            nrow(All.Resp.sub), ncol(All.Resp.sub)))

write.csv(All.Resp.sub, file.path(harmonized_dir, "harmonized_raw_responses.csv"), row.names = FALSE)

# --- Validate and export Landsat data ---
if (exists("Landsat.RR") && is.data.frame(Landsat.RR) && nrow(Landsat.RR) > 0) {
  cat(sprintf("  harmonized_landsat_rr.csv: %d rows x %d cols\n",
              nrow(Landsat.RR), ncol(Landsat.RR)))
  write.csv(Landsat.RR, file.path(harmonized_dir, "harmonized_landsat_rr.csv"), row.names = FALSE)
} else {
  warning("Landsat.RR not available; harmonized_landsat_rr.csv not exported.")
}

# --- Validate and export site metadata ---
site_required_cols <- c("CA_MPA_Name_Short", "MPA_Start")
stopifnot("Expected columns missing from Site" =
            all(site_required_cols %in% names(Site)))
cat(sprintf("  harmonized_site_metadata.csv: %d rows x %d cols\n",
            nrow(Site), ncol(Site)))

write.csv(Site, file.path(harmonized_dir, "harmonized_site_metadata.csv"), row.names = FALSE)

# --- Generate and report checksums ---
csv_files <- list.files(harmonized_dir, pattern = "\\.csv$", full.names = TRUE)
checksums <- vapply(csv_files, tools::md5sum, character(1))
writeLines(paste(checksums, basename(names(checksums))),
           file.path(harmonized_dir, "checksums.md5"))

cat("\nHarmonized CSVs exported to:", harmonized_dir, "\n")
for (f in csv_files) {
  size_kb <- round(file.size(f) / 1024, 1)
  cat(sprintf("  %-40s %s KB  [MD5: %s]\n", basename(f), size_kb, checksums[f]))
}
cat("Checksums written to:", file.path(harmonized_dir, "checksums.md5"), "\n")
