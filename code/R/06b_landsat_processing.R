# =============================================================================
# 06b_landsat_processing.R
# =============================================================================
#
# PURPOSE:
#   Process satellite-derived kelp canopy biomass data from Landsat imagery
#   for the MPA pBACIPS analysis.
#
# WHAT THIS SCRIPT DOES:
#   1. Imports Landsat kelp canopy area/biomass data
#   2. Transforms data to long format for analysis
#   3. Calculates annual mean biomass by MPA and status (MPA vs reference)
#   4. Converts to proportions of time series maximum
#   5. Calculates log response ratios
#   6. Assigns time since MPA and Before/After labels
#
# DATA SOURCE:
#   Landsat-derived kelp canopy area from remote sensing analysis.
#   This provides an independent measure of kelp abundance that:
#   - Covers a broader spatial extent than diver surveys
#   - Has consistent methodology back to 1984
#   - Measures canopy area (proxy for biomass) from above
#
# ADVANTAGES OF SATELLITE DATA:
#   - No diver bias or sampling error
#   - Continuous time series without gaps
#   - Can detect large-scale patterns
#
# LIMITATIONS:
#   - Only measures surface canopy (misses understory kelp)
#   - Affected by water clarity and atmospheric conditions
#   - Lower taxonomic resolution (just "kelp", not species)
#
# INPUTS:
#   - data/LANDSAT/MPA_Runs_new.csv
#   - Site object (from 03_data_import.R)
#   - sites.short.edit object
#
# OUTPUTS:
#   - Landsat.RR: Kelp biomass response ratios with metadata
#   - Landsat.resp: Raw kelp biomass data (long format)
#
# DEPENDENCIES:
#   Requires 00-03 scripts to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

cat("Processing Landsat satellite kelp data...\n")

# Source utility functions
source(here::here("code", "R", "01_utils.R"))

# ===========================================================================
# Section 1: Import Landsat kelp biomass data
# ===========================================================================
# Landsat satellite imagery allows us to estimate kelp canopy area remotely.
# The data is organized with one row per MPA-status-replicate combination and
# one column per year (wide format). Each annual value is the MAX (peak) kelp
# canopy biomass observed across the year, in kg wet weight, derived from
# Landsat reflectance via the Bell/SBC LTER kelp-watch product. Annual MAX
# (rather than mean) is used because kelp canopy is highly seasonal and the
# peak reflects the realized standing biomass each year.

landsat_data_path <- here::here("data", "LANDSAT", "MPA_Runs_new.csv")

# Check if file exists before trying to read it
# file.exists() returns TRUE/FALSE so the script handles missing data cleanly
if (!file.exists(landsat_data_path)) {
  warning("Landsat data file not found: ", landsat_data_path,
          "\nLandsat.RR and Landsat.resp will not be created.")
} else {

  fn <- read.csv(landsat_data_path)

  # Validate expected columns in Landsat data
  required_landsat_cols <- c("MPA", "status", "rep", "lat", "lon")
  missing_landsat_cols <- setdiff(required_landsat_cols, names(fn))
  if (length(missing_landsat_cols) > 0) {
    stop("Landsat CSV is missing required columns: ",
         paste(missing_landsat_cols, collapse = ", "),
         "\nCheck that the input file has not changed format.")
  }

  # Verify that year columns exist (should be named X1984, X1985, ..., X2023 etc.)
  year_cols <- grep("^X\\d{4}$", names(fn), value = TRUE)
  if (length(year_cols) == 0) {
    stop("Landsat CSV has no year columns (expected pattern X1984, X1985, ...).")
  }

  cat("  Landsat data loaded:", nrow(fn), "rows,", length(year_cols),
      "year columns (", gsub("^X", "", year_cols[1]), "-",
      gsub("^X", "", year_cols[length(year_cols)]), ")\n")

  # Validate that the Site object exists (needed for BA assignment)
  if (!exists("Site") || !is.data.frame(Site)) {
    stop("'Site' object not found. Ensure 03_data_import.R has been sourced first.")
  }

  # ===========================================================================
  # Section 2: Transform to long format
  # ===========================================================================
  # The raw data is in "wide" format with one column per year.
  # We need "long" format (one row per observation) for analysis.
  # gather() from tidyr converts wide to long: creates Year and Biomass columns.
  # The minus signs mean "don't gather these columns" (keep them as-is).

  bio.summarize.long <- gather(fn, Year, Biomass, -MPA, -status, -rep, -lat, -lon)

  # Derive numeric year directly from gathered column names (e.g., "X1984" -> 1984).
  # This avoids brittle assumptions about gather() row ordering.
  bio.summarize.long$year <- suppressWarnings(as.integer(gsub("^X", "", bio.summarize.long$Year)))
  bad_year <- sum(is.na(bio.summarize.long$year))
  if (bad_year > 0) {
    warning("Landsat year parsing: ", bad_year, " rows have non-year 'Year' values.")
  }
  bio.summarize.long <- bio.summarize.long %>% dplyr::arrange(MPA)  # Sort by MPA name

  stopifnot("Landsat long-format data must have > 0 rows" = nrow(bio.summarize.long) > 0)
  cat("  Landsat long format:", nrow(bio.summarize.long), "rows,",
      length(unique(bio.summarize.long$MPA)), "MPAs\n")

  # ===========================================================================
  # Section 3: Calculate annual mean biomass by MPA and status
  # ===========================================================================
  # Average across replicates within each MPA-status-year combination.
  # This produces a single mean biomass estimate per MPA per year for both
  # MPA (inside) and reference (outside) sites.

  bio.ave.mpa <- bio.summarize.long %>%
    dplyr::group_by(MPA, status, year) %>%
    dplyr::summarise_at(c("Biomass"), mean, na.rm = TRUE) %>%
    dplyr::ungroup()

  # Spread to wide format (one column per status) and remove incomplete cases
  # This creates mpa and reference columns needed for log response ratio calculation
  bio.sum.max.short <- bio.ave.mpa %>% spread(status, Biomass)
  bio.sum.max.short <- bio.sum.max.short[complete.cases(bio.sum.max.short), ]

  cat("  Landsat MPA-year averages:", nrow(bio.sum.max.short), "complete MPA-year pairs\n")

  # ===========================================================================
  # Section 4: Convert to proportions of maximum
  # ===========================================================================
  # Normalize each MPA's time series to proportion of its maximum observed value.
  # This puts all MPAs on a 0-1 scale for comparable response ratio calculation.

  bio.sum.max.long <- gather(bio.sum.max.short, status, meanBio, -MPA, -year)

  # Add required columns for calculate_proportions() utility function
  bio.sum.max.long$CA_MPA_Name_Short <- bio.sum.max.long$MPA
  bio.sum.max.long$y <- "Macrocystis pyrifera"

  # Use the standardized calculate_proportions function from 01_utils.R
  # This ensures consistent proportion calculation and zero-correction across all data sources
  bio.sum.max.long <- calculate_proportions(
    bio.sum.max.long,
    value_col = "meanBio",
    correction_method = "adaptive"  # Use adaptive instead of fixed 0.01 for consistency
  )

  # ===========================================================================
  # Section 5: Calculate log response ratios (MPA vs reference)
  # ===========================================================================
  # lnRR = ln(proportion_MPA / proportion_reference)
  # Positive lnRR means higher biomass inside MPA relative to outside.

  bio.sum.max.long$taxon_name <- "Macrocystis pyrifera"
  # Select columns needed for spread: MPA, year, taxon_name, status, PropCorr
  bio.sum.max.long.sub <- bio.sum.max.long[, c("MPA", "year", "taxon_name", "status", "PropCorr")]
  bio.sum.max.LANDSAT.diff <- bio.sum.max.long.sub %>%
    spread(status, PropCorr)
  bio.sum.max.LANDSAT.diff <- calculate_log_response_ratio(bio.sum.max.LANDSAT.diff)

  cat("  Landsat log response ratios:", nrow(bio.sum.max.LANDSAT.diff), "rows\n")

  # ===========================================================================
  # Section 6: Assign Before/After and time since MPA implementation
  # ===========================================================================

  # FIX [I19]: Include `Diff`, `mpa`, and `reference` columns produced by
  # calculate_log_response_ratio() (via spread + calculate_log_response_ratio).
  # These are expected by downstream RR_STANDARD_COLS in 07_combine_and_export.R.
  bio.sum.max.LANDSAT.diff <- assign_ba_from_site_table(
    data.frame(CA_MPA_Name_Short = bio.sum.max.LANDSAT.diff$MPA,
               year = bio.sum.max.LANDSAT.diff$year,
               lnDiff = bio.sum.max.LANDSAT.diff$lnDiff,
               Diff = bio.sum.max.LANDSAT.diff$Diff,
               mpa = bio.sum.max.LANDSAT.diff$mpa,
               reference = bio.sum.max.LANDSAT.diff$reference),
    Site
  )
  bio.sum.max.LANDSAT.diff <- assign_time_from_site_table(bio.sum.max.LANDSAT.diff, Site)

  # ===========================================================================
  # Section 7: Build Landsat.RR output (matching format of LTER/PISCO/KFM RR datasets)
  # ===========================================================================

  # Reconstruct full response ratio dataframe with columns matching other modules
  Landsat.RR <- bio.sum.max.LANDSAT.diff
  Landsat.RR$y <- "Macrocystis pyrifera"
  Landsat.RR$resp <- "Bio"  # Landsat only provides biomass (canopy area proxy), not density
  Landsat.RR$source <- "Landsat"

  # Merge with site metadata for type, Location, Hectares
  Landsat.RR <- merge(Landsat.RR, sites.short.edit,
                       by.x = "CA_MPA_Name_Short",
                       by.y = "CA_MPA_Name_Short",
                       all.x = TRUE)

  stopifnot("Landsat.RR must have > 0 rows" = nrow(Landsat.RR) > 0)

  # ===========================================================================
  # Section 8: Build Landsat.resp output (raw response data in long format)
  # ===========================================================================

  # Use the wide-format biomass (mpa vs reference columns) before log-ratio transformation
  Landsat.resp.wide <- bio.sum.max.short
  colnames(Landsat.resp.wide)[colnames(Landsat.resp.wide) == "MPA"] <- "CA_MPA_Name_Short"
  Landsat.resp.wide$taxon_name <- "Macrocystis pyrifera"
  Landsat.resp.wide$source <- "Landsat"

  Landsat.resp <- gather(Landsat.resp.wide, status, value, -CA_MPA_Name_Short, -year,
                         -taxon_name, -source)
  Landsat.resp$resp <- "Bio"

  # ===========================================================================
  # Section 9: Clean up intermediate objects
  # ===========================================================================

  rm(fn, bio.summarize.long, bio.ave.mpa, bio.sum.max.short,
     bio.sum.max.long, bio.sum.max.LANDSAT.diff, Landsat.resp.wide)

}

# ===========================================================================
# Section 10: Summary of outputs available for downstream scripts
# ===========================================================================
# Landsat.RR   - Kelp biomass log response ratios with time, BA, site info
#                Columns: CA_MPA_Name_Short, year, lnDiff, BA, time, y, resp, source,
#                         type, Location, Hectares
# Landsat.resp - Raw kelp biomass response data (mpa vs reference, long format)
#                Columns: CA_MPA_Name_Short, year, taxon_name, source, status, value, resp

if (exists("Landsat.RR")) {
  cat("Landsat processing complete.\n")
  cat("  Outputs: Landsat.RR (", nrow(Landsat.RR), " rows), ",
      "Landsat.resp (", nrow(Landsat.resp), " rows)\n", sep = "")
} else {
  cat("Landsat processing skipped (data file not found).\n")
}
