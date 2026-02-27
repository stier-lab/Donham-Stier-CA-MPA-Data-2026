# =============================================================================
# 04_pisco_processing.R
# =============================================================================
#
# PURPOSE:
#   Process PISCO (Partnership for Interdisciplinary Studies of Coastal Oceans)
#   kelp forest monitoring data for the MPA pBACIPS analysis.
#
# WHAT THIS SCRIPT DOES:
#   1. Imports PISCO swath data (urchins, kelp, lobsters)
#   2. Imports PISCO fish survey data (sheephead)
#   3. Calculates density (count per m2) for all taxa
#   4. Estimates biomass via:
#      - Direct conversion (lobster carapace length -> biomass)
#      - Bootstrap resampling from size frequency data (urchins)
#      - Stipe density conversion (kelp)
#   5. Calculates log response ratios (MPA vs reference sites)
#   6. Assigns time since MPA implementation and Before/After labels
#
# DATA SOURCES:
#   PISCO surveys kelp forests at permanent sites using:
#   - Swath transects: Counts of benthic invertebrates and kelp
#   - Fish transects: Visual census of reef fish
#   Two research groups contribute data:
#   - UCSB (UC Santa Barbara): Started 1999, sizes lobsters directly
#   - VRG (Vantuna Research Group): Started 2003, uses size frequency
#
# INPUTS:
#   - data/PISCO/MLPA_kelpforest_swath_2024.csv
#   - data/PISCO/MLPA_kelpforest_fish_2024.csv
#   - data/PISCO/master_site_table_Emilyedit.csv
#   - data/PISCO/spp_attribute_table_downloaded_9-13-22_SHSPUL.csv
#   - SizeFreq (from 03_data_import.R)
#
# OUTPUTS:
#   - Swath.join.sub: Response ratio data with time, BA, MPA features
#   - PISCO.resp: Raw density/biomass data (MPA vs reference, long format)
#   - All_PISCO.mpa: Averaged PISCO data by MPA
#
# DEPENDENCIES:
#   Requires 00_libraries.R, 01_utils.R, 03_data_import.R to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

cat("Processing PISCO data...\n")

# ===========================================================================
# Section 1: Import and Wrangle PISCO Swath Data
# ===========================================================================
# PISCO swath surveys count invertebrates and kelp along transects.

swath_path <- here::here("data", "PISCO", "MLPA_kelpforest_swath_2024.csv")
Swath <- tryCatch(
  read.csv(swath_path),
  error = function(e) stop("Failed to read PISCO swath data: ", swath_path,
                           "\n  ", e$message, call. = FALSE)
)

# Validate expected columns before proceeding
swath_required_cols <- c("campus", "zone", "site", "year", "month",
                         "transect", "classcode", "count", "size")
swath_missing_cols <- setdiff(swath_required_cols, names(Swath))
if (length(swath_missing_cols) > 0) {
  stop("PISCO swath data missing expected columns: ",
       paste(swath_missing_cols, collapse = ", "), call. = FALSE)
}

cat("  Loaded swath data:", nrow(Swath), "rows\n")

# Remove deep sites (only surveyed by VRG, not always kelp forest)
# subset() filters rows based on a logical condition
# The | operator means "OR" and & means "AND"
# This keeps: all VRG data, OR UCSB data that isn't from DEEP/SPECIALO zones
Swath.edit <- subset(Swath,
                     campus == "VRG" | (campus == "UCSB" & zone != "DEEP" & zone != "SPECIALO"))

# Import site table (edited so PISCO references run correctly with this script)
# This table maps PISCO site codes to MPA names and contains metadata about each site
site_table_path <- here::here("data", "PISCO", "master_site_table_Emilyedit.csv")
sites <- tryCatch(
  read.csv(site_table_path),
  error = function(e) stop("Failed to read PISCO site table: ", site_table_path,
                           "\n  ", e$message, call. = FALSE)
)

# Remove duplicate site-MPA combinations to get one row per unique site
# %>% is the "pipe" operator - it passes the result of the left side as input to the right
# filter() keeps rows that meet a condition; !duplicated() keeps only first occurrence
sites.short <- sites %>%
  dplyr::filter(!duplicated(cbind(site, CA_MPA_Name_Short)))

# Transform Macrocystis: count * size = total stipes (size = stipes per individual)
# ifelse() is a vectorized if-else: condition, value if TRUE, value if FALSE
# For giant kelp (MACPYRAD), size = stipes per individual, so count * size = total stipes
Swath.edit$count <- ifelse(Swath.edit$classcode == "MACPYRAD",
                           Swath.edit$count * Swath.edit$size,
                           Swath.edit$count)

# ===========================================================================
# Section 2: Create Complete Species Roster and Fill Zeros
# ===========================================================================
# IMPORTANT: Survey data only records species that were SEEN. If a species wasn't observed,
# there's no row for it. We need to create explicit zeros for "absence" data.
# This is critical for calculating accurate means and response ratios.

# unique() returns distinct combinations - here, all unique survey events
Site.yr <- unique(Swath.edit[, c("site", "year", "month", "zone", "transect", "campus")])
Classcodes <- unique(Swath.edit$classcode)

# Cross-join (by=NULL) creates every possible combination of surveys x species
# This gives us a row for every species at every survey event
Site.All.Classes <- merge(x = Site.yr, y = Classcodes, by = NULL)

# merge() combines two dataframes like a SQL JOIN
# all=TRUE makes this an "outer join" - keeps all rows from both tables
# Surveys without a species get NA for count (which we'll convert to 0)
Swath.site.PISCO <- merge(Site.All.Classes, Swath.edit,
                          by.x = c("site", "year", "month", "zone", "transect", "y", "campus"),
                          by.y = c("site", "year", "month", "zone", "transect", "classcode", "campus"),
                          all = TRUE)

# Zero-fill absent-species count: species not observed at a survey event have NA count
# after the outer join — replace with 0 (true absence, not missing data)
Swath.site.PISCO$count[is.na(Swath.site.PISCO$count)] <- 0

# Subset summer months only (May-October) for consistency across survey methods
# PISCO surveys in summer; other programs also focus on summer, so this keeps comparisons fair
Swath.site.PISCO.sub <- subset(Swath.site.PISCO,
                                month >= SURVEY_MONTH_START & month <= SURVEY_MONTH_END)

# ===========================================================================
# Section 3: Aggregate Transect -> Zone -> Site Annual Means
# ===========================================================================
# PISCO data has a hierarchical structure: transects nested in zones nested in sites.
# We aggregate from the finest scale (transect) up to the site level to get annual means.
# This is the standard approach for nested survey designs.

# group_by() defines groups; summarise_at() calculates summary stats for specified columns
# sum aggregates within groups (e.g., sum lobster counts across size classes per transect)
# ungroup() removes grouping so subsequent operations treat all rows equally
Swath.site.PISCO.sum.int <- Swath.site.PISCO.sub %>%
  dplyr::group_by(year, month, site, zone, transect, y) %>%
  dplyr::summarise_at(c("count"), sum) %>%
  dplyr::ungroup()

# Average across months within each year (some years have multiple survey visits)
# mean() calculates the arithmetic average
Swath.site.PISCO.sum <- Swath.site.PISCO.sum.int %>%
  dplyr::group_by(year, site, zone, transect, y) %>%
  dplyr::summarize_at(c("count"), mean)

# Average transects within zones to get zone-level means
Swath.ave.zone <- Swath.site.PISCO.sum %>%
  dplyr::group_by(year, site, y, zone) %>%
  dplyr::summarise_at(c("count"), mean)

# Average zones within sites to get site-level annual means
# This is our final aggregation level for analysis
Swath.ave <- Swath.ave.zone %>%
  dplyr::group_by(year, site, y) %>%
  dplyr::summarise_at(c("count"), mean)

# ===========================================================================
# Section 4: Merge with Site Table, Filter Southern CA, Subset Target Taxa
# ===========================================================================
# Now we join the aggregated data with site metadata to get MPA names and other info.
# merge() combines dataframes by matching values in specified columns (like VLOOKUP in Excel).
# by.x and by.y specify which columns to match; all.x=TRUE keeps all rows from first table.

Swath.ave.site <- merge(Swath.ave, sites.short, by.x = c("site"), by.y = c("site"),
                        all.x = TRUE)

# Remove sites not associated with an MPA (reference sites outside study area)
# complete.cases() returns TRUE for rows with no NA values in specified column
Swath.ave.site <- Swath.ave.site[complete.cases(Swath.ave.site$CA_MPA_Name_Short), ]

# Keep only needed columns using dplyr::select for robust column selection
# FIXED (2026-02-06): Converted from fragile column indices to explicit column names
Swath.ave.site.sub <- Swath.ave.site %>%
  dplyr::select(site, year, y, count, CA_MPA_Name_Short, site_designation, site_status, BaselineRegion)

# Only Southern California (our study region)
# BaselineRegion categorizes sites by MLPA planning region
Swath.ave.site.PISCO <- subset(Swath.ave.site.sub, BaselineRegion == "SOUTH")

# Only target taxa for this analysis:
# STRPURAD = purple urchin (Strongylocentrotus purpuratus)
# MESFRAAD = red urchin (Mesocentrotus franciscanus)
# MACPYRAD = giant kelp (Macrocystis pyrifera)
# PANINT = California spiny lobster (Panulirus interruptus)
# %in% checks if values are in a set (like "is y one of these?")
Swath.ave.site.PISCO.subset <- subset(Swath.ave.site.PISCO,
                                      y %in% c("STRPURAD", "MESFRAAD", "MACPYRAD", "PANINT"))

# ===========================================================================
# Section 5: Calculate Density (count per m2)
# ===========================================================================
# PISCO swath transects are 60 m2 (2m wide x 30m long).
# Density = count / area gives individuals per square meter for standardization.

Swath.ave.site.PISCO.subset$Density <- Swath.ave.site.PISCO.subset$count / PISCO_SWATH_AREA_M2
Swath.ave.site.PISCO.subset$biomass <- NA  # Placeholder - will be filled per-taxon below

cat("  Swath taxa subset:", nrow(Swath.ave.site.PISCO.subset), "rows\n")

# ===========================================================================
# Section 6a: Process Lobster Biomass -- UCSB (Direct CL-to-Biomass)
# ===========================================================================
# UCSB researchers measure lobster carapace length (CL) directly in the field.
# We can convert CL to biomass using allometric equations (bio_lobster() function).
# This section handles UCSB data separately because they have individual size data.

cat("  Processing UCSB lobster biomass...\n")

# Go back to raw transect-level data (before aggregation) to get size information
# We need individual sizes to calculate biomass accurately
Swath.site.PISCO.PANINT <- subset(Swath.site.PISCO.sub, y == "PANINT")

# Join lobster data with site metadata
PANINT <- merge(Swath.site.PISCO.PANINT, sites.short, by.x = c("site"), by.y = c("site"),
                all.x = TRUE)

# UCSB started measuring carapace length in UCSB_LOBSTER_SIZE_START_YEAR (2010).
# Before that, they only counted lobsters without measuring size.
PANINT.UCSB <- subset(PANINT, year >= UCSB_LOBSTER_SIZE_START_YEAR & campus.x == "UCSB")

# Convert CL (cm in data) to biomass using bio_lobster() from 01_utils.R
# bio_lobster() expects mm, so multiply by 10 to convert cm -> mm
# biomass = weight per individual * number of individuals
PANINT.UCSB$biomass <- bio_lobster(PANINT.UCSB$size * 10) * PANINT.UCSB$count

# Select only the columns we need for further processing
panint_ucsb_needed <- c(
  "year", "month", "site", "zone", "transect", "y",
  "count", "biomass",
  "CA_MPA_Name_Short", "site_designation", "site_status", "BaselineRegion"
)
panint_ucsb_missing <- setdiff(panint_ucsb_needed, names(PANINT.UCSB))
if (length(panint_ucsb_missing) > 0) {
  stop("PISCO PANINT.UCSB: missing expected columns: ",
       paste(panint_ucsb_missing, collapse = ", "),
       call. = FALSE)
}
PANINT.UCSB.sub <- PANINT.UCSB[, panint_ucsb_needed]

# DATA QUALITY CHECK: Find site-years where lobsters were counted but not sized.
# count > 0 means lobsters were seen; biomass == 0 means no size was recorded.
# We can't estimate accurate biomass without sizes, so remove entire site-years.
PANINT.UCSB.null <- subset(PANINT.UCSB.sub, count > 0 & biomass == 0)

# PERFORMANCE: Use anti-join instead of iterative subsetting to avoid O(n^2) copies.
# Build a key of site-year combos to exclude, then filter in one pass.
if (nrow(PANINT.UCSB.null) > 0) {
  null_keys <- unique(PANINT.UCSB.null[, c("site", "year")])
  null_keys$.exclude <- TRUE
  PANINT.UCSB.sub.sub <- merge(PANINT.UCSB.sub, null_keys,
                                 by = c("site", "year"), all.x = TRUE)
  PANINT.UCSB.sub.sub <- PANINT.UCSB.sub.sub[is.na(PANINT.UCSB.sub.sub$.exclude), ]
  PANINT.UCSB.sub.sub$.exclude <- NULL
  cat("    Removed", nrow(null_keys), "site-years with unsized lobsters (UCSB)\n")
} else {
  PANINT.UCSB.sub.sub <- PANINT.UCSB.sub
}

# Aggregate: sum within transect, average across months, average across zones
PANINT.All.sub.zone.sum <- PANINT.UCSB.sub.sub %>%
  dplyr::group_by(year, month, site, CA_MPA_Name_Short, site_designation, site_status,
           BaselineRegion, zone, transect, y) %>%
  dplyr::summarise_at(c("count", "biomass"), sum)

PANINT.All.sub.zone.mean <- PANINT.All.sub.zone.sum %>%
  dplyr::group_by(year, site, CA_MPA_Name_Short, site_designation, site_status,
           BaselineRegion, zone, y) %>%
  dplyr::summarise_at(c("count", "biomass"), mean)

PANINT.All.sub.mean <- PANINT.All.sub.zone.mean %>%
  dplyr::group_by(year, site, CA_MPA_Name_Short, site_designation, site_status,
           BaselineRegion, y) %>%
  dplyr::summarise_at(c("count", "biomass"), mean)

PANINT.All.sub.mean$Density <- PANINT.All.sub.mean$count / PISCO_SWATH_AREA_M2
PANINT.All.sub.mean$y <- "PANINT"
PANINT.All.sub.mean$campus <- "UCSB"

cat("    UCSB lobster: ", nrow(PANINT.All.sub.mean), " site-year records\n", sep = "")

# ===========================================================================
# Section 6b: Process Lobster Biomass -- VRG (Bootstrap Resampling)
# ===========================================================================
# VRG (Vantuna Research Group) doesn't measure every individual lobster.
# Instead, they measure a sample and record size frequencies (distribution of sizes).
# We use BOOTSTRAP RESAMPLING to estimate biomass:
#   1. For each site-year, get the count of lobsters seen
#   2. Randomly sample that many sizes from the size frequency distribution
#   3. Convert sizes to biomass
#   4. Repeat 1000 times and take the average
# This accounts for uncertainty in size distributions.

cat("  Processing VRG lobster biomass (bootstrap resampling)...\n")

# VRG started measuring sizes in VRG_LOBSTER_SIZE_START_YEAR (2011),
# stored in a separate size frequency table
PANINT.VRG <- subset(PANINT, year >= VRG_LOBSTER_SIZE_START_YEAR & campus.x == "VRG")

# Merge size frequency data with site metadata
# SizeFreq contains measured sizes from a subsample of individuals
SizeFreq.PANINT.VRG <- merge(SizeFreq, sites.short,
                              by.x = c("site", "site_new", "campus"),
                              by.y = c("site", "site_new", "campus"),
                              all.x = TRUE)
SizeFreq.PANINT.VRG <- subset(SizeFreq.PANINT.VRG, classcode == "PANINT")

# CACHING: Bootstrap resampling is slow (1000 iterations per site-year).
# We save results to an RDS file so we only run it once.
# .rds files store R objects; readRDS() loads them back into memory.
# If the cache file exists and FORCE_BOOTSTRAP isn't set, use cached results.
.cache_vrg_lob <- here::here("data", "cache", "vrg_panint_bootstrap.rds")
force_boot <- should_force_bootstrap()
if (file.exists(.cache_vrg_lob) && !force_boot) {
  cat("    Loading cached VRG lobster bootstrap results...\n")
  VRG.PANINT.site <- safe_read_rds(.cache_vrg_lob)
  if (is.null(VRG.PANINT.site)) {
    cat("    Cache invalid, will recompute...\n")
  }
}

if (!exists("VRG.PANINT.site") || is.null(VRG.PANINT.site)) {
  # Get unique zone-level records to iterate over
  # Each row = one site-year-zone combination that needs bootstrap estimation
  u <- unique(PANINT.VRG[, c("year", "site", "CA_MPA_Name_Short", "site_designation",
                              "site_status", "BaselineRegion", "zone")])

  # Pre-allocate list for results (PERFORMANCE: avoids O(n²) rbind copies)
  results_list <- vector("list", nrow(u))

  # Set random seed for reproducibility of bootstrap resampling
  set.seed(SEED_PISCO_LOBSTER)
  for (i in seq_len(nrow(u))) {
    t <- which(PANINT.VRG$site == u$site[i] &
               PANINT.VRG$year == u$year[i] &
               PANINT.VRG$zone == u$zone[i])
    n <- sum(PANINT.VRG$count[t])
    t2 <- which(SizeFreq.PANINT.VRG$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] &
                SizeFreq.PANINT.VRG$year == u$year[i] &
                SizeFreq.PANINT.VRG$site_status == u$site_status[i])

    # Wrap bootstrap in tryCatch for robustness
    result <- tryCatch(
      bootstrap_biomass(count = n,
                        size_freq_indices = t2,
                        size_freq_table = SizeFreq.PANINT.VRG,
                        biomass_fun = bio_lobster),
      error = function(e) {
        warning("Bootstrap failed for site ", u$site[i], " year ", u$year[i], ": ", e$message)
        list(biomass = NA, se = NA, count = n)
      }
    )

    results_list[[i]] <- data.frame(site = u$site[i],
                                     year = u$year[i],
                                     zone = u$zone[i],
                                     transect = length(t),
                                     biomass = result$biomass,
                                     count = result$count)
  }

  # Combine all results at once (efficient)
  VRG.PANINT.site <- do.call(rbind, results_list)

  # Save cache
  dir.create(dirname(.cache_vrg_lob), recursive = TRUE, showWarnings = FALSE)
  saveRDS(VRG.PANINT.site, .cache_vrg_lob)
  cat("    Cached VRG lobster bootstrap results.\n")
}

# Merge with site table
VRG.PANINT.site.merge <- merge(VRG.PANINT.site, sites.short,
                                by.x = c("site"), by.y = c("site"),
                                all.x = TRUE)

vrg_panint_needed <- c(
  "year", "site", "CA_MPA_Name_Short",
  "site_designation", "site_status", "BaselineRegion",
  "zone", "transect",
  "count", "biomass"
)
vrg_panint_missing <- setdiff(vrg_panint_needed, names(VRG.PANINT.site.merge))
if (length(vrg_panint_missing) > 0) {
  stop("PISCO VRG PANINT merge: missing expected columns: ",
       paste(vrg_panint_missing, collapse = ", "),
       call. = FALSE)
}
VRG.PANINT.site.merge <- VRG.PANINT.site.merge[, vrg_panint_needed]
# Zero-fill only numeric columns; preserve NAs in character/factor columns
num_cols <- sapply(VRG.PANINT.site.merge, is.numeric)
VRG.PANINT.site.merge[num_cols][is.na(VRG.PANINT.site.merge[num_cols])] <- 0

# Remove site-years with lobsters seen but not measured (same logic as UCSB above)
PANINT.VRG.null <- subset(VRG.PANINT.site.merge, count > 0 & biomass == 0)

# PERFORMANCE: Anti-join instead of iterative subsetting
if (nrow(PANINT.VRG.null) > 0) {
  null_keys <- unique(PANINT.VRG.null[, c("site", "year")])
  null_keys$.exclude <- TRUE
  VRG.PANINT.site.merge.sub <- merge(VRG.PANINT.site.merge, null_keys,
                                      by = c("site", "year"), all.x = TRUE)
  VRG.PANINT.site.merge.sub <- VRG.PANINT.site.merge.sub[is.na(VRG.PANINT.site.merge.sub$.exclude), ]
  VRG.PANINT.site.merge.sub$.exclude <- NULL
  cat("    Removed", nrow(null_keys), "site-years with unsized lobsters (VRG)\n")
} else {
  VRG.PANINT.site.merge.sub <- VRG.PANINT.site.merge
}

# Sum across transects, normalize by number of transects
VRG.PANINT.site.all.TransectSum <- VRG.PANINT.site.merge.sub %>%
  dplyr::group_by(year, site, CA_MPA_Name_Short, site_designation, site_status,
           BaselineRegion, zone, transect) %>%
  dplyr::summarise_at(c("count", "biomass"), sum, na.rm = TRUE)

VRG.PANINT.site.all.TransectSum$count <- VRG.PANINT.site.all.TransectSum$count /
  VRG.PANINT.site.all.TransectSum$transect
VRG.PANINT.site.all.TransectSum$biomass <- VRG.PANINT.site.all.TransectSum$biomass /
  VRG.PANINT.site.all.TransectSum$transect

# Average across zones
VRG.PANINT.site.all <- VRG.PANINT.site.all.TransectSum %>%
  dplyr::group_by(year, site, CA_MPA_Name_Short, site_designation, site_status, BaselineRegion) %>%
  dplyr::summarise_at(c("count", "biomass"), mean)

VRG.PANINT.site.all$Density <- VRG.PANINT.site.all$count / PISCO_SWATH_AREA_M2
VRG.PANINT.site.all$y <- "PANINT"

VRG.PANINT <- VRG.PANINT.site.all
# Reorder columns to match PANINT.All.sub.mean structure for rbind
# FIXED: Replaced fragile column indices c(1:6, 10, 7:9) with explicit column names
VRG.PANINT <- dplyr::select(VRG.PANINT,
  year, site, CA_MPA_Name_Short, site_designation, site_status, BaselineRegion,
  y, count, biomass, Density)
VRG.PANINT$campus <- "VRG"

cat("    VRG lobster: ", nrow(VRG.PANINT), " site-year records\n", sep = "")

# ===========================================================================
# Section 6c: Merge Lobster Datasets (UCSB + VRG)
# ===========================================================================

PANINT.all <- rbind(VRG.PANINT, PANINT.All.sub.mean)
# Reorder columns for downstream rbind with Swath.ave.site.PISCO.subset
# FIXED: Replaced fragile column indices c(2, 1, 7, 8, 11, 3:6, 10, 9) with explicit column names
# FIX [C5]: Removed `campus` — Swath.ave.site.PISCO.subset does not have it,
# so including it here would cause an rbind column mismatch at line 416.
# `campus` is not needed downstream (dropped at line 718 column selection).
PANINT.all <- dplyr::select(PANINT.all,
  site, year, y, count, CA_MPA_Name_Short, site_designation,
  site_status, BaselineRegion, Density, biomass)

Swath.PISCO.subset <- rbind(PANINT.all, Swath.ave.site.PISCO.subset)
Swath.PISCO <- Swath.PISCO.subset[!is.na(Swath.PISCO.subset$count), ]

cat("  Combined swath + lobster:", nrow(Swath.PISCO), "rows\n")

# ===========================================================================
# Section 7: Process Urchin Biomass via Bootstrap Resampling
# ===========================================================================

cat("  Processing urchin biomass (bootstrap resampling)...\n")

Swath.Urchin.Bio <- subset(Swath.site.PISCO.sum.int, y == "STRPURAD" | y == "MESFRAAD")
URCHINS <- merge(Swath.Urchin.Bio, sites.short, by.x = c("site"), by.y = c("site"),
                 all.x = TRUE)

# Prepare urchin size frequency data
SizeFreq.Urch <- subset(SizeFreq, classcode == "MESFRA" | classcode == "STRPUR")
SizeFreq.Urch <- dplyr::mutate(SizeFreq.Urch, classcode = dplyr::case_when(
  classcode == "STRPUR"  ~ "STRPURAD",
  classcode == "MESFRA"  ~ "MESFRAAD",
  TRUE ~ classcode))

SizeFreq.Urch <- merge(SizeFreq.Urch, sites.short,
                        by.x = c("site", "site_new", "campus"),
                        by.y = c("site", "site_new", "campus"),
                        all.x = TRUE)

# SIZE FILTER: Only urchins >= 25mm (PISCO minimum size threshold)
#
# RATIONALE: PISCO field protocols only record urchins >= 25mm test diameter
# (PISCO_URCHIN_MIN_SIZE_MM = 25, defined in 00c_analysis_constants.R).
# Smaller urchins are difficult to see and identify in the field, and their
# counts would be inconsistent across observers. This 25mm threshold is a
# PISCO-specific methodological standard.
#
# CROSS-PROGRAM NOTE: This 25mm filter is applied ONLY to PISCO's own
# urchin biomass estimation (in this script). KFM (05_kfm_processing.R)
# and LTER (06_lter_processing.R) do NOT apply this size filter because:
#   (1) KFM and LTER do not collect individual urchin sizes — they only
#       record total counts, which may include urchins < 25mm.
#   (2) To estimate KFM/LTER urchin biomass, we bootstrap-resample from
#       the UNFILTERED PISCO size-frequency data (SizeFreq.Urch.OG, saved
#       below). Using the unfiltered data ensures the bootstrap distribution
#       includes small urchins that KFM/LTER field crews would have counted.
#   (3) The filtered data (SizeFreq.Urch, after the subset below) is only
#       used for PISCO's own biomass estimation in this script.
#
# This creates a minor inconsistency: PISCO biomass is estimated from
# urchins >= 25mm only, while KFM/LTER biomass is estimated from the full
# size distribution. In practice, urchins < 25mm contribute negligible
# biomass (the allometric relationship is approximately cubic in diameter),
# so this has minimal effect on biomass estimates.

# Preserve unfiltered size frequency data for KFM and LTER BEFORE applying
# the 25mm PISCO-specific filter
SizeFreq.Urch.OG <- SizeFreq.Urch

SizeFreq.Urch <- subset(SizeFreq.Urch, size >= PISCO_URCHIN_MIN_SIZE_MM)

# Cache bootstrap results to avoid re-running the slow loop
.cache_urchin <- here::here("data", "cache", "pisco_urchin_bootstrap.rds")
if (file.exists(.cache_urchin) && !force_boot) {
  cat("    Loading cached urchin bootstrap results...\n")
  Urchin.site <- safe_read_rds(.cache_urchin)
  if (is.null(Urchin.site)) {
    cat("    Cache invalid, will recompute...\n")
  }
}

if (!exists("Urchin.site") || is.null(Urchin.site)) {
  # Unique zone-level records to iterate over
  u <- unique(URCHINS[, c("year", "site", "site_new", "CA_MPA_Name_Short",
                           "site_designation", "site_status", "BaselineRegion", "zone", "y")])

  # Pre-allocate list for results (PERFORMANCE: avoids O(n²) rbind copies)
  results_list <- vector("list", nrow(u))

  # Set random seed for reproducibility of bootstrap resampling
  set.seed(SEED_PISCO_URCHIN)
  for (i in seq_len(nrow(u))) {
    t <- which(URCHINS$site == u$site[i] &
               URCHINS$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] &
               URCHINS$year == u$year[i] &
               URCHINS$zone == u$zone[i] &
               URCHINS$y == u$y[i])
    n <- sum(URCHINS$count[t])
    t2 <- which(SizeFreq.Urch$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] &
                SizeFreq.Urch$site_status == u$site_status[i] &
                SizeFreq.Urch$year == u$year[i] &
                SizeFreq.Urch$classcode == u$y[i])

    # Choose biomass function based on species
    bio_fun <- if (u$y[i] == "MESFRAAD") bio_redurch else bio_purpurch

    # Wrap bootstrap in tryCatch for robustness
    result <- tryCatch(
      bootstrap_biomass(count = n,
                        size_freq_indices = t2,
                        size_freq_table = SizeFreq.Urch,
                        biomass_fun = bio_fun),
      error = function(e) {
        warning("Bootstrap failed for site ", u$site[i], " year ", u$year[i],
                " species ", u$y[i], ": ", e$message)
        list(biomass = NA, se = NA, count = n)
      }
    )

    results_list[[i]] <- data.frame(site = u$site[i],
                                     CA_MPA_Name_Short = u$CA_MPA_Name_Short[i],
                                     site_status = u$site_status[i],
                                     year = u$year[i],
                                     zone = u$zone[i],
                                     transect = length(t),
                                     count = result$count,
                                     biomass = result$biomass,
                                     y = u$y[i],
                                     site_new = u$site_new[i])
  }

  # Combine all results at once (efficient)
  Urchin.site <- do.call(rbind, results_list)

  # Save cache
  dir.create(dirname(.cache_urchin), recursive = TRUE, showWarnings = FALSE)
  saveRDS(Urchin.site, .cache_urchin)
  cat("    Cached urchin bootstrap results.\n")
}

# Merge with site table for additional columns
PISCO.Urchin.site.merge <- merge(Urchin.site, sites.short,
                                  by.x = c("site", "CA_MPA_Name_Short", "site_new"),
                                  by.y = c("site", "CA_MPA_Name_Short", "site_new"),
                                  all.x = TRUE)

# Select and reorder columns from merge result
# FIXED: Replaced fragile column indices c(2, 1, 3, 17, 5:10, 16, 20) with explicit column names
# site_status.y comes from sites.short (the authoritative source); site_designation and
# BaselineRegion also come from sites.short via the merge
PISCO.Urchin.site.merge <- dplyr::select(PISCO.Urchin.site.merge,
  CA_MPA_Name_Short, site, site_new, site_status.y,
  year, zone, transect, count, biomass, y,
  site_designation, BaselineRegion)
# Zero-fill only numeric columns; preserve NAs in character/factor columns
num_cols <- sapply(PISCO.Urchin.site.merge, is.numeric)
PISCO.Urchin.site.merge[num_cols][is.na(PISCO.Urchin.site.merge[num_cols])] <- 0

# Remove site-year-species combos where urchins were seen but not measured
PISCO.Urchin.null <- subset(PISCO.Urchin.site.merge, count > 0 & biomass == 0)

# PERFORMANCE: Anti-join instead of iterative subsetting (includes species key
# because both STRPURAD and MESFRAAD are processed together)
if (nrow(PISCO.Urchin.null) > 0) {
  null_keys <- unique(PISCO.Urchin.null[, c("site", "year", "y")])
  null_keys$.exclude <- TRUE
  PISCO.Urchin.site.merge.sub <- merge(PISCO.Urchin.site.merge, null_keys,
                                        by = c("site", "year", "y"), all.x = TRUE)
  PISCO.Urchin.site.merge.sub <- PISCO.Urchin.site.merge.sub[
    is.na(PISCO.Urchin.site.merge.sub$.exclude), ]
  PISCO.Urchin.site.merge.sub$.exclude <- NULL
  cat("    Removed", nrow(null_keys), "site-year-species combos with unsized urchins\n")
} else {
  PISCO.Urchin.site.merge.sub <- PISCO.Urchin.site.merge
}

names(PISCO.Urchin.site.merge.sub)[names(PISCO.Urchin.site.merge.sub) == "site_status.y"] <- "site_status"

# Sum across transects, normalize, average across zones
PISCO.Urchin.site.transectSum <- PISCO.Urchin.site.merge.sub %>%
  dplyr::group_by(year, site, CA_MPA_Name_Short, site_designation, site_status,
           BaselineRegion, zone, transect, y) %>%
  dplyr::summarise_at(c("count", "biomass"), sum, na.rm = TRUE)

PISCO.Urchin.site.transectSum$count <- PISCO.Urchin.site.transectSum$count /
  PISCO.Urchin.site.transectSum$transect
PISCO.Urchin.site.transectSum$biomass <- PISCO.Urchin.site.transectSum$biomass /
  PISCO.Urchin.site.transectSum$transect

PISCO.Urchin.site.all <- PISCO.Urchin.site.transectSum %>%
  dplyr::group_by(year, site, CA_MPA_Name_Short, site_designation, site_status,
           BaselineRegion, y) %>%
  dplyr::summarise_at(c("count", "biomass"), mean)

PISCO.Urchin.site.all$Density <- PISCO.Urchin.site.all$count / PISCO_SWATH_AREA_M2

# Merge urchin biomass back into main swath data
Swath.PISCO.allbio <- merge(Swath.PISCO, PISCO.Urchin.site.all,
                            by.x = c("site", "year", "y", "CA_MPA_Name_Short"),
                            by.y = c("site", "year", "y", "CA_MPA_Name_Short"),
                            all.x = TRUE)

# Fill in biomass from urchin estimates where swath biomass was NA
Swath.PISCO.allbio$biomass.x <- ifelse(is.na(Swath.PISCO.allbio$biomass.x),
                                       Swath.PISCO.allbio$biomass.y,
                                       Swath.PISCO.allbio$biomass.x)

# Keep merge keys + left-side (.x) columns, rename to drop .x suffix
# Also keep count.y (urchin count) for completeness
# FIXED: Replaced fragile column indices c(1:11, 15) with explicit column names and renames
Swath.PISCO <- dplyr::select(Swath.PISCO.allbio,
  site, year, y, CA_MPA_Name_Short,
  count = count.x,
  site_designation = site_designation.x,
  site_status = site_status.x,
  BaselineRegion = BaselineRegion.x,
  Density = Density.x,
  biomass = biomass.x,
  count.y)

# ===========================================================================
# Section 8: Import and Process PISCO Fish Data (Sheephead)
# ===========================================================================

cat("  Processing PISCO fish data (sheephead)...\n")

fish_att_path <- here::here("data", "PISCO", "spp_attribute_table_downloaded_9-13-22_SHSPUL.csv")
FishAtt <- tryCatch(
  read.csv(fish_att_path),
  error = function(e) stop("Failed to read fish attribute table: ", fish_att_path,
                           "\n  ", e$message, call. = FALSE)
)

fish_path <- here::here("data", "PISCO", "MLPA_kelpforest_fish_2024.csv")
Fish <- tryCatch(
  read.csv(fish_path),
  error = function(e) stop("Failed to read PISCO fish data: ", fish_path,
                           "\n  ", e$message, call. = FALSE)
)

# Validate expected columns
fish_required_cols <- c("campus", "classcode", "year", "month", "site",
                        "zone", "level", "transect", "count", "fish_tl")
fish_missing_cols <- setdiff(fish_required_cols, names(Fish))
if (length(fish_missing_cols) > 0) {
  stop("PISCO fish data missing expected columns: ",
       paste(fish_missing_cols, collapse = ", "), call. = FALSE)
}

cat("  Loaded fish data:", nrow(Fish), "rows\n")
Fish <- subset(Fish, campus == "UCSB" | campus == "VRG")

# Merge fish with species attributes for length-weight parameters
Fish.merged <- merge(Fish, FishAtt, by.x = "classcode", by.y = "pisco_classcode",
                     all.x = TRUE)
Fish.merged <- Fish.merged[complete.cases(Fish.merged$WL_input_length), ]

# ===========================================================================
# Section 9: Calculate Fish Biomass from Length-Weight Relationships
# ===========================================================================
# Fish biomass is estimated using allometric relationships: biomass = a * L^b
# where L is body length (total length TL or standard length SL).
# Different species require different length conversions (TYPICAL vs REVERSE)
# to get from observed TL to the length unit expected by the allometric equation.

Fish.merged$biomass <- NA
for (i in seq_len(nrow(Fish.merged))) {
  if (Fish.merged$WL_input_length[i] == "TL") {
    # bio = N * (a * TL^b)
    Fish.merged$biomass[i] <- Fish.merged$count[i] *
      (Fish.merged$WL_a[i] * (Fish.merged$fish_tl[i] ^ Fish.merged$WL_b[i]))
  } else if (Fish.merged$WL_input_length[i] == "SL" &
             Fish.merged$LC_type_for_WL[i] == "TYPICAL") {
    # bio = N * a * (TL * LCa + LCb)^b
    Fish.merged$biomass[i] <- Fish.merged$count[i] *
      (Fish.merged$WL_a[i] *
       ((Fish.merged$fish_tl[i] * Fish.merged$LC.a._for_WL[i] +
         Fish.merged$LC.b._for_WL[i]) ^ Fish.merged$WL_b[i]))
  } else if (Fish.merged$WL_input_length[i] == "SL" &
             Fish.merged$LC_type_for_WL[i] == "REVERSE") {
    # bio = N * a * ((TL - LCb)/LCa)^b
    Fish.merged$biomass[i] <- Fish.merged$count[i] *
      (Fish.merged$WL_a[i] *
       (((Fish.merged$fish_tl[i] - Fish.merged$LC.b._for_WL[i]) /
         Fish.merged$LC.a._for_WL[i]) ^ Fish.merged$WL_b[i]))
  } else {
    Fish.merged$biomass[i] <- NA
  }
}

# Subset to needed columns
Fish.sub <- dplyr::select(Fish.merged, classcode, year, month, site, zone, level,
                   transect, count, fish_tl, TrophicGroup, BroadTrophic, Targeted, biomass)

# Join with site table
Fish.sub.site <- merge(Fish.sub, sites.short, by.x = c("site"), by.y = c("site"),
                       all.x = TRUE)
# Keep original Fish.sub columns + key site metadata from sites.short
# FIXED: Replaced fragile column indices c(1:13, 20:22, 25) with explicit column names
Fish.sub.site <- dplyr::select(Fish.sub.site,
  site, classcode, year, month, zone, level, transect, count, fish_tl,
  TrophicGroup, BroadTrophic, Targeted, biomass,
  CA_MPA_Name_Short, site_designation, site_status, BaselineRegion)

# Remove canopy level (irrelevant for sheephead) and filter to study region
Fish.sub.site.NOcanopy <- subset(Fish.sub.site, level == "BOT")
Fish.sub.site.NOcanopy <- subset(Fish.sub.site.NOcanopy, BaselineRegion == "SOUTH")

# Sum biomass/counts across levels within each site-year-zone-transect
Fish.sum <- Fish.sub.site.NOcanopy %>%
  dplyr::group_by(CA_MPA_Name_Short, site, year, zone, transect, classcode,
           site_status, site_designation) %>%
  dplyr::summarise_at(c("biomass", "count"), sum)

# Create complete species roster for fish (fill zeros)
Site.yr.fish <- unique(Fish.sum[, c("site", "year", "zone", "transect",
                                     "CA_MPA_Name_Short", "site_status", "site_designation")])
Classcodes.fish <- unique(Fish.sum[, c("classcode")])
Site.All.Classes.fish <- merge(x = Site.yr.fish, y = Classcodes.fish, by = NULL)

PISCO.fish <- merge(Site.All.Classes.fish, Fish.sum,
                    by.x = c("site", "year", "zone", "transect", "CA_MPA_Name_Short",
                             "site_status", "classcode", "site_designation"),
                    by.y = c("site", "year", "zone", "transect", "CA_MPA_Name_Short",
                             "site_status", "classcode", "site_designation"),
                    all = TRUE)
# Zero-fill only numeric columns; preserve NAs in character/factor columns
num_cols <- sapply(PISCO.fish, is.numeric)
PISCO.fish[num_cols][is.na(PISCO.fish[num_cols])] <- 0

# Average: transect -> zone -> site
Fish.mean <- PISCO.fish %>%
  dplyr::group_by(CA_MPA_Name_Short, site, year, zone, classcode, site_status, site_designation) %>%
  dplyr::summarise_at(c("biomass", "count"), mean) %>%
  dplyr::ungroup()

Fish.mean.final <- Fish.mean %>%
  dplyr::group_by(CA_MPA_Name_Short, site, year, classcode, site_status, site_designation) %>%
  dplyr::summarise_at(c("biomass", "count"), mean) %>%
  dplyr::ungroup()

# Subset for sheephead (California sheephead, Semicossyphus pulcher)
Fish.SPUL.All <- subset(Fish.mean.final, classcode == "SPUL")

# TRANSECT AREA NOTE: Sheephead density is divided by PISCO_SWATH_AREA_M2
# (60 m^2), the same area used for benthic invertebrate swath transects.
# PISCO fish transects use the same 2m x 30m swath dimensions as benthic
# transects at most sites (see PISCO protocol documentation). The fish are
# counted along the same transect corridors as invertebrates, just recorded
# in a separate dataset with additional size/level information. The "level"
# column distinguishes bottom (BOT), midwater (MID), and canopy (CAN)
# observations; we filter to BOT only (line 738 above) for consistency with
# benthic survey methodology. Using PISCO_SWATH_AREA_M2 here is appropriate
# because PISCO fish and swath surveys share the same transect footprint.
Fish.SPUL.All$density <- Fish.SPUL.All$count / PISCO_SWATH_AREA_M2
# Rename 'classcode' column to 'y' to match swath data species column name
names(Fish.SPUL.All)[names(Fish.SPUL.All) == "classcode"] <- "y"

cat("  Sheephead records:", nrow(Fish.SPUL.All), "\n")

# ===========================================================================
# Section 10: Combine All PISCO Data and Filter for Minimum Data Coverage
# ===========================================================================

# Select and reorder columns to match Fish.SPUL.All structure for rbind
# Column order: CA_MPA_Name_Short, site, year, y, site_status, site_designation, biomass, count
# NOTE: Column indices depend on merge order - verify if upstream changes occur
swath_pisco_needed <- c(
  "CA_MPA_Name_Short", "site", "year", "y",
  "site_status", "site_designation",
  "biomass", "count"
)
swath_pisco_missing <- setdiff(swath_pisco_needed, names(Swath.PISCO))
if (length(swath_pisco_missing) > 0) {
  stop("PISCO Swath.PISCO: missing expected columns: ",
       paste(swath_pisco_missing, collapse = ", "),
       call. = FALSE)
}
Swath.PISCO.sub <- Swath.PISCO[, swath_pisco_needed]
Swath.PISCO.sub$density <- Swath.PISCO.sub$count / PISCO_SWATH_AREA_M2

All_PISCO <- rbind(Swath.PISCO.sub, Fish.SPUL.All)

cat("  Combined all PISCO taxa:", nrow(All_PISCO), "rows\n")

# Remove sites with fewer than MINIMUM_YEARS_OF_DATA years of data.
# Sites with too few years provide unreliable trend estimates for pBACIPS.
site_years <- All_PISCO %>%
  dplyr::distinct(site, year) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(n_years = dplyr::n()) %>%
  dplyr::filter(n_years >= MINIMUM_YEARS_OF_DATA)

All_PISCO.short <- All_PISCO[All_PISCO$site %in% site_years$site, ]

n_sites_dropped <- length(unique(All_PISCO$site)) - length(unique(All_PISCO.short$site))
cat("  Dropped", n_sites_dropped, "sites with <", MINIMUM_YEARS_OF_DATA,
    "years; kept", length(unique(All_PISCO.short$site)), "sites\n")

# ===========================================================================
# Section 11: Remove Excluded Sites and MPAs
# ===========================================================================
# EXCLUDED_REFERENCE_SITES and EXCLUDED_MPAS are defined in 01_utils.R.
# Reference site exclusions come from expert review (Scott Hamilton).
# MPA exclusions are based on data coverage or study design criteria.

All_PISCO.short.sub <- All_PISCO.short

# Remove poor reference sites (per Scott Hamilton)
All_PISCO.short.sub <- All_PISCO.short.sub[
  !All_PISCO.short.sub$site %in% EXCLUDED_REFERENCE_SITES, ]

# Remove MPAs that do not meet inclusion criteria
All_PISCO.short.sub <- All_PISCO.short.sub[
  !All_PISCO.short.sub$CA_MPA_Name_Short %in% EXCLUDED_MPAS, ]

cat("  After exclusions:", nrow(All_PISCO.short.sub), "rows,",
    length(unique(All_PISCO.short.sub$CA_MPA_Name_Short)), "MPAs\n")

# ===========================================================================
# Section 12: Assign Time Since MPA Implementation
# ===========================================================================

# NOTE: Re-reading Site from CSV overwrites the merged version from 03_data_import.R.
# This is intentional for PISCO-specific processing. The merged version (with type, Location,
# Hectares) is recreated via sites.short.edit below.
site_list_path <- here::here("data", "Site_List_All.csv")
Site <- tryCatch(
  read.csv(site_list_path),
  error = function(e) stop("Failed to read site list: ", site_list_path,
                           "\n  ", e$message, call. = FALSE)
)
stopifnot("Site table must contain 'CA_MPA_Name_Short' and 'MPA_Start' columns" =
            all(c("CA_MPA_Name_Short", "MPA_Start") %in% names(Site)))

# Use utility function instead of hundreds of lines of if-else year mapping
All_PISCO.short.sub <- assign_time_from_site_table(All_PISCO.short.sub, Site)

# ===========================================================================
# Section 13: Calculate Macrocystis Biomass from Stipe Density
# ===========================================================================
# Giant kelp (Macrocystis pyrifera) biomass is estimated from stipe density
# using an empirical conversion factor (bio_macro() in 01_utils.R).
# This is applied after time assignment because it only applies to MACPYRAD.

# bio_macro() from utils: stipe_density * ave_slope * 1000
All_PISCO.short.sub$biomass <- ifelse(All_PISCO.short.sub$y == "MACPYRAD",
                                      bio_macro(All_PISCO.short.sub$density),
                                      All_PISCO.short.sub$biomass)

# ===========================================================================
# Section 14: Average Across MPA/Reference Sites
# ===========================================================================
# Collapse individual sites into MPA-level means for response ratio calculation.
# Each MPA may have multiple survey sites inside and outside the reserve.

All_PISCO.mpa <- All_PISCO.short.sub %>%
  dplyr::group_by(CA_MPA_Name_Short, year, y, site_status, site_designation) %>%
  dplyr::summarise_at(c("biomass", "count", "density"), mean, na.rm = TRUE) %>%
  dplyr::ungroup()

# ===========================================================================
# Section 15: Calculate Proportions and Log Response Ratios
# ===========================================================================
# LOG RESPONSE RATIO (lnRR) is our primary metric for comparing MPA vs reference sites.
#
# The calculation has two steps:
# 1. PROPORTIONS: Convert raw values to proportion of time series maximum.
#    This standardizes across MPAs with different baseline abundances.
#    We add 0.01 (PropCorr) to avoid log(0) when proportions are very small.
#
# 2. LOG RESPONSE RATIO: lnRR = ln(MPA / Reference)
#    - Positive values: MPA site has higher abundance than reference
#    - Negative values: MPA site has lower abundance than reference
#    - Zero: No difference between MPA and reference

cat("  Calculating response ratios...\n")

## --- Biomass response ratio ---
# select() picks specific columns from a dataframe
All.bio <- dplyr::select(All_PISCO.mpa, CA_MPA_Name_Short, year, y, site_status, site_designation, biomass)
All.bio <- subset(All.bio, site_status != "")           # Remove empty status
All.bio <- subset(All.bio, site_designation != "N/A")   # Remove invalid designations
All.bio <- All.bio[!is.na(All.bio$biomass) & !is.nan(All.bio$biomass), ]  # Remove NA/NaN

# Convert to proportions of time series maximum using utility function from 01_utils.R
# This standardizes values so different MPAs are comparable
All.bio <- calculate_proportions(All.bio, "biomass")

# Create wide format (one column for mpa, one for reference) to calculate ratio
# spread() converts from long to wide format: key column becomes column names
All.bio.sub <- All.bio[, c("CA_MPA_Name_Short", "year", "y", "site_designation", "site_status", "PropCorr")]
Short.bio <- All.bio.sub %>%
  spread(site_status, PropCorr)

# Calculate log response ratio using utility function from 01_utils.R
Short.bio.diff <- calculate_log_response_ratio(Short.bio)
All.bio.diff.final <- Short.bio.diff
All.bio.diff.final$resp <- "Bio"  # Label this as biomass response

## --- Density response ratio ---
All.den <- dplyr::select(All_PISCO.mpa, CA_MPA_Name_Short, year, y, site_status, site_designation, density)
All.den <- subset(All.den, site_status != "")

# Convert to proportions of time series maximum
All.den <- calculate_proportions(All.den, "density")

All.den.sub <- All.den[, c("CA_MPA_Name_Short", "year", "y", "site_designation", "site_status", "PropCorr")]
Short.den <- All.den.sub %>%
  spread(site_status, PropCorr)
Short.den.diff <- calculate_log_response_ratio(Short.den)
All.den.diff.final <- Short.den.diff
All.den.diff.final$resp <- "Den"

# ===========================================================================
# Section 16: Create PISCO.resp Dataframe with Raw Density/Biomass
# ===========================================================================

# Density in wide format
All_PISCO.edit <- All_PISCO.mpa[, c("site_status", "CA_MPA_Name_Short", "year", "density", "y")]
All_PISCO.edit <- All_PISCO.edit[complete.cases(All_PISCO.edit), ]
All_PISCO.edit <- subset(All_PISCO.edit, site_status != "")
PISCO.den <- All_PISCO.edit %>%
  spread(site_status, density)
PISCO.den$source <- "PISCO"
PISCO.den <- PISCO.den[complete.cases(PISCO.den), ]
PISCO.den.long <- gather(PISCO.den, status, value, mpa, reference)
PISCO.den.long$resp <- "Den"

# Biomass: adjust units. Macrocystis biomass is already per m2 (computed from
# stipe density). Other taxa have biomass per transect, so divide by transect area.
All_PISCO.mpa$Bio <- ifelse(All_PISCO.mpa$y == "MACPYRAD",
                            All_PISCO.mpa$biomass,
                            All_PISCO.mpa$biomass / PISCO_SWATH_AREA_M2)

All_PISCO.edit.bio <- All_PISCO.mpa[, c("site_status", "CA_MPA_Name_Short", "year", "Bio", "y")]
All_PISCO.edit.bio <- All_PISCO.edit.bio[complete.cases(All_PISCO.edit.bio), ]
All_PISCO.edit.bio <- subset(All_PISCO.edit.bio, site_status != "")
PISCO.bio <- All_PISCO.edit.bio %>%
  spread(site_status, Bio)
PISCO.bio$source <- "PISCO"
PISCO.bio <- PISCO.bio[complete.cases(PISCO.bio), ]
PISCO.bio.long <- gather(PISCO.bio, status, value, mpa, reference)
PISCO.bio.long$resp <- "Bio"

PISCO.resp <- rbind(PISCO.den.long, PISCO.bio.long)

cat("  Biomass RR:", nrow(All.bio.diff.final), "rows;",
    "Density RR:", nrow(All.den.diff.final), "rows\n")
cat("  PISCO.resp:", nrow(PISCO.resp), "rows\n")

# ===========================================================================
# Section 17: Combine Response Ratios and Assign Time/BA Labels
# ===========================================================================

Swath.join <- rbind(All.bio.diff.final, All.den.diff.final)
Swath.join.sub <- Swath.join %>% dplyr::arrange(CA_MPA_Name_Short, year)
Swath.join.sub <- data.frame(Swath.join.sub)

# Assign time since MPA implementation using utility function
Swath.join.sub$time <- 0
Swath.join.sub <- assign_time_from_site_table(Swath.join.sub, Site)

# Merge with MPA feature data (size, type, location)
mpa_features_path <- here::here("data", "MPAfeatures_subset.csv")
Site.size <- tryCatch(
  read.csv(mpa_features_path),
  error = function(e) stop("Failed to read MPA features: ", mpa_features_path,
                           "\n  ", e$message, call. = FALSE)
)
Site.merged <- merge(Site, Site.size, by.x = c("Site"), by.y = c("NAME"), all = TRUE)
Site.merged <- Site.merged[!is.na(Site.merged$Lat), ]
sites.short.edit <- subset(Site.merged,
                           select = c(CA_MPA_Name_Short, type, Location, Hectares))

Swath.join.sub <- merge(Swath.join.sub, sites.short.edit,
                        by.x = c("CA_MPA_Name_Short"),
                        by.y = c("CA_MPA_Name_Short"),
                        all.x = TRUE)

Swath.join.sub <- Swath.join.sub[complete.cases(Swath.join.sub$year), ]
Swath.join.sub <- subset(Swath.join.sub,
                         !(CA_MPA_Name_Short %in% c("Anacapa Island SMR 1978",
                                                     "Anacapa Island SMCA",
                                                     "Point Dume SMR")))

Swath.join.sub$source <- "PISCO"

# NOTE [M4]: BA assignment here uses hardcoded logic (time==0 → Before, else After)
# rather than assign_ba_from_site_table(). This is functionally equivalent because
# `time` was already set by assign_time_from_site_table() at line 871, which uses
# the same year <= MPA_Start convention. Kept as-is for consistency with the
# safety check below (which validates BA against the time column).
Swath.join.sub$BA <- ifelse(Swath.join.sub$time == 0, "Before", "After")

# Safety check: flag MPAs where ALL observations are "Before" (potential lookup failure)
ba_check <- tapply(Swath.join.sub$BA, Swath.join.sub$CA_MPA_Name_Short, function(x) all(x == "Before"))
all_before_mpas <- names(ba_check)[ba_check]
if (length(all_before_mpas) > 0) {
  warning("PISCO BA assignment: ", length(all_before_mpas), " MPA(s) have ONLY 'Before' observations ",
          "(possible site table lookup failure): ", paste(all_before_mpas, collapse = ", "))
}

# ===========================================================================
# Section 18: Cleanup and Final Outputs
# ===========================================================================
# Final outputs:
#   Swath.join.sub  -- response ratio data with time, BA, MPA features
#   PISCO.resp      -- raw density/biomass inside/outside MPAs
#   All_PISCO.mpa   -- averaged PISCO data by MPA

# Clean up intermediate objects to free memory
# Keep only final outputs: Swath.join.sub, PISCO.resp, All_PISCO.mpa
rm(list = intersect(ls(), c(
  # Data import intermediates
  "Swath", "Swath.edit", "Swath.site.PISCO", "Swath.site.PISCO.sub",
  # Lobster processing intermediates
  "PANINT", "PANINT.UCSB", "PANINT.VRG", "PANINT.UCSB.null", "PANINT.VRG.null",
  "PANINT.UCSB.sub.sub", "VRG.PANINT.site", "VRG.PANINT.site.merge",
  # Urchin processing intermediates
  "URCHINS", "Urchin.site", "PISCO.Urchin.site.merge",
  # Size frequency intermediates (keep SizeFreq.Urch.OG for other scripts)
  "SizeFreq.PANINT.VRG",
  # Proportion/RR calculation intermediates
  "All.bio", "All.den", "All.bio.sub", "All.den.sub",
  "Short.bio", "Short.den", "Short.bio.diff", "Short.den.diff",
  # Bootstrap result lists
  "results_list", "u"
)))

# Force garbage collection after heavy processing
gc(verbose = FALSE)

cat("PISCO processing complete.\n")
cat("  Output objects: Swath.join.sub (", nrow(Swath.join.sub), " rows), ",
    "PISCO.resp (", nrow(PISCO.resp), " rows), ",
    "All_PISCO.mpa (", nrow(All_PISCO.mpa), " rows)\n", sep = "")
