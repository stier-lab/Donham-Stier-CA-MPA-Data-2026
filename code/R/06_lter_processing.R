# =============================================================================
# 06_lter_processing.R
# =============================================================================
#
# PURPOSE:
#   Process SBC LTER (Santa Barbara Coastal Long-Term Ecological Research)
#   data for the MPA pBACIPS analysis.
#
# WHAT THIS SCRIPT DOES:
#   1. Processes LTER urchin data (S. purpuratus, M. franciscanus)
#   2. Processes LTER lobster abundance and biomass data
#   3. Processes LTER Macrocystis kelp frond data
#   4. Processes LTER sheephead fish data
#   5. Calculates biomass via bootstrap resampling (urchins) or
#      allometric equations (lobster, fish)
#   6. Calculates log response ratios (MPA vs reference sites)
#   7. Assigns time since MPA and Before/After labels
#
# DATA SOURCES:
#   SBC LTER (https://sbclter.msi.ucsb.edu/) has monitored kelp forests
#   since 2000. Key datasets used:
#   - Annual_Kelp_All_Years: Frond counts (1984-present)
#   - Annual_fish_comb: Fish abundance by size (1984-present)
#   - Lobster_Abundance_All_Years: Lobster size-abundance (2012-present)
#
# KEY FEATURES OF LTER DATA:
#   - Only covers Campus Point SMCA and Naples SMCA (mainland MPAs)
#   - Has complementary datasets for same sites (fish, kelp, inverts)
#   - Standardized protocols since program inception
#   - Lobster data has individual sizes (no bootstrap needed)
#
# INPUTS:
#   - lter object from 05_kfm_processing.R (MBON subset)
#   - data/LTER/Annual_Kelp_All_Years_20240305.csv
#   - data/LTER/Annual_fish_comb_20240307.csv
#   - data/LTER/Lobster_Abundance_All_Years_20230831.csv
#   - SizeFreq.Urch.OG (from 04_pisco_processing.R)
#
# OUTPUTS:
#   - LTER.join.ave: Urchin density & biomass response ratios
#   - LTER.lob: Lobster density & biomass response ratios
#   - Short.lter.macro: Kelp frond density & biomass response ratios
#   - LTER.fish: Sheephead density & biomass response ratios
#   - LTER.resp, LTER.lob.resp, LTER.macro.resp, LTER.fish.resp: Raw data
#
# DEPENDENCIES:
#   Requires 00-05 scripts to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

cat("Processing LTER data...\n")

# Source utility functions
source(here::here("code", "R", "01_utils.R"))

# Respect FORCE_BOOTSTRAP only when it is explicitly TRUE.
# This prevents accidental cache-busting when the variable exists but is FALSE/0/"FALSE".
force_boot <- should_force_bootstrap()

# ===========================================================================
# Section 1: Process LTER urchin/macro data from MBON integrated dataset
# ===========================================================================
# The 'lter' object was created in 05_kfm_processing.R by subsetting MBON data.
# LTER = Santa Barbara Coastal Long-Term Ecological Research program.
# They monitor kelp forest sites along the mainland coast (not Channel Islands).

# Validate that the 'lter' object exists from 05_kfm_processing.R
if (!exists("lter") || !is.data.frame(lter) || nrow(lter) == 0) {
  stop("'lter' object not found or empty. Ensure 05_kfm_processing.R has been sourced first.")
}
if (!exists("Sites2") || !is.data.frame(Sites2)) {
  stop("'Sites2' object not found. Ensure 05_kfm_processing.R has been sourced first.")
}

cat("  LTER input from MBON:", nrow(lter), "rows\n")

# Merge LTER survey data with site metadata
# merge() joins two tables like a database JOIN or Excel VLOOKUP
lter.site <- validated_merge(lter, Sites2, by = c("site_id"), warn_threshold = 0)

# Remove sites with no MPA association and excluded LTER sites
# EXCLUDED_LTER_SITES (defined in 01_utils.R): "a-l-02" = Arroyo Honda,
# excluded due to data quality issues per Bob Miller (LTER PI)
lter.site <- subset(lter.site, CA_MPA_Name_Short != "" & !(site_id %in% EXCLUDED_LTER_SITES))

# Extract year from date column using lubridate's year() function
lter.site$year <- year(lter.site$date)

# Assign time since MPA implementation using utility function from 01_utils.R
# This calculates how many years after MPA start each observation occurred
# (replaces hundreds of if-else lines with a single function call)
lter.site <- assign_time_from_site_table(lter.site, Sites2)

# Subset to only the two LTER MPAs we're analyzing
# These are the mainland MPAs monitored by LTER with sufficient before/after data
lter.site <- subset(lter.site, CA_MPA_Name_Short == "Campus Point SMCA" |
                      CA_MPA_Name_Short == "Naples SMCA")

stopifnot("LTER site data must have > 0 rows after filtering" = nrow(lter.site) > 0)
cat("  LTER after site/MPA filtering:", nrow(lter.site), "rows\n")

# Subset target taxa: urchins from quads, macrocystis from swath
# Different taxa are surveyed using different methods:
# - Urchins: counted in quadrats (small fixed areas, typically 1 m^2)
# - Kelp: counted along swath transects (longer strips, typically 40 m^2)
# The | operator means "OR" - keep any row that matches any condition
lter.sub <- subset(lter.site,
                   (taxon_name == "Strongylocentrotus purpuratus" & sample_method == "quad") |
                     (taxon_name == "Mesocentrotus franciscanus" & sample_method == "quad") |
                     (taxon_name == "Macrocystis pyrifera" & sample_method == "swath"))

# Calculate density: count divided by survey area gives individuals per square meter
lter.sub$den <- safe_divide(lter.sub$count, lter.sub$area, context = "LTER den = count/area")

# Calculate site means using two-step aggregation (nested sampling design)
# Step 1: Average across replicates within each transect
# group_by() defines the grouping structure; summarise_at() calculates mean for 'den'
lter.ave <- lter.sub %>%
  group_by(site_id, transect_id, taxon_name, status, CA_MPA_Name_Short,
           ChannelIsland, MPA_Start, time, year) %>%
  summarise_at(c("den"), mean)

# Step 2: Average across transects to get site-level annual means
lter.ave <- lter.ave %>%
  group_by(site_id, taxon_name, status, CA_MPA_Name_Short,
           ChannelIsland, MPA_Start, time, year) %>%
  summarise_at(c("den"), mean)

cat("  LTER target taxa density averages:", nrow(lter.ave), "rows\n")

# ===========================================================================
# Section 2: Bootstrap urchin biomass from size frequency distributions
# ===========================================================================
# LTER counts urchins but doesn't measure individual sizes. To estimate biomass,
# we bootstrap-resample from the PISCO size-frequency distribution (SizeFreq.Urch.OG)
# collected from the same MPA/status/year, then apply species-specific allometric
# equations (bio_redurch for M. franciscanus, bio_purpurch for S. purpuratus).
# This is the same approach used in 05_kfm_processing.R for KFM urchins.
#
# CROSS-PROGRAM BOOTSTRAP ASSUMPTION:
#   LTER does not collect urchin size data — only total counts per quadrat.
#   To estimate biomass, we borrow size-frequency distributions from PISCO
#   (SizeFreq.Urch.OG), which measures individual urchin test diameters.
#
#   The bootstrap matches by: CA_MPA_Name_Short, site_status (mpa/reference),
#   year, and species (classcode). For each LTER site-year, we draw `n`
#   individuals (with replacement) from the matching PISCO size-frequency
#   pool and convert to biomass via allometric equations.
#
#   Key assumption: Urchin size structure at LTER mainland sites is similar
#   to PISCO sites within the same MPA and year. This is a necessary
#   approximation because LTER's monitoring protocol does not include size
#   measurements for urchins. When no matching PISCO size data exists for a
#   given MPA/year/status combination, biomass is set to NA and that
#   observation is excluded.
#
#   Note: SizeFreq.Urch.OG is the UNFILTERED size-frequency data (before
#   the PISCO-specific 25mm minimum filter applied in 04_pisco_processing.R).
#   See 05_kfm_processing.R Section 6 for the identical approach and
#   further discussion of this assumption.

URCHINS <- subset(lter.site,
                  taxon_name == "Strongylocentrotus purpuratus" |
                    taxon_name == "Mesocentrotus franciscanus")
URCHINS <- mutate(URCHINS, taxon_name = case_when(
  taxon_name == "Strongylocentrotus purpuratus" ~ "STRPURAD",
  taxon_name == "Mesocentrotus franciscanus"    ~ "MESFRAAD",
  TRUE ~ taxon_name))
URCHINS$replicate_id <- as.numeric(URCHINS$replicate_id)

# Get unique site-year-transect-taxon combinations
u_cols_required <- c("year", "site_id", "CA_MPA_Name_Short", "status",
                     "area", "taxon_name", "transect_id")
missing_u_cols <- setdiff(u_cols_required, names(URCHINS))
if (length(missing_u_cols) > 0) {
  stop("URCHINS is missing required columns for LTER urchin bootstrap: ",
       paste(missing_u_cols, collapse = ", "))
}
u_cols <- u_cols_required
if ("site_designation" %in% names(URCHINS)) u_cols <- c(u_cols, "site_designation")
if ("BaselineRegion" %in% names(URCHINS)) u_cols <- c(u_cols, "BaselineRegion")
u <- unique(URCHINS[, u_cols])

cat("  LTER urchin bootstrap:", nrow(u), "site-year-transect-species combinations\n")

# Prepare size frequency data
SizeFreq.Urch <- SizeFreq.Urch.OG
SizeFreq.Urch <- SizeFreq.Urch[complete.cases(SizeFreq.Urch$size), ]

# Cache bootstrap results to avoid re-running the slow loop
.cache_lter_urch <- here::here("data", "cache", "lter_urchin_bootstrap.rds")
cached <- NULL
if (file.exists(.cache_lter_urch) && !force_boot) {
  cat("    Loading cached LTER urchin bootstrap results...\n")
  cached <- safe_read_rds(.cache_lter_urch, default = NULL)
  if (is.null(cached)) {
    cat("    Cache invalid, recomputing...\n")
  }
}

if (!is.null(cached)) {
  # Upgrade old cache schemas: earlier caches may not include these metadata columns.
  if (!"site_designation" %in% names(cached)) cached$site_designation <- NA_character_
  if (!"BaselineRegion" %in% names(cached)) cached$BaselineRegion <- NA_character_

  required_core <- c(
    "year", "site", "site_status", "CA_MPA_Name_Short",
    "transect", "quad", "area", "y", "count", "biomass"
  )
  missing_core <- setdiff(required_core, names(cached))
  if (length(missing_core) > 0) {
    cat("    Cached LTER urchin bootstrap missing core columns (",
        paste(missing_core, collapse = ", "),
        "); recomputing...\n", sep = "")
    cached <- NULL
  }
}

if (!is.null(cached)) {
  Urchin.site <- cached
} else {
  # FIX [M5]: Replaced O(n^2) rbind-in-loop with pre-allocated list accumulation.
  # Previously, each iteration called rbind() which copies the entire growing
  # data.frame on every iteration. With ~1000+ iterations, this is very slow.
  # Now we store results in a pre-allocated list and bind once at the end.
  results_list <- vector("list", nrow(u))

  set.seed(SEED_LTER_URCHIN)  # Reproducibility for bootstrap resampling (from 00c_analysis_constants.R)
  for (i in seq_len(nrow(u))) {
    t <- which(URCHINS$site_id == u$site_id[i] &
                 URCHINS$year == u$year[i] &
                 URCHINS$transect_id == u$transect_id[i] &
                 URCHINS$taxon_name == u$taxon_name[i])
    n <- sum(URCHINS$count[t], na.rm = TRUE)
    t2 <- which(SizeFreq.Urch$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] &
                  SizeFreq.Urch$site_status == u$status[i] &
                  SizeFreq.Urch$year == u$year[i] &
                  SizeFreq.Urch$classcode == u$taxon_name[i])

    # Select appropriate biomass function based on species:
    #   bio_redurch: M. franciscanus (red urchin) allometric equation
    #   bio_purpurch: S. purpuratus (purple urchin) allometric equation
    bio_fun <- if (u$taxon_name[i] == "MESFRAAD") bio_redurch else bio_purpurch

    # Use bootstrap_biomass utility function
    result <- bootstrap_biomass(count = n,
                                size_freq_indices = t2,
                                size_freq_table = SizeFreq.Urch,
                                biomass_fun = bio_fun,
                                n_resamples = N_BOOTSTRAP_RESAMPLES)

    results_list[[i]] <- data.frame(
      site             = u$site_id[i],
      CA_MPA_Name_Short = u$CA_MPA_Name_Short[i],
      site_status      = u$status[i],
      site_designation = if ("site_designation" %in% names(u)) as.character(u$site_designation[i]) else NA_character_,
      BaselineRegion   = if ("BaselineRegion" %in% names(u)) as.character(u$BaselineRegion[i]) else NA_character_,
      year             = u$year[i],
      transect         = u$transect_id[i],
      quad             = length(t),
      area             = u$area[i],
      y                = u$taxon_name[i],
      biomass          = result$biomass,
      count            = result$count,
      stringsAsFactors = FALSE
    )
  }

  # Combine all results at once — O(n) instead of O(n^2)
  Urchin.site <- dplyr::bind_rows(results_list)

  # Save cache
  dir.create(dirname(.cache_lter_urch), recursive = TRUE, showWarnings = FALSE)
  saveRDS(Urchin.site, .cache_lter_urch)
  cat("    Cached LTER urchin bootstrap results.\n")
}

stopifnot("Urchin.site must have > 0 rows" = nrow(Urchin.site) > 0)
cat("  LTER urchin bootstrap complete:", nrow(Urchin.site), "rows\n")

# Bootstrapped biomass output already carries site metadata (from Sites2 via URCHINS).
# Avoid merging with sites.short (PISCO/KFM site table), which can silently mismatch and
# makes downstream column-index slicing brittle.
LTER.Urchin.site.merge <- Urchin.site
need_urch_cols <- c("year", "site", "site_status", "CA_MPA_Name_Short",
                    "site_designation", "BaselineRegion",
                    "transect", "quad", "area", "y", "count", "biomass")
missing_urch_cols <- setdiff(need_urch_cols, names(LTER.Urchin.site.merge))
if (length(missing_urch_cols) > 0) {
  stop("LTER.Urchin.site.merge missing required columns: ",
       paste(missing_urch_cols, collapse = ", "))
}
LTER.Urchin.site.merge <- LTER.Urchin.site.merge[, need_urch_cols]
# Zero-fill only numeric columns; preserve NAs in character/factor columns
# Log how many NAs are being replaced so silent data loss is detectable
num_cols <- sapply(LTER.Urchin.site.merge, is.numeric)
n_na_before_zerofill <- sum(is.na(LTER.Urchin.site.merge[num_cols]))
LTER.Urchin.site.merge[num_cols][is.na(LTER.Urchin.site.merge[num_cols])] <- 0
if (n_na_before_zerofill > 0) {
  cat(sprintf("  LTER urchin zero-fill: replaced %d NAs with 0 across %d numeric columns\n",
              n_na_before_zerofill, sum(num_cols)))
}

# Remove site-years where urchins were counted but not sized (can't estimate biomass)
LTER.Urchin.null <- subset(LTER.Urchin.site.merge, count > 0 & biomass == 0)
LTER.Urchin.site.merge.sub <- LTER.Urchin.site.merge
for (i in seq_len(nrow(LTER.Urchin.null))) {
  LTER.Urchin.site.merge.sub <- subset(LTER.Urchin.site.merge.sub,
                                       site != LTER.Urchin.null$site[i] |
                                         year != LTER.Urchin.null$year[i] |
                                         y != LTER.Urchin.null$y[i])
}

if (nrow(LTER.Urchin.null) > 0) {
  cat("  Removed", nrow(LTER.Urchin.null), "urchin rows with count > 0 but no biomass\n")
}

# Normalize by number of quads per transect
LTER.Urchin.site.transectSum <- LTER.Urchin.site.merge.sub
LTER.Urchin.site.transectSum$count <- safe_divide(
  LTER.Urchin.site.transectSum$count,
  LTER.Urchin.site.transectSum$quad,
  context = "LTER urchin count / quad"
)
LTER.Urchin.site.transectSum$biomass <- safe_divide(
  LTER.Urchin.site.transectSum$biomass,
  LTER.Urchin.site.transectSum$quad,
  context = "LTER urchin biomass / quad"
)

# Average across transects
LTER.Urchin.site.all <- LTER.Urchin.site.transectSum %>%
  group_by(year, site, CA_MPA_Name_Short, site_designation, site_status,
           BaselineRegion, y, area) %>%
  summarise_at(c("count", "biomass"), mean)
LTER.Urchin.site.all$Density <- safe_divide(
  LTER.Urchin.site.all$count,
  LTER.Urchin.site.all$area,
  context = "LTER urchin Density = count/area"
)

# Recode taxa in lter.ave to match abbreviated names used in bootstrap output
lter.ave <- mutate(lter.ave, taxon_name = case_when(
  taxon_name == "Strongylocentrotus purpuratus" ~ "STRPURAD",
  taxon_name == "Mesocentrotus franciscanus"    ~ "MESFRAAD",
  TRUE ~ taxon_name))

# Merge density and biomass data
LTER.allbio <- merge(lter.ave, LTER.Urchin.site.all,
                     by.x = c("site_id", "status", "year", "taxon_name", "CA_MPA_Name_Short"),
                     by.y = c("site", "site_status", "year", "y", "CA_MPA_Name_Short"),
                     all = TRUE)

Swath.LTER <- LTER.allbio
keep_swath_lter <- c("site_id", "status", "CA_MPA_Name_Short", "ChannelIsland", "MPA_Start",
                     "Density", "time", "year", "taxon_name", "den", "biomass")
missing_swath_lter <- setdiff(keep_swath_lter, names(Swath.LTER))
if (length(missing_swath_lter) > 0) {
  stop("Swath.LTER missing required columns after merge: ",
       paste(missing_swath_lter, collapse = ", "))
}
Swath.LTER <- Swath.LTER[, keep_swath_lter]
names(Swath.LTER)[names(Swath.LTER) == "Density"]   <- "count.y"
names(Swath.LTER)[names(Swath.LTER) == "status"]    <- "site_status"
names(Swath.LTER)[names(Swath.LTER) == "den"]       <- "Density"
names(Swath.LTER)[names(Swath.LTER) == "site_id"]   <- "site"
names(Swath.LTER)[names(Swath.LTER) == "taxon_name"] <- "y"

# ===========================================================================
# Section 3: Calculate proportions and log response ratios for urchins
# ===========================================================================

# Average across sites
lter.ave.ave <- Swath.LTER %>%
  group_by(y, site_status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, time, year) %>%
  summarise_at(vars("Density", "biomass"), mean)

# Restore full species names (reverse abbreviation from Section 2)
lter.ave.ave <- mutate(lter.ave.ave, y = case_when(
  y == "STRPURAD" ~ "Strongylocentrotus purpuratus",
  y == "MESFRAAD" ~ "Mesocentrotus franciscanus",
  TRUE ~ y))
names(lter.ave.ave)[names(lter.ave.ave) == "Density"] <- "den"

## --- LTER density proportions and log response ratios ---
# Exclude Macrocystis here — it is processed independently in Section 5 via standalone
# LTER kelp frond data. Including it here would create duplicate density response ratios
# with different source labels ("LTER" vs "LTER macro surveys"), bypassing deduplication.
All.den <- lter.ave.ave
All.den <- subset(All.den, y != "Macrocystis pyrifera")
All.den <- All.den[complete.cases(All.den$den), ]
All.den <- calculate_proportions(All.den, "den")

All.den.sub <- All.den[, c("CA_MPA_Name_Short", "year", "y", "site_status", "PropCorr")]
Short.den <- All.den.sub %>% spread(site_status, PropCorr)
Short.den.diff <- calculate_log_response_ratio(Short.den)
Short.den.diff$resp <- "Den"

## --- LTER biomass proportions and log response ratios ---
# Exclude Macrocystis (processed in Section 5); NA biomass from the urchin merge
# would be removed by complete.cases anyway, but filtering explicitly for clarity.
All.bio <- lter.ave.ave
All.bio <- subset(All.bio, y != "Macrocystis pyrifera")
All.bio <- All.bio[complete.cases(All.bio$biomass), ]
All.bio <- calculate_proportions(All.bio, "biomass")

All.bio.sub <- All.bio[, c("CA_MPA_Name_Short", "year", "y", "site_status", "PropCorr")]
Short.bio <- All.bio.sub %>% spread(site_status, PropCorr)
Short.bio.diff <- calculate_log_response_ratio(Short.bio)
Short.bio.diff$resp <- "Bio"

# Combine density and biomass response ratios
LTER.join.ave <- rbind(Short.den.diff, Short.bio.diff)

cat("  LTER urchin response ratios:", nrow(LTER.join.ave), "rows (density + biomass)\n")

## --- Create LTER urchin response dataframes (density/biomass, long format) ---
# These dataframes store the raw density/biomass values for MPA vs reference sites.
# They're used for plotting raw data time series and as input to other analyses.

# Select relevant columns for density response dataframe
lter.edit.den <- lter.ave.ave[, c("site_status", "CA_MPA_Name_Short", "MPA_Start", "year", "den", "y")]
lter.edit.den <- lter.edit.den[complete.cases(lter.edit.den), ]  # Remove rows with NAs

# spread() converts from long to wide format: site_status values become columns
LTER.den <- lter.edit.den %>% spread(site_status, den)
LTER.den <- subset(LTER.den, y != "Macrocystis pyrifera")  # Kelp processed separately in Section 5
LTER.den$source <- "LTER"  # Add source label for tracking
LTER.den <- LTER.den[complete.cases(LTER.den), ]

# gather() converts from wide back to long format for storage
# Columns mpa, reference are gathered into key-value pairs
# key = "site_status" (mpa or reference), value = "value" (the density)
LTER.den.long <- gather(LTER.den, site_status, value, mpa, reference)
LTER.den.long$resp <- "Den"  # Label as density response

# Same process for biomass
lter.edit.bio <- lter.ave.ave[, c("site_status", "CA_MPA_Name_Short", "MPA_Start", "year", "biomass", "y")]
lter.edit.bio <- lter.edit.bio[complete.cases(lter.edit.bio), ]
LTER.bio <- lter.edit.bio %>% spread(site_status, biomass)
LTER.bio$source <- "LTER"
LTER.bio <- LTER.bio[complete.cases(LTER.bio), ]
LTER.bio.long <- gather(LTER.bio, site_status, value, mpa, reference)
LTER.bio.long$resp <- "Bio"

# Combine density and biomass into one response dataframe
# rbind() stacks dataframes vertically (row bind)
LTER.resp <- rbind(LTER.bio.long, LTER.den.long)

## --- Assign time and BA labels to LTER.join.ave ---
LTER.join.ave <- assign_time_from_site_table(LTER.join.ave, Sites2)
LTER.join.ave <- merge(LTER.join.ave, sites.short.edit,
                       by.x = c("CA_MPA_Name_Short"),
                       by.y = c("CA_MPA_Name_Short"),
                       all.x = TRUE)
LTER.join.ave <- LTER.join.ave[complete.cases(LTER.join.ave$year), ]
LTER.join.ave$source <- "LTER"
LTER.join.ave <- assign_ba_from_site_table(LTER.join.ave, Sites2)

# ===========================================================================
# Section 4: Import and process LTER lobster data
# ===========================================================================
# LTER has a separate lobster monitoring dataset with individual size measurements.
# Unlike PISCO/KFM, we don't need bootstrap resampling because LTER measures
# each individual lobster's carapace length directly.
# Lobster (Panulirus interruptus) is a key urchin predator in the trophic cascade.
#
# BEFORE-PERIOD DATA LIMITATION:
#   SBC LTER began lobster abundance surveys in 2012 (the data file covers
#   2012-2023). The LTER MPAs (Campus Point SMCA and Naples SMCA) were
#   also established in 2012 (MLPA South Coast implementation).
#
#   Under the pipeline convention (year <= MPA_Start is "Before"), only
#   observations from 2012 itself are classified as "Before" period. This
#   means the before-period sample for LTER lobster is at most 1 year,
#   which is very limited for establishing a pre-MPA baseline.
#
#   In contrast, other LTER taxa (urchins, kelp, sheephead) have data
#   going back much further (to the 1980s-2000s via MBON integration),
#   providing robust before-period baselines.
#
#   The PISCO lobster data (04_pisco_processing.R) has a similar limitation:
#   PISCO size measurements started in 2010-2011 (UCSB/VRG respectively),
#   giving only 1-2 years of before-period lobster biomass data.
#
#   This edge case does not cause errors in the pipeline — lnRR is still
#   calculated for years where both MPA and reference site data exist. But
#   the limited before-period means lobster effect sizes rely primarily on
#   the after-period trajectory rather than before-after contrasts.

lter.lob <- read.csv(here::here("data", "LTER", "Lobster_Abundance_All_Years_20230831.csv"))

# Validate expected columns
required_lob_cols <- c("SITE", "YEAR", "TRANSECT", "REPLICATE", "COUNT", "SIZE_MM")
missing_lob_cols <- setdiff(required_lob_cols, names(lter.lob))
if (length(missing_lob_cols) > 0) {
  stop("LTER lobster CSV is missing required columns: ",
       paste(missing_lob_cols, collapse = ", "))
}

cat("  LTER lobster data loaded:", nrow(lter.lob), "rows\n")

# Merge with site metadata using SITE column as the key
lter.lob.site <- validated_merge(lter.lob, Sites2, by.x = c("SITE"), by.y = c("site_id"), warn_threshold = 0)
lter.lob.site <- subset(lter.lob.site, CA_MPA_Name_Short != "")

# Calculate biomass from size using allometric equation
# bio_lobster() from 01_utils.R converts carapace length (mm) to weight (g)
# Formula: biomass = 0.001352821 * carapace_length^2.913963113
# Total biomass = count * weight per individual
lter.lob.site$biomass <- lter.lob.site$COUNT * bio_lobster(lter.lob.site$SIZE_MM)
n_na_lob_bio <- sum(is.na(lter.lob.site$biomass))
lter.lob.site$biomass[is.na(lter.lob.site$biomass)] <- 0
if (n_na_lob_bio > 0) {
  cat(sprintf("  LTER lobster: zero-filled %d NA biomass values (NA sizes)\n", n_na_lob_bio))
}

# Summarise: replicate -> transect -> site -> MPA level
# This hierarchical averaging accounts for the nested sampling design
lter.lob.site.sum <- lter.lob.site %>%
  group_by(YEAR, SITE, TRANSECT, REPLICATE, status, CA_MPA_Name_Short,
           ChannelIsland, MPA_Start) %>%
  summarise_at(c("COUNT", "biomass"), sum)

lter.lob.site.ave <- lter.lob.site.sum %>%
  group_by(YEAR, SITE, TRANSECT, status, CA_MPA_Name_Short,
           ChannelIsland, MPA_Start) %>%
  summarise_at(c("COUNT", "biomass"), mean)

lter.lob.site.ave <- lter.lob.site.ave %>%
  group_by(YEAR, SITE, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>%
  summarise_at(c("COUNT", "biomass"), mean)

lter.lob.site.max <- lter.lob.site.ave %>%
  group_by(status, CA_MPA_Name_Short, MPA_Start, YEAR) %>%
  summarise_at(c("COUNT", "biomass"), mean)

## --- Lobster biomass proportions and log response ratios ---
# NOTE: Using adaptive correction (min/2) for consistency across all taxa/programs.
# This is statistically preferred over fixed +0.01 (Aitchison 1986) as it scales
# appropriately with the data and avoids inflating effect sizes for rare species.
lter.lob.site.max$taxon_name <- "Panulirus interruptus"
lter.lob.site.max$y <- lter.lob.site.max$taxon_name  # Required for calculate_proportions()

# Rename YEAR to year for consistency with other data sources
names(lter.lob.site.max)[names(lter.lob.site.max) == "YEAR"] <- "year"

# Biomass proportions - use standardized function with adaptive correction
lter.lob.site.max <- calculate_proportions(lter.lob.site.max, "biomass")

# Use column names instead of fragile indices
All.bio.panLTERall.sub <- lter.lob.site.max[, c("CA_MPA_Name_Short", "year", "y", "status", "PropCorr")]
Short.bio.panLTERall <- All.bio.panLTERall.sub %>% spread(status, PropCorr)
Short.bio.panLTERall.diff <- calculate_log_response_ratio(Short.bio.panLTERall)
Short.bio.panLTERall.diff$resp <- "Bio"

# Count/density proportions - use standardized function with adaptive correction
lter.lob.site.max <- calculate_proportions(lter.lob.site.max, "COUNT")

# Use column names instead of fragile indices
All.den.panLTERall.sub <- lter.lob.site.max[, c("CA_MPA_Name_Short", "year", "y", "status", "PropCorr")]
Short.den.panLTERall <- All.den.panLTERall.sub %>% spread(status, PropCorr)
Short.den.panLTERall.diff <- calculate_log_response_ratio(Short.den.panLTERall)
Short.den.panLTERall.diff$resp <- "Den"

# Combine lobster results
LTER.lob <- rbind(Short.bio.panLTERall.diff, Short.den.panLTERall.diff)
# Defensive check: ensure expected column names are present
stopifnot("year" %in% names(LTER.lob), "y" %in% names(LTER.lob),
          "CA_MPA_Name_Short" %in% names(LTER.lob))

cat("  LTER lobster response ratios:", nrow(LTER.lob), "rows\n")

# Assign time and BA labels using utility functions
LTER.lob <- assign_time_from_site_table(LTER.lob, Sites2)
LTER.lob <- merge(LTER.lob, sites.short.edit,
                  by.x = c("CA_MPA_Name_Short"),
                  by.y = c("CA_MPA_Name_Short"),
                  all.x = TRUE)
LTER.lob <- LTER.lob[complete.cases(LTER.lob$year), ]
LTER.lob$source <- "LTER lob surveys"
LTER.lob$BA <- "N/A"
LTER.lob <- assign_ba_from_site_table(LTER.lob, Sites2)

# NOTE: LTER lobster monitoring started in 2012 and most LTER MPAs were also
# implemented in 2012, so the Before-period typically has only 1 year (2012 itself).
# This is a known limitation that weakens the before-after contrast for LTER lobster.
n_before_lob <- sum(LTER.lob$BA == "Before", na.rm = TRUE)
n_after_lob <- sum(LTER.lob$BA == "After", na.rm = TRUE)
cat(sprintf("  LTER lobster BA: %d Before, %d After (limited before-period expected)\n",
            n_before_lob, n_after_lob))

## --- Create LTER lobster response dataframes ---
# Convert raw counts to per-m^2 density using LTER lobster plot area
# LTER_LOBSTER_PLOT_AREA_M2 = 300 m^2 (defined in 00c_analysis_constants.R)
lter.lob.site.max$Den <- lter.lob.site.max$COUNT / LTER_LOBSTER_PLOT_AREA_M2
lter.lob.site.max$Bio <- lter.lob.site.max$biomass / LTER_LOBSTER_PLOT_AREA_M2

# Use column names instead of fragile indices
lter.lob.edit.den <- lter.lob.site.max[, c("status", "CA_MPA_Name_Short", "MPA_Start", "year", "Den", "y")]
LTER.lob.den <- lter.lob.edit.den %>% spread(status, Den)
LTER.lob.den$source <- "LTER"
LTER.lob.den.sub <- LTER.lob.den[complete.cases(LTER.lob.den), ]
LTER.lob.den.long <- gather(LTER.lob.den.sub, status, value, mpa, reference)
LTER.lob.den.long$resp <- "Den"

lter.lob.edit.bio <- lter.lob.site.max[, c("status", "CA_MPA_Name_Short", "MPA_Start", "year", "Bio", "y")]
LTER.lob.bio <- lter.lob.edit.bio %>% spread(status, Bio)
LTER.lob.bio$source <- "LTER"
LTER.lob.bio.sub <- LTER.lob.bio[complete.cases(LTER.lob.bio), ]
LTER.lob.bio.long <- gather(LTER.lob.bio.sub, status, value, mpa, reference)
LTER.lob.bio.long$resp <- "Bio"

LTER.lob.resp <- rbind(LTER.lob.den.long, LTER.lob.bio.long)

# ===========================================================================
# Section 5: Import and process LTER kelp frond data
# ===========================================================================
# LTER monitors giant kelp (Macrocystis pyrifera) via frond counts.
# Important methodological note: LTER counts FRONDS (blade-bearing branches),
# while KFM counts STIPES (main stem structures). A single stipe can
# support multiple fronds, so these measurements are not directly comparable.
# See TODO [D5] below regarding the bio_macro() conversion.

lter.macro <- read.csv(here::here("data", "LTER", "Annual_Kelp_All_Years_20240305.csv"))

# Validate expected columns
required_macro_cols <- c("SITE", "YEAR", "TRANSECT", "QUAD", "SIDE",
                         "AREA", "FRONDS")
missing_macro_cols <- setdiff(required_macro_cols, names(lter.macro))
if (length(missing_macro_cols) > 0) {
  stop("LTER kelp CSV is missing required columns: ",
       paste(missing_macro_cols, collapse = ", "))
}

cat("  LTER kelp frond data loaded:", nrow(lter.macro), "rows\n")

lter.macro.site <- validated_merge(lter.macro, Sites2, by.x = c("SITE"), by.y = c("site_id"), warn_threshold = 0)

# Filter: remove empty MPA names, excluded LTER sites, invalid frond counts, non-target MPAs
# FRONDS = -99999 is LTER's placeholder for missing/invalid data
lter.macro.site <- subset(lter.macro.site, CA_MPA_Name_Short != "" & !(SITE %in% EXCLUDED_LTER_SITE_CODES))
lter.macro.site <- subset(lter.macro.site, FRONDS != -99999)  # -99999 = LTER missing data sentinel
lter.macro.site <- subset(lter.macro.site, CA_MPA_Name_Short == "Campus Point SMCA" |
                            CA_MPA_Name_Short == "Naples SMCA")

# Calculate frond density and biomass using bio_macro() utility function
lter.macro.site$frondDen <- safe_divide(lter.macro.site$FRONDS, lter.macro.site$AREA,
                                        context = "LTER macro frondDen = FRONDS/AREA")

# ===========================================================================
# WARNING [D5]: STIPE vs FROND MISMATCH — DOMAIN EXPERT REVIEW REQUIRED
# ===========================================================================
#
# bio_macro() was calibrated on STIPE counts (see 01_utils.R: "Giant kelp
# biomass from stipe count"), using an average slope from site-specific
# stipe-to-biomass regressions. However, LTER data records FROND counts
# (the FRONDS column in Annual_Kelp_All_Years), not stipe counts.
#
# Macrocystis pyrifera anatomy:
#   - A "stipe" is the main stem-like structure emerging from the holdfast.
#   - A "frond" is a blade-bearing lateral branch growing from the stipe.
#   - One stipe supports multiple fronds (typically 3-10+ depending on
#     plant age and environmental conditions; Reed et al. 2008 Ecology).
#
# Consequence: Passing frond density into bio_macro() treats each frond
# as if it were a stipe, which likely OVERESTIMATES biomass because frond
# counts per unit area are higher than stipe counts per unit area. The
# magnitude of overestimation depends on the fronds-per-stipe ratio.
#
# Comparison across programs:
#   - PISCO (04_pisco_processing.R): Uses count * size, where "size" is
#     the number of stipes per individual plant. So count * size = total
#     stipes. bio_macro() then converts stipe density to biomass. CORRECT.
#   - KFM (05_kfm_processing.R): Uses raw stipe counts from the KFM
#     Macrocystis stipe count dataset (KFM_Macrocystis_RawData). The
#     bootstrap resamples per-individual stipe counts. CORRECT.
#   - LTER (this script): Uses FRONDS column from Annual_Kelp_All_Years.
#     Fronds != stipes. POTENTIAL MISMATCH.
#
# NOTE: Because lnRR is computed on PROPORTIONS (proportion of time-series
# max), not raw biomass, and proportions are calculated within each
# MPA x taxon x source, the mismatch cancels out in the response ratio
# as long as the frond-to-stipe ratio is approximately constant over time
# and space within each MPA. The absolute biomass values in LTER.macro.resp
# are affected, but the lnRR (our primary metric) is robust to this
# multiplicative bias.
#
# DO NOT CHANGE this calculation without domain expert confirmation from
# the LTER kelp ecology group (e.g., Dan Reed, Tom Bell) regarding whether
# frond counts should be converted to stipe equivalents before applying
# the allometric relationship.
# ===========================================================================
lter.macro.site$biomass  <- bio_macro(lter.macro.site$frondDen)
cat("  NOTE [D5]: LTER kelp biomass uses bio_macro() on FROND density (not stipe density).\n",
    "  bio_macro() was calibrated on stipe counts — see WARNING [D5] above.\n",
    "  This is acceptable for lnRR because proportions cancel the bias.\n")

# Summarise: quad -> transect -> site -> MPA level
# This hierarchical averaging accounts for the nested sampling design
lter.macro.site.sum <- lter.macro.site %>%
  group_by(YEAR, SITE, TRANSECT, QUAD, SIDE, AREA, status,
           CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>%
  summarise_at(c("FRONDS", "biomass"), sum)
lter.macro.site.sum$den <- safe_divide(lter.macro.site.sum$FRONDS, lter.macro.site.sum$AREA,
                                       context = "LTER macro den = FRONDS/AREA")

lter.macro.site.ave <- lter.macro.site.sum %>%
  group_by(YEAR, TRANSECT, SITE, status, CA_MPA_Name_Short,
           ChannelIsland, MPA_Start) %>%
  summarise_at(c("den", "biomass"), mean)

lter.macro.site.ave <- lter.macro.site.ave %>%
  group_by(YEAR, SITE, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>%
  summarise_at(c("den", "biomass"), mean)

lter.macro.site.max <- lter.macro.site.ave %>%
  group_by(status, CA_MPA_Name_Short, MPA_Start, YEAR) %>%
  summarise_at(c("den", "biomass"), mean)

## --- Macrocystis frond density proportions and log response ratios ---
# NOTE: Using adaptive correction (min/2) for consistency across all taxa/programs.
# This is statistically preferred over fixed +0.01 (Aitchison 1986) as it scales
# appropriately with the data and avoids inflating effect sizes for rare species.
lter.macro.site.max$taxon_name <- "Macrocystis pyrifera"
lter.macro.site.max$y <- lter.macro.site.max$taxon_name  # Required for calculate_proportions()

# Rename YEAR to year for consistency with other data sources
names(lter.macro.site.max)[names(lter.macro.site.max) == "YEAR"] <- "year"

# Frond density proportions - use standardized function with adaptive correction
lter.macro.site.max <- calculate_proportions(lter.macro.site.max, "den")

# Use column names instead of fragile indices
All.lter.macro.frond.sub <- lter.macro.site.max[, c("CA_MPA_Name_Short", "year", "y", "status", "PropCorr")]
Short.lter.macro.frond <- All.lter.macro.frond.sub %>% spread(status, PropCorr)
Short.lter.macro.frond.diff <- calculate_log_response_ratio(Short.lter.macro.frond)
Short.lter.macro.frond.diff$resp <- "Den"

# Biomass proportions - use standardized function with adaptive correction
lter.macro.site.max <- calculate_proportions(lter.macro.site.max, "biomass")

# Use column names instead of fragile indices
All.lter.macro.bio.sub <- lter.macro.site.max[, c("CA_MPA_Name_Short", "year", "y", "status", "PropCorr")]
Short.lter.macro.bio <- All.lter.macro.bio.sub %>% spread(status, PropCorr)
Short.lter.macro.bio.diff <- calculate_log_response_ratio(Short.lter.macro.bio)
Short.lter.macro.bio.diff$resp <- "Bio"

# Combine macrocystis results
Short.lter.macro <- rbind(Short.lter.macro.frond.diff, Short.lter.macro.bio.diff)

cat("  LTER Macrocystis response ratios:", nrow(Short.lter.macro), "rows\n")

# Assign time and BA labels using utility functions
Short.lter.macro <- assign_time_from_site_table(Short.lter.macro, Sites2)
Short.lter.macro <- merge(Short.lter.macro, sites.short.edit,
                          by.x = c("CA_MPA_Name_Short"),
                          by.y = c("CA_MPA_Name_Short"),
                          all.x = TRUE)
Short.lter.macro <- Short.lter.macro[complete.cases(Short.lter.macro$year), ]
Short.lter.macro$source <- "LTER macro surveys"
Short.lter.macro <- assign_ba_from_site_table(Short.lter.macro, Sites2)

## --- Create LTER macrocystis response dataframes ---
# Use column names instead of fragile indices
lter.macro.edit.den <- lter.macro.site.max[, c("status", "CA_MPA_Name_Short", "MPA_Start", "year", "den", "y")]
LTER.macro.den <- lter.macro.edit.den %>% spread(status, den)
LTER.macro.den$source <- "LTER"
LTER.macro.den.sub <- LTER.macro.den[complete.cases(LTER.macro.den), ]
LTER.macro.den.long <- gather(LTER.macro.den.sub, status, value, mpa, reference)
LTER.macro.den.long$resp <- "Den"

lter.macro.edit.bio <- lter.macro.site.max[, c("status", "CA_MPA_Name_Short", "MPA_Start", "year", "biomass", "y")]
LTER.macro.bio <- lter.macro.edit.bio %>% spread(status, biomass)
LTER.macro.bio$source <- "LTER"
LTER.macro.bio.sub <- LTER.macro.bio[complete.cases(LTER.macro.bio), ]
LTER.macro.bio.long <- gather(LTER.macro.bio.sub, status, value, mpa, reference)
LTER.macro.bio.long$resp <- "Bio"

LTER.macro.resp <- rbind(LTER.macro.den.long, LTER.macro.bio.long)

# ===========================================================================
# Section 6: Import and process LTER fish data (sheephead only)
# ===========================================================================
# LTER conducts visual fish surveys along with their benthic monitoring.
# We focus on California sheephead (Semicossyphus pulcher), a key urchin predator.
# Sheephead are heavily fished outside MPAs, making them a good indicator species
# for MPA effectiveness — they should recover faster inside protected areas.

lter.fish <- read.csv(here::here("data", "LTER", "Annual_fish_comb_20240307.csv"))

# Validate expected columns
required_fish_cols <- c("SITE", "YEAR", "TRANSECT", "SCIENTIFIC_NAME",
                        "COUNT", "SIZE")
missing_fish_cols <- setdiff(required_fish_cols, names(lter.fish))
if (length(missing_fish_cols) > 0) {
  stop("LTER fish CSV is missing required columns: ",
       paste(missing_fish_cols, collapse = ", "))
}

cat("  LTER fish data loaded:", nrow(lter.fish), "rows\n")

# Subset to California sheephead only (Semicossyphus pulcher)
lter.fish.sub <- subset(lter.fish, SCIENTIFIC_NAME == "Semicossyphus pulcher")

# Merge with site metadata and filter to our study MPAs
lter.fish.sub.site <- validated_merge(lter.fish.sub, Sites2, by.x = c("SITE"), by.y = c("site_id"), warn_threshold = 0)
lter.fish.sub.site <- subset(lter.fish.sub.site, CA_MPA_Name_Short != "" & !(SITE %in% EXCLUDED_LTER_SITE_CODES))
lter.fish.sub.site <- subset(lter.fish.sub.site, CA_MPA_Name_Short == "Campus Point SMCA" |
                                CA_MPA_Name_Short == "Naples SMCA")

# Remove invalid count values (COUNT = -99999 is LTER's placeholder for missing data)
lter.fish.sub.site <- subset(lter.fish.sub.site, COUNT != -99999)  # -99999 = LTER missing data sentinel

# Calculate sheephead biomass using bio_sheephead() from 01_utils.R
# Allometric W = a * L^b where L = total length (cm)
# Coefficients from Cowen (1990), defined in 00c_analysis_constants.R
lter.fish.sub.site$biomass <- bio_sheephead(lter.fish.sub.site$SIZE)
n_na_fish_bio <- sum(is.na(lter.fish.sub.site$biomass))
lter.fish.sub.site$biomass[is.na(lter.fish.sub.site$biomass)] <- 0
if (n_na_fish_bio > 0) {
  cat(sprintf("  LTER sheephead: zero-filled %d NA biomass values (NA sizes)\n", n_na_fish_bio))
}

# Summarise: transect -> site -> MPA level (hierarchical averaging)
lter.fish.sum <- lter.fish.sub.site %>%
  group_by(YEAR, SITE, TRANSECT, status, CA_MPA_Name_Short,
           ChannelIsland, MPA_Start) %>%
  summarise_at(c("COUNT", "biomass"), sum)

lter.fish.ave <- lter.fish.sum %>%
  group_by(YEAR, SITE, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>%
  summarise_at(c("COUNT", "biomass"), mean)

lter.fish.max <- lter.fish.ave %>%
  group_by(status, CA_MPA_Name_Short, MPA_Start, YEAR) %>%
  summarise_at(c("COUNT", "biomass"), mean)

## --- Fish count proportions and log response ratios ---
# NOTE: Using adaptive correction (min/2) for consistency across all taxa/programs.
# This is statistically preferred over fixed +0.01 (Aitchison 1986) as it scales
# appropriately with the data and avoids inflating effect sizes for rare species.
lter.fish.max$taxon_name <- "SPUL"
lter.fish.max$y <- lter.fish.max$taxon_name  # Required for calculate_proportions()

# Rename YEAR to year for consistency with other data sources
names(lter.fish.max)[names(lter.fish.max) == "YEAR"] <- "year"

# Count proportions - use standardized function with adaptive correction
lter.fish.max <- calculate_proportions(lter.fish.max, "COUNT")

# Use column names instead of fragile indices
All.lter.fish.count.sub <- lter.fish.max[, c("CA_MPA_Name_Short", "year", "y", "status", "PropCorr")]
Short.lter.fish.count <- All.lter.fish.count.sub %>% spread(status, PropCorr)
Short.lter.fish.count.diff <- calculate_log_response_ratio(Short.lter.fish.count)
Short.lter.fish.count.diff$resp <- "Den"

# Biomass proportions - use standardized function with adaptive correction
lter.fish.max <- calculate_proportions(lter.fish.max, "biomass")

# Use column names instead of fragile indices
All.lter.spul.bio.sub <- lter.fish.max[, c("CA_MPA_Name_Short", "year", "y", "status", "PropCorr")]
Short.lter.spul.bio <- All.lter.spul.bio.sub %>% spread(status, PropCorr)
Short.lter.spul.bio.diff <- calculate_log_response_ratio(Short.lter.spul.bio)
Short.lter.spul.bio.diff$resp <- "Bio"

# Combine fish results
LTER.fish <- rbind(Short.lter.fish.count.diff, Short.lter.spul.bio.diff)

cat("  LTER fish response ratios:", nrow(LTER.fish), "rows\n")

# Assign time and BA labels using utility functions
LTER.fish <- assign_time_from_site_table(LTER.fish, Sites2)
LTER.fish <- merge(LTER.fish, sites.short.edit,
                   by.x = c("CA_MPA_Name_Short"),
                   by.y = c("CA_MPA_Name_Short"),
                   all.x = TRUE)
LTER.fish <- LTER.fish[complete.cases(LTER.fish$year), ]
LTER.fish$source <- "LTER"
LTER.fish <- assign_ba_from_site_table(LTER.fish, Sites2)

## --- Create LTER fish response dataframes ---
# Convert raw counts to per-m^2 density using LTER fish survey area
# LTER_FISH_SURVEY_AREA_M2 = 80 m^2 (defined in 00c_analysis_constants.R)
lter.fish.max$Bio <- lter.fish.max$biomass / LTER_FISH_SURVEY_AREA_M2
lter.fish.max$Den <- lter.fish.max$COUNT / LTER_FISH_SURVEY_AREA_M2

# Use column names instead of fragile indices
lter.fish.edit.den <- lter.fish.max[, c("status", "CA_MPA_Name_Short", "MPA_Start", "year", "Den", "y")]
LTER.fish.den <- lter.fish.edit.den %>% spread(status, Den)
LTER.fish.den$source <- "LTER"
LTER.fish.den.sub <- LTER.fish.den[complete.cases(LTER.fish.den), ]
LTER.fish.den.long <- gather(LTER.fish.den.sub, status, value, mpa, reference)
LTER.fish.den.long$resp <- "Den"

lter.fish.edit.bio <- lter.fish.max[, c("status", "CA_MPA_Name_Short", "MPA_Start", "year", "Bio", "y")]
LTER.fish.bio <- lter.fish.edit.bio %>% spread(status, Bio)
LTER.fish.bio$source <- "LTER"
LTER.fish.bio.sub <- LTER.fish.bio[complete.cases(LTER.fish.bio), ]
LTER.fish.bio.long <- gather(LTER.fish.bio.sub, status, value, mpa, reference)
LTER.fish.bio.long$resp <- "Bio"

LTER.fish.resp <- rbind(LTER.fish.den.long, LTER.fish.bio.long)

# ===========================================================================
# Section 7: Summary of outputs available for downstream scripts
# ===========================================================================
# LTER.join.ave    - Urchin density & biomass log response ratios with time, BA, site info
# LTER.lob         - Lobster density & biomass log response ratios with time, BA, site info
# Short.lter.macro - Macrocystis frond density & biomass log response ratios with time, BA
# LTER.fish        - Sheephead density & biomass log response ratios with time, BA
# LTER.resp        - Urchin response (den & bio, mpa vs reference, long format)
# LTER.lob.resp    - Lobster response (den & bio, mpa vs reference, long format)
# LTER.macro.resp  - Macrocystis response (den & bio, mpa vs reference, long format)
# LTER.fish.resp   - Fish response (den & bio, mpa vs reference, long format)

cat("LTER processing complete.\n")
cat("  Outputs: LTER.join.ave (", nrow(LTER.join.ave), " rows), ",
    "LTER.lob (", nrow(LTER.lob), " rows), ",
    "Short.lter.macro (", nrow(Short.lter.macro), " rows), ",
    "LTER.fish (", nrow(LTER.fish), " rows)\n", sep = "")
