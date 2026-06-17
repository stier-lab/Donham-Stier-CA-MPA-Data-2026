# =============================================================================
# 05_kfm_processing.R
# =============================================================================
#
# PURPOSE:
#   Process KFM (Kelp Forest Monitoring) data from the National Park Service
#   Channel Islands monitoring program for the MPA pBACIPS analysis.
#
# WHAT THIS SCRIPT DOES:
#   1. Imports MBON-integrated quad and swath data for KFM sites
#   2. Processes urchin data (S. purpuratus, M. franciscanus)
#   3. Processes Macrocystis kelp stipe data
#   4. Processes sheephead fish data (visual fish surveys)
#   5. Calculates biomass via bootstrap resampling
#   6. Calculates log response ratios (MPA vs reference sites)
#   7. Assigns time since MPA and Before/After labels
#
# DATA SOURCES:
#   KFM = NPS Kelp Forest Monitoring, Channel Islands National Park
#   Data comes via the MBON (Marine Biodiversity Observation Network) synthesis
#   which integrates multiple long-term monitoring programs.
#
# KEY DIFFERENCES FROM PISCO:
#   - KFM uses different transect sizes (urchin quads 1 m^2 in 1982-1984,
#     2 m^2 from 1985 onward; Macrocystis quads 1, 2, or 10 m^2; PISCO swath 60 m^2)
#   - KFM has longer time series (since 1982 for some sites)
#   - All KFM sites are in the Channel Islands
#   - Sheephead surveyed via visual fish counts and roving diver fish counts
#
# INPUTS:
#   - data/MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv
#   - data/MBON/SBCMBON_kelp_forest_integrated_fish_20231022.csv
#   - data/MBON/KFM_Macrocystis_RawData_1984-2023.csv
#   - data/MBON/SBCMBON_kelp_forest_site_geolocation_*.csv
#   - SizeFreq.Urch.OG (from 04_pisco_processing.R - unfiltered urchin sizes)
#
# OUTPUTS:
#   - KFM.join.ave: Response ratios for urchins and kelp
#   - KFM.resp: Raw density/biomass data (MPA vs reference)
#   - kfm.fish: Sheephead response ratios
#   - KFM.fish.den.long: Raw fish density data
#   - lter: LTER subset of MBON (passed to 06_lter_processing.R)
#
# DEPENDENCIES:
#   Requires 00-04 scripts to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

cat("Processing KFM/MBON data...\n")

# ===========================================================================
# Section 1: Import MBON quad/swath data and site geolocation
# ===========================================================================
# MBON (Marine Biodiversity Observation Network) integrates data from multiple programs.
# This file contains both KFM (Kelp Forest Monitoring - Channel Islands NPS) and
# LTER (Long-Term Ecological Research - mainland) data.
# We'll split them apart and process KFM here, LTER in 06_lter_processing.R.

MBON <- read.csv(here::here("data", "MBON", "SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv"))

# Validate that expected columns exist in the raw MBON data
required_mbon_cols <- c("data_source", "sample_method", "site_id", "subsite_id",
                        "transect_id", "replicate_id", "proj_taxon_id",
                        "auth_taxon_id", "auth_name", "taxon_name",
                        "site_name", "subsite_name", "area", "count",
                        "latitude", "longitude", "date")
missing_mbon_cols <- setdiff(required_mbon_cols, names(MBON))
if (length(missing_mbon_cols) > 0) {
  stop("MBON quad/swath CSV is missing required columns: ",
       paste(missing_mbon_cols, collapse = ", "),
       "\nCheck that the input file has not changed format.")
}

cat("  MBON quad/swath data loaded:", nrow(MBON), "rows,", ncol(MBON), "columns\n")

# Sites2 contains site coordinates and MPA information
# This table links site IDs to MPA names and geographic metadata
Sites2 <- read.csv(here::here("data", "MBON", "SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv"))
Sites2$site_id <- as.factor(Sites2$site_id)  # Convert to factor for proper merging

# Treat FORCE_BOOTSTRAP consistently across scripts: only force recompute when TRUE.
force_boot <- should_force_bootstrap()

# ===========================================================================
# Section 2: Process KFM subset: type conversions, date formatting, missing values
# ===========================================================================
# Data imported from CSV often needs type conversions.
# R may read text as "character" strings when we want "factor" (categorical variables)
# or numeric. This section standardizes column types.

# Factor conversions - factors are R's way of representing categorical variables
# They're more memory-efficient and have defined levels for grouping/plotting
factor_cols <- c("data_source", "sample_method", "site_id", "subsite_id",
                 "transect_id", "replicate_id", "proj_taxon_id",
                 "auth_taxon_id", "auth_name", "taxon_name",
                 "site_name", "subsite_name")

# Loop through columns and convert to factor if not already
# [[col]] accesses a column by name (like $col but works with variables)
for (col in factor_cols) {
  if (!is.factor(MBON[[col]])) MBON[[col]] <- as.factor(MBON[[col]])
}

# Numeric conversions
numeric_cols <- c("area", "count", "latitude", "longitude")
for (col in numeric_cols) {
  if (is.factor(MBON[[col]])) {
    MBON[[col]] <- as.numeric(levels(MBON[[col]]))[as.integer(MBON[[col]])]
  } else if (is.character(MBON[[col]])) {
    MBON[[col]] <- as.numeric(MBON[[col]])
  }
}

# Date conversion
tmpDateFormat <- "%Y-%m-%d"
tmp1date <- as.Date(MBON$date, format = tmpDateFormat)
if (length(tmp1date) == length(tmp1date[!is.na(tmp1date)])) {
  MBON$date <- tmp1date
} else {
  stop("Date conversion failed for MBON$date. ",
       sum(is.na(tmp1date)), " of ", length(tmp1date),
       " dates could not be parsed with format '", tmpDateFormat, "'.")
}
rm(tmpDateFormat, tmp1date)

# Convert missing value placeholders to NA
# The MBON CSV uses "." as a missing value indicator in some columns
MBON$count <- ifelse(trimws(as.character(MBON$count)) == trimws("."), NA, MBON$count)
suppressWarnings(MBON$count <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(MBON$count)) == as.character(as.numeric("."))), NA, MBON$count))
MBON$auth_taxon_id <- as.factor(ifelse(trimws(as.character(MBON$auth_taxon_id)) == trimws("."), NA, as.character(MBON$auth_taxon_id)))
MBON$auth_name <- as.factor(ifelse(trimws(as.character(MBON$auth_name)) == trimws("."), NA, as.character(MBON$auth_name)))
MBON$subsite_name <- as.factor(ifelse(trimws(as.character(MBON$subsite_name)) == trimws("."), NA, as.character(MBON$subsite_name)))
MBON$latitude <- ifelse(trimws(as.character(MBON$latitude)) == trimws("."), NA, MBON$latitude)
suppressWarnings(MBON$latitude <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(MBON$latitude)) == as.character(as.numeric("."))), NA, MBON$latitude))
MBON$longitude <- ifelse(trimws(as.character(MBON$longitude)) == trimws("."), NA, MBON$longitude)
suppressWarnings(MBON$longitude <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(MBON$longitude)) == as.character(as.numeric("."))), NA, MBON$longitude))

# Split into KFM and LTER subsets based on data_source column
# KFM data is processed here; LTER data is passed to 06_lter_processing.R
kfm <- subset(MBON, data_source == "kfm")
lter <- subset(MBON, data_source == "lter")  # Used by 06_lter_processing.R

stopifnot("KFM subset must have > 0 rows" = nrow(kfm) > 0)
stopifnot("LTER subset must have > 0 rows" = nrow(lter) > 0)
cat("  Split MBON: KFM =", nrow(kfm), "rows, LTER =", nrow(lter), "rows\n")

# ===========================================================================
# Section 3: Remove excluded KFM sites and join to site table
# ===========================================================================
# merge() joins the survey data with site metadata using site_id as the key
# This adds columns like CA_MPA_Name_Short, status (mpa/reference), etc.
kfm.site <- validated_merge(kfm, Sites2, by = c("site_id"), warn_threshold = 0)

# Extract year from date column using lubridate's year() function
kfm.site$year <- year(kfm.site$date)

# Remove empty MPA names (sites not associated with any MPA)
kfm.site <- subset(kfm.site, CA_MPA_Name_Short != "")

# Remove excluded KFM sites using constant from 00c_analysis_constants.R / 01_utils.R
# EXCLUDED_KFM_SITES contains site IDs that don't meet data quality standards
# (e.g., insufficient temporal coverage, inconsistent sampling protocols)
# The %in% operator checks membership; ! negates to keep sites NOT in the list
kfm.site <- kfm.site[!(kfm.site$site_id %in% EXCLUDED_KFM_SITES), ]

cat("  KFM after site filtering:", nrow(kfm.site), "rows\n")

# ===========================================================================
# Section 4: Assign time since MPA implementation
# ===========================================================================

kfm.site <- assign_time_from_site_table(kfm.site, Sites2)

# ===========================================================================
# Section 5: Subset target taxa
# ===========================================================================
# Focus on the 4 focal taxa in the trophic cascade:
#   - Panulirus interruptus (spiny lobster, predator)
#   - Strongylocentrotus purpuratus (purple urchin, herbivore)
#   - Mesocentrotus franciscanus (red urchin, herbivore)
#   - Macrocystis pyrifera (giant kelp, primary producer)
# Note: Sheephead (Semicossyphus pulcher) is processed separately from
# the MBON fish file in Section 13 below.

kfm.site <- subset(kfm.site, taxon_name == "Strongylocentrotus purpuratus" |
                     taxon_name == "Mesocentrotus franciscanus" |
                     taxon_name == "Macrocystis pyrifera" |
                     taxon_name == "Panulirus interruptus")
# KFM LOBSTER NOTE (audit 2026-05-04): KFM lobsters (Panulirus interruptus)
# are surveyed on belt transects whose dimensions changed in 1985, per
# Davis et al. 1997 KFM Handbook Vol 1, p.39:
#   - 1983-1984: 10 transects per site, each 2 m x 20 m = 40 m^2.
#   - 1985-present: 12 transects per site, each 3 m x 20 m = 60 m^2.
# Mechanism: two divers each lay a 10 m perpendicular tape from the main
# transect, then count organisms within 1.5 m on each side (each diver
# covers 3 m x 10 m = 30 m^2; CountA + CountB = 3 m x 20 m = 60 m^2).
# Density is computed at L206 below via safe_divide(count, area) so the
# row-level `area` is used and the temporal transect-size change is
# handled correctly. KFM does NOT collect lobster sizes - there is no
# bootstrap step for KFM lobster, so KFM lobster contributes ONLY to
# density (no biomass) in the harmonized outputs. PISCO and LTER lobster
# biomass come through their own size-aware pipelines (UCSB direct CL
# measurement, VRG roving-diver bootstrap, LTER direct CL measurement).

# Remove juvenile giant kelp (proj_taxon_id "t-k-002" = juvenile Macrocystis).
# PISCO and LTER record only adult Macrocystis; excluding KFM juveniles keeps
# the cross-program comparison on adult stipe counts.
kfm.site <- subset(kfm.site, taxon_name != "Macrocystis pyrifera" | proj_taxon_id != "t-k-002")

cat("  KFM target taxa rows:", nrow(kfm.site), "\n")

# Sum across the adult-Macrocystis proj_taxon_ids retained above to get the
# total stipe count per site-year-replicate. (Re: Emily's "where did holdfast
# come from?" review note) KFM records adult giant kelp under more than one
# proj_taxon_id; these are treated as non-overlapping tallies of stipes for the
# same replicate and summed. This assumes the IDs do not double-count the same
# plant. VERIFY against the KFM proj_taxon_id scheme if a definitive
# stipe-vs-holdfast distinction is needed; the bootstrap is otherwise robust to
# this because lnRR is proportion-based.
kfm.sum <- kfm.site %>%
  dplyr::group_by(site_id, taxon_name, status, CA_MPA_Name_Short, area, ChannelIsland, replicate_id, MPA_Start, time, year) %>%
  dplyr::summarise_at(c("count"), sum) %>%
  dplyr::ungroup()

# Transform to density since there are different survey methods
# Use safe_divide to handle any zero area values
kfm.sum$den <- safe_divide(kfm.sum$count, kfm.sum$area, context = "KFM density calculation")

# Average across sites and replicates
kfm.ave <- kfm.sum %>%
  dplyr::group_by(site_id, taxon_name, status, CA_MPA_Name_Short, area, ChannelIsland, MPA_Start, time, year) %>%
  dplyr::summarise_at(c("den"), mean) %>%
  dplyr::ungroup()

# ===========================================================================
# Section 6: Calculate urchin biomass via bootstrap resampling
# ===========================================================================
# KFM records urchin COUNTS per quadrat (no per-quadrat sizes), but it DOES run
# separate site-level urchin size-frequency surveys ("Natural Habitat", 1985-2023).
# Those size data are present in data/ALL_sizefreq_2024.csv (campus == "KFM") and
# are pulled into the combined urchin size pool by 03_data_import.R (recoded to
# STRPUR/MESFRA). We estimate biomass by bootstrap-resampling sizes from the
# combined pool SizeFreq.Urch.OG, then applying species-specific allometric
# equations (bio_redurch for M. franciscanus / red urchin, bio_purpurch for
# S. purpuratus / purple urchin).
#
# HISTORY (corrected 2026-06, Emily review): KFM urchin sizes were present in the
# pool but unusable -- 04_pisco_processing.R attached MPA metadata to the size
# rows via the PISCO-only site table, leaving KFM (and LTER) rows with
# CA_MPA_Name_Short = NA, so they never matched and KFM biomass fell back to PISCO
# sizes. (This bug was also present in Emily's V10, which merged the pool with the
# same PISCO-only table; she flagged the symptom in review.) 04 now backfills MPA
# metadata onto the KFM/LTER size rows, so the bootstrap below draws from the
# combined PISCO+KFM pool -- which at KFM's Channel-Island MPAs is overwhelmingly
# KFM's own sizes. This realizes Emily's intended design (review note #19).
#
#   MATCHING: the bootstrap draws sizes for each KFM count by
#   CA_MPA_Name_Short, site_status (mpa/reference), year, and species
#   (classcode) — see the loop below. When no size data exist for a given
#   MPA/year/status/species cell, biomass is set to NA and the observation is
#   excluded.
#
#   25 mm THRESHOLD: KFM uses SizeFreq.Urch.OG (the UNFILTERED pool); the 25 mm
#   PISCO cutoff is not applied to it, since KFM counts may include urchins
#   < 25 mm. Small urchins contribute negligible biomass (allometry ~cubic in
#   diameter), so this choice is immaterial to biomass.
#
#   LTER NOTE: 06_lter_processing.R uses the same combined pool (SizeFreq.Urch.OG)
#   and is fixed by the same 04 metadata backfill -- LTER urchins now draw from
#   LTER's own sizes (classcodes SFL/SFS/SPS/SPL) plus any PISCO sizes at the
#   same MPA/year.

URCHINS <- subset(kfm.site, taxon_name == "Strongylocentrotus purpuratus" |
                    taxon_name == "Mesocentrotus franciscanus")
URCHINS <- dplyr::mutate(URCHINS, taxon_name = dplyr::case_when(
  taxon_name == "Strongylocentrotus purpuratus" ~ "STRPURAD",
  taxon_name == "Mesocentrotus franciscanus" ~ "MESFRAAD",
  TRUE ~ taxon_name))

u <- unique(URCHINS[, c("year", "site_id", "CA_MPA_Name_Short", "status", "area", "taxon_name")])
cat("  KFM urchin bootstrap: ", nrow(u), " site-year-species combinations\n")

# Urchin size pool for the KFM bootstrap. SizeFreq.Urch.OG is the combined,
# all-program urchin size-frequency pool built in 04_pisco_processing.R -- now
# with MPA metadata attached to KFM/LTER rows as well (see the metadata-backfill
# step in 04). This realizes Emily's intended design (review note #19: the size
# pool is PISCO + KFM, mostly KFM, and is used for biomass across programs):
# each KFM count below draws sizes matched by MPA/status/year/species from the
# combined pool, which at KFM's Channel-Island MPAs is dominated by KFM's own
# size-frequency surveys. (Previously KFM/LTER size rows had CA_MPA_Name_Short =
# NA and never matched, so KFM biomass fell back to PISCO sizes -- the latent
# bug Emily flagged in her review.)
SizeFreq.Urch <- SizeFreq.Urch.OG

# Cache bootstrap results to avoid re-running the slow loop
.cache_kfm_urch <- here::here("data", "cache", "kfm_urchin_bootstrap.rds")
if (file.exists(.cache_kfm_urch) && !force_boot) {
  cat("    Loading cached KFM urchin bootstrap results...\n")
  Urchin.site <- safe_read_rds(.cache_kfm_urch)
  if (is.null(Urchin.site)) {
    cat("    Cache invalid, will recompute...\n")
  } else {
    # Cache schema can evolve; invalidate cache if required columns are missing.
    required_urch_cols <- c(
      "site", "CA_MPA_Name_Short", "site_status",
      "year", "area", "transect",
      "count", "biomass", "y"
    )
    missing_urch_cols <- setdiff(required_urch_cols, names(Urchin.site))
    if (length(missing_urch_cols) > 0) {
      cat("    Cached KFM urchin bootstrap missing columns (",
          paste(missing_urch_cols, collapse = ", "),
          "); will recompute...\n", sep = "")
      Urchin.site <- NULL
    }
  }
}

if (!exists("Urchin.site") || is.null(Urchin.site)) {
  # Pre-allocate list for results (PERFORMANCE: avoids O(n^2) rbind copies)
  results_list <- vector("list", nrow(u))

  # Set random seed for reproducibility of bootstrap resampling
  # SEED_KFM_URCHIN defined in 00c_analysis_constants.R
  set.seed(SEED_KFM_URCHIN)
  for (i in seq_len(nrow(u))) {
    t <- which(URCHINS$site_id == u$site_id[i] &
                 URCHINS$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] &
                 URCHINS$year == u$year[i] &
                 URCHINS$taxon_name == u$taxon_name[i])
    n <- sum(URCHINS$count[t], na.rm = TRUE)
    t2 <- which(SizeFreq.Urch$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] &
                  SizeFreq.Urch$site_status == u$status[i] &
                  SizeFreq.Urch$year == u$year[i] &
                  SizeFreq.Urch$classcode == u$taxon_name[i])

    # Select biomass function based on species:
    #   bio_redurch: M. franciscanus (red urchin) allometric equation
    #   bio_purpurch: S. purpuratus (purple urchin) allometric equation
    bio_fun <- if (u$taxon_name[i] == "MESFRAAD") bio_redurch else bio_purpurch

    # Wrap bootstrap in tryCatch to capture failures
    result <- tryCatch(
      bootstrap_biomass(
        count = n,
        size_freq_indices = t2,
        size_freq_table = SizeFreq.Urch,
        biomass_fun = bio_fun
      ),
      error = function(e) {
        warning("Bootstrap failed for site ", u$site_id[i], " year ", u$year[i],
                " species ", u$taxon_name[i], ": ", e$message)
        list(biomass = NA, se = NA, count = n)
      }
    )

    results_list[[i]] <- data.frame(
      site = u$site_id[i],
      CA_MPA_Name_Short = u$CA_MPA_Name_Short[i],
      site_status = u$status[i],
      year = u$year[i],
      area = u$area[i],
      transect = length(t),
      count = result$count,
      biomass = result$biomass,
      y = u$taxon_name[i]
    )
  }

  # Combine all results at once (efficient)
  Urchin.site <- do.call(rbind, results_list)

  # Save cache
  dir.create(dirname(.cache_kfm_urch), recursive = TRUE, showWarnings = FALSE)
  saveRDS(Urchin.site, .cache_kfm_urch)
  cat("    Cached KFM urchin bootstrap results.\n")
}

stopifnot("Urchin.site must have > 0 rows" = nrow(Urchin.site) > 0)
cat("  KFM urchin bootstrap complete:", nrow(Urchin.site), "rows\n")

# Urchin.site already contains KFM/MBON site_status and MPA assignment.
# Do NOT merge with PISCO's `sites.short` (which may exist in the session and has different site IDs).
KFM.Urchin.site.merge <- Urchin.site

# Remove Old Anacapa MPA and the SMCA since they don't meet criteria for inclusion
# (listed in EXCLUDED_MPAS in 00c_analysis_constants.R)
KFM.Urchin.site.merge <- subset(
  KFM.Urchin.site.merge,
  !(CA_MPA_Name_Short %in% c("Anacapa Island SMCA", "Anacapa Island SMR 1978"))
)

# Find sites where urchin was found on transect but no size was recorded
# These can't produce accurate biomass so they are excluded
KFM.Urchin.null <- subset(KFM.Urchin.site.merge, count > 0 & (is.na(biomass) | biomass == 0))
KFM.Urchin.site.merge.sub <- KFM.Urchin.site.merge

# Remove sites where urchins were seen but not measured (can't get accurate biomass)
for (i in seq_len(nrow(KFM.Urchin.null))) {
  KFM.Urchin.site.merge.sub <- subset(KFM.Urchin.site.merge.sub,
                                        site != KFM.Urchin.null$site[i] |
                                          year != KFM.Urchin.null$year[i] |
                                          y != KFM.Urchin.null$y[i])
}

if (nrow(KFM.Urchin.null) > 0) {
  cat("  Removed", nrow(KFM.Urchin.null), "urchin rows with count > 0 but no biomass\n")
}

# Account for different number of transects at sites
KFM.Urchin.site.merge.sub$count.ave <- safe_divide(
  KFM.Urchin.site.merge.sub$count,
  KFM.Urchin.site.merge.sub$transect,
  context = "KFM urchin count/transect"
)
KFM.Urchin.site.merge.sub$biomass.ave <- safe_divide(
  KFM.Urchin.site.merge.sub$biomass,
  KFM.Urchin.site.merge.sub$transect,
  context = "KFM urchin biomass/transect"
)

KFM.Urchin.site.all <- KFM.Urchin.site.merge.sub

# Convert per-quadrat counts/biomass to per-m^2 densities.
# KFM urchin quadrat sizes changed over time: 1 m^2 in 1982-1984, then 2 m^2 from
# 1985 onward (per data audit 2026-05-03). Use the row-level `area` column rather
# than a single constant. The `u` table at L254 deduplicates on `area`, so each
# bootstrap row carries the correct quad size for its year. Earlier code divided
# by a fixed KFM_QUAD_AREA_M2 = 2, which under-estimated densities by 50% for
# 1982-1984 observations.
KFM.Urchin.site.all$Density <- safe_divide(
  KFM.Urchin.site.all$count.ave, KFM.Urchin.site.all$area,
  context = "KFM urchin density per m^2"
)
KFM.Urchin.site.all$bio.m2 <- safe_divide(
  KFM.Urchin.site.all$biomass.ave, KFM.Urchin.site.all$area,
  context = "KFM urchin biomass per m^2"
)

# Restore full species names (reverse the abbreviation from Section 6)
KFM.Urchin.site.all <- dplyr::mutate(KFM.Urchin.site.all, y = dplyr::case_when(
  y == "STRPURAD" ~ "Strongylocentrotus purpuratus",
  y == "MESFRAAD" ~ "Mesocentrotus franciscanus",
  TRUE ~ y))

# ===========================================================================
# Section 7: Process KFM Macrocystis stipe data from raw stipe counts file
# ===========================================================================
# KFM divers count stipes on ~100 individual kelp plants per site per survey
# (the KFM_Macrocystis_RawData file). For our biomass conversion, "stipes"
# and "fronds" refer to the same structural unit. Li (LTER) confirmed the
# LTER and KFM raw counts are the same quantity despite the differing
# terminology in each program's metadata. We therefore treat stipes = fronds
# throughout and apply the bio_macro() conversion (01_utils.R) which uses an
# averaged stipe-to-biomass slope from published allometric relationships.

kfm.stipe <- read.csv(here::here("data", "MBON", "KFM_Macrocystis_RawData_1984-2023.csv"))

# Validate expected columns in stipe data
required_stipe_cols <- c("SiteNumber", "SurveyDate", "ScientificName", "Stipe_Count")
missing_stipe_cols <- setdiff(required_stipe_cols, names(kfm.stipe))
if (length(missing_stipe_cols) > 0) {
  stop("KFM Macrocystis raw data CSV is missing required columns: ",
       paste(missing_stipe_cols, collapse = ", "))
}

kfm.stipe <- kfm.stipe[!duplicated(kfm.stipe), ] # Data issues caused duplicates
cat("  KFM stipe data loaded:", nrow(kfm.stipe), "rows (after deduplication)\n")

# Create site_id matching site table format
# KFM site IDs follow the pattern "a-k-NN" with zero-padded two-digit numbers
kfm.stipe$site_id <- ifelse(kfm.stipe$SiteNumber > 9,
                              paste0("a-k-", kfm.stipe$SiteNumber),
                              paste0("a-k-0", kfm.stipe$SiteNumber))
kfm.stipe.site <- validated_merge(kfm.stipe, Sites2, by = c("site_id"), warn_threshold = 0)
kfm.stipe.site <- kfm.stipe.site %>%
  dplyr::mutate(date = as.Date(SurveyDate, format = "%m/%d/%Y"))
kfm.stipe.site$year <- year(kfm.stipe.site$date)

# Remove excluded KFM sites (from EXCLUDED_KFM_SITES constant)
kfm.stipe.site <- kfm.stipe.site[!(kfm.stipe.site$site_id %in% EXCLUDED_KFM_SITES), ]

# Calculate average number of stipes per site
kfm.stipe.mean <- kfm.stipe.site %>%
  dplyr::group_by(site_id, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>%
  dplyr::summarise_at(c("Stipe_Count"), mean)

# Calculate average number of stipes per MPA/Reference
kfm.stipe.max <- kfm.stipe.mean %>%
  dplyr::group_by(status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>%
  dplyr::summarise_at(c("Stipe_Count"), mean)

kfm.stipe.ave <- subset(kfm.stipe.max, CA_MPA_Name_Short != "")

# Assign time since MPA using utility function
kfm.stipe.ave <- assign_time_from_site_table(kfm.stipe.ave, Sites2)

# ===========================================================================
# Section 8: Calculate Macrocystis biomass via bootstrap resampling
# ===========================================================================
# For Macrocystis, the bootstrap resamples stipe counts per individual from the
# raw stipe count data, then applies the bio_macro() allometric function to
# convert stipe counts to biomass (grams wet weight).

MACRO <- subset(kfm.sum, taxon_name == "Macrocystis pyrifera")
# Remove Old Anacapa MPA and the SMCA (listed in EXCLUDED_MPAS)
MACRO <- subset(MACRO, CA_MPA_Name_Short != "Anacapa Island SMCA" &
                  CA_MPA_Name_Short != "Anacapa Island SMR 1978")

u <- unique(MACRO[, c("year", "site_id", "CA_MPA_Name_Short", "status", "area", "taxon_name")])
cat("  KFM Macrocystis bootstrap:", nrow(u), "site-year combinations\n")

kfm.stipe.OG <- kfm.stipe.site

# Cache bootstrap results to avoid re-running the slow loop
.cache_kfm_macro <- here::here("data", "cache", "kfm_macro_bootstrap.rds")
if (file.exists(.cache_kfm_macro) && !force_boot) {
  cat("    Loading cached KFM Macrocystis bootstrap results...\n")
  Macro.site <- safe_read_rds(.cache_kfm_macro)
  if (is.null(Macro.site)) {
    cat("    Cache invalid, will recompute...\n")
  } else {
    required_macro_cols <- c(
      "site", "CA_MPA_Name_Short", "site_status",
      "year", "transect", "area",
      "y", "count", "biomass", "ind"
    )
    missing_macro_cols <- setdiff(required_macro_cols, names(Macro.site))
    if (length(missing_macro_cols) > 0) {
      cat("    Cached KFM Macrocystis bootstrap missing columns (",
          paste(missing_macro_cols, collapse = ", "),
          "); will recompute...\n", sep = "")
      Macro.site <- NULL
    }
  }
}

if (!exists("Macro.site") || is.null(Macro.site)) {
  # Pre-allocate list for results (PERFORMANCE: avoids O(n^2) rbind copies)
  results_list <- vector("list", nrow(u))

  # N_BOOTSTRAP_RESAMPLES defined in 00c_analysis_constants.R (default 1000)
  n_boot <- N_BOOTSTRAP_RESAMPLES

  # Set random seed for reproducibility of bootstrap resampling
  # SEED_KFM_MACRO defined in 00c_analysis_constants.R
  set.seed(SEED_KFM_MACRO)
  for (i in seq_len(nrow(u))) {
    t <- which(MACRO$site_id == u$site_id[i] &
                 MACRO$year == u$year[i] &
                 MACRO$area == u$area[i] &
                 MACRO$taxon_name == u$taxon_name[i])
    n <- sum(MACRO$count[t], na.rm = TRUE)
    # FIX [C4]: Use derived `year` column (line 335, from SurveyDate) to match
    # the numeric year in `u` (from MACRO). `SurveyYear` is the raw CSV column
    # and could differ in type. `ScientificName` is correct (raw CSV column,
    # contains full binomial e.g. "Macrocystis pyrifera").
    t2 <- which(kfm.stipe.site$site_id == u$site_id[i] &
                  kfm.stipe.site$year == u$year[i] &
                  kfm.stipe.site$ScientificName == u$taxon_name[i])

    # Wrap bootstrap in tryCatch to capture failures
    result_row <- tryCatch({
      if (n != 0 & length(t2) != 0) {
        # Kelp found on transect and in stipe data: bootstrap stipe counts
        a <- kfm.stipe.site[t2, ]
        s <- data.frame(matrix(NA, nrow = n_boot, ncol = n))
        for (k in seq_len(n_boot)) {
          s[k, ] <- sample(a$Stipe_Count, n, replace = TRUE)
        }
        # Apply bio_macro() allometric function to convert stipe counts to biomass
        s_bio <- s %>%
          dplyr::mutate(dplyr::across(starts_with("X"), ~ bio_macro(.)))
        s_bio <- s_bio %>%
          dplyr::mutate(sum = rowSums(.))
        s <- s %>%
          dplyr::mutate(sum = rowSums(.))
        aveB <- mean(s_bio$sum)
        aveD <- mean(s$sum)

        data.frame(site = u$site_id[i], CA_MPA_Name_Short = u$CA_MPA_Name_Short[i],
                   site_status = u$status[i],
                   year = u$year[i], transect = length(t), area = u$area[i],
                   y = u$taxon_name[i], count = aveD, biomass = aveB, ind = n)

      } else if (n == 0 & length(t2) == 0) {
        # No kelp on transect or in stipe data: true zero
        data.frame(site = u$site_id[i], CA_MPA_Name_Short = u$CA_MPA_Name_Short[i],
                   site_status = u$status[i],
                   year = u$year[i], transect = length(t), area = u$area[i],
                   y = u$taxon_name[i], count = 0, biomass = 0, ind = 0)

      } else if (n != 0 & length(t2) == 0) {
        # Kelp on transect but no stipe data: can't estimate biomass, mark as NA
        data.frame(site = u$site_id[i], CA_MPA_Name_Short = u$CA_MPA_Name_Short[i],
                   site_status = u$status[i],
                   year = u$year[i], transect = length(t), area = u$area[i],
                   y = u$taxon_name[i], count = NA, biomass = NA, ind = n)

      } else {
        # n == 0 but stipe data exists: no individuals on transect
        data.frame(site = u$site_id[i], CA_MPA_Name_Short = u$CA_MPA_Name_Short[i],
                   site_status = u$status[i],
                   year = u$year[i], transect = length(t), area = u$area[i],
                   y = u$taxon_name[i], count = 0, biomass = 0, ind = n)
      }
    }, error = function(e) {
      warning("KFM Macro bootstrap failed for site ", u$site_id[i], " year ", u$year[i], ": ", e$message)
      data.frame(site = u$site_id[i], CA_MPA_Name_Short = u$CA_MPA_Name_Short[i],
                 site_status = u$status[i],
                 year = u$year[i], transect = NA, area = u$area[i],
                 y = u$taxon_name[i], count = NA, biomass = NA, ind = NA)
    })

    results_list[[i]] <- result_row
  }

  # Combine all results at once (efficient)
  Macro.site <- do.call(rbind, results_list)

  # Save cache
  dir.create(dirname(.cache_kfm_macro), recursive = TRUE, showWarnings = FALSE)
  saveRDS(Macro.site, .cache_kfm_macro)
  cat("    Cached KFM Macrocystis bootstrap results.\n")
}

stopifnot("Macro.site must have > 0 rows" = nrow(Macro.site) > 0)
cat("  KFM Macrocystis bootstrap complete:", nrow(Macro.site), "rows\n")

# Macro.site already contains KFM/MBON site_status and MPA assignment.
# Do NOT merge with PISCO's `sites.short` (which may exist in the session and has different site IDs).
KFM.Macro.site.merge <- Macro.site

# Convert to per m^2 by dividing by transects and area
# Use safe_divide to handle any zero area/transect values
area_transect <- KFM.Macro.site.merge$area * KFM.Macro.site.merge$transect
KFM.Macro.site.merge$biomassCorr <- safe_divide(KFM.Macro.site.merge$biomass, area_transect, context = "KFM macro biomass/m2")
KFM.Macro.site.merge$densityCorr <- safe_divide(KFM.Macro.site.merge$count, area_transect, context = "KFM macro density/m2")
KFM.Macro.site.merge$densityIndCorr <- safe_divide(KFM.Macro.site.merge$ind, area_transect, context = "KFM macro ind density/m2")

# ===========================================================================
# Section 9: Combine urchin and Macrocystis data, calculate density per m^2
# ===========================================================================

# Wrangle macro columns
KFM.Macro.site.all.sub <- KFM.Macro.site.merge[, c(
  "site", "site_status", "year", "y", "CA_MPA_Name_Short", "area", "transect",
  "count", "biomass", "ind",
  "biomassCorr", "densityCorr", "densityIndCorr"
)]
names(KFM.Macro.site.all.sub)[names(KFM.Macro.site.all.sub) == "biomass"] <- "biomassRaw"
names(KFM.Macro.site.all.sub)[names(KFM.Macro.site.all.sub) == "biomassCorr"] <- "biomass"
names(KFM.Macro.site.all.sub)[names(KFM.Macro.site.all.sub) == "densityCorr"] <- "density"

# Wrangle urchin columns
KFM.Urchin.site.all.sub <- KFM.Urchin.site.all[, c(
  "site", "CA_MPA_Name_Short", "site_status", "year", "area", "transect",
  "count", "biomass", "y",
  "Density", "bio.m2"
)]
names(KFM.Urchin.site.all.sub)[names(KFM.Urchin.site.all.sub) == "Density"] <- "density"
names(KFM.Urchin.site.all.sub)[names(KFM.Urchin.site.all.sub) == "biomass"] <- "biomassRaw"
names(KFM.Urchin.site.all.sub)[names(KFM.Urchin.site.all.sub) == "bio.m2"] <- "biomass"
KFM.Urchin.site.all.sub$ind <- NA_real_
KFM.Urchin.site.all.sub$densityIndCorr <- NA

# Ensure identical column sets + ordering before rbind (prevents 'numbers of columns' errors).
common_cols <- c(
  "site", "site_status", "year", "y", "CA_MPA_Name_Short", "area", "transect",
  "count", "biomassRaw", "ind", "biomass", "density", "densityIndCorr"
)
KFM.Macro.site.all.sub <- KFM.Macro.site.all.sub[, common_cols]
KFM.Urchin.site.all.sub <- KFM.Urchin.site.all.sub[, common_cols]

# Combine urchin + Macrocystis data
KFM.bio.site.all <- rbind(KFM.Macro.site.all.sub, KFM.Urchin.site.all.sub)

# Merge with average density data
# After merging: den = density of individuals; for kelp, density = stipe density
KFM.allbio <- merge(kfm.ave, KFM.bio.site.all,
                     by.x = c("site_id", "status", "year", "taxon_name", "CA_MPA_Name_Short", "area"),
                     by.y = c("site", "site_status", "year", "y", "CA_MPA_Name_Short", "area"),
                     all.x = TRUE)

Swath.KFM <- KFM.allbio
names(Swath.KFM)[names(Swath.KFM) == "taxon_name"] <- "y"
names(Swath.KFM)[names(Swath.KFM) == "count"] <- "countRaw"
names(Swath.KFM)[names(Swath.KFM) == "density"] <- "den.y"
names(Swath.KFM)[names(Swath.KFM) == "site_id"] <- "site"

Swath.KFM.Corr <- Swath.KFM
# For Macrocystis, use stipe density (den.y from biomass data) rather than individual
# density (den from count data). This makes kelp density comparable across programs
# since KFM counts stipes directly while other programs use different units.
Swath.KFM.Corr$den <- ifelse(Swath.KFM.Corr$y == "Macrocystis pyrifera",
                               Swath.KFM.Corr$den.y, Swath.KFM.Corr$den)

# Calculate the mean across sites
kfm.ave.ave <- Swath.KFM.Corr %>%
  dplyr::group_by(y, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, time, year, area) %>%
  dplyr::summarise_at(c("den", "biomass"), mean)

# Remove MPAs that don't meet criteria for inclusion
# (listed in EXCLUDED_MPAS in 00c_analysis_constants.R)
kfm.ave.ave <- subset(kfm.ave.ave,
                        CA_MPA_Name_Short != "Anacapa Island SMCA" &
                          CA_MPA_Name_Short != "Anacapa Island SMR 1978" &
                          CA_MPA_Name_Short != "Carrington Point SMR" &
                          CA_MPA_Name_Short != "Painted Cave SMCA")

cat("  KFM combined taxa data:", nrow(kfm.ave.ave), "rows\n")

# ===========================================================================
# Section 10: Calculate proportions and log response ratios
# ===========================================================================
# This section converts raw density/biomass values to LOG RESPONSE RATIOS (lnRR).
# lnRR compares MPA sites to reference sites, standardized across different MPAs.
#
# Steps:
# 1. calculate_proportions(): Convert to proportion of maximum (0-1 scale + adaptive correction)
# 2. spread(): Reshape from long (one row per status) to wide (mpa/reference columns)
# 3. calculate_log_response_ratio(): Calculate lnRR = ln(mpa / reference)

# -- KFM Density proportions and response ratios --
All.den <- kfm.ave.ave
All.den <- All.den[!is.na(All.den$den), ]  # Remove invalid values

# calculate_proportions() from 01_utils.R normalizes by max value within each MPA x taxa
All.den <- calculate_proportions(All.den, "den")

# Select columns needed for reshaping and calculate response ratio
All.den.sub <- All.den[, c("CA_MPA_Name_Short", "year", "y", "status", "area", "PropCorr")]

# spread() from tidyr: converts long format to wide format
# status column values (mpa, reference) become new columns with PropCorr values
Short.den <- All.den.sub %>%
  spread(status, PropCorr)

# calculate_log_response_ratio() from 01_utils.R: lnRR = ln(mpa / reference)
Short.den.diff <- calculate_log_response_ratio(Short.den)
Short.den.diff$resp <- "Den"  # Label as density response

# -- KFM Biomass proportions and response ratios --
# Same process as density, but for biomass values
All.bio <- kfm.ave.ave
All.bio <- All.bio[!is.na(All.bio$biomass), ]
All.bio <- calculate_proportions(All.bio, "biomass")

All.bio.sub <- All.bio[, c("CA_MPA_Name_Short", "year", "y", "status", "area", "PropCorr")]
Short.bio <- All.bio.sub %>%
  spread(status, PropCorr)
Short.bio.diff <- calculate_log_response_ratio(Short.bio)
Short.bio.diff$resp <- "Bio"  # Label as biomass response

# Combine density and biomass response ratios into one dataframe
# rbind() stacks dataframes vertically (row bind)
KFM.join.ave <- rbind(Short.den.diff, Short.bio.diff)

cat("  KFM response ratios:", nrow(KFM.join.ave), "rows (density + biomass)\n")

# ===========================================================================
# Section 11: Create KFM.resp dataframe with raw density/biomass
# ===========================================================================

# Filter Macrocystis to the modern (1996+) 5 m quadrat protocol.
#
# AREA == 10 EXPLAINED (audit 2026-05-04):
# In the MBON-integrated CSV each Macrocystis row has area = 1, 2, or 10 m^2.
# These correspond to KFM's three quadrat eras (per Davis et al. 1997
# Handbook Vol 1, pp.13, 35, 38):
#   - area = 1 m^2 (1982-1984): legacy 1 m quadrat protocol, 30-40 quads/site.
#   - area = 2 m^2 (1985+):     1 m quadrat protocol with paired-diver
#                                 redesign (two divers each sample adjacent
#                                 1 m^2 quads on opposite sides; entered as
#                                 the 2 m^2 paired total).
#   - area = 10 m^2 (1996+):    Macrocystis-specific 5 m x 1 m quadrat
#                                 protocol added in 1996. Each diver samples
#                                 a 5 m^2 half-quad; the paired total
#                                 (CountA + CountB) covers 10 m^2. MBON
#                                 stores this paired total as area = 10.
# We keep only the 10 m^2 records (the modern adult-Macrocystis-targeted
# protocol). This drops the 1 m^2 era entirely and the 2 m^2 records, which
# would mix sampling units across years. See METHODS_REFERENCES.md for the
# full provenance.
kfm.edit.den <- subset(kfm.ave.ave, y != "Macrocystis pyrifera" | area != 2)
kfm.edit.den <- subset(kfm.edit.den, y != "Macrocystis pyrifera" | area != 1)

kfm.edit.den <- kfm.edit.den[, c("status", "CA_MPA_Name_Short", "MPA_Start", "year", "den", "y")]
kfm.edit.den <- kfm.edit.den[complete.cases(kfm.edit.den), ]
KFM.den <- kfm.edit.den %>%
  spread(status, den)
KFM.den$source <- "KFM"
KFM.den <- KFM.den[complete.cases(KFM.den), ]
# FIXED: Replaced positional column indices 5:6 with explicit column names
KFM.den.long <- gather(KFM.den, status, value, mpa, reference)
KFM.den.long$resp <- "Den"

kfm.edit.bio <- subset(kfm.ave.ave, y != "Macrocystis pyrifera" | area != 2)
kfm.edit.bio <- subset(kfm.edit.bio, y != "Macrocystis pyrifera" | area != 1)

kfm.edit.bio <- kfm.edit.bio[, c("status", "CA_MPA_Name_Short", "MPA_Start", "year", "biomass", "y")]
kfm.edit.bio <- kfm.edit.bio[complete.cases(kfm.edit.bio), ]
KFM.bio <- kfm.edit.bio %>%
  spread(status, biomass)
KFM.bio$source <- "KFM"
KFM.bio <- KFM.bio[complete.cases(KFM.bio), ]
# FIXED: Replaced positional column indices 5:6 with explicit column names
KFM.bio.long <- gather(KFM.bio, status, value, mpa, reference)
KFM.bio.long$resp <- "Bio"

KFM.resp <- rbind(KFM.den.long, KFM.bio.long)

# ===========================================================================
# Section 12: Assign time and Before/After labels to KFM.join.ave
# ===========================================================================

# Assign time since MPA
KFM.join.ave <- assign_time_from_site_table(KFM.join.ave, Sites2)

# Merge with site metadata (type, location, hectares)
KFM.join.ave <- merge(KFM.join.ave, sites.short.edit,
                        by.x = c("CA_MPA_Name_Short"),
                        by.y = c("CA_MPA_Name_Short"),
                        all.x = TRUE)
KFM.join.ave <- KFM.join.ave[complete.cases(KFM.join.ave$year), ]
KFM.join.ave$source <- "KFM"

# Assign Before/After labels using utility function
KFM.join.ave <- assign_ba_from_site_table(KFM.join.ave, Sites2)

# ===========================================================================
# Section 13: Process KFM fish data (sheephead from MBON fish file)
# ===========================================================================
# KFM sheephead (Semicossyphus pulcher) data comes from a separate MBON fish file.
# Unlike PISCO, KFM does NOT record individual fish sizes, only counts per
# transect, so the KFM sheephead pipeline produces density only. Biomass is
# unavailable for KFM sheephead and is excluded from biomass-side analyses.
# KFM uses two fish survey methods (per Davis et al. 1997 KFM Handbook
# Vol 1, p.45-48):
#   - "visualfish" (Visual Fish Transect): since 1985. Transect dimensions
#     changed in 1996 from 3 m x 2 m x 100 m to 3 m x 2 m x 50 m (four
#     transects per site throughout; transects 1+2 from 1996+ sum to old
#     transect 1, transects 3+4 sum to old transect 2 - comparable across
#     the change).
#   - "rdfc" (Roving Diver Fish Count): implemented in 1996
#     (KFM_RDFC_SURVEY_START_YEAR). 30-minute roving counts, abundance
#     scored on a 10-point time-of-encounter scale.
# We process both methods separately and combine their response ratios.

mbon.fish <- read.csv(here::here("data", "MBON", "SBCMBON_kelp_forest_integrated_fish_20231022.csv"))

# Validate expected columns
required_fish_cols <- c("data_source", "site_id", "proj_taxon_id", "date",
                        "sample_method", "count", "area")
missing_fish_cols <- setdiff(required_fish_cols, names(mbon.fish))
if (length(missing_fish_cols) > 0) {
  stop("MBON fish CSV is missing required columns: ",
       paste(missing_fish_cols, collapse = ", "))
}

mbon.fish <- subset(mbon.fish, data_source == "kfm")
cat("  KFM fish data loaded:", nrow(mbon.fish), "rows\n")

# Subset sheephead (Semicossyphus pulcher) using KFM project taxon IDs:
#   t-k-127 = S. pulcher from visual fish counts
#   t-k-233 = S. pulcher from roving diver fish counts
mbon.fish.spul <- subset(mbon.fish, proj_taxon_id == "t-k-127" | proj_taxon_id == "t-k-233")

# Date conversion
tmpDateFormat <- "%Y-%m-%d"
tmp1date <- as.Date(mbon.fish.spul$date, format = tmpDateFormat)
if (length(tmp1date) == length(tmp1date[!is.na(tmp1date)])) {
  mbon.fish.spul$date <- tmp1date
} else {
  stop("Date conversion failed for mbon.fish.spul$date. ",
       sum(is.na(tmp1date)), " of ", length(tmp1date), " dates failed.")
}
rm(tmpDateFormat, tmp1date)
mbon.fish.spul$year <- year(mbon.fish.spul$date)

# Subset relevant survey methods:
#   - RDFC (roving diver fish count) only from 2003 onward when protocol was standardized
#     (KFM_RDFC_SURVEY_START_YEAR from 00c_analysis_constants.R)
#   - Visual fish surveys included for all years
mbon.fish.spul.sub <- subset(mbon.fish.spul,
                               (sample_method == "rdfc" & year >= KFM_RDFC_SURVEY_START_YEAR) |
                                 (sample_method == "visualfish"))

# Join to site table
mbon.spul.site <- validated_merge(mbon.fish.spul.sub, Sites2, by = c("site_id"), warn_threshold = 0)
mbon.spul.site$count <- as.numeric(mbon.spul.site$count)
mbon.spul.site <- subset(mbon.spul.site, CA_MPA_Name_Short != "")
mbon.spul.site$den <- mbon.spul.site$count / mbon.spul.site$area

# Calculate sheephead density per MPA/Reference - averaging hierarchy:
# replicate -> transect -> site -> MPA level
mbon.spul.sum <- mbon.spul.site %>%
  dplyr::group_by(site_id, sample_method, transect_id, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>%
  dplyr::summarise_at(c("den"), mean)

mbon.spul.ave <- mbon.spul.sum %>%
  dplyr::group_by(site_id, sample_method, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>%
  dplyr::summarise_at(c("den"), mean)

mbon.spul.max <- mbon.spul.ave %>%
  dplyr::group_by(status, sample_method, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>%
  dplyr::summarise_at(c("den"), mean)

# Remove old Anacapa SMR (pre-MLPA, listed in EXCLUDED_MPAS)
mbon.spul.max <- subset(mbon.spul.max, CA_MPA_Name_Short != "Anacapa Island SMR 1978")
colnames(mbon.spul.max)[colnames(mbon.spul.max) == "den"] <- "count"

mbon.spul.max$Prop <- 0
mbon.spul.max$PropCorr <- 0
mbon.spul.max$taxon_name <- "SPUL"

mbon.spul.max.rd <- subset(mbon.spul.max, sample_method == "rdfc")
mbon.spul.max.vf <- subset(mbon.spul.max, sample_method == "visualfish")

# -- KFM fish raw density dataframe (visual fish method) --
mbon.fish.edit <- mbon.spul.max.vf[, c("status", "CA_MPA_Name_Short", "MPA_Start", "year", "count", "taxon_name")]
KFM.fish.den <- mbon.fish.edit %>%
  spread(status, count)
KFM.fish.den$source <- "KFM"
KFM.fish.den.sub <- KFM.fish.den[complete.cases(KFM.fish.den), ]
# FIXED: Replaced positional column indices 5:6 with explicit column names
KFM.fish.den.long <- gather(KFM.fish.den.sub, status, value, mpa, reference)
KFM.fish.den.long$resp <- "Den"

# -- Proportions for RDFC method --
# Use standardized calculate_proportions() function from 01_utils.R
# Add 'y' column required by the function (taxon_name already set to "SPUL")
# NOTE: Using adaptive correction (min/2) for consistency across all taxa/programs.
# This is statistically preferred over fixed +0.01 (Aitchison 1986) as it scales
# appropriately with the data and avoids inflating effect sizes for rare species.
mbon.spul.max.rd$y <- mbon.spul.max.rd$taxon_name
mbon.spul.max.rd <- calculate_proportions(mbon.spul.max.rd, "count")

All.den.spul.kfm.sub.rd <- mbon.spul.max.rd[, c("CA_MPA_Name_Short", "sample_method", "year", "taxon_name", "status", "PropCorr")]
Short.den.spul.kfm.sub.rd <- All.den.spul.kfm.sub.rd %>%
  spread(status, PropCorr)
Short.den.spul.kfm.sub.rd.diff <- calculate_log_response_ratio(Short.den.spul.kfm.sub.rd)
names(Short.den.spul.kfm.sub.rd.diff)[names(Short.den.spul.kfm.sub.rd.diff) == "taxon_name"] <- "y"
Short.den.spul.kfm.sub.rd.diff$resp <- "RD"

# -- Proportions for visual fish method --
# Use standardized calculate_proportions() function from 01_utils.R
# NOTE: Using adaptive correction (min/2) for consistency - see RDFC method comment above.
mbon.spul.max.vf$y <- mbon.spul.max.vf$taxon_name
mbon.spul.max.vf <- calculate_proportions(mbon.spul.max.vf, "count")

All.den.spul.kfm.sub.vf <- mbon.spul.max.vf[, c("CA_MPA_Name_Short", "sample_method", "year", "taxon_name", "status", "PropCorr")]
Short.den.spul.kfm.sub.vf <- All.den.spul.kfm.sub.vf %>%
  spread(status, PropCorr)
Short.den.spul.kfm.sub.vf.diff <- calculate_log_response_ratio(Short.den.spul.kfm.sub.vf)
names(Short.den.spul.kfm.sub.vf.diff)[names(Short.den.spul.kfm.sub.vf.diff) == "taxon_name"] <- "y"
Short.den.spul.kfm.sub.vf.diff$resp <- "Den"

# Combine RDFC and visual fish response ratios
Short.den.spul.kfm.sub.diff <- rbind(Short.den.spul.kfm.sub.vf.diff, Short.den.spul.kfm.sub.rd.diff)

# Assign time since MPA
Short.den.spul.kfm.sub.diff <- assign_time_from_site_table(Short.den.spul.kfm.sub.diff, Sites2)

# Merge with site metadata
kfm.fish <- merge(Short.den.spul.kfm.sub.diff, sites.short.edit,
                    by.x = c("CA_MPA_Name_Short"),
                    by.y = c("CA_MPA_Name_Short"),
                    all.x = TRUE)
kfm.fish <- kfm.fish[complete.cases(kfm.fish$year), ]
kfm.fish$source <- "KFM"

# Assign Before/After labels
kfm.fish <- assign_ba_from_site_table(kfm.fish, Sites2)

cat("  KFM fish response ratios:", nrow(kfm.fish), "rows\n")

# ===========================================================================
# Section 14: Final outputs summary
# ===========================================================================
#   KFM.join.ave  - Response ratio data for urchins, macro (density + biomass)
#   KFM.resp      - Raw density/biomass data inside/outside MPAs
#   kfm.fish      - Sheephead fish response ratio data
#   KFM.fish.den.long - Raw fish density data
#   Sites2        - Site geolocation table (used by downstream scripts)
#   lter          - LTER subset of MBON (used by 06_lter_processing.R)
#   kfm.stipe.site, kfm.stipe.OG - Stipe data (used downstream)

cat("KFM/MBON processing complete.\n")
cat("  Outputs: KFM.join.ave (", nrow(KFM.join.ave), " rows), ",
    "KFM.resp (", nrow(KFM.resp), " rows), ",
    "kfm.fish (", nrow(kfm.fish), " rows)\n", sep = "")
