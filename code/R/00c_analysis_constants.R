# =============================================================================
# 00c_analysis_constants.R
# =============================================================================
#
# PURPOSE:
#   Define named constants for all hardcoded values used across the analysis
#   pipeline. Centralizing these values keeps assumptions explicit and makes
#   maintenance easier.
#
# USAGE:
#   source(here::here("code", "R", "00c_analysis_constants.R"))
#
# NOTE ON SHARED CONSTANTS:
#   This file is duplicated in both the data-processing repo
#   (Donham-Stier-CA-MPA-Data-2026) and the analysis repo
#   (sbc-2026-donham-kelp-mpa-cascade). Constants used by both repos are marked
#   with "SHARED" comments. Changes to shared constants must be mirrored
#   in both copies.
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================


# ===========================================================================
# Survey Area Standards (m^2 per transect by program)
# ===========================================================================
# Each monitoring program uses different sampling units. These constants
# convert raw counts to per-square-meter densities.

#' PISCO swath transect area (2m wide x 30m long = 60 m^2)
#' Used in: 04_pisco_processing.R (density normalization for PISCO swath counts)
#' Rationale: Standard PISCO swath protocol; see Carr et al. (2004) Ecological Applications
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
PISCO_SWATH_AREA_M2 <- 60

#' KFM quadrat area (default = 2 m^2 for urchins from 1985 onward)
#' Used in: 05_kfm_processing.R (legacy default; current code uses row-level `area`)
#' Rationale: KFM quadrat protocol changed through time:
#'   - Urchins: 1 m^2 in 1982-1984, then 2 m^2 from 1985 onward
#'   - Macrocystis: 1 m^2 (1982-1984), 2 m^2 (1985+), 10 m^2 (1996+); analysis
#'     uses the 10 m^2 subset only (filtered in 05_kfm_processing.R Section 11)
#' Density normalization in 05_kfm_processing.R now divides by the row-level
#' `area` column rather than this constant, so quad-size changes are handled
#' correctly. Constant retained for documentation/back-compat.
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
KFM_QUAD_AREA_M2 <- 2

#' LTER lobster plot area (300 m^2)
#' Used in: 06_lter_processing.R (density normalization for LTER lobster transects)
#' Rationale: SBC LTER lobster visual transect protocol (60m x 5m benthic area)
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
LTER_LOBSTER_PLOT_AREA_M2 <- 300

#' LTER fish survey transect area (80 m^2) for the "BIG" fish protocol
#' Used in: 06_lter_processing.R (density normalization for LTER sheephead)
#' Rationale: SBC LTER fish survey protocol uses TWO transect sizes:
#'   - 80 m^2 (40 m x 2 m) for BIG mobile fish (sheephead and other
#'     large-bodied species, 60,304 obs in Annual_fish_comb_20240307.csv)
#'   - 20 m^2 for CRYPTIC small benthic fish (gobies, blennies, sculpins,
#'     etc.: 32,877 obs)
#'   - 40 m^2 for an intermediate subset (630 obs)
#' Sheephead (the only LTER fish species used in this analysis) is uniformly
#' surveyed at 80 m^2, verified against AREA column on 2026-05-04. A
#' defensive check in 06_lter_processing.R Section 6 will warn if a future
#' LTER data release ever changes this.
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
LTER_FISH_SURVEY_AREA_M2 <- 80


# ===========================================================================
# Size Thresholds (minimum sizes for adult classification)
# ===========================================================================
# Minimum sizes for including organisms in analyses, based on field protocols.

#' PISCO minimum urchin size (mm test diameter)
#' Used in: 04_pisco_processing.R (size filtering for PISCO urchin records)
#' Rationale: PISCO field protocols only record urchins >= 25mm because
#' smaller urchins are cryptic and difficult to identify consistently in
#' the field. Applying a matching threshold ensures comparable counts
#' across programs.
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
PISCO_URCHIN_MIN_SIZE_MM <- 25


# ===========================================================================
# Temporal Parameters (year ranges and time-based filters)
# ===========================================================================
# Year ranges and time-based filters for analysis.

#' Survey season months (May-October)
#' Used in: 04_pisco_processing.R (seasonal filtering of PISCO surveys)
#' Rationale: Restricting to summer months (May-Oct) ensures consistency
#' across survey programs, which vary in winter effort. Most kelp forest
#' surveys are conducted during calm-water summer months.
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
SURVEY_MONTH_START <- 5   # May
SURVEY_MONTH_END <- 10    # October

#' Year lobster size measurements began by campus
#' Used in: 04_pisco_processing.R (filtering lobster records to years with size data)
#' Rationale: Carapace length measurements started at different times at each
#' PISCO campus. Before these years, only counts were recorded, so biomass
#' estimation via bootstrap is not possible.
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
UCSB_LOBSTER_SIZE_START_YEAR <- 2010
VRG_LOBSTER_SIZE_START_YEAR <- 2011

#' KFM RDFC (roving diver fish count) survey methodology start year
#' Used in: 05_kfm_processing.R (filtering KFM fish counts to standardized protocol)
#' Rationale: KFM implemented the Roving Diver Fish Count protocol in 1996
#' (Davis et al. 1997, Kelp Forest Monitoring Handbook Vol 1, p.48). The
#' earlier Visual Fish Transect protocol started in 1985 with a transect-
#' length change in 1996 (3 m x 2 m x 100 m -> 3 m x 2 m x 50 m). The two
#' methods coexist post-1996. We previously had this constant set to 2003
#' (incorrect); corrected 2026-05-04 after a metadata audit against the
#' NPS protocol manual at https://irma.nps.gov/DataStore/DownloadFile/485444.
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
KFM_RDFC_SURVEY_START_YEAR <- 1996

#' Minimum years of data required for site inclusion
#' Used in: 07_combine_and_export.R (filtering sites with insufficient data)
#' Rationale: Sites with fewer than 5 years of data provide unreliable
#' trend estimates for the pBACIPS methodology; this threshold balances
#' inclusion (more MPAs) against statistical power.
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
MINIMUM_YEARS_OF_DATA <- 5

#' Standardized time point for effect size extraction (years post-MPA implementation)
#' Used in: 08_effect_sizes.R (analysis repo) for NLS effect size prediction
#' Rationale: Using a fixed time point (t=11) rather than the maximum observed time
#' allows for comparable effect sizes across MPAs with different establishment dates.
#' We chose 11 years because it corresponds to the age of our youngest MPA in 2023
#' (MLPA South Coast MPAs implemented in 2012). This controls for differences in
#' effect size that could arise simply from longer protection duration.
#' See: Thiault et al. (2017) Methods in Ecology and Evolution for pBACIPS methodology.
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
EFFECT_SIZE_TIME_YEARS <- 11


# ===========================================================================
# Random Seeds (for reproducible bootstrap sampling)
# ===========================================================================
# Parameters for bootstrap resampling of biomass estimates.
# Seeds use sequential numbering (12345-12349) for easy identification.

#' Number of bootstrap resamples for biomass estimation
#' Used in: 04_pisco_processing.R, 05_kfm_processing.R, 06_lter_processing.R
#' Rationale: 1000 resamples balances precision and computation time;
#' bootstrap SEs stabilize well before 1000 iterations.
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
N_BOOTSTRAP_RESAMPLES <- 1000

#' Random seed for PISCO lobster bootstrap biomass estimation
#' Used in: 04_pisco_processing.R (set.seed before lobster bootstrap loop)
SEED_PISCO_LOBSTER <- 12345

#' Random seed for PISCO urchin bootstrap biomass estimation
#' Used in: 04_pisco_processing.R (set.seed before urchin bootstrap loop)
SEED_PISCO_URCHIN <- 12346

#' Random seed for KFM urchin bootstrap biomass estimation
#' Used in: 05_kfm_processing.R (set.seed before urchin bootstrap loop)
SEED_KFM_URCHIN <- 12347

#' Random seed for KFM Macrocystis bootstrap biomass estimation
#' Used in: 05_kfm_processing.R (set.seed before kelp bootstrap loop)
SEED_KFM_MACRO <- 12348

#' Random seed for LTER urchin bootstrap biomass estimation
#' Used in: 06_lter_processing.R (set.seed before urchin bootstrap loop)
SEED_LTER_URCHIN <- 12349


# ===========================================================================
# Statistical Thresholds (significance levels and diagnostic cutoffs)
# ===========================================================================
# Significance levels and diagnostic thresholds used primarily in the
# analysis repo (08_effect_sizes.R, 09_meta_analysis.R). Defined here for
# consistency across both repos.

#' Standard significance level (alpha = 0.05)
#' Used in: 12_results_summary.R (analysis repo), reporting significant effects
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
SIGNIFICANCE_ALPHA <- 0.05

#' Z-critical value for 95% confidence intervals
#' Used in: 08_effect_sizes.R (analysis repo), constructing CIs for effect sizes
#' Rationale: qnorm(0.975) = 1.959964... rounded to 1.96 by convention
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
Z_CRITICAL_95_CI <- 1.96

#' Shapiro-Wilk test sample size bounds
#' Used in: 08_effect_sizes.R (analysis repo), NLS residual normality testing
#' Rationale: R's shapiro.test() requires 3-5000 observations
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
MIN_OBS_SHAPIRO_TEST <- 3
MAX_OBS_SHAPIRO_TEST <- 5000

#' DHARMa simulation count for residual diagnostics
#' Used in: 08_effect_sizes.R (analysis repo), model diagnostics via DHARMa
#' Rationale: 250 simulations provides adequate residual diagnostic power
#' while keeping computation tractable for many models
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
DHARMA_N_SIMULATIONS <- 250

#' Heteroscedasticity correlation thresholds
#' Used in: 08_effect_sizes.R (analysis repo), flagging models with non-constant variance
#' Rationale: LM (linear models) are more strict because they assume homoscedasticity;
#' NLS (non-linear models) get a more lenient threshold due to inherent complexity
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
MAX_HETEROSCEDASTICITY_CORRELATION_LM <- 0.5
MAX_HETEROSCEDASTICITY_CORRELATION_NLS <- 0.6

#' Maximum outliers allowed in NLS diagnostics
#' Used in: 08_effect_sizes.R (analysis repo), flagging NLS models with excessive outliers
#' Rationale: A single outlier is tolerated; more than 1 suggests model misspecification
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
MAX_OUTLIERS_NLS <- 1

#' NLS normality threshold (more lenient than default 0.05)
#' Used in: 08_effect_sizes.R (analysis repo), Shapiro-Wilk p-value cutoff for NLS residuals
#' Rationale: NLS residuals are expected to be somewhat non-normal; using 0.01
#' instead of 0.05 avoids excessive false-positive flagging
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
NLS_NORMALITY_THRESHOLD <- 0.01

#' Cook's distance threshold numerator for meta-analysis outlier detection
#' Used in: 09_meta_analysis.R, 13_additional_analyses.R (analysis repo)
#' Rationale: Standard rule. Outliers have Cook's D > 4/n (Bollen & Jackman 1990)
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
COOKS_DISTANCE_NUMERATOR <- 4


# ===========================================================================
# Allometric Coefficients (length-weight relationships)
# ===========================================================================
# Length-weight relationships from published literature for converting
# size measurements to biomass.

#' California sheephead (Semicossyphus pulcher) length-weight relationship
#' biomass = a * length^b (where length is in cm and biomass in grams)
#' Used in: 06_lter_processing.R (converting LTER sheephead length to biomass)
#' Source: Cowen (1990) Environmental Biology of Fishes
#' SHARED: Also defined in analysis repo 00c_analysis_constants.R
SHEEPHEAD_BIOMASS_ALLOMETRIC_A <- 0.0144
SHEEPHEAD_BIOMASS_ALLOMETRIC_B <- 3.04


# ===========================================================================
# Trophic Assignments (species-to-trophic-level mapping)
# ===========================================================================
# Maps species names (both full and abbreviated forms) to trophic levels.
# Used by temporal analysis (10_temporal_analysis.R), figures (11_figures.R),
# and other scripts in the analysis repo.
#
# Trophic cascade order (top predators first):
#   Predators: P. interruptus (spiny lobster), S. pulcher (sheephead)
#   Urchins:   S. purpuratus (purple urchin), M. franciscanus (red urchin)
#   Kelp:      M. pyrifera (giant kelp)
#
# SHARED: Also defined in analysis repo 00c_analysis_constants.R

trophic_assignment <- c(
  "Panulirus interruptus" = "Predators",
  "Semicossyphus pulcher" = "Predators",
  "P. interruptus" = "Predators",
  "S. pulcher" = "Predators",
  "Strongylocentrotus purpuratus" = "Urchins",
  "Mesocentrotus franciscanus" = "Urchins",
  "S. purpuratus" = "Urchins",
  "M. franciscanus" = "Urchins",
  "Macrocystis pyrifera" = "Kelp",
  "M. pyrifera" = "Kelp"
)


# =============================================================================
# Confirmation message
# =============================================================================
cat("Analysis constants loaded.\n")
