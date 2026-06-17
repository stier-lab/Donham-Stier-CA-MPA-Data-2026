# =============================================================================
# 01_utils.R  (Data Processing Repo)
# =============================================================================
#
# PURPOSE:
#   Utility functions and constants for data processing (scripts 03-07).
#   This is a subset of the analysis repo's utils. Only functions needed
#   for data import, biomass conversion, bootstrap resampling, response
#   ratio calculation, and site exclusion logic.
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest. Data Processing
# =============================================================================


# =============================================================================
# SECTION 1: TIME AND PERIOD ASSIGNMENT FUNCTIONS
# =============================================================================

#' @title Calculate time since MPA implementation
#' @description Computes the number of years elapsed since an MPA was
#'   implemented. Years at or before implementation return 0 (the "Before"
#'   period). Years after implementation return years - mpa_start.
#' @param years Numeric vector of observation years.
#' @param mpa_start Numeric scalar; the year the MPA was implemented.
#'   Convention: the implementation year itself is treated as "Before" (time=0).
#' @return Numeric vector of the same length as \code{years}, with non-negative
#'   integer values representing years since MPA establishment.
#' @examples
#'   calculate_time_since_mpa(c(2008, 2012, 2015), mpa_start = 2012)
#'   # Returns: c(0, 0, 3)
calculate_time_since_mpa <- function(years, mpa_start) {
  ifelse(years <= mpa_start, 0, years - mpa_start)
}

#' @title Assign Before/After period labels
#' @description Labels each observation year as "Before" or "After" relative to
#'   MPA implementation. The implementation year itself is classified as
#'   "Before" (consistent with pBACIPS convention).
#' @param years Numeric vector of observation years.
#' @param mpa_start Numeric scalar; the year the MPA was implemented.
#' @return Character vector of "Before" or "After" labels, same length as
#'   \code{years}.
#' @examples
#'   assign_ba_period(c(2008, 2012, 2015), mpa_start = 2012)
#'   # Returns: c("Before", "Before", "After")
assign_ba_period <- function(years, mpa_start) {
  ifelse(years <= mpa_start, "Before", "After")
}

#' @title Assign time since MPA to a dataframe via site lookup table
#' @description Iterates over unique MPAs in the dataframe, looks up each MPA's
#'   implementation year from the site table, and applies
#'   \code{calculate_time_since_mpa()} to compute the \code{time} column.
#'   MPAs not found in the site table default to time=0 with a warning.
#' @param df Data frame containing an MPA identifier column and a \code{year}
#'   column.
#' @param site_table Data frame with MPA metadata, must contain a column
#'   matching \code{mpa_col} and a \code{MPA_Start} column.
#' @param mpa_col Character; name of the MPA identifier column in both
#'   \code{df} and \code{site_table}. Defaults to \code{"CA_MPA_Name_Short"}.
#' @return The input data frame with a \code{time} column added or overwritten.
assign_time_from_site_table <- function(df, site_table, mpa_col = "CA_MPA_Name_Short") {
  df$time <- 0
  mpas <- unique(df[[mpa_col]])
  matched_mpas <- character(0)

  for (mpa in mpas) {
    idx <- which(df[[mpa_col]] == mpa)
    j <- which(site_table[[mpa_col]] == mpa)
    if (length(j) > 0) {
      mpa_start <- site_table$MPA_Start[j[1]]
      df$time[idx] <- calculate_time_since_mpa(df$year[idx], mpa_start)
      matched_mpas <- c(matched_mpas, mpa)
    }
  }

  unmatched <- setdiff(mpas, matched_mpas)
  if (length(unmatched) > 0) {
    warning("assign_time_from_site_table: ", length(unmatched),
            " MPA(s) not found in site table (time defaults to 0): ",
            paste(unmatched, collapse = ", "), call. = FALSE)
  }
  df
}

#' @title Assign Before/After labels to a dataframe via site lookup table
#' @description Iterates over unique MPAs in the dataframe, looks up each MPA's
#'   implementation year from the site table, and applies
#'   \code{assign_ba_period()} to compute the \code{BA} column.
#'   MPAs not found in the site table default to BA=NA with a warning.
#' @param df Data frame containing an MPA identifier column and a \code{year}
#'   column.
#' @param site_table Data frame with MPA metadata, must contain a column
#'   matching \code{mpa_col} and a \code{MPA_Start} column.
#' @param mpa_col Character; name of the MPA identifier column in both
#'   \code{df} and \code{site_table}. Defaults to \code{"CA_MPA_Name_Short"}.
#' @return The input data frame with a \code{BA} column added or overwritten.
assign_ba_from_site_table <- function(df, site_table, mpa_col = "CA_MPA_Name_Short") {
  df$BA <- NA
  mpas <- unique(df[[mpa_col]])

  for (mpa in mpas) {
    idx <- which(df[[mpa_col]] == mpa)
    j <- which(site_table[[mpa_col]] == mpa)
    if (length(j) > 0) {
      mpa_start <- site_table$MPA_Start[j[1]]
      df$BA[idx] <- assign_ba_period(df$year[idx], mpa_start)
    }
  }

  unmatched <- unique(df[[mpa_col]][is.na(df$BA)])
  if (length(unmatched) > 0) {
    warning("assign_ba_from_site_table: ", length(unmatched),
            " MPA(s) not found in site table (defaulting to BA=NA): ",
            paste(unmatched, collapse = ", "), call. = FALSE)
  }
  df
}


# =============================================================================
# SECTION 2: BIOMASS CONVERSION FUNCTIONS
# =============================================================================
# Allometric length-weight (or stipe-weight) relationships used to convert
# size-frequency data into biomass estimates during bootstrap resampling.
# Each function takes a size measurement and returns biomass in grams.
# Coefficients are from published allometric relationships for each species.

# NOTE ON ALLOMETRIC SOURCES (per Emily Donham, 2026-04 review):
#   The length-weight / size-biomass relationships used below come from the
#   SBC LTER's calibration work (compiled and provided by the LTER team), not
#   from the originally cited textbooks/papers. The original textbook
#   citations (Barsky 2001, Ebert 2010, Cowen 1990) are retained here only as
#   the published basis the LTER calibrations build on.

#' @title Lobster biomass from carapace length
#' @description Allometric conversion of California spiny lobster
#'   (\emph{Panulirus interruptus}) carapace length to biomass (grams)
#'   using a power-law relationship: W = a * L^b.
#' @param size Numeric vector of carapace lengths (mm).
#'   NOTE: PISCO raw data stores carapace length in cm; the conversion to mm
#'   happens upstream in 04_pisco_processing.R (size * 10).
#' @return Numeric vector of biomass estimates (grams).
#' @source SBC LTER allometric calibration (cf. Barsky 2001, CDFG Lobster Assessment).
bio_lobster <- function(size) {
  0.001352821 * (size) ^ 2.913963113
}

#' @title Red urchin biomass from test diameter
#' @description Allometric conversion of red sea urchin
#'   (\emph{Mesocentrotus franciscanus}) test diameter to biomass (grams)
#'   using a power-law relationship: W = a * D^b.
#' @param size Numeric vector of test diameters (mm).
#' @return Numeric vector of biomass estimates (grams).
#' @source SBC LTER allometric calibration (cf. Ebert 2010, Sea Urchins).
bio_redurch <- function(size) {
  0.00059 * (size) ^ 2.917
}

#' @title Purple urchin biomass from test diameter
#' @description Allometric conversion of purple sea urchin
#'   (\emph{Strongylocentrotus purpuratus}) test diameter to biomass (grams)
#'   using a power-law relationship: W = a * D^b.
#' @param size Numeric vector of test diameters (mm).
#' @return Numeric vector of biomass estimates (grams).
#' @source SBC LTER allometric calibration (cf. Ebert 2010, Sea Urchins).
bio_purpurch <- function(size) {
  0.00059 * (size) ^ 2.870
}

#' @title Sheephead biomass from total length
#' @description Allometric conversion of California sheephead
#'   (\emph{Semicossyphus pulcher}) total length to biomass (grams)
#'   using a power-law relationship: W = a * L^b.
#' @param size_cm Numeric vector of total lengths (cm).
#' @return Numeric vector of biomass estimates (grams).
#' @source SBC LTER allometric calibration (cf. Cowen 1990, Env. Biol. Fish.).
bio_sheephead <- function(size_cm) {
  SHEEPHEAD_BIOMASS_ALLOMETRIC_A * (size_cm) ^ SHEEPHEAD_BIOMASS_ALLOMETRIC_B
}

# =============================================================================
# Macrocystis stipe-to-biomass slope
# =============================================================================
# DECISION (2026-06, Emily review): use Reed et al. (2009)'s published pooled
# stipe-density-to-standing-crop slope of 0.08.
#
#   Reed, D.C., Rassweiler, A. & Arkema, K.K. (2009). Density derived estimates
#   of standing crop and net primary production in the giant kelp Macrocystis
#   pyrifera. Marine Biology 156: 2077-2083. doi:10.1007/s00227-009-1238-6
#   (open access, PMC3873067). Fig. 2b: standing crop = 0.08 x stipe density +
#   0.01 (r^2 = 0.79), an SBC-LTER calibration (Jan 2003 - Dec 2008).
#
# This replaces an earlier value, mean(c(0.10386, 0.10103, 0.09267, 0.09204,
# 0.08054, 0.08505)) = 0.0925, carried over from the archived pipeline
# (pBACIPS_PISCO_V10.R L674) where it was labelled "slopes May-Oct from LTER"
# but could not be traced to any published source (audited 2026-05-05 against
# Reed 2009 and Rassweiler et al. 2018, Ecology 99(9):2132 — the six values
# appear in neither). Reed 2009's pooled 0.08 is the canonical, open-access,
# already-cited SBC LTER calibration, so it is the defensible project standard.
#
# RESULT IMPACT: the analysis uses proportion-based lnRR, which is invariant
# under multiplicative rescaling of biomass, so meta-analysis effect sizes,
# CIs, and weights are UNCHANGED. Only the absolute Macrocystis biomass values
# in output/harmonized/harmonized_raw_responses.csv shift (~14% lower than the
# previous 0.0925).
# =============================================================================
MACRO_AVE_SLOPE <- 0.08  # Reed et al. 2009, Mar Biol 156:2077-2083

#' @title Giant kelp biomass from stipe count
#' @description Converts \emph{Macrocystis pyrifera} stipe counts to biomass
#'   (grams) using MACRO_AVE_SLOPE (Reed et al. 2009 pooled slope = 0.08;
#'   see note above the constant definition), multiplied by 1000 to convert
#'   kg to g.
#' @param stipe Numeric vector of stipe counts (or stipe density per m^2).
#'   PISCO records stipes per individual; KFM records per-plant stipe counts
#'   in the KFM_Macrocystis_RawData file. LTER stores its counts in a column
#'   named FRONDS, but per Li (LTER) these are the same structural unit
#'   (stipes = fronds for this purpose), so all three programs feed compatible
#'   counts into bio_macro().
#' @note MACRO_AVE_SLOPE = 0.08 (Reed et al. 2009). Because the analysis uses
#'   proportion-based lnRR, the choice of slope does not affect meta-analysis
#'   effect sizes; it only rescales absolute biomass.
#' @return Numeric vector of biomass estimates (grams).
bio_macro <- function(stipe) {
  stipe * MACRO_AVE_SLOPE * 1000
}


# =============================================================================
# SECTION 3: BOOTSTRAP BIOMASS ESTIMATION
# =============================================================================

#' @title Bootstrap biomass estimation from size-frequency data
#' @description Estimates total biomass and its standard error for a given
#'   observation by bootstrap resampling from a size-frequency distribution.
#'   For each resample, \code{count} individuals are drawn (with replacement)
#'   from the available size classes, converted to biomass via
#'   \code{biomass_fun}, and summed. The mean and SD across resamples give
#'   the biomass estimate and its SE.
#'
#'   VARIANCE PROPAGATION NOTE:
#'   The bootstrap SE returned here quantifies uncertainty in the biomass
#'   estimate arising from unknown individual sizes (i.e., we know how many
#'   urchins/lobsters were counted, but not their exact sizes). This SE is
#'   stored in the cached bootstrap results but is NOT propagated through
#'   to the final lnRR variance used in the meta-analysis. Instead, the
#'   meta-analysis (in the analysis repo, 08_effect_sizes.R) derives lnRR
#'   variance analytically from the NLS/linear model fits (via delta method
#'   SEs for NLS models, or standard regression SEs for linear models).
#'   This means that uncertainty from the biomass bootstrap is not reflected
#'   in the meta-analytic confidence intervals or weights. In practice, this
#'   is a minor limitation because (a) the bootstrap SE is typically small
#'   relative to the temporal variance captured by the NLS/lmer models, and
#'   (b) the proportion-based lnRR partially absorbs this uncertainty since
#'   biomass estimation error affects both MPA and reference sites similarly.
#'
#' @param count Integer; number of individuals observed at this site-year.
#'   If NA, returns NA for biomass and SE.
#' @param size_freq_indices Integer vector; row indices into
#'   \code{size_freq_table} identifying the size-frequency records that match
#'   this observation (same site, species, time period).
#' @param size_freq_table Data frame with at least \code{size} (numeric) and
#'   \code{count} (integer) columns representing the size-frequency
#'   distribution to resample from.
#' @param biomass_fun Function; an allometric conversion function (e.g.,
#'   \code{bio_lobster}, \code{bio_redurch}) that takes a size vector and
#'   returns a biomass vector.
#' @param n_resamples Integer; number of bootstrap resamples (default 1000).
#' @return A named list with three elements:
#'   \describe{
#'     \item{biomass}{Numeric; mean total biomass across bootstrap resamples.}
#'     \item{se}{Numeric; standard deviation of total biomass across resamples
#'       (bootstrap SE). See VARIANCE PROPAGATION NOTE above. This SE is not
#'       carried through to the final meta-analysis variance.}
#'     \item{count}{Integer; the input count value.}
#'   }
bootstrap_biomass <- function(count, size_freq_indices, size_freq_table,
                              biomass_fun, n_resamples = 1000) {
  n <- count
  t2 <- size_freq_indices

  # Case 1: count is NA, cannot estimate biomass
  if (is.na(n)) {
    return(list(biomass = NA, se = NA, count = n))
  }

  # Case 2: count is zero, biomass is definitionally zero
  if (n == 0) {
    return(list(biomass = 0, se = 0, count = 0))
  }

  # Case 3: count > 0 but no size-frequency data available
  if (length(t2) == 0) {
    return(list(biomass = NA, se = NA, count = n))
  }

  # Case 4: count > 0 and size-frequency data exists, run bootstrap
  a <- numeric(0)
  for (j in seq_along(t2)) {
    reps <- rep(size_freq_table$size[t2[j]], size_freq_table$count[t2[j]])
    a <- c(a, reps)
  }
  a <- as.numeric(na.omit(a))

  if (length(a) == 0) {
    return(list(biomass = NA, se = NA, count = n))
  }

  s <- matrix(NA, nrow = n_resamples, ncol = n)
  for (k in 1:n_resamples) {
    s[k, ] <- sample(a, n, replace = TRUE)
  }

  # Direct vectorized call. biomass_fun is already vectorized, so applying
  # it to the full matrix is equivalent to apply(s, c(1,2), biomass_fun)
  # but substantially faster.
  s_bio <- biomass_fun(s)
  bootstrap_sums <- rowSums(s_bio)
  ave_biomass <- mean(bootstrap_sums)
  se_biomass <- sd(bootstrap_sums)

  return(list(biomass = ave_biomass, se = se_biomass, count = n))
}


# =============================================================================
# SECTION 4: RESPONSE RATIO CALCULATIONS
# =============================================================================

#' @title Calculate proportions with zero-correction for log response ratios
#' @description For each MPA-taxon combination, computes the proportion of the
#'   observed value relative to the maximum value across all years. Then applies
#'   a zero-correction to avoid log(0) when computing log response ratios.
#'   Three correction methods are supported:
#'   \describe{
#'     \item{adaptive}{Adds half the minimum non-zero proportion (default).}
#'     \item{fixed}{Adds a constant 0.01.}
#'     \item{none}{No correction applied.}
#'   }
#'
#' METHODOLOGICAL NOTE: PROPORTION-BASED lnRR:
#'   This project computes log response ratios on PROPORTIONS of time-series
#'   maxima rather than on raw density or biomass values. This is non-standard
#'   compared to typical meta-analysis practice (e.g., Hedges et al. 1999),
#'   where lnRR = ln(mean_treatment / mean_control) uses raw means.
#'
#'   Rationale for the proportion approach:
#'   (1) Standardization across monitoring programs. PISCO, KFM, and LTER
#'       use different transect areas, survey methods, and spatial extents.
#'       Raw density/biomass values are not directly comparable across
#'       programs. Proportions normalize each MPA-taxon time series to a
#'       common 0-1 scale, making cross-program effect sizes comparable.
#'   (2) Focus on relative change. The scientific question is whether
#'       MPA sites changed more (or less) relative to their own baseline
#'       maximum than reference sites. Proportions capture this relative
#'       change while removing absolute scale differences.
#'   (3) Precedent. Proportion-based approaches are used in BACI designs
#'       where absolute values differ across sites; see Thiault et al.
#'       (2017) Methods in Ecology and Evolution for the pBACIPS framework.
#'
#'   Implication: Because proportions are calculated within each
#'   MPA x taxon x status (mpa/reference) combination, any multiplicative
#'   bias in the underlying measurements (e.g., the stipe-vs-frond issue
#'   in LTER kelp data) cancels out in the proportion calculation, provided
#'   the bias is constant across years within each group.
#'
#'   See also: Manuscript methods section and calculate_log_response_ratio()
#'   in this file.
#'
#' @param df Data frame with columns \code{CA_MPA_Name_Short}, \code{y}
#'   (taxon), and a numeric value column specified by \code{value_col}.
#' @param value_col Character; name of the column containing raw values
#'   (density or biomass) to convert to proportions.
#' @param correction_method Character; one of \code{"adaptive"} (default),
#'   \code{"fixed"}, or \code{"none"}.
#' @return The input data frame with two new columns:
#'   \describe{
#'     \item{Prop}{Raw proportion (value / max value within MPA-taxon).}
#'     \item{PropCorr}{Zero-corrected proportion for safe log transformation.}
#'   }
calculate_proportions <- function(df, value_col, correction_method = "adaptive") {
  df$Prop <- NA_real_
  df$PropCorr <- NA_real_

  mpas <- unique(df$CA_MPA_Name_Short)
  taxa <- unique(df$y)

  for (mpa in mpas) {
    for (taxon in taxa) {
      idx <- which(df$CA_MPA_Name_Short == mpa & df$y == taxon)
      if (length(idx) > 0) {
        max_val <- max(df[[value_col]][idx], na.rm = TRUE)
        if (is.finite(max_val) && max_val > 0) {
          df$Prop[idx] <- df[[value_col]][idx] / max_val

          if (correction_method == "adaptive") {
            non_zero_props <- df$Prop[idx][df$Prop[idx] > 0]
            if (length(non_zero_props) > 0) {
              correction <- min(non_zero_props, na.rm = TRUE) / 2
            } else {
              warning("calculate_proportions(): MPA '", mpa, "' taxon '", taxon,
                      "' has all-zero proportions (no variation in data). ",
                      "Using fixed correction = 0.01.", call. = FALSE)
              correction <- 0.01
            }
            df$PropCorr[idx] <- df$Prop[idx] + correction
          } else if (correction_method == "fixed") {
            df$PropCorr[idx] <- df$Prop[idx] + 0.01
          } else if (correction_method == "none") {
            df$PropCorr[idx] <- df$Prop[idx]
          } else {
            warning("Unknown correction_method '", correction_method, "'. Using 'adaptive'.")
            non_zero_props <- df$Prop[idx][df$Prop[idx] > 0]
            correction <- ifelse(length(non_zero_props) > 0,
                                  min(non_zero_props, na.rm = TRUE) / 2, 0.01)
            df$PropCorr[idx] <- df$Prop[idx] + correction
          }
        } else {
          # max_val is non-finite or <= 0: all values are NA, zero, or negative
          warning("calculate_proportions(): MPA '", mpa, "' taxon '", taxon,
                  "' has max_val=", max_val, ", setting Prop and PropCorr to NA.",
                  call. = FALSE)
          df$Prop[idx] <- NA_real_
          df$PropCorr[idx] <- NA_real_
        }
      }
    }
  }
  df
}

#' @title Calculate log response ratio (lnRR)
#' @description Computes the log response ratio as ln(MPA / Reference) for
#'   each observation. Rows with missing MPA or reference values are dropped.
#'   Zero or negative ratios are set to NA to avoid undefined log values.
#'   Uses \code{safe_divide()} to handle division-by-zero gracefully.
#'
#'   NOTE: In this pipeline, the \code{mpa} and \code{reference} columns
#'   contain zero-corrected PROPORTIONS (PropCorr), not raw density or
#'   biomass values. Proportions are generated by \code{calculate_proportions()}
#'   upstream. See the methodological note in that function for rationale.
#'
#' @param df Data frame with columns \code{mpa} (numeric; MPA site value) and
#'   \code{reference} (numeric; reference site value). May contain additional
#'   columns which are preserved in the output.
#' @return The input data frame (minus incomplete rows) with two new columns:
#'   \describe{
#'     \item{Diff}{Numeric; the ratio mpa/reference.}
#'     \item{lnDiff}{Numeric; the natural log of Diff (i.e., lnRR).}
#'   }
calculate_log_response_ratio <- function(df) {
  if (!all(c("mpa", "reference") %in% names(df))) {
    stop("calculate_log_response_ratio(): expected columns 'mpa' and 'reference'.")
  }

  keep <- complete.cases(df[, c("mpa", "reference")])
  df <- df[keep, , drop = FALSE]

  df$Diff <- safe_divide(df$mpa, df$reference, replace_zero = TRUE, context = "mpa/reference")
  df$Diff[!is.finite(df$Diff) | df$Diff <= 0] <- NA_real_
  df$lnDiff <- log(df$Diff)
  df
}


# =============================================================================
# SECTION 5: SPECIES NAME STANDARDIZATION
# =============================================================================

#' @title Standardize species codes to full scientific names
#' @description Converts abbreviated PISCO-style species codes to full
#'   binomial scientific names. Names that don't match any known code are
#'   returned unchanged (they may already be full names from KFM/LTER).
#'
#'   Mapping:
#'   \itemize{
#'     \item PANINT -> Panulirus interruptus (California spiny lobster)
#'     \item SPUL -> Semicossyphus pulcher (California sheephead)
#'     \item MACPYRAD -> Macrocystis pyrifera (giant kelp)
#'     \item MESFRAAD -> Mesocentrotus franciscanus (red sea urchin)
#'     \item STRPURAD -> Strongylocentrotus purpuratus (purple sea urchin)
#'   }
#' @param names Character vector of species codes or names.
#' @return Character vector with codes replaced by full scientific names.
#' @examples
#'   standardize_species_names(c("PANINT", "SPUL", "Macrocystis pyrifera"))
#'   # Returns: c("Panulirus interruptus", "Semicossyphus pulcher", "Macrocystis pyrifera")
standardize_species_names <- function(names) {
  names[names == "PANINT"] <- "Panulirus interruptus"
  names[names == "SPUL"] <- "Semicossyphus pulcher"
  names[names == "MACPYRAD"] <- "Macrocystis pyrifera"
  names[names == "MESFRAAD"] <- "Mesocentrotus franciscanus"
  names[names == "STRPURAD"] <- "Strongylocentrotus purpuratus"
  names
}


# =============================================================================
# SECTION 6: SITE AND MPA EXCLUSION CONSTANTS
# =============================================================================
#
# These constants control which sites, MPAs, and KFM transects are included
# or excluded at various stages of the data processing pipeline. They are
# defined here in 01_utils.R so all downstream scripts (03-07) share the
# same exclusion criteria.
#
# EXCLUDED_REFERENCE_SITES
#   PISCO reference (control) site codes to exclude from the analysis.
#   These sites are dropped during PISCO processing (04_pisco_processing.R)
#   because they lack matching MPA-reference pairs or have data quality issues.
#
# EXCLUDED_MPAS
#   MPA names to exclude entirely from all analyses. Reasons include:
#   insufficient temporal coverage, no matching reference site, sites outside
#   the study region, or MPAs with pre-existing protection (e.g., Anacapa
#   Island SMR 1978). Used in 04_pisco_processing.R and 07_combine_and_export.R.
#
# EXCLUDED_KFM_SITES
#   KFM (Kelp Forest Monitoring) transect codes to exclude. These are
#   Anacapa Island transects that fall within the 1978 SMR boundary (pre-MLPA
#   protection) and would confound the MLPA-era before/after analysis.
#   Used in 05_kfm_processing.R.
#
# SHEEPHEAD_ONLY_MPAS
#   MPAs that contribute data only for sheephead (Semicossyphus pulcher),
#   typically from PISCO fish surveys. These are excluded from the main
#   cross-taxa analysis in 07_combine_and_export.R to prevent sheephead
#   from being over-represented, then selectively re-included for
#   sheephead-specific analyses via reinclude_sheephead_mpas().
#
# SHEEPHEAD_REINCLUSION_MPAS
#   Subset of SHEEPHEAD_ONLY_MPAS that are re-included specifically for
#   sheephead analysis. This is a curated list (not all sheephead-only
#   MPAs qualify). Santa Barbara Island SMR has special handling: only
#   PISCO sheephead data and all KFM data are re-included.
#   Used by reinclude_sheephead_mpas() in 07_combine_and_export.R.

EXCLUDED_REFERENCE_SITES <- c(
  "ANACAPA_BLACK_SEA_BASS",
  "SMI_BAY_POINT",
  "SCAI_SHIP_ROCK",
  "SMI_CROOK_POINT_E",
  "SMI_CROOK_POINT_W",
  "SBI_SUTIL",
  "SCI_YELLOWBANKS_CEN",
  "SCI_YELLOWBANKS_W",
  "SRI_FORD_POINT",
  "SMI_PRINCE_ISLAND_CEN",
  "SMI_PRINCE_ISLAND_N",
  "SMI_HARE_ROCK",
  "SMI_TYLER_BIGHT_E",
  "SMI_TYLER_BIGHT_W",
  "SBI_GRAVEYARD_CANYON",
  "SBI_GRAVEYARD_CANYON_N"
)

EXCLUDED_MPAS <- c(
  "Carrington Point SMR",
  "N/A",
  "Arrow Point to Lion Head Point SMCA",
  "Crystal Cove SMCA",
  "Laguna Beach SMR",
  "South La Jolla SMR",
  "Vandenberg SMR",
  "Point Conception SMR",
  "Anacapa Island SMR 1978",
  "Painted Cave SMCA",
  "Anacapa Island SMCA"
)

EXCLUDED_KFM_SITES <- c(
  "a-k-21", "a-k-05", "a-k-07", "a-k-26", "a-k-27", "a-k-28",
  "a-k-29", "a-k-30", "a-k-36", "a-k-37", "a-k-12", "a-k-13",
  "a-k-31", "a-k-35"
)

# NOTE: "Painted Cave SMCA" and "Anacapa Island SMCA" appear in both EXCLUDED_MPAS
# and SHEEPHEAD_ONLY_MPAS. This is intentional belt-and-suspenders: they are excluded
# during per-program processing (via EXCLUDED_MPAS), but the SHEEPHEAD_ONLY_MPAS
# filter in script 07 provides a safety net in case any rows survive. They are NOT
# in SHEEPHEAD_REINCLUSION_MPAS, so they stay permanently excluded.
SHEEPHEAD_ONLY_MPAS <- c(
  "Blue Cavern Onshore SMCA",
  "Painted Cave SMCA",
  "Dana Point SMCA",
  "Farnsworth Onshore SMCA",
  "Point Dume SMR",
  "Cat Harbor SMCA",
  "Swamis SMCA",
  "Anacapa Island SMCA",
  "Long Point SMR",
  "Point Dume SMCA",
  "Santa Barbara Island SMR"
)

# Curated subset of SHEEPHEAD_ONLY_MPAS re-included for sheephead-specific
# analyses. Inclusion criterion: regulations protect SHEEPHEAD specifically
# (i.e., no take of finfish at the SMCA), so the MPA functions like an SMR
# from sheephead's perspective.
#
# Removed 2026-04 (per Emily review):
#   - "Dana Point SMCA": permits take of lobster, finfish, and urchins,
#     so not a sheephead refuge in any meaningful sense.
#   - "Cat Harbor SMCA": permits commercial take of lobster and urchin and
#     recreational take of finfish (including sheephead).
#
# The take rules and reinclusion decision for every sheephead-only MPA are
# tabulated in data/mpa_take_regulations.csv (a provenance/decision record;
# cells marked VERIFY still need confirmation against CDFW regulations).
SHEEPHEAD_REINCLUSION_MPAS <- c(
  "Blue Cavern Onshore SMCA",
  "Farnsworth Onshore SMCA",
  "Swamis SMCA",
  "Long Point SMR"
)

#' @title Re-include sheephead-only MPAs for sheephead analysis
#' @description After sheephead-only MPAs are removed from the main dataset
#'   (to prevent over-representation), this function extracts the rows that
#'   should be added back for sheephead-specific analyses. It handles
#'   Santa Barbara Island SMR as a special case: only PISCO sheephead data
#'   and all KFM data (regardless of taxon) are re-included from that site.
#' @param full_data Data frame; the complete (pre-exclusion) dataset.
#' @param mpa_col Character; name of the MPA identifier column.
#'   Defaults to \code{"CA_MPA_Name_Short"}.
#' @param taxon_col Character; name of the taxon identifier column.
#'   Defaults to \code{"y"} (abbreviated codes). Use
#'   \code{"taxon_name"} for raw response data with full scientific names.
#' @param sheephead_value Character; the value in \code{taxon_col} that
#'   identifies sheephead records. Defaults to \code{"SPUL"} for response
#'   ratio data; use \code{"Semicossyphus pulcher"} for raw response data.
#' @return A data frame containing only the rows to be re-included.
#'   This should be rbind()-ed to the main (exclusion-filtered) dataset.
reinclude_sheephead_mpas <- function(full_data, mpa_col = "CA_MPA_Name_Short",
                                     taxon_col = "y",
                                     sheephead_value = "SPUL") {
  standard_rows <- full_data[[mpa_col]] %in% SHEEPHEAD_REINCLUSION_MPAS &
                   full_data[[taxon_col]] == sheephead_value
  sbi_pisco <- full_data[[mpa_col]] == "Santa Barbara Island SMR" &
               full_data[[taxon_col]] == sheephead_value &
               full_data$source == "PISCO"
  sbi_kfm <- full_data[[mpa_col]] == "Santa Barbara Island SMR" &
              full_data$source == "KFM"
  full_data[standard_rows | sbi_pisco | sbi_kfm, ]
}


# =============================================================================
# SECTION 7: FILE I/O UTILITIES
# =============================================================================

#' @title Read CSV with existence check
#' @description Wrapper around \code{read.csv()} that stops with an informative
#'   error message if the file does not exist, rather than producing a cryptic
#'   R error.
#' @param file_path Character; path to the CSV file.
#' @param ... Additional arguments passed to \code{read.csv()}.
#' @return A data frame read from the CSV file.
safe_read_csv <- function(file_path, ...) {
  if (!file.exists(file_path)) {
    stop("Required data file not found: ", file_path,
         "\nSee README.md for data setup instructions.")
  }
  read.csv(file_path, ...)
}

#' @title Read RDS with graceful fallback
#' @description Attempts to read an RDS file. Returns \code{default} if the
#'   file does not exist or cannot be read (e.g., corrupted cache).
#' @param file_path Character; path to the RDS file.
#' @param default Value to return if the file is missing or unreadable.
#'   Defaults to \code{NULL}.
#' @return The deserialized R object, or \code{default} on failure.
safe_read_rds <- function(file_path, default = NULL) {
  if (!file.exists(file_path)) return(default)
  tryCatch(
    readRDS(file_path),
    error = function(e) {
      warning("Failed to read cache file ", file_path, ": ", e$message)
      default
    }
  )
}

#' @title Run bootstrap computation with RDS caching
#' @description Executes a bootstrap computation and caches the result as an
#'   RDS file under \code{data/cache/}. On subsequent calls, loads the cached
#'   result instead of recomputing (unless forced). This avoids re-running
#'   expensive bootstrap resampling during iterative development.
#' @param cache_name Character; base name for the cache file (without path
#'   or extension). The file is saved as
#'   \code{data/cache/<cache_name>.rds}.
#' @param seed Integer; random seed set before calling \code{compute_fn} to
#'   ensure reproducible bootstrap results.
#' @param compute_fn Function (no arguments); the bootstrap computation to
#'   execute. Should return the result to be cached.
#' @param force Logical; if \code{TRUE}, recompute even if cache exists.
#'   Also respects a global \code{FORCE_BOOTSTRAP} flag in the global
#'   environment.
#' @return The (possibly cached) result of \code{compute_fn()}.
run_cached_bootstrap <- function(cache_name, seed, compute_fn, force = FALSE) {
  cache_path <- here::here("data", "cache", paste0(cache_name, ".rds"))
  force_global <- should_force_bootstrap()

  if (!force && !force_global && file.exists(cache_path)) {
    cat("    Loading cached", cache_name, "results...\n")
    result <- safe_read_rds(cache_path)
    if (!is.null(result)) return(result)
    cat("    Cache invalid, recomputing...\n")
  }

  cat("    Computing", cache_name, "(this may take a while)...\n")
  set.seed(seed)
  result <- compute_fn()

  cache_dir <- dirname(cache_path)
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  saveRDS(result, cache_path)
  cat("    Cached", cache_name, "results to", cache_path, "\n")
  result
}

#' @title Merge with row-count validation
#' @description Wrapper around \code{merge()} that warns if the merge drops
#'   more rows than expected. Useful for catching silent data loss from
#'   mismatched join keys.
#' @param x Data frame; left table.
#' @param y Data frame; right table.
#' @param by Character vector; column names to join on.
#' @param all.x Logical; if \code{TRUE}, keep all rows from \code{x}
#'   (left join). Defaults to \code{TRUE}.
#' @param all.y Logical; if \code{TRUE}, keep all rows from \code{y}
#'   (right join). Defaults to \code{FALSE}.
#' @param warn_threshold Integer; number of dropped rows allowed before
#'   triggering a warning. Defaults to 0 (warn on any row loss).
#' @param ... Additional arguments passed to \code{merge()}.
#' @return The merged data frame.
validated_merge <- function(x, y, by = NULL, by.x = by, by.y = by,
                            all.x = TRUE, all.y = FALSE,
                            warn_threshold = 0, ...) {
  n_before <- nrow(x)
  result <- merge(x, y, by.x = by.x, by.y = by.y, all.x = all.x, all.y = all.y, ...)
  n_after <- nrow(result)

  if (n_after < n_before && (n_before - n_after) > warn_threshold) {
    warning("Merge dropped ", n_before - n_after, " rows from ",
            deparse(substitute(x)), " (",
            round(100 * (n_before - n_after) / n_before, 1), "%)",
            call. = FALSE)
  }
  if (n_after > n_before) {
    warning("Merge INCREASED row count from ", n_before, " to ", n_after,
            " (+", n_after - n_before, " rows). This may indicate a many-to-many join ",
            "due to non-unique keys in ", deparse(substitute(y)), ".",
            call. = FALSE)
  }
  result
}

#' @title Safe division with zero-handling
#' @description Divides numerator by divisor, replacing zero denominators with
#'   NA (and issuing a warning) to avoid Inf results. Used primarily in
#'   \code{calculate_log_response_ratio()} where reference-site values can
#'   occasionally be zero.
#' @param numerator Numeric vector; the dividend.
#' @param divisor Numeric vector; the divisor. Zeros are replaced with NA
#'   if \code{replace_zero = TRUE}.
#' @param replace_zero Logical; if \code{TRUE} (default), replace zero
#'   divisors with NA before division.
#' @param context Character; optional label for the warning message to help
#'   identify which calculation triggered the zero-division warning.
#' @return Numeric vector of quotients (same length as inputs).
safe_divide <- function(numerator, divisor, replace_zero = TRUE, context = "") {
  n_zero <- sum(divisor == 0, na.rm = TRUE)
  if (n_zero > 0) {
    msg <- paste0("Division by zero: ", n_zero, " zero values in divisor")
    if (nchar(context) > 0) msg <- paste0(msg, " (", context, ")")
    warning(msg)
    if (replace_zero) divisor[divisor == 0] <- NA
  }
  numerator / divisor
}


# =============================================================================
# SECTION 8: FORCE_BOOTSTRAP UTILITY
# =============================================================================

#' @title Check whether bootstrap recomputation is requested
#' @description Centralized check for the FORCE_BOOTSTRAP global flag.
#'   Returns TRUE only when the flag exists in the global environment AND
#'   is exactly TRUE. Warns if the flag is set to a non-logical value
#'   (e.g., "yes", 1, "TRUE") to catch common mistakes.
#' @return Logical scalar.
should_force_bootstrap <- function() {
  if (!exists("FORCE_BOOTSTRAP", envir = .GlobalEnv)) {
    return(FALSE)
  }
  force_flag <- get("FORCE_BOOTSTRAP", envir = .GlobalEnv)
  if (isTRUE(force_flag)) {
    return(TRUE)
  }
  if (!isFALSE(force_flag)) {
    warning("FORCE_BOOTSTRAP must be TRUE or FALSE, got: ",
            paste(class(force_flag), collapse = "/"), " = ", deparse(force_flag),
            ". Treating as FALSE.", call. = FALSE)
  }
  return(FALSE)
}


# =============================================================================
# SECTION 9: COLUMN VALIDATION UTILITY
# =============================================================================

#' @title Validate required columns and subset a data frame
#' @description Checks that all required columns exist in a data frame.
#'   Stops with an informative error if any are missing. Optionally warns
#'   about extra columns being dropped.
#' @param df Data frame to validate and subset.
#' @param required_cols Character vector of required column names.
#' @param context Character; label for error/warning messages (e.g., "Swath.join.sub").
#' @param drop_extra Logical; if TRUE (default), return only required_cols.
#'   If FALSE, return all columns but still validate that required ones exist.
#' @return The data frame, subset to required_cols if drop_extra = TRUE.
validate_and_subset_columns <- function(df, required_cols, context = "", drop_extra = TRUE) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns", if (nchar(context) > 0) paste0(" in ", context),
         ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (drop_extra) {
    df[, required_cols, drop = FALSE]
  } else {
    df
  }
}


# =============================================================================
# SECTION 10: EXCLUDED LTER SITES
# =============================================================================
# LTER sites excluded from analysis. These are maintained as a constant
# for consistency with how PISCO and KFM exclusions are handled.

#' LTER site exclusions
#' Arroyo Hondo (AHND / a-l-02) is excluded because the site does not meet
#' the kelp-forest habitat criteria for this analysis: it is dominated by
#' surfgrass / sand with sparse, ephemeral kelp coverage and lacks a
#' persistent giant-kelp forest community, so it is not a valid kelp-forest
#' control or treatment site for the trophic-cascade question. (Updated
#' 2026-04 from earlier "data quality" framing, per Emily review.)
#' Two IDs because MBON uses "a-l-02" while standalone LTER CSVs use "AHND".
EXCLUDED_LTER_SITES <- c("a-l-02")
EXCLUDED_LTER_SITE_CODES <- c("AHND")


# =============================================================================
# Confirmation message
# =============================================================================
cat("Utility functions and constants loaded.\n")
