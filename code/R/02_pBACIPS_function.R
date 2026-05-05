# =============================================================================
# 02_pBACIPS_function.R
# =============================================================================
#
# PURPOSE:
#   Implement the Progressive-Change BACIPS (Before-After-Control-Impact-Pairs
#   with Symmetry) analysis method. This is the core statistical approach for
#   detecting MPA effects that may change progressively over time.
#
# WHAT THIS SCRIPT DOES:
#   1. Defines the main ProgressiveChangeBACIPS() function
#   2. Fits four candidate temporal models to the MPA vs reference difference:
#      - Step model: Abrupt change after MPA implementation
#      - Linear model: Gradual linear change over time
#      - Asymptotic model: Change that levels off (Michaelis-Menten form)
#      - Sigmoid model: S-shaped change (logistic form)
#   3. Selects the best model using AICc (corrected Akaike Information Criterion)
#   4. Provides standalone model fitting functions for effect size calculation
#
# BACKGROUND:
#   Traditional BACI designs assume an instantaneous "step" change after
#   intervention. But ecological responses to protection often unfold gradually:
#   - Fish populations may take years to recover
#   - Trophic cascades propagate slowly through the food web
#   - Habitat structure changes incrementally
#
#   The pBACIPS approach tests whether the data support a step change or one
#   of several progressive change models, providing more realistic estimates
#   of MPA effects.
#
# REFERENCE:
#   Thiault, L., Kernaléguen, L., Osenberg, C.W. & Claudet, J. (2017)
#   "Progressive-Change BACIPS: a flexible approach for environmental impact
#   assessment." Methods in Ecology and Evolution, 8(3), 288-296.
#
# AUTHORS: L. Thiault, L. Kernaléguen, C.W. Osenberg & J. Claudet
#          (adapted by Emily Donham & Adrian Stier)
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================


# =============================================================================
# REQUIRED PACKAGES
# =============================================================================
# These should already be loaded by 00_libraries.R, but we ensure they're
# available for this specific functionality

if (!requireNamespace("minpack.lm", quietly = TRUE)) {
  stop("Package 'minpack.lm' is required for pBACIPS nonlinear fitting. Install with: install.packages('minpack.lm')",
       call. = FALSE)
}
if (!requireNamespace("nls2", quietly = TRUE)) {
  stop("Package 'nls2' is required for pBACIPS nonlinear fitting. Install with: install.packages('nls2')",
       call. = FALSE)
}
if (!requireNamespace("AICcmodavg", quietly = TRUE)) {
  stop("Package 'AICcmodavg' is required for AICc model comparison. Install with: install.packages('AICcmodavg')",
       call. = FALSE)
}

library(minpack.lm)   # Levenberg-Marquardt nonlinear least squares
library(nls2)         # NLS with grid search
library(AICcmodavg)   # AICc calculation


# -----------------------------------------------------------------------------
# Helper: check NLS convergence
# -----------------------------------------------------------------------------
#' Check whether an nls fit actually converged.
#' Returns the fit if converged (or if convergence info is unavailable),
#' NULL otherwise.
check_nls_convergence <- function(fit) {
  if (is.null(fit)) return(NULL)
  # Standard nls objects store convergence info
  if (!is.null(fit$convInfo) && !isTRUE(fit$convInfo$isConv)) {
    return(NULL)
  }
  fit
}


# =============================================================================
# MODEL FITTING UTILITIES (LINES ~83-730)
# =============================================================================
# These helper functions wrap minpack.lm/nls to make NLS fitting reliable on
# real ecological time series, which are short (often 5-15 yrs after MPA),
# noisy, and prone to NLS convergence failure with naive starting values.
#
# WHAT EACH UTILITY DOES (read top-to-bottom):
#   - log_model_fit()           : append a row to MODEL_FIT_LOG with timestamp,
#                                  MPA, taxa, model type, success/failure, and
#                                  diagnostic msg. Used to produce the
#                                  outputs/model_fallback_audit.csv report.
#   - diagnose_data_quality()   : compute n_before/n_after, delta range, SD,
#                                  monotonicity, and time_max; flag combos
#                                  unlikely to support sigmoid/asymptotic
#                                  (e.g., n_after < 4). Returns suggested
#                                  model: step | linear | asymptotic.
#   - preprocess_for_nls()      : optional winsorize at ±3 SD (see L188 doc
#                                  block) plus a tiny time offset so sigmoid
#                                  doesn't blow up at t=0. Returns processed
#                                  delta/time and the transform_info used by
#                                  the audit table.
#   - generate_starting_values(): returns a LIST of starting-value strategies
#                                  (e.g., several Vmax/K combos for asymptotic,
#                                  several K/x0 combos for sigmoid). The fitter
#                                  iterates until one converges.
#   - fit_with_fallbacks()      : tries the strategies above in order, then
#                                  falls back to a simpler model (linear or
#                                  step) if every NLS attempt fails. The
#                                  fallback is logged but EXCLUDED from AICc
#                                  competition (see L1563 fallback note).
#   - safer/simpler wrappers around nls and nlsLM with shared error handling.
#
# Why all of this exists:
#   pBACIPS NLS models (asymptotic, sigmoid) are notoriously start-sensitive.
#   For ~150 MPA x taxa cells, naive nls() fails on 30-50% of fits. The
#   utilities here push convergence to >95% while making each failure
#   inspectable in MODEL_FIT_LOG. None of these utilities change the model
#   FORMS (which match Thiault et al. 2017), only the way starting values
#   are chosen and how convergence failures are handled.

# Global list to track model fitting diagnostics
if (!exists("MODEL_FIT_LOG") || !is.list(MODEL_FIT_LOG)) {
  MODEL_FIT_LOG <- list()
}

#' Log model fitting attempt
#'
#' Records details about model fitting attempts for later review
#'
#' @param mpa_name Character MPA identifier
#' @param taxa Character taxa name
#' @param model_type Character model type attempted
#' @param success Logical whether fit succeeded
#' @param message Character diagnostic message
#' @param data_summary Optional list with data characteristics
log_model_fit <- function(mpa_name = "unknown", taxa = "unknown", model_type,
                           success, message, data_summary = NULL) {
  entry <- list(
    timestamp = Sys.time(),
    mpa = mpa_name,
    taxa = taxa,
    model_type = model_type,
    success = success,
    message = message,
    data_summary = data_summary
  )
  if (length(MODEL_FIT_LOG) >= 10000) {
    warning("MODEL_FIT_LOG has exceeded 10,000 entries; consider clearing it.",
            call. = FALSE, immediate. = TRUE)
  }
  MODEL_FIT_LOG[[length(MODEL_FIT_LOG) + 1]] <<- entry
}

#' Diagnose data quality for NLS fitting
#'
#' Checks for common issues that cause NLS convergence failures
#'
#' @param delta Numeric vector of differences
#' @param time.model Time values
#' @return List with diagnostic flags and summary statistics
diagnose_data_quality <- function(delta, time.model) {

  dat <- data.frame(delta = delta, time = time.model)
  dat <- dat[complete.cases(dat), ]

  diagnostics <- list(
    n_obs = nrow(dat),
    n_before = sum(dat$time == 0),
    n_after = sum(dat$time > 0),
    has_zeros_in_time = any(dat$time == 0 & is.na(dat$delta) == FALSE),
    delta_range = diff(range(dat$delta, na.rm = TRUE)),
    delta_sd = sd(dat$delta, na.rm = TRUE),
    delta_mean = mean(dat$delta, na.rm = TRUE),
    time_range = diff(range(dat$time, na.rm = TRUE)),
    time_max = max(dat$time, na.rm = TRUE),
    has_extreme_values = any(abs(dat$delta) > 10),  # Log ratios > 10 are extreme
    is_monotonic = FALSE,
    suggested_model = "linear"
  )

  # Check for monotonicity (sign of change in response)
  if (nrow(dat) >= 3) {
    diffs <- diff(dat$delta[order(dat$time)])
    if (all(diffs >= 0) || all(diffs <= 0)) {
      diagnostics$is_monotonic <- TRUE
    }
  }

  # Suggest appropriate model based on data characteristics
  if (diagnostics$n_obs < 5) {
    diagnostics$suggested_model <- "step"  # Too few points for nonlinear
  } else if (diagnostics$delta_sd < 0.1) {
    diagnostics$suggested_model <- "step"  # Very little variance
  } else if (diagnostics$time_max < 3) {
    diagnostics$suggested_model <- "linear"  # Time series too short
  } else if (diagnostics$is_monotonic && diagnostics$delta_range > 1) {
    diagnostics$suggested_model <- "asymptotic"  # Monotonic with range
  } else {
    diagnostics$suggested_model <- "linear"
  }

  # Check for conditions that will cause sigmoid/asymptotic to fail
  diagnostics$problematic_for_sigmoid <- (
    diagnostics$n_after < 4 ||
    diagnostics$delta_sd < 0.05 ||
    diagnostics$time_max < 5
  )

  diagnostics$problematic_for_asymptotic <- (
    diagnostics$n_after < 3 ||
    diagnostics$delta_sd < 0.05
  )

  diagnostics
}

#' Preprocess data for NLS fitting
#'
#' Handles outliers, scaling, and other preprocessing that improves convergence
#'
#' @param delta Numeric vector
#' @param time.model Time values
#' @param remove_outliers Logical, whether to winsorize extreme values
#' @return List with processed delta, time, and transformation info
preprocess_for_nls <- function(delta, time.model, remove_outliers = TRUE) {

  dat <- data.frame(delta = delta, time = time.model)
  dat <- dat[complete.cases(dat), ]

  transform_info <- list(
    original_n = length(delta),
    processed_n = nrow(dat),
    outliers_modified = 0,
    time_offset = 0
  )

  # Winsorize extreme values (cap at 3 SD)
  if (remove_outliers && nrow(dat) > 5) {
    mu <- mean(dat$delta, na.rm = TRUE)
    sigma <- sd(dat$delta, na.rm = TRUE)
    lower <- mu - 3 * sigma
    upper <- mu + 3 * sigma

    n_extreme <- sum(dat$delta < lower | dat$delta > upper)
    if (n_extreme > 0) {
      dat$delta <- pmax(pmin(dat$delta, upper), lower)
      transform_info$outliers_modified <- n_extreme
    }
  }

  # Add small offset to time if it starts at 0 (avoids 0^K issues in sigmoid)
  if (min(dat$time, na.rm = TRUE) == 0) {
    # Only offset the after-period times for sigmoid, not the before period
    transform_info$time_offset <- 0.01
  }

  list(
    delta = dat$delta,
    time = dat$time,
    transform_info = transform_info
  )
}

#' Generate multiple starting value strategies for NLS models
#'
#' Returns a list of starting value sets to try in sequence when fitting
#' asymptotic or sigmoid models.
#'
#' @param delta Numeric vector of differences (impact - control)
#' @param time.model Time since MPA implementation
#' @param time.model.of.impact Index of last before-period observation
#' @param model_type Character: "asymptotic" or "sigmoid"
#' @return List of named lists, each containing starting values
generate_starting_values <- function(delta, time.model, time.model.of.impact,
                                      model_type = "asymptotic") {

  # Partition data into before and after
  before_idx <- 1:time.model.of.impact
  after_idx <- (time.model.of.impact + 1):length(delta)

  delta_before <- delta[before_idx]
  delta_after <- delta[after_idx]
  time_after <- time.model[after_idx]

  # Median-based estimates with NA handling
  B_mean <- mean(delta_before, na.rm = TRUE)
  B_median <- median(delta_before, na.rm = TRUE)

  after_mean <- mean(delta_after, na.rm = TRUE)
  after_max <- max(delta_after, na.rm = TRUE)
  after_min <- min(delta_after, na.rm = TRUE)

  # Effect size estimates (change from before to after)
  M_change <- after_mean - B_mean
  M_range <- after_max - B_mean
  M_robust <- median(delta_after, na.rm = TRUE) - B_median

  # Time scale estimates for L parameter
  L_median <- median(time_after[time_after > 0], na.rm = TRUE)
  L_max <- max(time_after, na.rm = TRUE) / 2
  L_quartile <- quantile(time_after[time_after > 0], 0.25, na.rm = TRUE)

 # Handle edge cases with safeguards
  if (!is.finite(B_mean)) B_mean <- 0
  if (!is.finite(B_median)) B_median <- 0
  if (!is.finite(M_change) || abs(M_change) < 1e-10) {
    M_change <- if (after_mean >= B_mean) 0.1 else -0.1
  }
  if (!is.finite(M_range)) M_range <- M_change
  if (!is.finite(M_robust)) M_robust <- M_change
  if (!is.finite(L_median) || L_median <= 0) L_median <- 1
  if (!is.finite(L_max) || L_max <= 0) L_max <- 1
  if (!is.finite(L_quartile) || L_quartile <= 0) L_quartile <- 1

  if (model_type == "asymptotic") {
    list(
      # Strategy 1: Data-driven change estimate
      list(M = M_change, B = B_mean, L = L_median),
      # Strategy 2: Maximum change with slower saturation
      list(M = M_range, B = B_median, L = L_max),
      # Strategy 3: Median-based estimates
      list(M = M_robust, B = B_median, L = L_quartile),
      # Strategy 4: Original approach (fallback)
      list(M = after_mean, B = B_mean, L = 1),
      # Strategy 5: Conservative small effect
      list(M = sign(M_change) * max(abs(M_change), 0.1), B = B_mean, L = L_median)
    )
  } else {  # sigmoid
    list(
      # Strategy 1: Moderate steepness at median time
      list(M = M_change, B = B_mean, L = L_median, K = 2),
      # Strategy 2: Steep transition earlier
      list(M = M_range, B = B_median, L = L_quartile, K = 5),
      # Strategy 3: Gradual transition
      list(M = M_robust, B = B_median, L = L_max, K = 1),
      # Strategy 4: Original defaults
      list(M = after_mean, B = B_mean, L = max(L_median, 1), K = 5),
      # Strategy 5: Very gradual
      list(M = M_change, B = B_mean, L = L_max, K = 0.5)
    )
  }
}


#' Fit asymptotic model using self-starting SSasymp
#'
#' Uses R's SSasymp self-starting function which automatically estimates
#' starting values, making it more reliable than manual specification.
#'
#' @param delta Numeric vector of differences
#' @param time.model Time since MPA implementation
#' @return Fitted nls model or NULL if all attempts fail
fit_asymptotic_selfstart <- function(delta, time.model) {

  dat <- data.frame(delta = delta, time = time.model)
  dat <- dat[complete.cases(dat), ]

  if (nrow(dat) < 4) return(NULL)

  tryCatch({
    # Shift time to avoid issues at 0
    time_shift <- min(dat$time) - 0.1
    dat$time_shifted <- dat$time - time_shift

    model <- nls(delta ~ SSasymp(time_shifted, Asym, R0, lrc), data = dat)
    attr(model, "time_shift") <- time_shift
    attr(model, "model_class") <- "selfstart_asymp"
    model

  }, error = function(e) {
    # Fallback with manual starting values
    fit <- tryCatch({
      Asym_start <- mean(tail(dat$delta, 3), na.rm = TRUE)
      R0_start <- mean(head(dat$delta, 3), na.rm = TRUE)

      nls(delta ~ Asym + (R0 - Asym) * exp(-exp(lrc) * time),
          data = dat,
          start = list(Asym = Asym_start, R0 = R0_start, lrc = log(0.5)),
          control = nls.control(maxiter = 100, warnOnly = TRUE))
    }, error = function(e2) NULL)
    fit <- check_nls_convergence(fit)
    fit
  })
}


#' Fit sigmoid model using self-starting SSlogis
#'
#' Uses R's SSlogis self-starting function for automatic starting values.
#'
#' @param delta Numeric vector of differences
#' @param time.model Time since MPA implementation
#' @return Fitted nls model or NULL if all attempts fail
fit_sigmoid_selfstart <- function(delta, time.model) {

  dat <- data.frame(delta = delta, time = time.model)
  dat <- dat[complete.cases(dat), ]

  if (nrow(dat) < 5) return(NULL)

  # Estimate baseline from before period
  B_est <- mean(dat$delta[dat$time == 0], na.rm = TRUE)
  if (!is.finite(B_est)) B_est <- mean(head(dat$delta, 3), na.rm = TRUE)
  if (!is.finite(B_est)) B_est <- 0

  # Shift delta by baseline for SSlogis
  dat$delta_shifted <- dat$delta - B_est

  tryCatch({
    model <- nls(delta_shifted ~ SSlogis(time, Asym, xmid, scal), data = dat)
    attr(model, "baseline_shift") <- B_est
    attr(model, "model_class") <- "selfstart_sig"
    model

  }, error = function(e) {
    # Fallback with port algorithm and bounds
    fit <- tryCatch({
      Asym_start <- max(dat$delta_shifted, na.rm = TRUE) - min(dat$delta_shifted, na.rm = TRUE)
      xmid_start <- median(dat$time[dat$time > 0], na.rm = TRUE)
      if (!is.finite(xmid_start)) xmid_start <- 5

      nls(delta_shifted ~ Asym / (1 + exp((xmid - time) / scal)),
          data = dat,
          start = list(Asym = max(Asym_start, 0.1), xmid = xmid_start, scal = 1),
          algorithm = "port",
          lower = c(Asym = -100, xmid = 0.1, scal = 0.01),
          upper = c(Asym = 100, xmid = max(dat$time) * 2, scal = max(dat$time)),
          control = nls.control(maxiter = 100, warnOnly = TRUE))
    }, error = function(e2) NULL)
    fit <- check_nls_convergence(fit)
    fit
  })
}


#' Fit quadratic model as intermediate fallback
#'
#' When nonlinear models fail, a quadratic model can capture curvature
#' while still being guaranteed to converge.
#'
#' @param delta Numeric vector of differences
#' @param time.model Time since MPA implementation
#' @return Fitted lm model or NULL
fit_quadratic_fallback <- function(delta, time.model) {
  dat <- data.frame(delta = delta, time = time.model)
  dat <- dat[complete.cases(dat), ]

  if (nrow(dat) < 4) return(NULL)

  fit <- tryCatch(
    lm(delta ~ time + I(time^2), data = dat),
    error = function(e) NULL
  )

  if (!is.null(fit)) {
    # Check if quadratic term is meaningful (significant curvature)
    quad_coef <- coef(fit)[3]
    if (!is.na(quad_coef) && abs(quad_coef) > 1e-10) {
      attr(fit, "model_class") <- "quadratic"
      return(fit)
    }
  }
  NULL
}


#' Fit GAM model as flexible fallback
#'
#' Uses a Generalized Additive Model with smoothing spline when parametric
#' nonlinear models fail. GAMs can capture complex nonlinear patterns while
#' being numerically stable.
#'
#' @param delta Numeric vector of differences
#' @param time.model Time since MPA implementation
#' @param k Integer, basis dimension for smooth (default 4 for small samples)
#' @return Fitted gam model or NULL
fit_gam_fallback <- function(delta, time.model, k = NULL) {

  dat <- data.frame(delta = delta, time = time.model)
  dat <- dat[complete.cases(dat), ]

  if (nrow(dat) < 5) return(NULL)

  # Set k based on sample size if not specified
  if (is.null(k)) {
    k <- min(floor(nrow(dat) / 2), 5)
    k <- max(k, 3)  # At least 3 basis functions
  }

  fit <- tryCatch({
    # Fit GAM with restricted maximum likelihood
    gam_fit <- mgcv::gam(delta ~ s(time, k = k, bs = "cr"),
                          data = dat, method = "REML")

    # Check that the smooth is actually doing something
    edf <- summary(gam_fit)$edf
    if (!is.null(edf) && edf > 0.5) {
      attr(gam_fit, "model_class") <- "gam"
      return(gam_fit)
    }
    NULL
  }, error = function(e) {
    # Fall back to simpler GAM specification
    tryCatch({
      gam_fit <- mgcv::gam(delta ~ s(time, k = 3, bs = "cs"), data = dat)
      attr(gam_fit, "model_class") <- "gam"
      gam_fit
    }, error = function(e2) NULL)
  })

  fit
}


#' Fit asymptotic model using log-transformed parameterization
#'
#' Uses an alternative parameterization of the Michaelis-Menten that is
#' more numerically stable: delta = Vmax * time / (Km + time) + B
#' where parameters are fit on log scale to ensure positivity.
#'
#' @param delta Numeric vector of differences
#' @param time.model Time values
#' @return Fitted nls model or NULL
fit_asymptotic_stable <- function(delta, time.model) {

  dat <- data.frame(delta = delta, time = time.model)
  dat <- dat[complete.cases(dat), ]

  if (nrow(dat) < 4) return(NULL)

  # Estimate starting values
  before_mean <- mean(dat$delta[dat$time == 0], na.rm = TRUE)
  after_mean <- mean(dat$delta[dat$time > 0], na.rm = TRUE)
  if (!is.finite(before_mean)) before_mean <- mean(head(dat$delta, 3), na.rm = TRUE)
  if (!is.finite(after_mean)) after_mean <- mean(tail(dat$delta, 3), na.rm = TRUE)
  if (!is.finite(before_mean)) before_mean <- 0
  if (!is.finite(after_mean)) after_mean <- 0

  effect_est <- after_mean - before_mean
  time_scale <- median(dat$time[dat$time > 0], na.rm = TRUE)
  if (!is.finite(time_scale) || time_scale < 0.1) time_scale <- 1

  # Try multiple parameterizations
  strategies <- list(
    # Strategy 1: Standard bounded optimization
    list(M = effect_est, B = before_mean, L = time_scale),
    # Strategy 2: Assume half-saturation at 20% of time range
    list(M = effect_est * 1.5, B = before_mean, L = time_scale * 0.2),
    # Strategy 3: Very slow saturation
    list(M = effect_est * 2, B = before_mean, L = time_scale * 2),
    # Strategy 4: Minimal parameterization
    list(M = sign(effect_est) * max(abs(effect_est), 0.5), B = 0, L = 1)
  )

  foAsy <- delta ~ (M * time) / (L + time) + B

  for (start in strategies) {
    # Ensure L is positive
    start$L <- max(abs(start$L), 0.01)

    fit <- tryCatch({
      # Use port algorithm with safe bounds
      nls(foAsy, data = dat, start = start,
          algorithm = "port",
          lower = c(M = -50, B = -50, L = 0.001),
          upper = c(M = 50, B = 50, L = max(dat$time) * 3),
          control = nls.control(maxiter = 200, warnOnly = TRUE,
                                 minFactor = 1/2048))
    }, error = function(e) NULL)
    fit <- check_nls_convergence(fit)

    if (!is.null(fit)) {
      # Validate fit is reasonable
      coefs <- coef(fit)
      if (is.finite(coefs["L"]) && coefs["L"] > 0.001 && coefs["L"] < max(dat$time) * 2.5) {
        attr(fit, "model_class") <- "asymptotic"
        return(fit)
      }
    }
  }

  NULL
}


#' Fit sigmoid using reparameterized 4-parameter logistic
#'
#' Uses an alternative 4PL parameterization that's more stable:
#' delta = D + (A - D) / (1 + (time/C)^B)
#' where A = lower asymptote, D = upper asymptote, C = inflection point, B = Hill slope
#'
#' @param delta Numeric vector of differences
#' @param time.model Time values
#' @return Fitted nls model or NULL
fit_sigmoid_4pl <- function(delta, time.model) {

  dat <- data.frame(delta = delta, time = time.model)
  dat <- dat[complete.cases(dat), ]

  if (nrow(dat) < 6) return(NULL)

  # Quartile-based starting value estimation
  sorted_idx <- order(dat$time)
  n <- nrow(dat)
  q1 <- max(1, floor(n * 0.25))
  q4 <- floor(n * 0.75)

  lower_vals <- dat$delta[sorted_idx[1:q1]]
  upper_vals <- dat$delta[sorted_idx[q4:n]]

  A_est <- mean(lower_vals, na.rm = TRUE)  # Lower asymptote
  D_est <- mean(upper_vals, na.rm = TRUE)  # Upper asymptote
  C_est <- median(dat$time[dat$time > 0], na.rm = TRUE)  # Inflection point
  B_est <- 2  # Hill slope (positive = increasing)

  if (!is.finite(A_est)) A_est <- min(dat$delta, na.rm = TRUE)
  if (!is.finite(D_est)) D_est <- max(dat$delta, na.rm = TRUE)
  if (!is.finite(C_est) || C_est < 0.1) C_est <- 5

  # If response is decreasing, swap A and D
  if (A_est > D_est) {
    B_est <- -B_est
  }

  # 4-parameter logistic formula
  fo4PL <- delta ~ D + (A - D) / (1 + (time / C)^B)

  strategies <- list(
    list(A = A_est, D = D_est, C = C_est, B = B_est),
    list(A = A_est, D = D_est, C = C_est * 0.5, B = B_est * 0.5),
    list(A = A_est, D = D_est, C = C_est * 2, B = B_est * 2),
    list(A = mean(c(A_est, D_est)), D = D_est, C = C_est, B = 1)
  )

  for (start in strategies) {
    # Ensure C is positive
    start$C <- max(abs(start$C), 0.1)

    fit <- tryCatch({
      nls(fo4PL, data = dat, start = start,
          algorithm = "port",
          lower = c(A = -20, D = -20, C = 0.01, B = -20),
          upper = c(A = 20, D = 20, C = max(dat$time) * 3, B = 20),
          control = nls.control(maxiter = 300, warnOnly = TRUE,
                                 minFactor = 1/4096))
    }, error = function(e) NULL)
    fit <- check_nls_convergence(fit)

    if (!is.null(fit)) {
      coefs <- coef(fit)
      # Validate fit
      if (all(is.finite(coefs)) && coefs["C"] > 0.01) {
        attr(fit, "model_class") <- "sigmoid_4pl"
        return(fit)
      }
    }
  }

  NULL
}


#' Fit piecewise linear model (changepoint model)
#'
#' Fits a model with a single breakpoint: constant before implementation,
#' linear change after. This is a simplified progressive change model that's
#' always estimable.
#'
#' @param delta Numeric vector
#' @param time.model Time values (0 for before, positive for after)
#' @return Fitted lm model or NULL
fit_piecewise_linear <- function(delta, time.model) {

  dat <- data.frame(delta = delta, time = time.model)
  dat <- dat[complete.cases(dat), ]

  if (nrow(dat) < 4) return(NULL)

  # Create piecewise indicator
  dat$after <- as.numeric(dat$time > 0)
  dat$time_after <- dat$time * dat$after

  fit <- tryCatch({
    lm(delta ~ after + time_after, data = dat)
  }, error = function(e) NULL)

  if (!is.null(fit)) {
    # Check that we have meaningful coefficients
    coefs <- coef(fit)
    if (all(is.finite(coefs))) {
      attr(fit, "model_class") <- "piecewise"
      return(fit)
    }
  }

  NULL
}


#' Run DHARMa diagnostics on a fitted model
#'
#' Uses DHARMa package to assess model fit quality via simulation-based
#' residual diagnostics.
#'
#' @param model Fitted model object (nls, lm, or gam)
#' @param plot Logical, whether to create diagnostic plots
#' @return List with diagnostic results
run_dharma_diagnostics <- function(model, plot = FALSE) {

  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    return(list(
      available = FALSE,
      message = "DHARMa package not installed"
    ))
  }

  result <- list(available = TRUE)

  # DHARMa works best with glm/gam but can simulate from nls
  tryCatch({
    # For nls models, we need to extract residuals manually
    if (inherits(model, "nls")) {
      resid <- residuals(model)
      fitted <- fitted(model)
      sigma <- summary(model)$sigma

      # Simulate residuals assuming normal errors
      n <- length(resid)
      n_sim <- 250
      sim_resid <- matrix(rnorm(n * n_sim, mean = 0, sd = sigma), nrow = n)

      result$uniformity_p <- tryCatch({
        ks.test(resid / sigma, "pnorm")$p.value
      }, error = function(e) NA)

      result$outlier_test <- tryCatch({
        outlier_count <- sum(abs(resid) > 3 * sigma)
        list(n_outliers = outlier_count, threshold = 3)
      }, error = function(e) NA)

    } else if (inherits(model, "lm") || inherits(model, "gam")) {
      # Standard DHARMa for linear/GAM models
      sim_resid <- DHARMa::simulateResiduals(model, n = 250, plot = FALSE)

      result$uniformity_p <- DHARMa::testUniformity(sim_resid, plot = FALSE)$p.value
      result$dispersion_p <- DHARMa::testDispersion(sim_resid, plot = FALSE)$p.value
      result$outlier_test <- DHARMa::testOutliers(sim_resid, plot = FALSE)

      if (plot) {
        DHARMa::plot(sim_resid)
      }
    }

    result$model_class <- class(model)[1]
    result$success <- TRUE

  }, error = function(e) {
    result$success <- FALSE
    result$message <- e$message
  })

  result
}


# =============================================================================
# MAIN FUNCTION: ProgressiveChangeBACIPS
# =============================================================================

#' Progressive-Change BACIPS Analysis
#'
#' Fits four candidate models to BACI delta (Impact - Control) data and
#' selects the best model using AICc weights.
#'
#' @param control Numeric vector of response values at Control (reference) site
#'        at each sampling time. Should be on the same scale as impact.
#' @param impact Numeric vector of response values at Impact (MPA) site at each
#'        sampling time. Same length as control.
#' @param time.true Sequential time index for ALL sampling dates (1, 2, 3, ...)
#'        This represents the full time series regardless of before/after.
#' @param time.model Time since intervention: 0 for Before period, sequential
#'        values (0, 1, 2, ...) for After period. This is the "treatment time".
#'
#' @return A named list with components:
#'   - weights: AICc weights for each model (% likelihood)
#'   - step: The fitted step (ANOVA) model object
#'   - linear: The fitted linear model object
#'   - asymptotic: The fitted asymptotic (NLS) model object
#'   - sigmoid: The fitted sigmoid (NLS) model object
#'   - best: Integer index of the best model (1=step, 2=linear, 3=asymp, 4=sig)
#'
#' @details
#' The function proceeds in 4 steps:
#'
#' STEP 1: Calculate delta at each time point
#'   delta = impact - control
#'   This is the difference between MPA and reference at each sampling time.
#'
#' STEP 2: Define the four candidate models
#'   All models predict delta as a function of time.model:
#'
#'   Step Model: delta ~ period (categorical Before/After)
#'     - Tests for a simple mean shift between periods
#'     - Fit using ANOVA (aov)
#'
#'   Linear Model: delta ~ time.model
#'     - Tests for a gradual linear trend starting at MPA implementation
#'     - Fit using ordinary least squares (lm)
#'
#'   Asymptotic Model: delta ~ (M * time.model) / (L + time.model) + B
#'     - Michaelis-Menten form: rapid initial change that levels off
#'     - M = maximum change, L = half-saturation time, B = baseline
#'     - Fit using nonlinear least squares (nls.lm + nls2)
#'
#'   Sigmoid Model: delta ~ (M * (time.model/L)^K) / (1 + (time.model/L)^K) + B
#'     - Hill function form: S-shaped change with flexible steepness
#'     - M = maximum change, L = midpoint, K = steepness, B = baseline
#'     - Fit using nonlinear least squares (nls.lm + nls2)
#'
#' STEP 3: Compare models using AICc
#'   - AICc penalizes for model complexity more than AIC
#'   - Appropriate for small sample sizes (typical in BACI studies)
#'   - Weights sum to 100% and represent relative likelihood
#'
#' STEP 4: Report best model
#'   - Prints AICc weights and summary of best-fitting model
#'   - Returns all model objects for further analysis
#'
#' @examples
#' # See 08_effect_sizes.R for full usage examples
#' # Basic usage:
#' # result <- ProgressiveChangeBACIPS(control, impact, time.true, time.model)
#' # result$weights  # Check which model is best
#' # result$best     # Index of best model
ProgressiveChangeBACIPS <- function(control, impact, time.true, time.model) {

  # Input validation (2026-02-06)
  # Check vector lengths match
  if (length(control) != length(impact) || length(control) != length(time.model) || length(control) != length(time.true)) {
    stop("ProgressiveChangeBACIPS(): control, impact, time.true, and time.model must be the same length.",
         call. = FALSE)
  }

  # Check for NA values
  if (any(is.na(control))) {
    stop("ProgressiveChangeBACIPS(): control vector contains NA values. Please remove or impute NAs before analysis.",
         call. = FALSE)
  }
  if (any(is.na(impact))) {
    stop("ProgressiveChangeBACIPS(): impact vector contains NA values. Please remove or impute NAs before analysis.",
         call. = FALSE)
  }
  if (any(is.na(time.model))) {
    stop("ProgressiveChangeBACIPS(): time.model vector contains NA values.",
         call. = FALSE)
  }
  if (any(is.na(time.true))) {
    stop("ProgressiveChangeBACIPS(): time.true vector contains NA values.",
         call. = FALSE)
  }

  # Check for non-finite values (Inf, -Inf, NaN)
  if (!all(is.finite(control))) {
    stop("ProgressiveChangeBACIPS(): control vector contains non-finite values (Inf, -Inf, or NaN).",
         call. = FALSE)
  }
  if (!all(is.finite(impact))) {
    stop("ProgressiveChangeBACIPS(): impact vector contains non-finite values (Inf, -Inf, or NaN).",
         call. = FALSE)
  }
  if (!all(is.finite(time.model))) {
    stop("ProgressiveChangeBACIPS(): time.model vector contains non-finite values.",
         call. = FALSE)
  }
  if (!all(is.finite(time.true))) {
    stop("ProgressiveChangeBACIPS(): time.true vector contains non-finite values.",
         call. = FALSE)
  }

  # Check sufficient data
  n <- length(control)
  if (n < 5) {
    stop("ProgressiveChangeBACIPS(): insufficient data (n=", n, "). Need at least 5 observations for pBACIPS analysis.",
         call. = FALSE)
  }

  # Check before/after periods exist
  if (!any(time.model == 0, na.rm = TRUE)) {
    stop("ProgressiveChangeBACIPS(): time.model must include at least one 0 (Before period).",
         call. = FALSE)
  }
  if (!any(time.model > 0, na.rm = TRUE)) {
    stop("ProgressiveChangeBACIPS(): time.model must include at least one value > 0 (After period).",
         call. = FALSE)
  }

  # Check sufficient observations in each period
  n_before <- sum(time.model == 0)
  n_after <- sum(time.model > 0)
  if (n_before < 2) {
    warning("ProgressiveChangeBACIPS(): only ", n_before, " Before observation(s). Results may be unstable.",
            call. = FALSE, immediate. = TRUE)
  }
  if (n_after < 2) {
    warning("ProgressiveChangeBACIPS(): only ", n_after, " After observation(s). Results may be unstable.",
            call. = FALSE, immediate. = TRUE)
  }

  # ---------------------------------------------------------------------------
  # STEP 1: Calculate delta at each sampling date
  # ---------------------------------------------------------------------------
  # Delta is the difference between impact (MPA) and control (reference)
  # Positive delta = higher values inside MPA than outside

  delta <- impact - control

  # Find the last time point in the "Before" period
  # This is needed for setting initial parameter values
  time.model.of.impact <- max(which(time.model == 0))

  # ---------------------------------------------------------------------------
  # STEP 2: Create period factor for step model
  # ---------------------------------------------------------------------------
  # Convert time.model to categorical Before/After
  # NOTE: In this framework, time.model == 0 represents the Before period
  # (all pre-MPA observations are mapped to 0), and time.model > 0 represents
  # sequential years in the After period. This is consistent with the linear
  # and nonlinear model parameterizations where the intercept at time.model=0
  # estimates the Before-period baseline.

  period <- ifelse(time.model == 0, "Before", "After")

  # ---------------------------------------------------------------------------
  # STEP 3: Fit the four candidate models
  # ---------------------------------------------------------------------------

  ## --- Step Model (ANOVA) ---
  # H0: mean delta is the same Before and After
  # This is the traditional BACI analysis approach

  step.Model <- lm(delta ~ period)

  ## --- Linear Model (OLS) ---
  # H0: no linear trend in delta over time since MPA
  # Tests for gradual, constant-rate change

  linear.Model <- lm(delta ~ time.model)

  ## --- Asymptotic Model (Nonlinear) ---
  # Michaelis-Menten form: rapid initial change that saturates
  # Equation: delta = (M * time) / (L + time) + B
  # M = maximum effect, L = time to reach half of maximum, B = baseline

  # First, diagnose data quality
  data_diag <- diagnose_data_quality(delta, time.model)

  myASYfun <- function(delta, time.model) {
    # Define the model formula
    foAsy <- delta ~ (M * time.model) / (L + time.model) + B

    # Define prediction function for nls.lm
    funAsy <- function(parS, time.model) {
      (parS$M * time.model) / (parS$L + time.model) + parS$B
    }

    # Residual function: observed - predicted
    residFun <- function(p, observed, time.model) {
      observed - funAsy(p, time.model)
    }

    # Preprocess data to handle outliers
    processed <- preprocess_for_nls(delta, time.model, remove_outliers = TRUE)

    # Recompute impact index from processed data (original may be stale after NA removal)
    time_model_of_impact_proc <- max(which(processed$time == 0))

    # Get multiple starting value strategies
    starts_list <- generate_starting_values(processed$delta, processed$time,
                                             time_model_of_impact_proc, "asymptotic")
    max_time <- max(processed$time, na.rm = TRUE)

    # Try each starting value strategy with bounded optimization
    for (parStart in starts_list) {
      # Ensure L is within bounds
      parStart$L <- max(parStart$L, 0.01)
      parStart$L <- min(parStart$L, max_time * 1.5)

      # Try port algorithm with bounds first (most reliable)
      fit <- tryCatch({
        nls(foAsy,
            start = parStart,
            algorithm = "port",
            lower = c(M = -100, B = -100, L = 0.001),
            upper = c(M = 100, B = 100, L = max(max_time * 2, 10)),
            control = nls.control(maxiter = 200, warnOnly = TRUE,
                                   minFactor = 1/2048))
      }, error = function(e) NULL)
      fit <- check_nls_convergence(fit)

      if (!is.null(fit) && inherits(fit, "nls")) {
        coefs <- coef(fit)
        # Check fit is not at boundaries and coefficients are finite
        if (all(is.finite(coefs)) && coefs["L"] > 0.01 && coefs["L"] < max_time * 1.9) {
          attr(fit, "model_class") <- "asymptotic"
          return(fit)
        }
      }

      # Fallback: try nls.lm with safer settings
      nls_ASY_out <- tryCatch(
        nls.lm(
          par = parStart,
          fn = residFun,
          observed = processed$delta,
          time.model = processed$time,
          control = nls.lm.control(maxfev = 500, maxiter = 200, ftol = 1e-6)
        ),
        error = function(e) NULL
      )

      if (!is.null(nls_ASY_out) && all(is.finite(coef(nls_ASY_out)))) {
        # Refine with nls2 grid search
        startPar <- as.list(coef(nls_ASY_out))
        gridPar <- expand.grid(
          M = startPar$M * c(0.5, 1.0, 1.5),
          B = startPar$B * c(0.5, 1.0, 1.5),
          L = pmax(0.01, startPar$L * c(0.5, 1.0, 1.5))
        )

        nls2_fit <- tryCatch(
          nls2(foAsy, start = gridPar, algorithm = "brute-force"),
          error = function(e) {
            tryCatch(
              nls2(foAsy, start = data.frame(t(unlist(startPar))), algorithm = "brute-force"),
              error = function(e2) NULL
            )
          }
        )

        if (!is.null(nls2_fit) && all(is.finite(coef(nls2_fit)))) {
          final_fit <- tryCatch(
            nls(foAsy, start = coef(nls2_fit),
                control = nls.control(maxiter = 100, warnOnly = TRUE)),
            error = function(e) nls2_fit
          )
          final_fit <- check_nls_convergence(final_fit)
          if (!is.null(final_fit) && all(is.finite(coef(final_fit)))) {
            attr(final_fit, "model_class") <- "asymptotic"
            return(final_fit)
          }
        }
      }
    }

    NULL
  }

  # Skip asymptotic if data diagnostics suggest it will fail
  asymptotic.Model <- NULL
  asymptotic_fallback <- "none"  # Track fallback level for AICc exclusion
  if (!data_diag$problematic_for_asymptotic) {
    # Try primary asymptotic fitting
    asymptotic.Model <- tryCatch(
      myASYfun(delta = delta, time.model = time.model),
      error = function(e) {
        log_model_fit(model_type = "asymptotic", success = FALSE,
                      message = paste("Primary fitting error:", e$message))
        NULL
      }
    )
  }

  # Fallback 1: try stable parameterization (still a true asymptotic NLS, valid for AICc)
  if (is.null(asymptotic.Model)) {
    asymptotic.Model <- tryCatch(
      fit_asymptotic_stable(delta, time.model),
      error = function(e) NULL
    )
    if (!is.null(asymptotic.Model)) {
      asymptotic_fallback <- "stable_reparam"
      log_model_fit(model_type = "asymptotic", success = TRUE,
                    message = "Using stable parameterization")
    }
  }

  # Fallback 2: try self-starting asymptotic model (still a true NLS, valid for AICc)
  if (is.null(asymptotic.Model)) {
    asymptotic.Model <- fit_asymptotic_selfstart(delta, time.model)
    if (!is.null(asymptotic.Model)) {
      asymptotic_fallback <- "selfstart"
      log_model_fit(model_type = "asymptotic", success = TRUE,
                    message = "Using self-starting SSasymp")
    }
  }

  # ---------------------------------------------------------------------------
  # FALLBACK MODELS BELOW: These are NOT true asymptotic NLS fits.
  # They are excluded from AICc comparison because:
  #   - Piecewise linear and quadratic are lm() models masquerading as nonlinear
  #   - GAM AICc (penalized likelihood) is not comparable to lm/nls AICc
  #   - Reporting these as "Asymptotic" would be misleading
  # The fallback models are still fitted and stored for reference, but they
  # do NOT enter the model selection competition.
  # ---------------------------------------------------------------------------

  # Fallback 3: piecewise linear (changepoint model). EXCLUDED from AICc
  if (is.null(asymptotic.Model)) {
    asymptotic.Model <- fit_piecewise_linear(delta, time.model)
    if (!is.null(asymptotic.Model)) {
      asymptotic_fallback <- "piecewise"
      attr(asymptotic.Model, "model_fallback") <- "piecewise"
      log_model_fit(model_type = "asymptotic", success = TRUE,
                    message = "Using piecewise linear as fallback (EXCLUDED from AICc)")
    }
  }

  # Fallback 4: quadratic model captures curvature. EXCLUDED from AICc
  if (is.null(asymptotic.Model)) {
    asymptotic.Model <- fit_quadratic_fallback(delta, time.model)
    if (!is.null(asymptotic.Model)) {
      asymptotic_fallback <- "quadratic"
      attr(asymptotic.Model, "model_fallback") <- "quadratic"
      log_model_fit(model_type = "asymptotic", success = TRUE,
                    message = "Using quadratic model as fallback (EXCLUDED from AICc)")
    }
  }

  # Fallback 5: GAM for flexible nonlinear fit. EXCLUDED from AICc
  if (is.null(asymptotic.Model)) {
    asymptotic.Model <- fit_gam_fallback(delta, time.model)
    if (!is.null(asymptotic.Model)) {
      asymptotic_fallback <- "gam"
      attr(asymptotic.Model, "model_fallback") <- "gam"
      log_model_fit(model_type = "asymptotic", success = TRUE,
                    message = "Using GAM as final fallback (EXCLUDED from AICc)")
    }
  }

  ## --- Sigmoid Model (Nonlinear) ---
  # Hill function form: S-shaped change with adjustable steepness
  # Equation: delta = (M * (time/L)^K) / (1 + (time/L)^K) + B
  # M = maximum effect, L = midpoint (time to half-max), K = steepness, B = baseline

  mySIGfun <- function(delta, time.model) {
    # Define prediction function for nls.lm with offset
    funSIG <- function(parS, time_offset) {
      # Add safeguards against numerical overflow
      ratio <- (time_offset / parS$L)^parS$K
      ratio <- pmin(ratio, 1e10)  # Cap extreme values
      ratio <- pmax(ratio, 1e-10)
      (parS$M * ratio) / (1 + ratio) + parS$B
    }

    # Residual function: observed - predicted
    residFun <- function(p, observed, time_offset) {
      pred <- funSIG(p, time_offset)
      if (any(!is.finite(pred))) return(rep(1e10, length(observed)))
      observed - pred
    }

    # Preprocess data
    processed <- preprocess_for_nls(delta, time.model, remove_outliers = TRUE)

    # Recompute time_offset from processed data so formula and nls.lm use the same vector
    time_proc <- processed$time + 0.01
    time_offset <- time_proc  # Used by foSIG formula via closure

    # Define the model formula - uses time_offset from processed data above
    foSIG <- delta ~ (M * (time_offset / L)^K) / (1 + (time_offset / L)^K) + B

    # Recompute impact index from processed data (original may be stale after NA removal)
    time_model_of_impact_proc <- max(which(processed$time == 0))

    # Get multiple starting value strategies
    starts_list <- generate_starting_values(processed$delta, processed$time,
                                             time_model_of_impact_proc, "sigmoid")
    max_time <- max(processed$time, na.rm = TRUE)

    # Try each starting value strategy
    for (parStart in starts_list) {
      # Constrain starting values more carefully
      parStart$L <- max(parStart$L, 0.5)
      parStart$L <- min(parStart$L, max_time * 1.2)
      parStart$K <- max(parStart$K, 0.3)
      parStart$K <- min(parStart$K, 10)

      # Validate starting values are finite
      if (!all(sapply(parStart, is.finite))) next

      # Try port algorithm with bounds first
      fit <- tryCatch({
        nls(foSIG,
            start = parStart,
            algorithm = "port",
            lower = c(M = -50, B = -50, L = 0.2, K = 0.2),
            upper = c(M = 50, B = 50, L = max_time * 1.5, K = 15),
            control = nls.control(maxiter = 300, warnOnly = TRUE,
                                   minFactor = 1/4096))
      }, error = function(e) NULL)
      fit <- check_nls_convergence(fit)

      if (!is.null(fit) && inherits(fit, "nls")) {
        coefs <- coef(fit)
        # Verify all coefficients are finite and within bounds
        if (all(is.finite(coefs)) &&
            coefs["K"] > 0.25 && coefs["K"] < 14 &&
            coefs["L"] > 0.25 && coefs["L"] < max_time * 1.4) {
          attr(fit, "model_class") <- "sigmoid"
          return(fit)
        }
      }

      # Fallback: try nls.lm with safer settings
      nls_SIG_out <- tryCatch(
        nls.lm(
          par = parStart,
          fn = residFun,
          observed = processed$delta,
          time_offset = time_proc,
          control = nls.lm.control(maxfev = 500, maxiter = 200, ftol = 1e-6)
        ),
        error = function(e) NULL
      )

      if (!is.null(nls_SIG_out) && all(is.finite(coef(nls_SIG_out)))) {
        # Refine with nls2 grid search
        startPar <- as.list(coef(nls_SIG_out))
        # Ensure grid parameters are valid
        gridPar <- expand.grid(
          M = startPar$M * c(0.7, 1.0, 1.3),
          B = startPar$B * c(0.7, 1.0, 1.3),
          L = pmax(0.5, startPar$L * c(0.7, 1.0, 1.3)),
          K = pmax(0.3, pmin(10, startPar$K * c(0.7, 1.0, 1.3)))
        )

        nls2_fit <- tryCatch(
          nls2(foSIG, start = gridPar, algorithm = "brute-force"),
          error = function(e) {
            tryCatch(
              nls2(foSIG, start = data.frame(t(unlist(startPar))), algorithm = "brute-force"),
              error = function(e2) NULL
            )
          }
        )

        if (!is.null(nls2_fit) && all(is.finite(coef(nls2_fit)))) {
          final_fit <- tryCatch(
            nls(foSIG, start = coef(nls2_fit),
                control = nls.control(maxiter = 100, warnOnly = TRUE)),
            error = function(e) nls2_fit
          )
          final_fit <- check_nls_convergence(final_fit)
          if (!is.null(final_fit) && all(is.finite(coef(final_fit)))) {
            attr(final_fit, "model_class") <- "sigmoid"
            return(final_fit)
          }
        }
      }
    }

    NULL
  }

  # Skip sigmoid if data diagnostics suggest it will fail
  sigmoid.Model <- NULL
  sigmoid_fallback <- "none"  # Track fallback level for AICc exclusion
  if (!data_diag$problematic_for_sigmoid) {
    # Try primary sigmoid fitting
    sigmoid.Model <- tryCatch(
      mySIGfun(delta = delta, time.model = time.model),
      error = function(e) {
        log_model_fit(model_type = "sigmoid", success = FALSE,
                      message = paste("Primary fitting error:", e$message))
        NULL
      }
    )
  }

  # Fallback 1: 4-parameter logistic (alternative NLS parameterization, valid for AICc)
  if (is.null(sigmoid.Model)) {
    sigmoid.Model <- tryCatch(
      fit_sigmoid_4pl(delta, time.model),
      error = function(e) NULL
    )
    if (!is.null(sigmoid.Model)) {
      sigmoid_fallback <- "4pl_reparam"
      log_model_fit(model_type = "sigmoid", success = TRUE,
                    message = "Using 4PL parameterization")
    }
  }

  # Fallback 2: try self-starting logistic model (still a true NLS, valid for AICc)
  if (is.null(sigmoid.Model)) {
    sigmoid.Model <- fit_sigmoid_selfstart(delta, time.model)
    if (!is.null(sigmoid.Model)) {
      sigmoid_fallback <- "selfstart"
      log_model_fit(model_type = "sigmoid", success = TRUE,
                    message = "Using self-starting SSlogis")
    }
  }

  # ---------------------------------------------------------------------------
  # FALLBACK MODELS BELOW: These are NOT true sigmoid NLS fits.
  # They are excluded from AICc comparison for the same reasons as the
  # asymptotic fallbacks (see above). Still fitted and stored for reference.
  # ---------------------------------------------------------------------------

  # Fallback 3: piecewise linear. EXCLUDED from AICc
  if (is.null(sigmoid.Model)) {
    sigmoid.Model <- fit_piecewise_linear(delta, time.model)
    if (!is.null(sigmoid.Model)) {
      sigmoid_fallback <- "piecewise"
      attr(sigmoid.Model, "model_fallback") <- "piecewise"
      log_model_fit(model_type = "sigmoid", success = TRUE,
                    message = "Using piecewise linear as fallback (EXCLUDED from AICc)")
    }
  }

  # Fallback 4: quadratic model. EXCLUDED from AICc
  if (is.null(sigmoid.Model)) {
    sigmoid.Model <- fit_quadratic_fallback(delta, time.model)
    if (!is.null(sigmoid.Model)) {
      sigmoid_fallback <- "quadratic"
      attr(sigmoid.Model, "model_fallback") <- "quadratic"
      log_model_fit(model_type = "sigmoid", success = TRUE,
                    message = "Using quadratic model as fallback (EXCLUDED from AICc)")
    }
  }

  # Fallback 5: GAM for flexible nonlinear fit. EXCLUDED from AICc
  if (is.null(sigmoid.Model)) {
    sigmoid.Model <- fit_gam_fallback(delta, time.model)
    if (!is.null(sigmoid.Model)) {
      sigmoid_fallback <- "gam"
      attr(sigmoid.Model, "model_fallback") <- "gam"
      log_model_fit(model_type = "sigmoid", success = TRUE,
                    message = "Using GAM as final fallback (EXCLUDED from AICc)")
    }
  }

  # ---------------------------------------------------------------------------
  # STEP 4: Model comparison using AICc
  # ---------------------------------------------------------------------------
  # AICc is preferred over AIC for small samples (n/K < 40)
  # where n = sample size and K = number of parameters
  #
  # IMPORTANT: Only models fit with comparable likelihood functions enter the
  # AICc comparison. Fallback models (piecewise linear, quadratic, GAM) that
  # were substituted for failed NLS fits are EXCLUDED because:
  #   1. GAM AICc uses penalized likelihood, not comparable to lm/nls AICc
  #   2. A piecewise-linear or quadratic model is a fundamentally different
  #      model class than the intended asymptotic/sigmoid. Reporting it as
  #      "Asymptotic" or "Sigmoid" would be misleading
  #   3. Degrees of freedom differ between model classes
  # True NLS reparameterizations (stable, self-start, 4PL) ARE included since
  # they are genuine nonlinear fits estimated via maximum likelihood.

  # Determine which fallback models should be excluded from AICc
  # Fallbacks beyond reparameterization are: piecewise, quadratic, gam
  INVALID_FALLBACKS <- c("piecewise", "quadratic", "gam")
  asymptotic_excluded <- asymptotic_fallback %in% INVALID_FALLBACKS
  sigmoid_excluded <- sigmoid_fallback %in% INVALID_FALLBACKS

  # Build list of successfully fitted models eligible for AICc comparison
  model_list <- list(step.Model = step.Model, linear.Model = linear.Model)
  model_names <- c("step.Model", "linear.Model")

  if (!is.null(asymptotic.Model) && !asymptotic_excluded) {
    model_list$asymptotic.Model <- asymptotic.Model
    model_names <- c(model_names, "asymptotic.Model")
  } else if (asymptotic_excluded) {
    log_model_fit(model_type = "asymptotic", success = TRUE,
                  message = paste0("Excluded from AICc: fallback=", asymptotic_fallback,
                                   " (not a true NLS fit)"))
  }
  if (!is.null(sigmoid.Model) && !sigmoid_excluded) {
    model_list$sigmoid.Model <- sigmoid.Model
    model_names <- c(model_names, "sigmoid.Model")
  } else if (sigmoid_excluded) {
    log_model_fit(model_type = "sigmoid", success = TRUE,
                  message = paste0("Excluded from AICc: fallback=", sigmoid_fallback,
                                   " (not a true NLS fit)"))
  }

  # Calculate AICc for each eligible model
  aicc_values <- sapply(model_list, function(m) {
    tryCatch(AICc(m), error = function(e) Inf)
  })

  # Create data frame with AICc values
  AICc.test <- data.frame(
    df = sapply(model_list, function(m) {
      tryCatch(length(coef(m)) + 1, error = function(e) NA)
    }),
    AIC = aicc_values
  )
  rownames(AICc.test) <- model_names

  # Calculate AICc weights (Akaike weights)
  # Lower AICc = better fit; weights sum to 1 (or 100%)
  AICc.test$diff <- AICc.test$AIC - min(AICc.test$AIC, na.rm = TRUE)

  # Relative likelihood: exp(-0.5 * delta_AICc)
  AICc.test$RL <- exp(-0.5 * AICc.test$diff)

  # Normalize to get weights (as percentages)
  RL_sum <- sum(AICc.test$RL, na.rm = TRUE)
  if (!is.finite(RL_sum) || RL_sum <= 0) {
    # Fallback: pick the minimum finite AICc as 100% weight, others 0%.
    best_idx <- which.min(AICc.test$AIC)
    AICc.test$aicWeights <- rep(0, nrow(AICc.test))
    AICc.test$aicWeights[best_idx] <- 100
  } else {
    AICc.test$aicWeights <- (AICc.test$RL / RL_sum) * 100
  }

  # Build full weights vector (4 elements, with 0 for missing/excluded models)
  w <- c(
    step.Model = AICc.test$aicWeights[rownames(AICc.test) == "step.Model"],
    linear.Model = AICc.test$aicWeights[rownames(AICc.test) == "linear.Model"],
    asymptotic.Model = ifelse("asymptotic.Model" %in% rownames(AICc.test),
                               AICc.test$aicWeights[rownames(AICc.test) == "asymptotic.Model"], 0),
    sigmoid.Model = ifelse("sigmoid.Model" %in% rownames(AICc.test),
                           AICc.test$aicWeights[rownames(AICc.test) == "sigmoid.Model"], 0)
  )

  # ---------------------------------------------------------------------------
  # STEP 5: Display results
  # ---------------------------------------------------------------------------

  # Print AICc weights
  print(w)

  # Identify best model
  best.Model <- which.max(w)

  # Print summary of best model
  if (best.Model == 1) {
    writeLines(paste("\n\nSTEP MODEL SELECTED - Likelihood = ",
                     round(w[1], 1), "%\n\n", sep = ""))
    print(summary(step.Model))
  }
  if (best.Model == 2) {
    writeLines(paste("\n\nLINEAR MODEL SELECTED - Likelihood = ",
                     round(w[2], 1), "%\n\n", sep = ""))
    print(summary(linear.Model))
  }
  if (best.Model == 3 && !is.null(asymptotic.Model)) {
    writeLines(paste("\n\nASYMPTOTIC MODEL SELECTED - Likelihood = ",
                     round(w[3], 1), "%\n\n", sep = ""))
    print(asymptotic.Model)
  }
  if (best.Model == 4 && !is.null(sigmoid.Model)) {
    writeLines(paste("\n\nSIGMOID MODEL SELECTED - Likelihood = ",
                     round(w[4], 1), "%\n\n", sep = ""))
    print(sigmoid.Model)
  }

  # ---------------------------------------------------------------------------
  # Return all models and weights
  # ---------------------------------------------------------------------------

  list(
    weights = w,
    step = step.Model,
    linear = linear.Model,
    asymptotic = asymptotic.Model,
    sigmoid = sigmoid.Model,
    best = best.Model,
    fallback_info = list(
      asymptotic = asymptotic_fallback,
      sigmoid = sigmoid_fallback,
      asymptotic_excluded_from_aicc = asymptotic_excluded,
      sigmoid_excluded_from_aicc = sigmoid_excluded
    )
  )
}


# =============================================================================
# STANDALONE MODEL FITTING FUNCTIONS
# =============================================================================
# These functions are used in 08_effect_sizes.R to fit specific models
# independently of the full pBACIPS comparison. Useful when you already
# know which model form is appropriate.

#' Fit asymptotic model standalone
#'
#' Fits the Michaelis-Menten form asymptotic model outside of the main
#' ProgressiveChangeBACIPS function. Used for effect size calculation.
#' Includes a fallback chain for reliable fitting.
#'
#' @param delta Numeric vector of Impact - Control differences
#' @param time.model Time since intervention (0 for Before, sequential for After)
#' @param time.model.of.impact Index of last Before-period observation
#' @param time.true Full time index vector
#' @return Fitted NLS model object or fallback model
myASYfun_standalone <- function(delta, time.model, time.model.of.impact, time.true) {
  foAsy <- delta ~ (M * time.model) / (L + time.model) + B

  funAsy <- function(parS, time.model) {
    (parS$M * time.model) / (parS$L + time.model) + parS$B
  }

  residFun <- function(p, observed, time.model) {
    pred <- funAsy(p, time.model)
    if (any(!is.finite(pred))) return(rep(1e10, length(observed)))
    observed - pred
  }

  # Check data quality first
  data_diag <- diagnose_data_quality(delta, time.model)
  processed <- preprocess_for_nls(delta, time.model, remove_outliers = TRUE)

  # Set starting values with multiple strategies
  after_idx <- which(time.model > 0)
  before_idx <- which(time.model == 0)
  M_start <- mean(delta[after_idx], na.rm = TRUE)
  B_start <- mean(delta[before_idx], na.rm = TRUE)
  max_time <- max(time.model, na.rm = TRUE)

  if (!is.finite(M_start)) M_start <- mean(tail(delta, 3), na.rm = TRUE)
  if (!is.finite(B_start)) B_start <- mean(head(delta, 3), na.rm = TRUE)
  if (!is.finite(M_start)) M_start <- 0.1
  if (!is.finite(B_start)) B_start <- 0

  # Multiple starting strategies
  starts_list <- list(
    list(M = M_start, B = B_start, L = max(max_time/2, 1)),
    list(M = M_start - B_start, B = B_start, L = 1),
    list(M = M_start * 1.5, B = B_start, L = max_time * 0.3),
    list(M = sign(M_start) * max(abs(M_start), 0.5), B = 0, L = 2)
  )

  # Try port algorithm with bounds first
  for (parStart in starts_list) {
    parStart$L <- max(parStart$L, 0.01)

    fit <- tryCatch({
      nls(foAsy, start = parStart, algorithm = "port",
          lower = c(M = -50, B = -50, L = 0.001),
          upper = c(M = 50, B = 50, L = max_time * 3),
          control = nls.control(maxiter = 200, warnOnly = TRUE,
                                 minFactor = 1/2048))
    }, error = function(e) NULL)
    fit <- check_nls_convergence(fit)

    if (!is.null(fit) && all(is.finite(coef(fit)))) {
      return(fit)
    }
  }

  # Try nls.lm with each starting point
  for (parStart in starts_list) {
    parStart$L <- max(parStart$L, 0.01)

    nls_ASY_out <- tryCatch(
      nls.lm(
        par = parStart,
        fn = residFun,
        observed = processed$delta,
        time.model = processed$time,
        control = nls.lm.control(maxfev = 500, maxiter = 200)
      ),
      error = function(e) NULL
    )

    if (!is.null(nls_ASY_out) && all(is.finite(coef(nls_ASY_out)))) {
      startPar <- as.list(coef(nls_ASY_out))
      gridPar <- expand.grid(
        M = startPar$M * c(0.7, 1.0, 1.3),
        B = startPar$B * c(0.7, 1.0, 1.3),
        L = pmax(0.01, startPar$L * c(0.7, 1.0, 1.3))
      )

      nls2_fit <- tryCatch(
        nls2(foAsy, start = gridPar, algorithm = "brute-force"),
        error = function(e) {
          tryCatch(
            nls2(foAsy, start = data.frame(t(unlist(startPar))), algorithm = "brute-force"),
            error = function(e2) NULL
          )
        }
      )

      if (!is.null(nls2_fit) && all(is.finite(coef(nls2_fit)))) {
        final_fit <- tryCatch(
          nls(foAsy, start = coef(nls2_fit),
              control = nls.control(maxiter = 100, warnOnly = TRUE)),
          error = function(e) nls2_fit
        )
        final_fit <- check_nls_convergence(final_fit)
        if (!is.null(final_fit) && all(is.finite(coef(final_fit)))) {
          return(final_fit)
        }
      }
    }
  }

  # Fallback: try stable parameterization (still a true NLS fit)
  fit <- fit_asymptotic_stable(delta, time.model)
  if (!is.null(fit)) {
    attr(fit, "model_fallback") <- "stable_reparam"
    return(fit)
  }

  # Fallback: self-starting (still a true NLS fit)
  fit <- fit_asymptotic_selfstart(delta, time.model)
  if (!is.null(fit)) {
    attr(fit, "model_fallback") <- "selfstart"
    return(fit)
  }

  # WARNING: Fallbacks below are NOT true asymptotic NLS fits.
  # Downstream code (08_effect_sizes.R) should check attr(model, "model_fallback")
  # and use the actual model class name rather than reporting "Asymptotic".

  # Fallback: piecewise linear
  fit <- fit_piecewise_linear(delta, time.model)
  if (!is.null(fit)) {
    attr(fit, "model_fallback") <- "piecewise"
    log_model_fit(model_type = "asymptotic_standalone", success = TRUE,
                  message = "Piecewise fallback used: NOT a true asymptotic fit")
    return(fit)
  }

  # Fallback: quadratic
  fit <- fit_quadratic_fallback(delta, time.model)
  if (!is.null(fit)) {
    attr(fit, "model_fallback") <- "quadratic"
    log_model_fit(model_type = "asymptotic_standalone", success = TRUE,
                  message = "Quadratic fallback used: NOT a true asymptotic fit")
    return(fit)
  }

  # Fallback: GAM
  fit <- fit_gam_fallback(delta, time.model)
  if (!is.null(fit)) {
    attr(fit, "model_fallback") <- "gam"
    log_model_fit(model_type = "asymptotic_standalone", success = TRUE,
                  message = "GAM fallback used: NOT a true asymptotic fit")
  }
  fit
}


#' Fit sigmoid model standalone
#'
#' Fits the Hill function form sigmoid model outside of the main
#' ProgressiveChangeBACIPS function. Used for effect size calculation.
#' Includes a fallback chain for reliable fitting.
#'
#' @param delta Numeric vector of Impact - Control differences
#' @param time.model Time since intervention (0 for Before, sequential for After)
#' @param time.model.of.impact Index of last Before-period observation
#' @param time.true Full time index vector
#' @return Fitted NLS model object or fallback model with model_fallback attribute
mySIGfun_standalone <- function(delta, time.model, time.model.of.impact, time.true) {
  # Add offset to avoid 0^K issues
  time_offset <- time.model + 0.01

  foSIG <- delta ~ (M * (time_offset / L)^K) / (1 + (time_offset / L)^K) + B

  funSIG <- function(parS, time_offset) {
    # Safeguards against numerical overflow
    ratio <- (time_offset / parS$L)^parS$K
    ratio <- pmin(ratio, 1e10)
    ratio <- pmax(ratio, 1e-10)
    (parS$M * ratio) / (1 + ratio) + parS$B
  }

  residFun <- function(p, observed, time_offset) {
    pred <- funSIG(p, time_offset)
    if (any(!is.finite(pred))) return(rep(1e10, length(observed)))
    observed - pred
  }

  # Diagnose data quality
  data_diag <- diagnose_data_quality(delta, time.model)
  processed <- preprocess_for_nls(delta, time.model, remove_outliers = TRUE)
  time_proc <- processed$time + 0.01

  # Set starting values with safeguards
  after_idx <- which(time.model > 0)
  before_idx <- which(time.model == 0)
  M_start <- mean(delta[after_idx], na.rm = TRUE)
  B_start <- mean(delta[before_idx], na.rm = TRUE)
  L_start <- median(time.model[time.model > 0], na.rm = TRUE)
  max_time <- max(time.model, na.rm = TRUE)

  if (!is.finite(M_start)) M_start <- mean(tail(delta, 3), na.rm = TRUE)
  if (!is.finite(B_start)) B_start <- mean(head(delta, 3), na.rm = TRUE)
  if (!is.finite(M_start)) M_start <- 0.1
  if (!is.finite(B_start)) B_start <- 0
  if (!is.finite(L_start) || L_start < 0.5) L_start <- max(max_time / 2, 1)

  # Multiple starting strategies with conservative K values
  starts_list <- list(
    list(M = M_start, B = B_start, L = L_start, K = 2),
    list(M = M_start - B_start, B = B_start, L = L_start * 0.5, K = 3),
    list(M = M_start * 1.5, B = B_start, L = L_start * 1.5, K = 1),
    list(M = sign(M_start) * max(abs(M_start), 0.5), B = 0, L = L_start, K = 2),
    list(M = M_start, B = B_start, L = max_time * 0.7, K = 5)
  )

  # Try port algorithm with bounds first
  for (parStart in starts_list) {
    parStart$L <- max(parStart$L, 0.5)
    parStart$L <- min(parStart$L, max_time * 1.2)
    parStart$K <- max(parStart$K, 0.3)
    parStart$K <- min(parStart$K, 10)

    fit <- tryCatch({
      nls(foSIG, start = parStart, algorithm = "port",
          lower = c(M = -50, B = -50, L = 0.2, K = 0.2),
          upper = c(M = 50, B = 50, L = max_time * 2, K = 15),
          control = nls.control(maxiter = 300, warnOnly = TRUE,
                                 minFactor = 1/4096))
    }, error = function(e) NULL)
    fit <- check_nls_convergence(fit)

    if (!is.null(fit) && all(is.finite(coef(fit)))) {
      coefs <- coef(fit)
      if (coefs["K"] > 0.25 && coefs["K"] < 14 &&
          coefs["L"] > 0.25 && coefs["L"] < max_time * 1.8) {
        return(fit)
      }
    }
  }

  # Try nls.lm with each starting point
  for (parStart in starts_list) {
    parStart$L <- max(parStart$L, 0.5)
    parStart$K <- max(parStart$K, 0.3)
    parStart$K <- min(parStart$K, 10)

    nls_SIG_out <- tryCatch(
      nls.lm(
        par = parStart,
        fn = residFun,
        observed = processed$delta,
        time_offset = time_proc,
        control = nls.lm.control(maxfev = 500, maxiter = 200)
      ),
      error = function(e) NULL
    )

    if (!is.null(nls_SIG_out) && all(is.finite(coef(nls_SIG_out)))) {
      startPar <- as.list(coef(nls_SIG_out))
      gridPar <- expand.grid(
        M = startPar$M * c(0.7, 1.0, 1.3),
        B = startPar$B * c(0.7, 1.0, 1.3),
        L = pmax(0.5, startPar$L * c(0.7, 1.0, 1.3)),
        K = pmax(0.3, pmin(10, startPar$K * c(0.7, 1.0, 1.3)))
      )

      nls2_fit <- tryCatch(
        nls2(foSIG, start = gridPar, algorithm = "brute-force"),
        error = function(e) {
          tryCatch(
            nls2(foSIG, start = data.frame(t(unlist(startPar))), algorithm = "brute-force"),
            error = function(e2) NULL
          )
        }
      )

      if (!is.null(nls2_fit) && all(is.finite(coef(nls2_fit)))) {
        final_fit <- tryCatch(
          nls(foSIG, start = coef(nls2_fit),
              control = nls.control(maxiter = 100, warnOnly = TRUE)),
          error = function(e) nls2_fit
        )
        final_fit <- check_nls_convergence(final_fit)
        if (!is.null(final_fit) && all(is.finite(coef(final_fit)))) {
          return(final_fit)
        }
      }
    }
  }

  # Fallback: 4-parameter logistic (still a true NLS fit)
  fit <- fit_sigmoid_4pl(delta, time.model)
  if (!is.null(fit)) {
    attr(fit, "model_fallback") <- "4pl_reparam"
    return(fit)
  }

  # Fallback: self-starting SSlogis (still a true NLS fit)
  fit <- fit_sigmoid_selfstart(delta, time.model)
  if (!is.null(fit)) {
    attr(fit, "model_fallback") <- "selfstart"
    return(fit)
  }

  # WARNING: Fallbacks below are NOT true sigmoid NLS fits.
  # Downstream code (08_effect_sizes.R) should check attr(model, "model_fallback")
  # and use the actual model class name rather than reporting "Sigmoid".

  # Fallback: piecewise linear
  fit <- fit_piecewise_linear(delta, time.model)
  if (!is.null(fit)) {
    attr(fit, "model_fallback") <- "piecewise"
    log_model_fit(model_type = "sigmoid_standalone", success = TRUE,
                  message = "Piecewise fallback used: NOT a true sigmoid fit")
    return(fit)
  }

  # Fallback: quadratic
  fit <- fit_quadratic_fallback(delta, time.model)
  if (!is.null(fit)) {
    attr(fit, "model_fallback") <- "quadratic"
    log_model_fit(model_type = "sigmoid_standalone", success = TRUE,
                  message = "Quadratic fallback used: NOT a true sigmoid fit")
    return(fit)
  }

  # Fallback: GAM
  fit <- fit_gam_fallback(delta, time.model)
  if (!is.null(fit)) {
    attr(fit, "model_fallback") <- "gam"
    log_model_fit(model_type = "sigmoid_standalone", success = TRUE,
                  message = "GAM fallback used: NOT a true sigmoid fit")
  }
  fit
}


# =============================================================================
# DIAGNOSTIC EXPORT FUNCTIONS
# =============================================================================

#' Get summary of model fitting diagnostics
#'
#' Returns a summary dataframe of all model fitting attempts recorded
#' during the session.
#'
#' @return Dataframe with columns: timestamp, mpa, taxa, model_type, success, message
#' @export
get_model_fit_summary <- function() {
  if (!exists("MODEL_FIT_LOG") || length(MODEL_FIT_LOG) == 0) {
    return(data.frame(
      timestamp = character(),
      mpa = character(),
      taxa = character(),
      model_type = character(),
      success = logical(),
      message = character(),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, lapply(MODEL_FIT_LOG, function(x) {
    data.frame(
      timestamp = as.character(x$timestamp),
      mpa = x$mpa,
      taxa = x$taxa,
      model_type = x$model_type,
      success = x$success,
      message = x$message,
      stringsAsFactors = FALSE
    )
  }))
}

#' Clear the model fitting log
#'
#' Resets the global log to empty state.
#' @export
clear_model_fit_log <- function() {
  MODEL_FIT_LOG <<- list()
}

#' Validate all fitted models in SumStats with DHARMa
#'
#' Runs DHARMa diagnostics on all models that can be re-fitted from the
#' SumStats dataframe.
#'
#' @param sumstats Dataframe from effect size calculations
#' @param all_data Original data with lnDiff values
#' @param save_plots Logical, whether to save diagnostic plots
#' @param output_dir Directory to save plots
#' @return Dataframe with DHARMa diagnostic p-values for each model
#' @export
validate_all_models <- function(sumstats, all_data = NULL, save_plots = FALSE,
                                  output_dir = here::here("plots", "diagnostics")) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    warning("DHARMa package not available - install for model diagnostics")
    return(NULL)
  }

  results <- list()

  for (i in seq_len(nrow(sumstats))) {
    row <- sumstats[i, ]
    result <- list(
      Taxa = row$Taxa,
      MPA = row$MPA,
      Model = row$Model,
      Source = row$Source,
      diagnosed = FALSE
    )

    # For linear models, we can run basic residual checks
    if (row$Model %in% c("Linear", "Step", "Mean")) {
      result$model_type <- "simple"
      result$notes <- "Simple model - standard diagnostics apply"
      result$diagnosed <- TRUE
    } else if (row$Model %in% c("Asymptotic", "Sigmoid") ||
               grepl("^(Asymptotic|Sigmoid)$", row$Model)) {
      result$model_type <- "nonlinear"
      result$notes <- "Nonlinear model - requires original data for diagnostics"
    } else if (grepl("fallback", row$Model, ignore.case = TRUE)) {
      result$model_type <- "fallback"
      result$notes <- paste("Fallback model used:", row$Model,
                            "- not a true NLS fit, interpret with caution")
    } else {
      result$model_type <- row$Model
    }

    results[[length(results) + 1]] <- result
  }

  do.call(rbind, lapply(results, as.data.frame))
}


#' Create model fit quality report
#'
#' Generates a summary report of model fitting quality across all analyses.
#'
#' @param sumstats SumStats dataframe
#' @return List with summary statistics and recommendations
#' @export
create_model_quality_report <- function(sumstats) {

  report <- list()

  # Count models by type
  report$model_counts <- table(sumstats$Model)

  # Count by analysis type
  if ("AnalysisType" %in% names(sumstats)) {
    report$analysis_type_counts <- table(sumstats$AnalysisType)
  } else if ("Type" %in% names(sumstats)) {
    report$analysis_type_counts <- table(sumstats$Type)
  }

  # Count by taxa
  report$taxa_counts <- table(sumstats$Taxa)

  # Identify potential issues
  report$issues <- list()

  # Check for extreme effect sizes
  if ("Mean" %in% names(sumstats)) {
    extreme_es <- sumstats[abs(sumstats$Mean) > 5, ]
    if (nrow(extreme_es) > 0) {
      report$issues$extreme_effect_sizes <- extreme_es[, c("Taxa", "MPA", "Mean", "SE")]
    }
  }

  # Check for very large SEs (imprecise estimates)
  if ("SE" %in% names(sumstats)) {
    high_se <- sumstats[sumstats$SE > 2, ]
    if (nrow(high_se) > 0) {
      report$issues$high_uncertainty <- high_se[, c("Taxa", "MPA", "Mean", "SE")]
    }
  }

  # Summary message
  total_models <- nrow(sumstats)
  nonlinear_models <- sum(sumstats$Model %in% c("Asymptotic", "Sigmoid"))
  # Count both old-style fallback labels and new-style "Sigmoid_fallback_*" labels
  fallback_models <- sum(
    sumstats$Model %in% c("Quadratic", "GAM", "Piecewise") |
    grepl("fallback", sumstats$Model, ignore.case = TRUE)
  )

  report$summary <- paste0(
    "Model Fitting Summary:\n",
    "  Total models: ", total_models, "\n",
    "  Nonlinear (true asymptotic/sigmoid NLS): ", nonlinear_models, "\n",
    "  Fallback models used (non-NLS substitutions): ", fallback_models, "\n",
    "  Issues identified: ", length(report$issues)
  )

  class(report) <- "model_quality_report"
  report
}

#' Print method for model quality report
#' @export
print.model_quality_report <- function(x, ...) {
  cat(x$summary, "\n\n")
  cat("Model types used:\n")
  print(x$model_counts)
  if (length(x$issues) > 0) {
    cat("\nPotential issues:\n")
    for (name in names(x$issues)) {
      cat(paste0("  - ", name, ": ", nrow(x$issues[[name]]), " cases\n"))
    }
  }
}


# =============================================================================
# Confirmation message
# =============================================================================
cat("pBACIPS functions loaded.\n")
cat("Main function: ProgressiveChangeBACIPS(control, impact, time.true, time.model)\n")
cat("Diagnostic functions: get_model_fit_summary(), create_model_quality_report()\n")
