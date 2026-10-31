# Functions for scripts 02 and 04. Sourcing this file only defines functions:
#   source("R/03_functions.R")
# Run R/01_data_overview.R once to prepare output/data/peak_aligned.rds.
# Forecasts are tracked in tracker units (log1p(counts) or counts); predictions,
# scores and expert weights use counts / count_scale.

# 1. SHARED SETTINGS -----------------------------------------------------
# Scripts 02 and 04 read these settings, so a change here applies to both.
multiqt_settings <- function() list(
  location_ids = c("US", "06", "48", "12", "36"),  # US, California, Texas, Florida, New York.
  horizons = -1:3,
  transform = "log1p",              # MultiQT tracks log1p(counts); "identity": raw counts.
  scale_method = "max",             # Scores/plots/weighting losses: "max", "median" or "none".
  expert_multipliers = c(0, .01, .015, .02, .03, .05, .075, .1, .15, .2, .3),
  warmup_rounds = 20L,              # First released rows that freeze the error scale S0.
  scale_floor = 1e-6,               # Guard on the error scale, not each rate.
  reference_multiplier = .05,       # Fixed MultiQT: eta = .05 * S0.
  H = 8L,                           # Lookback for learned lambda; sharing uses 1/(2H).
  lambda_start = 2.72,
  # No forecasts were issued 2025-10-25 to 2025-11-15, so a longer pre-peak window
  # misses weeks for locations with December/early-January peaks.
  weeks_before_peak = 3L,           # before = weeks -3..0, including the peak week.
  weeks_after_peak = 8L)            # after = weeks 1..8.

# 2. CASES, SCALES AND WINDOWS --------------------------------------------
# Prepared cases are named disease_location_horizon, e.g. covid_06_0.
get_case <- function(prepared, disease, location, horizon) {
  case <- prepared$cases[[paste(disease, location, horizon, sep = "_")]]
  if (is.null(case)) stop("No prepared case ", disease, " ", location, " h", horizon,
    ". Rerun Rscript R/01_data_overview.R.")
  case
}

to_tracker <- function(x, transform) {
  stopifnot(transform %in% c("log1p", "identity"))
  if (transform == "log1p") log1p(x) else x
}
from_tracker <- function(x, transform) {
  stopifnot(transform %in% c("log1p", "identity"))
  if (transform == "log1p") expm1(x) else x
}

# Count scale for scores, plots and weighting losses; MultiQT never tracks it.
# It uses only reports published before the location's first forecast.
choose_count_scale <- function(case, method = "max") {
  history <- case$scale_history[is.finite(case$scale_history)]
  scale <- 1
  used <- method
  if (method != "none" && length(history)) {
    scale <- if (method == "median") median(history) else max(history)
    if (scale <= 0) {
      scale <- max(history)
      used <- "max"
    }
  }
  if (!is.finite(scale) || scale <= 0 || !length(history)) {
    scale <- 1
    used <- "none"
  }
  list(scale = scale, used = used, method = method, reference_weeks = length(history))
}

# Frozen error scale S0: Q90 of the absolute base errors in tracker units, pooled
# over all quantile levels of the first `rounds` released rows. eta = m * S0.
# https://delphi-forecasting-reports.netlify.app/calibration_experiments/e01_eta_settings
# The whole release batch of the freeze week is consumed as warmup (no replay), but
# only its first `rounds` rows enter S0. With fewer releases every week is warmup.
error_scale <- function(case, rounds = 20L, floor = 1e-6, transform = "log1p") {
  feedback_by_week <- case$feedback_by_week
  week <- match(TRUE, cumsum(lengths(feedback_by_week)) >= rounds)
  initialized <- !is.na(week)
  warmup_weeks <- if (initialized) week else length(feedback_by_week)
  released <- Reduce(c, feedback_by_week[seq_len(warmup_weeks)], init = integer())
  rows <- head(released, rounds)
  raw_scale <- scale <- NA_real_
  if (initialized) {
    residuals <- sweep(to_tracker(case$base[rows, , drop = FALSE], transform), 1,
      to_tracker(case$feedback_observed[rows], transform), "-")
    raw_scale <- as.numeric(quantile(abs(residuals), .9, names = FALSE))
    scale <- max(floor, raw_scale)
  }
  list(initialized = initialized, warmup_weeks = warmup_weeks,
    date = if (initialized) case$reference_date[warmup_weeks] else as.Date(NA),
    rows = rows, consumed_rows = released, raw_scale = raw_scale, scale = scale,
    rounds = rounds, floor = floor)
}

# Peak-relative windows, as a logical vector over the case's weeks. The peak is
# the largest latest-reported count in the prepared winter (hindsight label only).
# before = weeks -3..0 and after = weeks 1..8 inside the winter; overall = all weeks.
phase_weeks <- function(case, phase, before = 3L, after = 8L) {
  week <- case$relative_week
  switch(phase,
    before = case$within_winter & week >= -before & week <= 0,
    after = case$within_winter & week >= 1 & week <= after,
    overall = rep(TRUE, length(week)),
    stop("Unknown phase ", phase))
}

# A maximum on the winter edge may be the tail of an earlier wave, and missing
# weeks shorten the before/after windows. Report both, but keep the location.
check_peak_window <- function(case, before = 3L, after = 8L) {
  label <- paste0(case$disease, " | ", case$location_name, " (", case$location,
    ") | horizon ", case$horizon)
  if (case$peak_boundary)
    warning(label, ": winter maximum on the edge of the winter window (",
      case$peak_date, "); location kept.", call. = FALSE, immediate. = TRUE)
  weeks <- seq.int(-before, after)
  i <- match(weeks, case$relative_week)
  ok <- !is.na(i)
  ok[ok] <- case$within_winter[i[ok]] & is.finite(case$truth_latest[i[ok]])
  if (!all(ok))
    warning(label, ": no scored target in peak-relative week(s) ",
      paste(weeks[!ok], collapse = ","), "; location kept.", call. = FALSE, immediate. = TRUE)
  invisible(all(ok) && !case$peak_boundary)
}

# 3. SCORING --------------------------------------------------------------
# Coverage and pinball loss of quantile forecasts q [week x level] against one
# truth per week. Pass counts with their count_scale to compare coverage on the
# count scale (exact integer ties stay ties); the loss is in counts / count_scale.
# The pinball loss penalizes underprediction by alpha, overprediction by 1 - alpha.
score_quantiles <- function(q, truth, alpha, count_scale = 1) {
  stopifnot(ncol(q) == length(alpha), length(truth) == nrow(q))
  error <- truth / count_scale - q / count_scale
  list(covered = truth <= q,
    loss = pmax(error, 0) * rep(alpha, each = nrow(q)) +
      pmax(-error, 0) * rep(1 - alpha, each = nrow(q)))
}

# Average scores over the given weeks. Calibration error (CE) is checked per level,
# |coverage - alpha|, and then averaged over levels; every level has equal weight.
summarize_scores <- function(scores, alpha, weeks) {
  if (!length(weeks)) {
    na <- rep(NA_real_, length(alpha))
    return(list(by_level = data.frame(alpha, coverage = na, coverage_error = na,
      quantile_loss = na), coverage_error = NA_real_, quantile_loss = NA_real_, n_weeks = 0L))
  }
  coverage <- colMeans(scores$covered[weeks, , drop = FALSE])
  loss <- colMeans(scores$loss[weeks, , drop = FALSE])
  by_level <- data.frame(alpha, coverage, coverage_error = abs(coverage - alpha),
    quantile_loss = loss)
  list(by_level = by_level, coverage_error = mean(by_level$coverage_error),
    quantile_loss = mean(loss), n_weeks = length(weeks))
}

# 4. STANDARD MULTIQT -----------------------------------------------------
# base: issued quantiles [week x level] in counts; feedback: first-use reported
# counts; feedback_by_week[[w]]: the earlier rows whose feedback is released at
# issue week w. eta: one rate, or one per week, in tracker units. The first
# warmup_weeks forecasts keep the base and their reports are never replayed.
# Optional leaky shift of all levels (leak = 0: standard MultiQT): each report adds
# leak * mean(alpha - covered); the shift is multiplied by decay (one value, or one
# per week) at every issue week and capped at +-leak_bound.
# Returns the adjusted forecasts in counts.
multiqt <- function(base, feedback, feedback_by_week, alpha, eta,
                    warmup_weeks = 0L, transform = "log1p",
                    leak = 0, decay = .5, leak_bound = .2) {
  n <- nrow(base)
  # Nothing to learn (eta = 0, or S0 never initialized): keep the exact base.
  if (warmup_weeks >= n || (isTRUE(all(eta == 0)) && leak == 0)) return(base)
  if (length(eta) == 1L) eta <- rep(eta, n)
  if (length(decay) == 1L) decay <- rep(decay, n)
  base_tracker <- to_tracker(base, transform)
  feedback_tracker <- to_tracker(feedback, transform)
  stopifnot(ncol(base) == length(alpha), length(feedback) == n,
    length(feedback_by_week) == n, length(eta) == n, all(is.finite(eta)), all(eta >= 0),
    all(is.finite(base_tracker)), all(is.finite(feedback_tracker[unlist(feedback_by_week)])))
  offsets <- numeric(length(alpha))
  shift <- 0
  q_tracker <- base_tracker
  for (week in seq_len(n)) {
    if (week <= warmup_weeks) next
    shift <- shift * decay[week]
    for (past_week in feedback_by_week[[week]]) {
      # Compare with the HISTORICAL ISSUED forecast in tracker units.
      # At alpha=.9, eta=.01, a miss adds .009; coverage subtracts .001.
      covered <- feedback_tracker[past_week] <= q_tracker[past_week, ]
      offsets <- offsets + eta[week] * (alpha - covered)
      shift <- shift + leak * mean(alpha - covered)
    }
    shift <- min(leak_bound, max(-leak_bound, shift))
    # Project only the issued forecasts, NOT the offsets ("lazy" MultiQT).
    # isoreg restores quantile order: e.g., (12,10) becomes (11,11).
    q_tracker[week, ] <- isoreg(base_tracker[week, ] + offsets + shift)$yf
  }
  # Do not clip negative predictions: that would change the issued coverage indicators.
  q <- from_tracker(q_tracker, transform)
  # Keep unchanged bases exact, so transformation roundoff cannot flip ties.
  unchanged <- q_tracker == base_tracker
  q[unchanged] <- base[unchanged]
  q
}

# 5. EXPERT CREATORS ------------------------------------------------------
# All creators return the same expert object:
#   forecasts  [week x expert x level] in counts / count_scale, NA before birth;
#   birth_week the first issue week of each expert;
#   experts    metadata: expert name, eta, start date (and multiplier).
# scale is the error_scale() result; eta = multiplier * S0.

# One MultiQT expert per rate multiplier. All experts issue from the first week,
# keep the base during the S0 warmup and learn from their own issued forecasts.
make_rate_experts <- function(case, multipliers, scale, count_scale = 1, transform = "log1p") {
  n <- nrow(case$base)
  eta <- multipliers * scale$scale
  eta[multipliers == 0] <- 0
  names <- paste0("m=", multipliers)
  forecasts <- array(NA_real_, c(n, length(multipliers), length(case$alpha)),
    dimnames = list(as.character(case$reference_date), names, as.character(case$alpha)))
  for (j in seq_along(multipliers))
    forecasts[, j, ] <- multiqt(case$base, case$feedback_observed, case$feedback_by_week,
      case$alpha, eta[j], scale$warmup_weeks, transform) / count_scale
  list(forecasts = forecasts, birth_week = rep(1L, length(multipliers)),
    experts = data.frame(expert = names, multiplier = multipliers, eta,
      start_date = case$reference_date[1]))
}

# Restart experts: after the S0 warmup, start one fresh MultiQT copy per multiplier
# at every issue week. A copy forecasts the base at birth and learns only from
# reports on its own forecasts; delayed reports on earlier rows are ignored.
make_restart_experts <- function(case, multipliers, scale, count_scale = 1, transform = "log1p") {
  n <- nrow(case$base)
  births <- if (scale$warmup_weeks < n) seq.int(scale$warmup_weeks + 1L, n) else integer()
  if (length(births)) stopifnot(is.finite(scale$scale), scale$scale > 0)
  birth_week <- rep(births, each = length(multipliers))
  multiplier <- rep(multipliers, times = length(births))
  names <- if (length(birth_week))
    paste0("start=", case$reference_date[birth_week], ";m=", multiplier) else character()
  forecasts <- array(NA_real_, c(n, length(birth_week), length(case$alpha)),
    dimnames = list(as.character(case$reference_date), names, as.character(case$alpha)))
  for (j in seq_along(birth_week)) {
    issued <- seq.int(birth_week[j], n)
    # Renumber the release schedule from the copy's birth; drop rows before it.
    schedule <- lapply(case$feedback_by_week[issued], function(rows)
      rows[rows >= issued[1]] - issued[1] + 1L)
    forecasts[issued, j, ] <- multiqt(case$base[issued, , drop = FALSE],
      case$feedback_observed[issued], schedule, case$alpha,
      multiplier[j] * scale$scale, transform = transform) / count_scale
  }
  list(forecasts = forecasts, birth_week = birth_week,
    experts = data.frame(expert = names, multiplier, eta = multiplier * scale$scale,
      start_date = case$reference_date[birth_week]))
}

# Diverse experts: a small bank whose experts differ in where their correction
# comes from. Every expert except the base runs a slow MultiQT (eta = slow, tracker
# units) on its own starting forecast, so each stays calibrated in the long run.
# All experts issue the base during the S0 warmup.
#   base     the original forecast;
#   slow     slow MultiQT on the base;
#   leaky    slow MultiQT plus a leaky shift of all levels (see multiqt()), halved
#            every week since the previous issue;
#   damp=k   start from the base minus k times its implied growth: base median minus
#            the latest released report, in tracker units, clipped to [-1, 1]
#            (k = 1 centres on the latest report, k < 0 extrapolates further);
#   width=w  start from the base quantiles moved w times as far from the base median.
make_diverse_experts <- function(case, scale, count_scale = 1, transform = "log1p",
                                 slow = .01, leak = .2, damping = c(-.5, .5, 1),
                                 widths = c(.75, 1.5)) {
  n <- nrow(case$base)
  base_tracker <- to_tracker(case$base, transform)
  mid <- base_tracker[, match(.5, case$alpha)]
  # Latest row whose first report is released by each issue week (none: no growth).
  latest <- cummax(vapply(case$feedback_by_week, function(rows) max(0L, rows), integer(1)))
  implied <- mid - to_tracker(case$feedback_observed[ifelse(latest > 0, latest, NA)], transform)
  implied <- pmin(1, pmax(-1, ifelse(is.finite(implied), implied, 0)))
  # Starting forecasts in counts, floored at zero counts.
  starts <- list(slow = case$base, leaky = case$base)
  for (k in damping) starts[[paste0("damp=", k)]] <-
    from_tracker(pmax(base_tracker - k * implied, 0), transform)
  for (w in widths) starts[[paste0("width=", w)]] <-
    from_tracker(pmax(mid + w * (base_tracker - mid), 0), transform)
  decay <- .5^c(1, pmax(1, as.numeric(diff(case$reference_date)) / 7))
  warmup <- seq_len(n) <= scale$warmup_weeks
  names <- c("base", names(starts))
  forecasts <- array(NA_real_, c(n, length(names), length(case$alpha)),
    dimnames = list(as.character(case$reference_date), names, as.character(case$alpha)))
  forecasts[, "base", ] <- case$base / count_scale
  for (name in names(starts)) {
    start <- starts[[name]]
    start[warmup, ] <- case$base[warmup, ]
    forecasts[, name, ] <- multiqt(start, case$feedback_observed, case$feedback_by_week,
      case$alpha, slow, scale$warmup_weeks, transform,
      leak = if (name == "leaky") leak else 0, decay = decay) / count_scale
  }
  list(forecasts = forecasts, birth_week = rep(1L, length(names)),
    experts = data.frame(expert = names, eta = c(0, rep(slow, length(starts))),
      start_date = case$reference_date[1]))
}

# Keep some experts of an expert object, e.g. the restart copies with one multiplier.
select_experts <- function(experts, keep) list(
  forecasts = experts$forecasts[, keep, , drop = FALSE],
  birth_week = experts$birth_week[keep],
  experts = experts$experts[keep, , drop = FALSE])

# 6. EXPONENTIAL WEIGHTING OF EXPERTS --------------------------------------
# Gibbs and Candes (2024), Algorithms 1-2, motivate exponential weighting and sharing.
# https://jmlr.org/papers/v25/22-1218.html
# https://github.com/isgibbs/DtACI/blob/main/DtACI.R
# experts: object from make_rate_experts() or make_restart_experts().
# base, feedback: counts / count_scale; the weights learn from pinball losses in
#   these units, released by feedback_by_week like the MultiQT updates.
# lambda: "learned" starts at lambda_start and, after H updates, uses
#   sqrt((log(2 K H) + 1) / sum of the previous H second moments); a number fixes
#   lambda; 0 gives the uniform average of the active experts.
# sharing: after each update, share (default 1/(2H)) of the mass is spread uniformly.
# per_quantile: one weight vector per level. FALSE: one weight vector for all levels,
#   learned from the mean pinball loss over levels, with scalar lambda and moment.
#   A common-weight average of ordered quantiles stays ordered, so it is not repaired.
# Later starts: an expert joins at its birth week with 1/K of the weight among the K
#   active experts; the others keep their relative weights. A released loss of row t
#   updates only experts that issued row t; the others keep their mass.
# The first warmup_weeks forecasts keep the base and their reports are not replayed.
# Returns the mixture q [week x level] in counts / count_scale, the issued weights
# [week x expert x level, or x 1 if common] and the lambda used at each issue week.
exp_weights <- function(experts, base, feedback, feedback_by_week, alpha,
                        lambda = "learned", sharing = TRUE, H = 8L, share = 1 / (2 * H),
                        per_quantile = TRUE, lambda_start = 2.72, warmup_weeks = 0L) {
  forecasts <- experts$forecasts
  birth <- experts$birth_week
  n <- nrow(base); k <- length(birth); a <- length(alpha)
  learned <- identical(lambda, "learned")
  stopifnot(identical(dim(forecasts), c(n, k, a)), length(feedback) == n,
    length(feedback_by_week) == n, learned || (is.numeric(lambda) &&
      length(lambda) == 1L && is.finite(lambda) && lambda >= 0),
    !sharing || (is.finite(share) && share > 0 && share <= 1))
  learn <- learned || lambda > 0
  columns <- if (per_quantile) a else 1L     # One weight vector per level, or one for all.
  learning_rate <- rep(if (learned) lambda_start else lambda, columns)
  p <- matrix(0, k, columns)                 # Weights; zero before an expert's birth.
  log_weights <- matrix(-Inf, k, columns)    # Kept across updates, even if exp() underflows.
  moments <- matrix(NA_real_, length(unlist(feedback_by_week)), columns)
  n_updates <- 0L
  weights <- array(0, c(n, k, columns), dimnames = list(dimnames(forecasts)[[1]],
    dimnames(forecasts)[[2]], if (per_quantile) dimnames(forecasts)[[3]] else "all levels"))
  lambda_used <- matrix(NA_real_, n, columns)
  log_mass <- function(x) apply(x, 2, function(v) max(v) + log(sum(exp(v - max(v)))))
  q <- base
  for (week in seq_len(n)) {
    admitted <- which(birth < week)
    if (learn && week > warmup_weeks) for (past in feedback_by_week[[week]]) {
      eligible <- admitted[birth[admitted] <= past]
      if (length(eligible) < 2L) next          # A single expert keeps all of its mass.
      sleeping <- setdiff(admitted, eligible)
      if (learned && n_updates >= H) {
        # Use the PREVIOUS H loss moments, excluding the feedback being processed now.
        recent <- seq.int(n_updates - H + 1L, n_updates)
        denominator <- colSums(moments[recent, , drop = FALSE])
        positive <- denominator > 0
        constant <- log(2 * length(eligible) * H) + 1
        learning_rate[positive] <- sqrt(constant / denominator[positive])
        # Tiny positive moments can overflow the ratio before its square root.
        overflow <- positive & !is.finite(learning_rate)
        learning_rate[overflow] <- sqrt(constant) / sqrt(denominator[overflow])
        # If recent losses are all zero, keep the previous finite learning rate.
      }
      error <- feedback[past] - matrix(forecasts[past, eligible, ], length(eligible), a)
      loss <- pmax(error, 0) * rep(alpha, each = length(eligible)) +
        pmax(-error, 0) * rep(1 - alpha, each = length(eligible))
      if (!per_quantile) loss <- matrix(rowMeans(loss), ncol = 1L)
      stopifnot(all(is.finite(loss)))
      # Experts without a forecast for this row sleep and keep their mass.
      mass <- 1 - colSums(p[sleeping, , drop = FALSE])
      # Paper's second moment: average squared expert losses using pre-update weights.
      n_updates <- n_updates + 1L
      moments[n_updates, ] <- colSums(sweep(p[eligible, , drop = FALSE], 2, mass, `/`) * loss^2)
      updated <- log_weights[eligible, , drop = FALSE] - sweep(loss, 2, learning_rate, `*`)
      if (length(sleeping)) {
        # Shift the eligible experts' log weights so that their total mass is unchanged.
        updated <- sweep(updated, 2, log_mass(log_weights[eligible, , drop = FALSE]) -
          log_mass(updated), `+`)
      }
      log_weights[eligible, ] <- updated
      # Exponential weights, stabilized and normalized separately for each column.
      log_weights[admitted, ] <- sweep(log_weights[admitted, , drop = FALSE], 2,
        apply(log_weights[admitted, , drop = FALSE], 2, max), `-`)
      p[admitted, ] <- exp(log_weights[admitted, , drop = FALSE])
      p[admitted, ] <- sweep(p[admitted, , drop = FALSE], 2,
        colSums(p[admitted, , drop = FALSE]), `/`)
      if (sharing) {
        p[eligible, ] <- (1 - share) * p[eligible, , drop = FALSE] +
          rep(share * mass / length(eligible), each = length(eligible))
        log_weights[admitted, ] <- log(p[admitted, , drop = FALSE])
      }
    }
    newborn <- which(birth == week)
    active <- c(admitted, newborn)
    if (!learn) {
      p[active, ] <- 1 / length(active)
    } else if (length(newborn)) {
      # Each newborn gets 1/K among K active experts, before issuing its first forecast.
      # Preserve the incumbents' relative weights after processing released losses.
      new_mass <- length(newborn) / length(active)
      p <- p * (1 - new_mass)
      p[newborn, ] <- new_mass / length(newborn)
      log_weights[active, ] <- log(p[active, , drop = FALSE])
    }
    if (length(active) && week > warmup_weeks) {
      stopifnot(all(is.finite(p)), max(abs(colSums(p) - 1)) < 1e-10)
      current <- matrix(forecasts[week, active, ], length(active), a)
      # Weight first, then repair crossings of per-level mixtures. Keep an all-base mixture exact.
      if (all(current == rep(base[week, ], each = length(active)))) {
        q[week, ] <- base[week, ]
      } else {
        mixture <- colSums(current * p[active, rep(seq_len(columns), length.out = a), drop = FALSE])
        q[week, ] <- if (per_quantile) isoreg(mixture)$yf else mixture
      }
    }
    weights[week, , ] <- p
    lambda_used[week, ] <- learning_rate
  }
  list(q = q, weights = weights, lambda = lambda_used)
}
