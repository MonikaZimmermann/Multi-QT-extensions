# Functions for scripts 02 and 04. Sourcing this file only defines functions:
#   source("R/03_functions.R")
# Run R/01_data_overview.R once to prepare output/data/peak_aligned.rds.
# Forecasts are tracked in tracker units (log1p(counts) or counts); predictions,
# scores and expert weights use counts / count_scale.

# 1. SHARED SETTINGS -----------------------------------------------------
# Scripts 02 and 04 read these settings, so a change here applies to both.
multiqt_settings <- function() list(
  location_ids = c("US", "06", "48", "36", "42", "17"),  # US, California, Texas, New York, Pennsylvania, Illinois.
  horizons = -1:3,
  transform = "identity",              # MultiQT tracks log1p(counts); "identity": raw counts.
  scale_method = "max",             # Scores/plots/weighting losses: "max", "median" or "none".
  expert_multipliers = c(0, .01, .015, .02, .03, .05, .075, .1, .15, .2, .3),
  warmup_rounds = 20L,              # First released rows that freeze the error scale S0.
  replay_warmup = TRUE,             # TRUE: learn from the first week on, S0 from those rows in hindsight.
  feedback_delay = 2L,              # Weeks after the first report before MultiQT learns from a forecast.
  score_target = "latest",          # Scores against the "latest" report or the "feedback" learned from.
  score_filled = TRUE,              # Score filled rounds (fill_gaps_up_to in R/01_data_overview.R).
  scale_floor = 1e-6,               # Guard on S0, sigma_t and the level scales, not each rate.
  reference_multiplier = .05,       # Fixed MultiQT: eta = .05 * S0.
  H = 8L,                           # Lookback for learned lambda; sharing uses 1/(2H).
  lambda_start = 2.72,
  # Window lengths in weeks. No forecasts were issued 2025-10-25 to 2025-11-15, so
  # COVID after-peak windows of September 2025 peaks miss some weeks.
  weeks_before_peak = 8L,           # before = weeks -7..0: 8 weeks, including the peak week.
  weeks_after_peak = 8L)            # after = weeks 1..8: 8 weeks.

# 2. CASES, SCALES AND WINDOWS --------------------------------------------
# Prepared cases are named disease_location_horizon, e.g. covid_06_0.
# feedback_delay selects the online feedback: each forecast is learned from the
# report available feedback_delay weeks after its first report (0: first report).
get_case <- function(prepared, disease, location, horizon, feedback_delay = 0L) {
  case <- prepared$cases[[paste(disease, location, horizon, sep = "_")]]
  if (is.null(case)) stop("No prepared case ", disease, " ", location, " h", horizon,
    ". Rerun Rscript R/01_data_overview.R.")
  feedback <- case$feedback_by_delay[[as.character(feedback_delay)]]
  if (is.null(feedback)) stop("No feedback with a ", feedback_delay, "-week delay. Add it to ",
    "feedback_delays in R/01_data_overview.R and rerun it.")
  case$feedback_delay <- feedback_delay
  case$feedback_observed <- feedback$observed
  case$feedback_date <- feedback$date
  case$feedback_by_week <- feedback$by_week
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
# replay = TRUE keeps S0 but has no warmup weeks: MultiQT learns from the first
# week on, so S0 uses reports released later (hindsight, acceptable here).
error_scale <- function(case, rounds = 20L, floor = 1e-6, transform = "log1p",
                        replay = FALSE) {
  feedback_by_week <- case$feedback_by_week
  week <- match(TRUE, cumsum(lengths(feedback_by_week)) >= rounds)
  initialized <- !is.na(week)
  freeze_week <- if (initialized) week else length(feedback_by_week)
  released <- Reduce(c, feedback_by_week[seq_len(freeze_week)], init = integer())
  rows <- head(released, rounds)
  raw_scale <- scale <- NA_real_
  if (initialized) {
    residuals <- sweep(to_tracker(case$base[rows, , drop = FALSE], transform), 1,
      to_tracker(case$feedback_observed[rows], transform), "-")
    raw_scale <- as.numeric(quantile(abs(residuals), .9, names = FALSE))
    scale <- max(floor, raw_scale)
  }
  warmup_weeks <- if (replay && initialized) 0L else freeze_week
  list(initialized = initialized, warmup_weeks = warmup_weeks,
    date = if (initialized) case$reference_date[freeze_week] else as.Date(NA),
    rows = rows, consumed_rows = if (warmup_weeks > 0L) released else integer(),
    raw_scale = raw_scale, scale = scale, rounds = rounds, floor = floor)
}

# Scale sigma_t of the median residuals y - b (tracker units) for the studentized
# median MultiQT: one estimate of their standard deviation per issue week.
#   "residuals" root mean square of the last `window` median residuals released by
#               issue week t (only released reports, no hindsight). NA until min_rows
#               are released. A small window adapts quickly but is noisier. Reports
#               arrive weeks after issue, so this scale still reacts to a wave late.
#   "base"      the sd implied by the issued base: the width between its quantile `levels`
#               divided by that of a standard normal, e.g. (q.9 - q.1) / 2.563 for the 80%
#               interval. No feedback lag; a point-mass base gives the floor.
# warmup_weeks: the weeks before the first sigma_t; the method issues the original base there.
median_sigma <- function(case, method = "residuals", window = 5L, min_rows = 5L,
                         floor = 1e-6, transform = "log1p", levels = c(.1, .9)) {
  stopifnot(window >= 1L, min_rows >= 1L)   # Otherwise a later report enters sigma_t.
  base_level <- function(alpha) to_tracker(case$base[, match(alpha, case$alpha)], transform)
  scale <- switch(method,
    residuals = {
      released <- unlist(case$feedback_by_week)              # Rows in release order.
      residual <- to_tracker(case$feedback_observed[released], transform) -
        base_level(.5)[released]
      n_released <- cumsum(lengths(case$feedback_by_week))   # By each issue week.
      vapply(n_released, function(k) if (k < min_rows) NA_real_ else
        max(floor, sqrt(mean(residual[max(1L, k - window + 1L):k]^2))), numeric(1))
    },
    base = {
      stopifnot(length(levels) == 2L, levels[1] < levels[2], all(levels %in% case$alpha))
      pmax(floor, (base_level(levels[2]) - base_level(levels[1])) / diff(qnorm(levels)))
    },
    stop("Unknown median sigma ", method))
  first <- match(TRUE, is.finite(scale))
  list(scale = scale, warmup_weeks = if (is.na(first)) length(scale) else first - 1L,
    method = method, window = window, min_rows = min_rows, floor = floor)
}

# Quantile-specific scales sigma_t^alpha (tracker units) for the studentized MultiQT on the
# full base: the normal sd implied by the base's central interval between levels alpha and
# 1 - alpha (at the median, the nearest pair). Symmetric in alpha and wider in heavier
# tails. Returns a matrix [week x level]. Needs levels symmetric about .5.
# Level scales often violate the research note's (sufficient) inward-flow condition, so
# MultiQT's long-run calibration guarantee covers a common scale, not these scales.
level_sigma <- function(case, floor = 1e-6, transform = "log1p") {
  q <- to_tracker(case$base, transform)
  m <- ncol(q)
  stopifnot(isTRUE(all.equal(case$alpha, 1 - rev(case$alpha))))
  lower <- pmin(seq_len(m), m + 1L - seq_len(m))     # Mirror pairs by index, not 1 - alpha.
  lower[lower == m + 1L - lower] <- lower[lower == m + 1L - lower] - 1L
  upper <- m + 1L - lower
  z <- qnorm(case$alpha)
  scale <- pmax((q[, upper, drop = FALSE] - q[, lower, drop = FALSE]) /
    rep(z[upper] - z[lower], each = nrow(q)), floor)
  dimnames(scale) <- dimnames(q)
  scale
}

# Quantile-specific scales sigma_t^alpha from recent residuals (tracker units): for every
# level, the RMS of y_s - b_s^alpha over the last `window` rows s released by issue week t
# (no hindsight). Rows are NA until min_rows are released, like median_sigma("residuals").
# Returns a matrix [week x level]. Like level_sigma(), no calibration guarantee.
level_residual_sigma <- function(case, window = 5L, min_rows = 5L, floor = 1e-6,
                                 transform = "log1p") {
  stopifnot(window >= 1L, min_rows >= 1L)   # Otherwise a later report enters sigma_t.
  released <- unlist(case$feedback_by_week)                # Rows in release order.
  residual <- to_tracker(case$feedback_observed[released], transform) -
    to_tracker(case$base[released, , drop = FALSE], transform)
  n_released <- cumsum(lengths(case$feedback_by_week))     # By each issue week.
  scale <- t(vapply(n_released, function(k) if (k < min_rows) rep(NA_real_, ncol(residual)) else
    pmax(sqrt(colMeans(residual[max(1L, k - window + 1L):k, , drop = FALSE]^2)), floor),
    numeric(ncol(residual))))
  dimnames(scale) <- dimnames(case$base)
  scale
}

# Peak-relative windows, as a logical vector over the case's weeks. The peak is
# the largest latest-reported count in the disease's peak window (hindsight label only).
# before and after are window lengths: before = weeks -(before - 1)..0, including the
# peak week, and after = weeks 1..after, inside the peak window; overall = all weeks.
phase_weeks <- function(case, phase, before = 8L, after = 8L) {
  week <- case$relative_week
  switch(phase,
    before = case$within_peak_window & week > -before & week <= 0,
    after = case$within_peak_window & week >= 1 & week <= after,
    overall = rep(TRUE, length(week)),
    stop("Unknown phase ", phase))
}

# A maximum on the window edge may be the tail of an earlier wave, and missing
# weeks shorten the before/after windows. Report both, but keep the location.
check_peak_window <- function(case, before = 8L, after = 8L, truth = case$truth_latest) {
  label <- paste0(case$disease, " | ", case$location_name, " (", case$location,
    ") | horizon ", case$horizon)
  if (case$peak_boundary)
    warning(label, ": maximum on the edge of the peak window (",
      case$peak_date, "); location kept.", call. = FALSE, immediate. = TRUE)
  weeks <- seq.int(1L - before, after)
  i <- match(weeks, case$relative_week)
  ok <- !is.na(i)
  ok[ok] <- case$within_peak_window[i[ok]] & is.finite(truth[i[ok]])
  if (!all(ok))
    warning(label, ": no scored target in peak-relative week(s) ",
      paste(weeks[!ok], collapse = ","), "; location kept.", call. = FALSE, immediate. = TRUE)
  invisible(all(ok) && !case$peak_boundary)
}

# 3. SCORING --------------------------------------------------------------
# Truth for coverage, coverage residuals and quantile losses. "latest": the latest
# stored report of each target week. "feedback": the report the online updates
# learned from at the case's feedback delay (NA for weeks never released).
# filled = FALSE leaves filled rounds unscored; MultiQT still issues and learns there.
score_truth <- function(case, target = "latest", filled = TRUE) {
  truth <- switch(target,
    latest = case$truth_latest,
    feedback = case$feedback_observed,
    stop("Unknown score target ", target))
  if (!filled) truth[case$filled] <- NA
  truth
}
score_truth_label <- function(target, delay = 0L) switch(target,
  latest = "latest reported counts",
  feedback = paste0("reports ", delay, " week(s) after the first report"),
  stop("Unknown score target ", target))

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

# 4. MULTIQT: STANDARD AND STUDENTIZED -----------------------------------
# Weighted isotonic regression (pool adjacent violators): the increasing q that
# minimizes sum(w * (q - v)^2). Equal weights give isoreg().
weighted_isoreg <- function(v, w) {
  value <- weight <- numeric(length(v))
  size <- integer(length(v))
  k <- 0L
  for (i in seq_along(v)) {
    k <- k + 1L
    value[k] <- v[i]; weight[k] <- w[i]; size[k] <- 1L
    while (k > 1L && value[k - 1L] > value[k]) {         # Pool the violating neighbours.
      total <- weight[k - 1L] + weight[k]
      value[k - 1L] <- (weight[k - 1L] * value[k - 1L] + weight[k] * value[k]) / total
      weight[k - 1L] <- total; size[k - 1L] <- size[k - 1L] + size[k]
      k <- k - 1L
    }
  }
  rep(value[seq_len(k)], size[seq_len(k)])
}

# base: issued quantiles [week x level] in counts; feedback: reported counts at
# the case's feedback delay; feedback_by_week[[w]]: the earlier rows whose feedback
# is released at issue week w. eta: one rate, or one per week, in tracker units. The first
# warmup_weeks forecasts keep the base and their reports are never replayed.
# Optional leaky shift of all levels (leak = 0: standard MultiQT): each report adds
# leak * mean(alpha - covered); the shift is multiplied by decay (one value, or one
# per week) at every issue week and capped at +-leak_bound.
# Optional studentized offsets (sigma = 1, start = 0: standard MultiQT): week w issues
# base + sigma[w] * offsets. sigma: one value, or one per week, in tracker units; then
# offsets and eta are in units of sigma, and a change of sigma rescales the whole
# accumulated offset (the leaky shift stays in tracker units). A matrix sigma
# [week x level] gives each level its own scale; the projection is then weighted
# isotonic regression with weights 1 / sigma^2 (Euclidean in units of sigma).
# start: the offsets before the first update, one value or one per level.
# Returns the adjusted forecasts in counts.
multiqt <- function(base, feedback, feedback_by_week, alpha, eta,
                    warmup_weeks = 0L, transform = "log1p",
                    leak = 0, decay = .5, leak_bound = .2, sigma = 1, start = 0) {
  n <- nrow(base)
  # Only warmup weeks (e.g. S0 never initialized), or nothing to learn and start = 0:
  # keep the exact base.
  if (warmup_weeks >= n || (isTRUE(all(eta == 0)) && leak == 0 && isTRUE(all(start == 0))))
    return(base)
  if (length(eta) == 1L) eta <- rep(eta, n)
  if (length(decay) == 1L) decay <- rep(decay, n)
  if (length(sigma) == 1L) sigma <- rep(sigma, n)
  by_level <- is.matrix(sigma)              # One scale per level: weighted projection.
  base_tracker <- to_tracker(base, transform)
  feedback_tracker <- to_tracker(feedback, transform)
  warmup <- seq_len(n) <= warmup_weeks      # sigma may be missing in warmup weeks.
  issued_sigma <- if (by_level) sigma[!warmup, , drop = FALSE] else sigma[!warmup]
  stopifnot(ncol(base) == length(alpha), length(feedback) == n,
    length(feedback_by_week) == n, length(eta) == n, all(is.finite(eta)), all(eta >= 0),
    NROW(sigma) == n, !by_level || ncol(sigma) == length(alpha),
    all(is.finite(issued_sigma) & issued_sigma > 0),
    length(start) %in% c(1L, length(alpha)), all(is.finite(start)),
    all(is.finite(base_tracker)), all(is.finite(feedback_tracker[unlist(feedback_by_week)])))
  offsets <- rep_len(start, length(alpha))
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
    # Weights (min sigma / sigma)^2: the same projection as 1 / sigma^2, without overflow.
    q_tracker[week, ] <- if (by_level) weighted_isoreg(base_tracker[week, ] +
      sigma[week, ] * offsets + shift, (min(sigma[week, ]) / sigma[week, ])^2) else
      isoreg(base_tracker[week, ] + sigma[week] * offsets + shift)$yf
  }
  # Do not clip negative predictions: that would change the issued coverage indicators.
  q <- from_tracker(q_tracker, transform)
  # Keep unchanged bases exact, so transformation roundoff cannot flip ties.
  unchanged <- q_tracker == base_tracker
  q[unchanged] <- base[unchanged]
  q
}

# Studentized median MultiQT: the base is reduced to its median b_t, and week t issues
# b_t + sigma_t * offsets at every level (tracker units; isoreg restores the order).
# Coverage and update are those of multiqt(): 1{y_s <= q_s} = 1{(y_s - b_s) / sigma_s <=
# issued offset_s}, so offsets and eta are in units of sigma. sigma: one value, or
# sigma_t per week (median_sigma()). start: the offsets before the first update, e.g.
# qnorm(alpha), a normal forecast with sd sigma_t, or qt(alpha, window) for sigma_t from
# the RMS of `window` residuals. The first warmup_weeks forecasts keep the original base,
# like every method; reports on them released later are compared with it.
# Returns the adjusted forecasts in counts.
median_multiqt <- function(case, eta, sigma, start, warmup_weeks = 0L, transform = "log1p") {
  median_base <- case$base
  median_base[] <- case$base[, match(.5, case$alpha)]   # Every level: the base median.
  warmup <- seq_len(nrow(median_base)) <= warmup_weeks
  median_base[warmup, ] <- case$base[warmup, ]
  multiqt(median_base, case$feedback_observed, case$feedback_by_week, case$alpha, eta,
    warmup_weeks, transform, sigma = sigma, start = start)
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
  # First reports at any feedback delay: the starting point uses what is known then.
  first <- case$feedback_by_delay[["0"]]
  latest <- cummax(vapply(first$by_week, function(rows) max(0L, rows), integer(1)))
  implied <- mid - to_tracker(first$observed[ifelse(latest > 0, latest, NA)], transform)
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
