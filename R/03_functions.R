# Functions for scripts 02 and 04. Sourcing this file only defines functions:
#   source("R/03_functions.R")
# Run R/01_data_overview.R once to prepare output/data/peak_aligned.rds.
# Forecasts are tracked in tracker units (log1p(counts) or counts; Delphi's production
# trackers always on the sqrt scale); predictions
# and scores use counts / count_scale (per 100k residents by default); the expert
# weights learn from losses relative to the base's (weighting_loss_scale()).

# 1. SHARED SETTINGS -----------------------------------------------------
# Scripts 02, 04 and 06 read these settings, so a change here applies to all three.
multiqt_settings <- function() list(
  location_ids = c("US", "06", "48", "36", "42"),  # US, California, Texas, New York, Pennsylvania.
  summary_exclude = "US",           # Left out of the pooled all-location summaries (02, 04, 06).
  horizons = -1:3,
  transform = "identity",           # MultiQT tracks raw counts; "log1p": log1p(counts).
  scale_method = "population",      # Scores/plots: per 100k residents ("population"), or "max", "median" or "none".
  expert_multipliers = c(0, .01, .015, .02, .03, .05, .075, .1, .15, .2, .3),
  warmup_rounds = 20L,              # First released rows that freeze the error scale S0.
  replay_warmup = TRUE,             # TRUE: learn from the first week on, S0 from those rows in hindsight.
  feedback_delay = 0L,              # Wait this many weeks after the first report, then learn once from the current one.
  feedback_revisions = FALSE,       # TRUE: re-run each week on the reports known then (Delphi).
  score_target = "latest",          # Scores against the "latest" report or the "feedback" learned from.
  score_filled = TRUE,              # Score filled rounds (fill_gaps_up_to in R/01_data_overview.R).
  scale_floor = 1e-6,               # Guard on S0, sigma_t and the level scales, not each rate.
  # The methods' own settings (rates, H, lambda) are in method_parameters() (section 7).
  # Window lengths in weeks. No forecasts were issued 2025-10-25 to 2025-11-15, so
  # COVID after-peak windows of September 2025 peaks miss some weeks.
  weeks_before_peak = 8L,           # before = weeks -7..0: 8 weeks, including the peak week.
  weeks_after_peak = 8L)            # after = weeks 1..8: 8 weeks.

# 2. DELAYED FEEDBACK, CASES, SCALES AND WINDOWS ----------------------------
# The delayed-feedback protocol of notes/notes.tex, Section 2.1.1 (wait for the delay, then
# learn once from the report current at the update). R/01_data_overview.R builds every
# case's feedback with feedback_schedule() and stores this protocol with the prepared data;
# R/05 and R/04 write its id into the tuned settings and the scores. get_case(),
# tuned_parameters() and R/06 stop on any other id, so results of different feedback
# definitions are never combined. Change the id when the protocol changes. The
# initialization is a separate convention (error_scale(), replay_warmup).
feedback_protocol <- function() list(
  id = "latest-at-use-v2",
  cutoff = paste("Issue cutoff of a round: the base's actual generation date where known",
    "(windowed_seasonal; CMU hub rounds that Delphi's production schedule ran late),",
    "otherwise the hub's Wednesday deadline (reference date - 3, also for filled rounds).",
    "Reports have dates only: a report dated on the cutoff day is available (same-day",
    "convention)."),
  vintage = paste("Delay d: each forecast is learned once, at the first later round whose",
    "cutoff is at least its first valid report + 7d days and whose current report (the",
    "latest by that cutoff) is valid, from that report; later revisions are ignored. An",
    "invalid current report waits for the next round; an older report is never used",
    "instead. Re-run: the same rows in the same rounds, each with its report at the re-run",
    "round's cutoff; invalid ones are left out."))
# Stop unless `found` (the protocol id of prepared data, tuned settings or scores, read
# from `what`) is the current one.
check_feedback_protocol <- function(found, what, rerun) {
  found <- unique(found)
  if (!length(found)) found <- "none (made before the protocol was versioned)"
  if (!identical(found, feedback_protocol()$id)) stop(what, " uses feedback protocol ",
    paste(found, collapse = ", "), ", not ", feedback_protocol()$id, ". Rerun ", rerun, ".",
    call. = FALSE)
  invisible(TRUE)
}
# A report is valid if it is a finite, nonnegative count; an NA in a stored version is an
# explicitly missing report.
valid_report <- function(x) is.finite(x) & x >= 0

# The report of each of `rows` that is current at `cutoff` (one date): the latest version
# published by then. versions: every stored report of a case, columns row, version and
# value, sorted by row and version. Value NA: no report yet, or a missing latest report
# (never replaced by an older valid one).
report_snapshot <- function(versions, rows, cutoff) {
  stopifnot(length(cutoff) == 1L, !is.na(cutoff))
  known <- versions[versions$version <= cutoff & versions$row %in% rows, , drop = FALSE]
  known <- known[!duplicated(known$row, fromLast = TRUE), , drop = FALSE]
  i <- match(rows, known$row)
  data.frame(row = rows, version = known$version[i], value = known$value[i])
}

# The delayed feedback of one case with a delay of `delay` weeks (notes/notes.tex, Section
# 2.1.1). versions: every stored report of each row's target week (row, version, value and
# optionally source; sorted by row and version); cutoff: the issue cutoff of each row's
# round, non-decreasing. For each row s:
#   first_report  nu_s, the date of its first valid report;
#   maturity      m_s = nu_s + 7 * delay days: calendar weeks after the first report,
#                 not forecast rounds, so neither a missing round nor a revision restarts
#                 this clock;
#   consumed_row  tau(s): the first later row t > s whose cutoff is at least m_s and whose
#                 current report there (the latest version on or before the cutoff, a
#                 report of the cutoff day included) is valid (NA: none, no feedback). An
#                 invalid current report waits for the next round (invalid_cutoffs counts
#                 these rounds); an older valid report is never used instead. The strict
#                 t > s keeps a forecast's own error out of it, also at horizon -1, whose
#                 first report is often known at its own issue;
#   selected      that current report, learned once at tau(s); later revisions are ignored.
# Returns observed (the selected values of the consumed rows, NA otherwise), date (the
# reference date of each row's consumption round), by_week (the rows consumed at each row,
# before its forecast is issued) and audit.
feedback_schedule <- function(versions, cutoff, delay, reference_date = NULL) {
  n <- length(cutoff)
  stopifnot(length(delay) == 1L, is.finite(delay), delay >= 0, !anyNA(cutoff),
    !is.unsorted(cutoff), all(versions$row %in% seq_len(n)), !anyNA(versions$version),
    identical(order(versions$row, versions$version), seq_len(nrow(versions))),
    !anyDuplicated(versions[c("row", "version")]),
    is.null(reference_date) || length(reference_date) == n)
  source <- if (is.null(versions$source)) rep("archive", nrow(versions)) else versions$source
  first_report <- maturity <- consumed_cutoff <- selected_version <- rep(as.Date(NA), n)
  selected_value <- rep(NA_real_, n)
  first_report_source <- selected_source <- rep(NA_character_, n)
  consumed_row <- invalid_cutoffs <- rep(NA_integer_, n)
  by_row <- split(seq_len(nrow(versions)), factor(versions$row, seq_len(n)))
  for (s in seq_len(n)) {
    k <- by_row[[s]]
    ok <- valid_report(versions$value[k])
    first <- match(TRUE, ok)
    if (is.na(first)) next                        # Never a valid report: no feedback.
    v <- versions$version[k]
    first_report[s] <- v[first]
    first_report_source[s] <- source[k[first]]
    maturity[s] <- v[first] + 7 * delay
    # Eligible rounds: later rows whose cutoff is at least the maturity, each with its
    # current report (the latest version by its cutoff; at least the first valid one).
    eligible <- which(seq_len(n) > s & cutoff >= maturity[s])
    current <- findInterval(as.numeric(cutoff[eligible]), as.numeric(v))
    use <- match(TRUE, ok[current])
    invalid_cutoffs[s] <- if (is.na(use)) length(eligible) else use - 1L
    if (is.na(use)) next                          # No later round with a valid report.
    consumed_row[s] <- eligible[use]
    consumed_cutoff[s] <- cutoff[eligible[use]]
    selected_version[s] <- v[current[use]]
    selected_value[s] <- versions$value[k[current[use]]]
    selected_source[s] <- source[k[current[use]]]
  }
  consumed <- which(!is.na(consumed_row))
  by_week <- lapply(seq_len(n), function(t) which(consumed_row == t))
  stopifnot(all(consumed_row[consumed] > consumed),
    all(maturity[consumed] <= consumed_cutoff[consumed]),
    all(selected_version[consumed] <= consumed_cutoff[consumed]),
    all(valid_report(selected_value[consumed])), !anyDuplicated(unlist(by_week)))
  observed <- rep(NA_real_, n)
  observed[consumed] <- selected_value[consumed]
  date <- if (is.null(reference_date)) rep(as.Date(NA), n) else reference_date[consumed_row]
  list(observed = observed, date = date, by_week = by_week,
    audit = data.frame(row = seq_len(n), issue_cutoff = cutoff, delay = as.integer(delay),
      first_report, first_report_source, maturity, consumed_row, consumed_date = date,
      consumed_cutoff, invalid_cutoffs, selected_version, selected_value, selected_source))
}
# Every case's feedback schedule in one table (one row per case row and delay), with the
# reference and target dates, for audits: feedback_audit(readRDS("output/data/peak_aligned.rds")).
feedback_audit <- function(prepared) do.call(rbind, lapply(prepared$cases, function(case)
  do.call(rbind, lapply(case$feedback_by_delay, function(feedback) data.frame(
    disease = case$disease, location = case$location, horizon = case$horizon,
    reference_date = case$reference_date, target_date = case$target_date,
    filled = case$filled, cutoff_source = case$cutoff_source, feedback$audit)))))

# Prepared cases are named disease_location_horizon, e.g. covid_06_0. feedback_delay
# selects the stateful feedback of feedback_schedule(): each forecast is learned once, at
# the first later round at least feedback_delay weeks after its first valid report whose
# current report is valid, from that report.
# revisions = TRUE: the weekly re-run design (Delphi's "exact" re-run). At every issue week
# each method is re-run from the first week, with the same rows in the same rounds but each
# with its report at that week's cutoff (rerun_methods(), case_as_of()).
get_case <- function(prepared, disease, location, horizon, feedback_delay = 0L,
                     revisions = FALSE) {
  check_feedback_protocol(prepared$settings$feedback_protocol$id,
    "output/data/peak_aligned.rds", "Rscript R/01_data_overview.R")
  case <- prepared$cases[[paste(disease, location, horizon, sep = "_")]]
  if (is.null(case)) stop("No prepared case ", disease, " ", location, " h", horizon,
    ". Rerun Rscript R/01_data_overview.R.")
  feedback <- case$feedback_by_delay[[as.character(feedback_delay)]]
  if (is.null(feedback)) stop("No feedback with a ", feedback_delay, "-week delay. Add it to ",
    "feedback_delays in R/01_data_overview.R and rerun it.")
  case$feedback_protocol <- prepared$settings$feedback_protocol$id
  case$feedback_delay <- feedback_delay
  case$feedback_revised <- revisions
  case$feedback_observed <- feedback$observed
  case$feedback_date <- feedback$date
  case$feedback_by_week <- feedback$by_week
  # Cross-horizon inputs for the soft-switch method. Keep the workflow's selected
  # base and its original availability/filling flags; no outcomes enter this panel.
  case$context_forecasts <- do.call(rbind, lapply(-1:3, function(h) {
    other <- prepared$cases[[paste(disease, location, h, sep = "_")]]
    if (is.null(other) || !identical(other$base_name, case$base_name)) return(NULL)
    mid <- match(.5, other$alpha)
    if (is.na(mid)) return(NULL)
    data.frame(reference_date = other$reference_date, target_date = other$target_date,
      as_of = other$as_of, horizon = h, median = other$base[, mid], filled = other$filled)
  }))
  case
}
# Label of a feedback design.
feedback_label <- function(delay, revisions = FALSE)
  paste0("delay ", delay, if (revisions) " + revisions (re-run)")

# The case as re-run at issue week `week` (re-run design): its first `week` rows. The rows
# consumed by then keep their rounds of the stateful schedule (every delay; later
# revisions never move them), but each is learned from its report current at this week's
# cutoff (case$as_of[week]); a row whose current report is invalid is left out of this
# replay (rerun_omitted), not replaced by an older or later report. The diverse experts'
# delay-0 feedback is replayed alike. S0 and the weighting loss scale stay frozen
# (frozen_s0, frozen_loss_scale: those of the stateful design with the same delay, see
# rerun_methods()).
case_as_of <- function(case, week, frozen_s0 = NULL, frozen_loss_scale = NULL) {
  rows <- seq_len(week)
  current <- report_snapshot(case$report_versions, rows, case$as_of[week])$value
  replay <- function(feedback) {
    by_week <- feedback$by_week[rows]
    released <- unlist(by_week)
    omitted <- released[!valid_report(current[released])]
    observed <- rep(NA_real_, week)
    kept <- setdiff(released, omitted)
    observed[kept] <- current[kept]
    list(observed = observed, date = feedback$date[rows],
      by_week = lapply(by_week, setdiff, omitted), omitted = omitted)
  }
  case$feedback_by_delay <- lapply(case$feedback_by_delay, replay)
  own <- case$feedback_by_delay[[as.character(case$feedback_delay)]]
  for (name in c("reference_date", "target_date", "truth_finalized", "filled", "filled_from",
    "truth_latest", "as_of", "cutoff_source", "season", "period", "peak_index", "peak_date",
    "relative_week")) case[[name]] <- case[[name]][rows]
  case$base <- case$base[rows, , drop = FALSE]
  case$feedback_observed <- own$observed
  case$feedback_date <- own$date
  case$feedback_by_week <- own$by_week
  case$rerun_omitted <- own$omitted
  case$feedback_revised <- FALSE
  case$frozen_s0 <- frozen_s0
  case$frozen_loss_scale <- frozen_loss_scale
  case
}

# Three prospective trend signals, in log1p(counts) for either tracker. Each
# lookup uses an exact reference/target date, excludes filled forecasts, and must
# have been available by that issue's as-of cutoff. Missing signals are neutral.
# The gate describes current issue context, not a peak probability or the target
# horizon's future state. In particular, h=-1 is available at reference t, not t-1.
# weights: of the path, nowcast and revision scores in the gate (equal: their mean;
# all 0: a neutral gate).
transition_gate <- function(case, strength = 1, width = .05, threshold = .05,
                            weights = c(1, 1, 1)) {
  stopifnot(length(strength) == 1L, is.finite(strength), strength >= 0, strength <= 1,
    length(width) == 1L, is.finite(width), width > 0,
    length(threshold) == 1L, is.finite(threshold), threshold >= 0,
    length(weights) == 3L, all(is.finite(weights)), all(weights >= 0))
  dates <- case$reference_date
  signals <- matrix(NA_real_, length(dates), 3L,
    dimnames = list(NULL, c("path", "nowcast", "revision")))
  history <- case$context_forecasts
  if (!is.null(history) && nrow(history)) {
    stopifnot(!anyDuplicated(history[c("reference_date", "horizon")]),
      all(history$target_date == history$reference_date + 7 * history$horizon))
    keys <- paste(history$reference_date, history$horizon)
    lookup <- function(reference, horizon, cutoff) {
      j <- match(paste(reference, horizon), keys)
      ok <- !is.na(j) & !is.na(history$as_of[j]) & !is.na(cutoff) &
        history$as_of[j] <= cutoff & !is.na(history$filled[j]) & !history$filled[j] &
        is.finite(history$median[j]) & history$median[j] >= 0
      value <- rep(NA_real_, length(j))
      value[ok] <- log1p(history$median[j[ok]])
      value
    }
    for (i in seq_along(dates)) {
      t <- dates[i]; cutoff <- case$as_of[i]
      current <- lookup(t, -1:3, cutoff)
      # Four current/future forecasts suffice if h=-1 is missing. Otherwise include
      # its edge to h=0 too, as in the descriptive five-point curve.
      if (all(is.finite(current[2:5]))) signals[i, "path"] <-
        max(-diff(current), na.rm = TRUE)
      previous_row <- match(paste(t - 7, -1), keys)
      previous_cutoff <- if (is.na(previous_row)) as.Date(NA) else history$as_of[previous_row]
      if (is.na(cutoff) || is.na(previous_cutoff) || previous_cutoff > cutoff) next
      previous <- lookup(t - 7, -1, previous_cutoff)
      before_previous <- lookup(t - 14, -1, previous_cutoff)
      growth <- current[1] - previous
      previous_growth <- previous - before_previous
      if (is.finite(growth) && is.finite(previous_growth))
        signals[i, "nowcast"] <- max(-growth, min(previous_growth, previous_growth - growth))
      # Compare each nowcast with older horizons j for its SAME target. The prior
      # revision is reconstructed at its own cutoff, not with information from t.
      revision <- current[1] - lookup(t - 7 * (1:4), 0:3, cutoff)
      previous_revision <- previous - lookup(t - 7 * (2:5), 0:3, previous_cutoff)
      common <- is.finite(revision) & is.finite(previous_revision)
      if (sum(common) >= 2L) signals[i, "revision"] <- min(-median(revision[common]),
        -median(revision[common] - previous_revision[common]))
    }
  }
  scores <- plogis((signals - log1p(threshold)) / log1p(width))
  scores[!is.finite(scores)] <- .5         # Unknown evidence pulls towards neutral.
  mixed <- if (sum(weights) == 0) rep(.5, nrow(scores)) else
    if (length(unique(weights)) == 1L) rowMeans(scores) else
      drop(scores %*% weights) / sum(weights)
  gate <- pmin(.95, pmax(.05, .5 + strength * (mixed - .5)))
  data.frame(reference_date = dates, signals, path_score = scores[, "path"],
    nowcast_score = scores[, "nowcast"], revision_score = scores[, "revision"],
    gate, n_available = rowSums(is.finite(signals)))
}

# Tracker units: "log1p", "identity" (counts) or "sqrt" (Delphi's REF-op; negative
# tracker values map back to 0, as in Delphi's calibration code).
to_tracker <- function(x, transform) {
  stopifnot(transform %in% c("log1p", "identity", "sqrt"))
  switch(transform, log1p = log1p(x), identity = x, sqrt = sqrt(x))
}
from_tracker <- function(x, transform) {
  stopifnot(transform %in% c("log1p", "identity", "sqrt"))
  switch(transform, log1p = expm1(x), identity = x, sqrt = pmax(x, 0)^2)
}

# Count scale for scores and plots; MultiQT never tracks it, and the expert weights
# divide it out (weighting_loss_scale()). "population": counts per 100k residents
# (2019 census, stored by R/01_data_overview.R). "max" and "median": of the reports
# published before the location's first forecast. "none": counts.
choose_count_scale <- function(case, method = "population") {
  if (method == "population") {
    if (is.null(case$population)) stop("No population in the prepared data. ",
      "Rerun Rscript R/01_data_overview.R.")
    return(list(scale = case$population / 1e5, used = method, method = method,
      reference_weeks = NA_integer_))
  }
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
count_scale_label <- function(method) switch(method,
  population = "counts per 100k residents",
  none = "counts",
  paste0("counts / historical ", method))
# The pooled all-location summaries of 02, 04 and 06 leave out settings$summary_exclude
# (the US: a sum over the states, which would count them twice). pooled_rows() keeps the
# other rows of a score table; pooled_label() names the pooled set in titles (the left-out
# US is not named there).
pooled_rows <- function(rows, settings) rows[!rows$location %in% settings$summary_exclude, ,
  drop = FALSE]
pooled_label <- function(n, settings) paste(n, "locations")

# Frozen error scale S0: Q90 of the absolute base errors in tracker units, pooled
# over all quantile levels of the first `rounds` released rows. eta = m * S0.
# https://delphi-forecasting-reports.netlify.app/calibration_experiments/e01_eta_settings
# The whole release batch of the freeze week is consumed as warmup (no replay), but
# only its first `rounds` rows enter S0. With fewer releases every week is warmup.
# replay = TRUE keeps S0 but has no warmup weeks: MultiQT learns from the first
# week on, so S0 uses reports released later (a retrospective initialization, declared
# separately from the feedback protocol: forecasts before the 20th release depend on it).
# The released rows and values are those of the case's stateful schedule (the reports
# learned from); the re-run design keeps this S0 of its stateful design (rerun_methods()).
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
    raw_scale <- abs_error_q90(case, rows, transform)
    scale <- max(floor, raw_scale)
  }
  warmup_weeks <- if (replay && initialized) 0L else freeze_week
  list(initialized = initialized, warmup_weeks = warmup_weeks,
    date = if (initialized) case$reference_date[freeze_week] else as.Date(NA),
    rows = rows, consumed_rows = if (warmup_weeks > 0L) released else integer(),
    raw_scale = raw_scale, scale = scale, rounds = rounds, floor = floor)
}
# Q90 of the absolute base errors in tracker units, pooled over all quantile levels of
# the given rows, against the reports learned from: S0 above, Delphi's adaptive rate below.
abs_error_q90 <- function(case, rows, transform) {
  residuals <- sweep(to_tracker(case$base[rows, , drop = FALSE], transform), 1,
    to_tracker(case$feedback_observed[rows], transform), "-")
  as.numeric(quantile(abs(residuals), .9, names = FALSE))
}

# Delphi's adaptive rate ("adaptive+" in their calibration code, used by REF-op): at each
# issue week, mult * abs_error_q90() of the last `rounds` rows learned by then (Inf: all
# of them), at least floor. The statistic of S0, but rolling with the feedback instead of
# frozen, so it needs no warmup and no hindsight. 0 before the first learned row (no
# update then). In the re-run design it uses the reports current at the issue week.
adaptive_rate <- function(case, mult, rounds = 20L, floor = 0, transform = "sqrt") {
  released <- integer()
  eta <- numeric(nrow(case$base))
  for (week in seq_along(eta)) {
    released <- c(released, case$feedback_by_week[[week]])
    if (length(released)) eta[week] <- max(floor, mult * abs_error_q90(case,
      if (is.finite(rounds)) tail(released, rounds) else released, transform))
  }
  eta
}

# Reference loss that makes the expert weights unit-free: the base's mean pinball loss
# over the S0 rows (error_scale()), against the reports learned from, in counts /
# count_scale like the weighting losses. exp_weights() divides its losses by it, so
# lambda_start and fixed lambdas act alike under every count scale and location.
# Like S0, it uses the S0 rows in hindsight when the warmup is replayed. 1 without S0
# (then every week is warmup and the weights never update).
weighting_loss_scale <- function(case, scale, count_scale = 1, floor = 1e-6) {
  if (!scale$initialized) return(1)
  rows <- scale$rows
  error <- (case$feedback_observed[rows] - case$base[rows, , drop = FALSE]) / count_scale
  loss <- pmax(error, 0) * rep(case$alpha, each = length(rows)) +
    pmax(-error, 0) * rep(1 - case$alpha, each = length(rows))
  max(floor, mean(loss))
}

# Scale of recent residuals r: their root mean square, or the median absolute residual
# divided by 0.6745 (a normal sd that one outlier cannot inflate).
residual_scale <- function(r, estimator = "rms") switch(estimator,
  rms = sqrt(mean(r^2)),
  median = median(abs(r)) / qnorm(.75),
  stop("Unknown residual scale estimator ", estimator))

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
# estimator (residuals): "rms", or "median": median |residual| / 0.6745 (the normal sd), which
# one outlying residual (e.g. from a broken base forecast) cannot inflate.
median_sigma <- function(case, method = "residuals", window = 5L, min_rows = 5L,
                         floor = 1e-6, transform = "log1p", levels = c(.1, .9),
                         estimator = "rms") {
  stopifnot(window >= 1L, min_rows >= 1L)   # Otherwise a later report enters sigma_t.
  base_level <- function(alpha) to_tracker(case$base[, match(alpha, case$alpha)], transform)
  scale <- switch(method,
    residuals = {
      released <- unlist(case$feedback_by_week)              # Rows in release order.
      residual <- to_tracker(case$feedback_observed[released], transform) -
        base_level(.5)[released]
      n_released <- cumsum(lengths(case$feedback_by_week))   # By each issue week.
      vapply(n_released, function(k) if (k < min_rows) NA_real_ else
        max(floor, residual_scale(residual[max(1L, k - window + 1L):k], estimator)), numeric(1))
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
# level, the RMS (or robust scale, estimator as in median_sigma()) of y_s - b_s^alpha over
# the last `window` rows s released by issue week t (no hindsight). Rows are NA until
# min_rows are released, like median_sigma("residuals").
# Returns a matrix [week x level]. Like level_sigma(), no calibration guarantee.
level_residual_sigma <- function(case, window = 5L, min_rows = 5L, floor = 1e-6,
                                 transform = "log1p", estimator = "rms") {
  stopifnot(window >= 1L, min_rows >= 1L)   # Otherwise a later report enters sigma_t.
  released <- unlist(case$feedback_by_week)                # Rows in release order.
  residual <- to_tracker(case$feedback_observed[released], transform) -
    to_tracker(case$base[released, , drop = FALSE], transform)
  n_released <- cumsum(lengths(case$feedback_by_week))     # By each issue week.
  scale <- t(vapply(n_released, function(k) if (k < min_rows) rep(NA_real_, ncol(residual)) else
    pmax(apply(residual[max(1L, k - window + 1L):k, , drop = FALSE], 2, residual_scale,
      estimator), floor),
    numeric(ncol(residual))))
  dimnames(scale) <- dimnames(case$base)
  scale
}

# Peak-relative windows, as a logical vector over the case's weeks, within the given
# period(s) of the target week: "evaluation", "tuning" or "other" (the design in
# R/01_data_overview.R). A peak is a local maximum of the latest-reported counts
# (find_peaks() in R/01_data_overview.R; hindsight label only); each week counts from
# its nearest peak. before and after are window lengths:
# before = weeks -(before - 1)..0, including the peak week, and after = weeks
# 1..after; overall = all weeks of the period(s).
phase_weeks <- function(case, phase, before = 8L, after = 8L, period = "evaluation") {
  week <- case$relative_week
  in_period <- case$period %in% period
  switch(phase,
    before = in_period & week > -before & week <= 0,
    after = in_period & week >= 1 & week <= after,
    overall = in_period,
    stop("Unknown phase ", phase))
}

# Missing weeks shorten the before/after windows: report them, but keep the peak.
check_peak_window <- function(case, before = 8L, after = 8L, truth = case$truth_latest,
                              period = "evaluation") {
  label <- paste0(case$disease, " | ", case$location_name, " (", case$location,
    ") | horizon ", case$horizon)
  weeks <- seq.int(1L - before, after)
  ok_all <- TRUE
  for (k in which(case$peaks$period %in% period)) {
    peak <- case$peaks[k, ]
    peak_label <- paste0(label, " | peak ", peak$peak_date)
    rows <- which(case$peak_index == k)
    i <- rows[match(weeks, case$relative_week[rows])]
    ok <- !is.na(i)
    ok[ok] <- case$period[i[ok]] %in% period & is.finite(truth[i[ok]])
    if (!all(ok))
      warning(peak_label, ": no scored target in peak-relative week(s) ",
        paste(weeks[!ok], collapse = ","), "; peak kept.", call. = FALSE, immediate. = TRUE)
    ok_all <- ok_all && all(ok)
  }
  invisible(ok_all)
}

# 3. SCORING --------------------------------------------------------------
# Truth for coverage, coverage residuals and quantile losses. "latest": the latest
# stored report of each target week (the archive's last version, prepared$settings$
# truth_version), the same truth for every feedback design. "feedback": the report the
# stateful updates learned from at the case's feedback delay, current at the update (NA
# for rows never consumed; in the re-run design, too, the stateful design's report).
# Calibration against these learning vintages is not calibration against the latest counts.
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
  feedback = paste0("the reports learned from (delay ", delay, ": current at the update)"),
  stop("Unknown score target ", target))

# Coverage and pinball loss of quantile forecasts q [week x level] against one
# truth per week. Pass counts with their count_scale to compare coverage on the
# count scale (exact integer ties stay ties); the loss is in counts / count_scale.
# The pinball loss penalizes underprediction by alpha, overprediction by 1 - alpha.
# score_scale "log": the loss of log(counts + 1) instead (q and truth in counts;
# count_scale unused), which scores relative errors and growth rates (Bosse et al.
# 2023, PLOS Comput Biol 19(8): e1011393). Negative quantiles count as 0 there.
# Quantiles commute with the monotone log, so the score stays proper; coverage is
# the same on both scales.
score_quantiles <- function(q, truth, alpha, count_scale = 1, score_scale = "natural") {
  stopifnot(ncol(q) == length(alpha), length(truth) == nrow(q))
  error <- switch(score_scale,
    natural = truth / count_scale - q / count_scale,
    log = log1p(truth) - log1p(pmax(q, 0)),
    stop("Unknown score scale ", score_scale))
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

# Plots: an NA row after every gap of more than max_days between consecutive dates, so
# lines and polygons break there (e.g. May to November 2024, between the seasons).
# Returns the dates and each further vector or matrix with the same rows.
break_gaps <- function(dates, ..., max_days = 28) {
  gap <- which(diff(as.numeric(dates)) > max_days)
  index <- sort(c(seq_along(dates), gap + .5))
  rows <- ifelse(index %% 1 == 0, index, NA)
  c(list(dates = dates[rows]), lapply(list(...), function(x)
    if (is.matrix(x)) x[rows, , drop = FALSE] else x[rows]))
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

# base: issued quantiles [week x level] in counts; feedback: the reported counts learned
# from (the case's feedback schedule); feedback_by_week[[w]]: the earlier rows consumed at
# issue week w, before its forecast is issued. Each consumed row is compared with ITS
# historical issued forecast, and the steps of a batch are summed (not averaged), all with
# eta[w], the rate of the update week. eta: one rate, or one per week, in tracker units.
# Only the issued forecasts are projected; the offsets stay unprojected. The first
# warmup_weeks forecasts keep the base and their reports are never replayed.
# Optional leaky shift of all levels on the log1p(counts) scale, for either tracker
# (leak = 0: standard MultiQT): each report adds leak * mean(alpha - covered); the shift
# is multiplied by decay (one value, or one per week) at every issue week and capped at
# +-leak_bound. The identity tracker shifts negative quantiles like 0.
# Optional studentized offsets (sigma = 1, start = 0: standard MultiQT): week w issues
# base + sigma[w] * offsets. sigma: one value, or one per week, in tracker units; then
# offsets and eta are in units of sigma, and a change of sigma rescales the whole
# accumulated offset (not the leaky shift). A matrix sigma
# [week x level] gives each level its own scale; the projection is then weighted
# isotonic regression with weights 1 / sigma^2 (Euclidean in units of sigma).
# start: the offsets before the first update, one value or one per level.
# Optional slow term and leaky offsets (Delphi's REF-op; eta_slow = 0, fast_decay = 0:
# standard MultiQT): week w issues base + sigma[w] * (offsets + slow), both terms take the
# same steps, the slow one with eta_slow (one rate, or one per week), and the offsets are
# multiplied by 1 - fast_decay at every issue week after the first, before its update;
# the slow term never decays.
# Returns the adjusted forecasts in counts.
multiqt <- function(base, feedback, feedback_by_week, alpha, eta,
                    warmup_weeks = 0L, transform = "log1p",
                    leak = 0, decay = .5, leak_bound = .2, sigma = 1, start = 0,
                    eta_slow = 0, fast_decay = 0) {
  n <- nrow(base)
  # Only warmup weeks (e.g. S0 never initialized), or nothing to learn and start = 0:
  # keep the exact base.
  if (warmup_weeks >= n || (isTRUE(all(eta == 0)) && isTRUE(all(eta_slow == 0)) &&
    leak == 0 && isTRUE(all(start == 0))))
    return(base)
  if (length(eta) == 1L) eta <- rep(eta, n)
  if (length(eta_slow) == 1L) eta_slow <- rep(eta_slow, n)
  if (length(decay) == 1L) decay <- rep(decay, n)
  if (length(sigma) == 1L) sigma <- rep(sigma, n)
  by_level <- is.matrix(sigma)              # One scale per level: weighted projection.
  base_tracker <- to_tracker(base, transform)
  feedback_tracker <- to_tracker(feedback, transform)
  warmup <- seq_len(n) <= warmup_weeks      # sigma may be missing in warmup weeks.
  issued_sigma <- if (by_level) sigma[!warmup, , drop = FALSE] else sigma[!warmup]
  stopifnot(ncol(base) == length(alpha), length(feedback) == n,
    length(feedback_by_week) == n, length(eta) == n, all(is.finite(eta)), all(eta >= 0),
    length(eta_slow) == n, all(is.finite(eta_slow)), all(eta_slow >= 0),
    length(fast_decay) == 1L, fast_decay >= 0, fast_decay < 1,
    NROW(sigma) == n, !by_level || ncol(sigma) == length(alpha),
    all(is.finite(issued_sigma) & issued_sigma > 0),
    length(start) %in% c(1L, length(alpha)), all(is.finite(start)),
    all(is.finite(base_tracker)), all(is.finite(feedback_tracker[unlist(feedback_by_week)])))
  offsets <- rep_len(start, length(alpha))
  slow <- rep(0, length(alpha))
  shift <- 0
  # The leaky shift on the log1p scale: inside the projection for log1p, after it for
  # identity (isoreg commutes with a constant shift, so log1p stays exactly as before).
  inner_shift <- function() if (transform == "log1p") shift else 0
  log_shift <- function(q) if (shift == 0) q else expm1(log1p(pmax(q, 0)) + shift) + pmin(q, 0)
  q_tracker <- base_tracker
  for (week in seq_len(n)) {
    if (week <= warmup_weeks) next
    shift <- shift * decay[week]
    if (fast_decay > 0 && week > warmup_weeks + 1L) offsets <- offsets * (1 - fast_decay)
    for (past_week in feedback_by_week[[week]]) {
      # Compare with the HISTORICAL ISSUED forecast in tracker units.
      # At alpha=.9, eta=.01, a miss adds .009; coverage subtracts .001.
      covered <- feedback_tracker[past_week] <= q_tracker[past_week, ]
      offsets <- offsets + eta[week] * (alpha - covered)
      slow <- slow + eta_slow[week] * (alpha - covered)
      shift <- shift + leak * mean(alpha - covered)
    }
    shift <- min(leak_bound, max(-leak_bound, shift))
    # Project only the issued forecasts, NOT the offsets ("lazy" MultiQT).
    # isoreg restores quantile order: e.g., (12,10) becomes (11,11).
    # Weights (min sigma / sigma)^2: the same projection as 1 / sigma^2, without overflow.
    q_tracker[week, ] <- if (by_level) weighted_isoreg(base_tracker[week, ] +
      sigma[week, ] * (offsets + slow) + inner_shift(),
      (min(sigma[week, ]) / sigma[week, ])^2) else
      isoreg(base_tracker[week, ] + sigma[week] * (offsets + slow) + inner_shift())$yf
    if (transform != "log1p" && shift != 0) q_tracker[week, ] <- to_tracker(log_shift(
      from_tracker(q_tracker[week, ], transform)), transform)
  }
  # Do not clip negative predictions: that would change the issued coverage indicators.
  # (The sqrt tracker maps negative values back to 0 counts, as Delphi's code does; its
  # coverage indicators above stay those of the unclipped tracker values.)
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
#   leaky    slow MultiQT plus a leaky shift of all levels on the log1p scale (see
#            multiqt()), halved every week since the previous issue;
#   damp=k   start from the base minus k times its implied growth on the log1p scale
#            (either tracker): log1p(base median) minus log1p(the latest released
#            report), clipped to [-1, 1] (k = 1 centres on the latest report, k < 0
#            extrapolates further);
#   width=w  start from the base quantiles moved w times as far from the base median
#            (tracker units).
make_diverse_experts <- function(case, scale, count_scale = 1, transform = "log1p",
                                 slow = .01, leak = .2, damping = c(-.5, .5, 1),
                                 widths = c(.75, 1.5)) {
  n <- nrow(case$base)
  base_tracker <- to_tracker(case$base, transform)
  mid <- base_tracker[, match(.5, case$alpha)]
  # Latest row released by each issue week under delay 0 (none: no growth), with its report
  # current at that release, at any feedback delay: the starting point uses what is known
  # then; in a re-run, case_as_of() replaces it by the current report.
  first <- case$feedback_by_delay[["0"]]
  latest <- cummax(vapply(first$by_week, function(rows) max(0L, rows), integer(1)))
  implied <- log1p(case$base[, match(.5, case$alpha)]) -
    log1p(first$observed[ifelse(latest > 0, latest, NA)])
  implied <- pmin(1, pmax(-1, ifelse(is.finite(implied), implied, 0)))
  # Starting forecasts in counts, floored at zero counts.
  starts <- list(slow = case$base, leaky = case$base)
  for (k in damping) starts[[paste0("damp=", k)]] <-
    expm1(pmax(log1p(case$base) - k * implied, 0))
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
#   these units divided by loss_scale (weighting_loss_scale(): the base's mean loss
#   over the S0 rows), released by feedback_by_week like the MultiQT updates. The
#   division makes lambda unit-free: the count scale does not change the weights.
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
# [week x expert x level, or x 1 if common], the lambda used at each issue week and the
# experts' forecasts it weighs (experts$forecasts).
# update_weight: nonnegative multiplier saved for each forecast at issuance; a
# released row uses update_weight[past], never the current context. It scales the
# loss (including its second moment) and the fixed-share exposure. Ones preserve
# ordinary weighting; zeros leave the bank untouched for that released row.
exp_weights <- function(experts, base, feedback, feedback_by_week, alpha,
                        lambda = "learned", sharing = TRUE, H = 8L, share = 1 / (2 * H),
                        per_quantile = TRUE, lambda_start = 2.72, warmup_weeks = 0L,
                        loss_scale = 1, update_weight = rep(1, nrow(base))) {
  forecasts <- experts$forecasts
  birth <- experts$birth_week
  n <- nrow(base); k <- length(birth); a <- length(alpha)
  learned <- identical(lambda, "learned")
  stopifnot(identical(dim(forecasts), c(n, k, a)), length(feedback) == n,
    length(feedback_by_week) == n, learned || (is.numeric(lambda) &&
      length(lambda) == 1L && is.finite(lambda) && lambda >= 0),
    !sharing || (is.finite(share) && share > 0 && share <= 1),
    length(loss_scale) == 1L, is.finite(loss_scale), loss_scale > 0,
    length(update_weight) == n, all(is.finite(update_weight)), all(update_weight >= 0))
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
      if (update_weight[past] == 0) next
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
      loss <- (pmax(error, 0) * rep(alpha, each = length(eligible)) +
        pmax(-error, 0) * rep(1 - alpha, each = length(eligible))) / loss_scale
      if (!per_quantile) loss <- matrix(rowMeans(loss), ncol = 1L)
      loss <- loss * update_weight[past]
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
        share_now <- if (update_weight[past] == 1) share else
          -expm1(update_weight[past] * log1p(-share))
        p[eligible, ] <- (1 - share_now) * p[eligible, , drop = FALSE] +
          rep(share_now * mass / length(eligible), each = length(eligible))
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
  list(q = q, weights = weights, lambda = lambda_used, experts = forecasts)
}

# Two context-conditioned weight banks for the SAME experts, learned separately
# for this location and forecast horizon. Twice the saved membership normalizes
# neutral context (.5) to the ordinary update, so switch_strength=0 is exactly the
# common-weight comparator. Both banks use the ordinary own-horizon delayed loss.
# Their learned lambda uses context-scaled moments: this is a heuristic extension
# of exponential weighting, not a new calibration or peak-detection guarantee.
soft_switch_weights <- function(experts, base, feedback, feedback_by_week, alpha,
                                gate, ...) {
  stopifnot(length(gate) == nrow(base), all(is.finite(gate)), all(gate >= 0 & gate <= 1))
  baseline <- exp_weights(experts, base, feedback, feedback_by_week, alpha,
    per_quantile = FALSE, update_weight = 2 * (1 - gate), ...)
  declining <- exp_weights(experts, base, feedback, feedback_by_week, alpha,
    per_quantile = FALSE, update_weight = 2 * gate, ...)
  weights <- sweep(baseline$weights, 1, 1 - gate, `*`) +
    sweep(declining$weights, 1, gate, `*`)
  q <- (1 - gate) * baseline$q + gate * declining$q
  # Keep identical forecasts, in particular the base during warmup, exact.
  same <- rowSums(baseline$q != declining$q) == 0
  q[same, ] <- baseline$q[same, , drop = FALSE]
  regimes <- array(c(baseline$weights, declining$weights),
    dim = c(dim(weights)[1:2], 2L),
    dimnames = c(dimnames(weights)[1:2], list(c("baseline", "declining"))))
  # lambda is an average diagnostic for the existing reporting interface. The
  # actual updates use the two separate learning rates returned in regime_lambda.
  list(q = q, weights = weights,
    lambda = (1 - gate) * baseline$lambda + gate * declining$lambda,
    regime_weights = regimes,
    regime_lambda = cbind(baseline = baseline$lambda[, 1], declining = declining$lambda[, 1]),
    experts = experts$forecasts)
}

# 7. THE RECALIBRATION METHODS OF SCRIPTS 04 AND 05 ----------------------------
# Settings of the methods. Rates are multiples m of S0 (eta = m * S0, so they act alike
# for both trackers and every location), except median_rate (in units of sigma_t), the
# diverse leak and damping (log1p units for either tracker) and REF-op's rates (Delphi's
# rolling adaptive rate, adaptive_rate()). by_method overrides a
# setting for one method, e.g. list(shared = list(H = 4L)); R/05_tuning.R tunes one
# setting per method (tuned_setting) on the tuning season.
method_parameters <- function() list(
  fixed_multiplier = .05,                # Fixed MultiQT: eta = m * S0.
  # Delphi's REF-op as in covid prod, cold (no burn-in season), with Delphi's constants,
  # not tuned (pipelines/covid_hosp_prod.R, cal_ref_args("op") in R/calibration/views.R of
  # cmu-delphi/delphi-hub-forecasting): MultiQT on sqrt(counts) for either tracker, with
  # the fast rate ref_op_multiplier * Q90 of the last ref_op_window learned rows (at least
  # ref_op_floor) and a slow term with ref_op_slow_multiplier * Q90 of all learned rows;
  # the fast offsets shrink by ref_op_fast_decay every week, the slow ones never.
  ref_op_transform = "sqrt", ref_op_multiplier = .03, ref_op_window = 20L,
  ref_op_floor = 1e-3, ref_op_slow_multiplier = .003, ref_op_fast_decay = .1,
  # Delphi's flu prod since 2026-10-07 (pipelines/flu_hosp_prod.R, "sqrt constant 0.018 per
  # 100k", cold), not tuned: Fixed MultiQT on sqrt(counts per 100k residents) for either
  # tracker, with this constant rate in those units (Delphi chose it on flu 2024/25-2025/26).
  sqrt_constant_rate = .018,
  # Slope rate: eta_t = min(cap * S0, min * S0 + |weekly slope of the base median|).
  slope_min_multiplier = .004, slope_cap_multiplier = .08,
  H = 8L, lambda_start = 2.72,           # Learned lambda (exp_weights()); sharing uses 1/(2H).
  outer_multiplier = .004,               # Outer MultiQT on a mixture's own forecasts.
  # Weekly restarts: a fresh copy per multiplier at every issue week, fixed meta-rate.
  restart_multiplier = .05, restart_multipliers = c(.01, .05, .2), restart_lambda = 2.72,
  # Diverse experts: slow MultiQT (eta = m * S0) on differently corrected starting forecasts.
  diverse_slow_multiplier = .008, diverse_leak = .2, diverse_damping = c(-.5, .5, 1),
  diverse_widths = c(.75, 1.5),
  # Soft switch: sigmoid scale/centre on log1p(counts), independent of tracker units.
  # strength=0 is the ordinary diverse common-weight mixture; 1 uses all the signal.
  switch_strength = 1, switch_width = .05, switch_threshold = .05,
  # Weights of the gate's path, nowcast and revision signals (1, 1, 1: their mean).
  switch_weight_path = 1, switch_weight_nowcast = 1, switch_weight_revision = 1,
  # Studentized MultiQT: sigma_t from the RMS of the last median_window released
  # residuals (the base until median_min_rows exist) or from the base width between
  # median_base_levels; median_start: "normal" or "t" start of the residual median methods.
  median_window = 5L, median_min_rows = 5L, median_base_levels = c(.1, .9),
  median_rate = .1, median_start = "normal",
  # Residual scales (studentized median, constant sigma, level RMS): "rms" or "median"
  # (robust to one outlying residual; see residual_scale()).
  residual_estimator = "rms",
  # Every studentized scale is clipped to these multiples of S0 (tracker units). Near 0
  # (e.g. a base of zero width) the offsets grow without effect for weeks; a jump (e.g.
  # a broken base forecast or report) then multiplies all of them. c(0, Inf): no limits.
  sigma_bounds = c(.05, 1),
  by_method = list())
method_parameter <- function(params, method, name) {
  value <- params$by_method[[method]][[name]]
  # Reuse the comparator's tuned expert bank unless explicitly overridden.
  if (is.null(value) && method == "diverse_soft_switch")
    value <- params$by_method$diverse_shared_common[[name]]
  if (is.null(value)) params[[name]] else value
}
# The values R/05_tuning.R chose for one disease, tracker transform and feedback design
# (delay, revisions) under a design of R/01_data_overview.R, as by_method overrides of
# params; the defaults where none were tuned (with a message). Settings tuned under another
# feedback protocol (feedback_protocol()) stop with an error instead of being used.
tuned_parameters <- function(disease, transform, params = method_parameters(),
                             design = "seasons", dir = "output/results/tuning",
                             feedback_delay = 0L, revisions = FALSE) {
  path <- file.path(dir, paste0("tuned_settings_", design, ".csv"))
  tuned <- if (file.exists(path)) read.csv(path) else data.frame(disease = character(),
    transform = character(), feedback_delay = integer(), feedback_revisions = logical(),
    feedback_protocol = character())
  if (file.exists(path)) check_feedback_protocol(tuned$feedback_protocol, path,
    "Rscript R/05_tuning.R")
  tuned <- tuned[tuned$disease == disease & tuned$transform == transform &
    tuned$feedback_delay == feedback_delay & tuned$feedback_revisions == revisions, ,
    drop = FALSE]
  if (!nrow(tuned)) message("No tuned settings for ", disease, " (", transform, ", feedback ",
    feedback_label(feedback_delay, revisions), ") in ", path,
    ": the defaults of method_parameters(). Run Rscript R/05_tuning.R.")
  for (i in seq_len(nrow(tuned)))
    params$by_method[[tuned$method[i]]][[tuned$setting[i]]] <- tuned$value[i]
  params
}
# The setting R/05_tuning.R tunes for each method; the other methods have none.
tuned_setting <- c(fixed = "fixed_multiplier", absolute_slope = "slope_min_multiplier",
  weighted = "H", weighted_common = "H", shared = "H", shared_common = "H",
  restart_weighted = "restart_lambda", restart_grid = "restart_lambda",
  diverse_uniform = "diverse_slow_multiplier", diverse_shared = "diverse_slow_multiplier",
  diverse_shared_common = "diverse_slow_multiplier",
  diverse_soft_switch = "switch_strength",
  studentized = "median_rate", studentized_base = "median_rate", median_constant = "median_rate",
  studentized_full = "median_rate", studentized_levels = "median_rate",
  studentized_level_rms = "median_rate",
  shared_outer = "outer_multiplier", shared_common_outer = "outer_multiplier",
  restart_uniform_outer = "outer_multiplier", restart_weighted_outer = "outer_multiplier",
  restart_grid_outer = "outer_multiplier")
method_names <- c("original", "fixed", "ref_op", "sqrt_constant", "absolute_slope", "studentized",
  "studentized_base",
  "median_constant", "studentized_full", "studentized_levels", "studentized_level_rms",
  "oracle_target", "uniform", "weighted", "weighted_common", "shared", "shared_common",
  "restart_uniform", "restart_weighted", "restart_grid", "diverse_uniform", "diverse_shared",
  "diverse_shared_common", "diverse_soft_switch", "restart_persistent", "shared_outer", "shared_common_outer",
  "restart_uniform_outer", "restart_weighted_outer", "restart_grid_outer")
# The mixture each outer MultiQT corrects.
outer_inner <- c(shared_outer = "shared", shared_common_outer = "shared_common",
  restart_uniform_outer = "restart_uniform", restart_weighted_outer = "restart_weighted",
  restart_grid_outer = "restart_grid")
# The method families of 04's plots (forecast pages, trade-off colours); their methods are
# also the candidates of 06's tree.
method_groups <- list(
  "Single trackers" = c("original", "fixed", "ref_op", "sqrt_constant", "absolute_slope",
    "oracle_target"),
  "Rate-grid experts" = c("original", "fixed", "uniform", "weighted", "weighted_common",
    "shared", "shared_common", "shared_outer", "shared_common_outer"),
  "Weekly starts: immediate admission" = c("original", "restart_persistent", "restart_uniform",
    "restart_weighted", "restart_grid", "restart_uniform_outer", "restart_weighted_outer",
    "restart_grid_outer"),
  "Diverse experts" = c("original", "shared", "shared_common", "diverse_uniform",
    "diverse_shared", "diverse_shared_common", "diverse_soft_switch"),
  "Studentized median" = c("original", "fixed", "studentized", "studentized_base",
    "median_constant"),
  "Studentized full base" = c("original", "fixed", "studentized_full", "studentized_levels",
    "studentized_level_rms"))
# Short method names without settings (06's tree; 04's labels add the tuned values).
method_short_labels <- c(original = "Original", fixed = "Fixed MultiQT",
  ref_op = "Delphi REF-op (covid prod)", sqrt_constant = "Delphi sqrt 0.018 (flu prod)",
  absolute_slope = "Slope rate", oracle_target = "Target oracle (hindsight)",
  uniform = "Rate grid: equal", weighted = "Rate grid: learned, no sharing",
  weighted_common = "Rate grid: learned, common weights",
  shared = "Rate grid: learned + sharing",
  shared_common = "Rate grid: learned + sharing, common weights",
  shared_outer = "Rate grid: sharing + outer",
  shared_common_outer = "Rate grid: sharing, common weights + outer",
  restart_persistent = "First start only", restart_uniform = "Weekly starts: equal",
  restart_weighted = "Weekly starts: learned", restart_grid = "Starts x rate grid: learned",
  restart_uniform_outer = "Starts equal + outer", restart_weighted_outer = "Starts learned + outer",
  restart_grid_outer = "Starts x rate grid + outer", diverse_uniform = "Diverse experts: equal",
  diverse_shared = "Diverse experts: learned + sharing",
  diverse_shared_common = "Diverse experts: learned + sharing, common weights",
  diverse_soft_switch = "Diverse experts: soft switching",
  studentized = "Studentized median (residual RMS)",
  studentized_base = "Studentized median (base width)",
  median_constant = "Median, constant sigma",
  studentized_full = "Studentized full base (base width)",
  studentized_levels = "Studentized full base (level scales)",
  studentized_level_rms = "Studentized full base (level RMS)")
# Colours of 04's trade-off plots and 06's tree: a method takes the colour family of the
# first group of method_groups that lists it, one shade per method; common-weight and outer
# variants share their method's shade. Every method group needs a palette: methods of a
# group without one get no colour and are not drawn.
group_palettes <- c("Single trackers" = "Blues 3", "Rate-grid experts" = "Oranges",
  "Weekly starts: immediate admission" = "Greens 3", "Diverse experts" = "Purples 3",
  "Studentized median" = "Burg", "Studentized full base" = "BluGrn")
home_group <- local({
  home <- character()
  for (group in names(method_groups)) for (method in method_groups[[group]])
    if (!method %in% c("original", names(home))) home[method] <- group
  home
})
tradeoff_colours <- local({
  base_method <- sub("_outer$", "", sub("_common", "", names(home_group)))
  colours <- c(original = "black")
  for (group in names(group_palettes)) {
    members <- home_group == group
    bases <- unique(base_method[members])
    # Sequential palettes run dark to light; drop the two lightest shades.
    shades <- hcl.colors(length(bases) + 2L, group_palettes[[group]])[seq_along(bases)]
    colours[names(home_group)[members]] <- shades[match(base_method[members], bases)]
  }
  colours
})

# Issue the forecasts of the given methods for one case, with every location fitted
# independently: offsets, expert weights, lambda and loss history start fresh. Returns
# the predictions in counts / count_scale (original first) and what 04 plots and saves.
# cache (R/05_tuning.R): an environment for one case and fixed settings apart from the
# transform. Scales, expert banks and mixtures computed for one grid value are reused for
# the next; each key holds every setting the stored object depends on. NULL: no reuse.
# Re-run design (get_case(revisions = TRUE)): rerun_methods() below; weeks: the issue
# weeks it re-runs (default: all), e.g. only the scored ones.
apply_methods <- function(case, settings, params = method_parameters(), methods = method_names,
                          cache = NULL, weeks = NULL) {
  if (isTRUE(case$feedback_revised))
    return(rerun_methods(case, settings, params, methods, cache, weeks))
  par <- function(method, name) method_parameter(params, method, name)
  methods <- union(methods, outer_inner[intersect(names(outer_inner), methods)])
  want <- function(...) any(c(...) %in% methods)
  transform <- settings$transform
  memo <- function(key, value) {
    if (is.null(cache)) return(value)
    key <- paste(transform, key)
    if (is.null(cache[[key]])) cache[[key]] <- value
    cache[[key]]
  }
  scaling <- choose_count_scale(case, settings$scale_method)
  count_scale <- scaling$scale
  s0 <- if (!is.null(case$frozen_s0)) case$frozen_s0 else memo("s0", error_scale(case,
    settings$warmup_rounds, settings$scale_floor, transform, settings$replay_warmup))
  # Keep stored predictions and the scoring/weighting inputs in counts / count_scale.
  base <- case$base / count_scale
  feedback <- case$feedback_observed / count_scale
  adjust <- function(counts, eta, warmup_weeks) multiqt(counts, case$feedback_observed,
    case$feedback_by_week, case$alpha, eta, warmup_weeks, transform) / count_scale

  # One tracker: fixed rate, slope rate or studentized.
  # Fixed MultiQT: eta = m * S0, using the same frozen scale as the experts.
  predictions <- list(original = base)
  if (want("fixed", "oracle_target"))
    predictions$fixed <- adjust(case$base, par("fixed", "fixed_multiplier") * s0$scale,
      s0$warmup_weeks)
  # Delphi's REF-op (method_parameters()): Fixed MultiQT's update on its own sqrt tracker,
  # with Delphi's rolling rates from the rows learned so far instead of S0 (re-run: as
  # known at the issue week, as Delphi's prod re-runs), the slow term and the leak; no
  # warmup and no hindsight.
  if (want("ref_op")) {
    tracker <- par("ref_op", "ref_op_transform")
    rate <- function(multiplier, rounds, floor) adaptive_rate(case,
      par("ref_op", multiplier), rounds, floor, tracker)
    predictions$ref_op <- multiqt(case$base, case$feedback_observed, case$feedback_by_week,
      case$alpha, rate("ref_op_multiplier", par("ref_op", "ref_op_window"),
        par("ref_op", "ref_op_floor")), 0L, tracker,
      eta_slow = rate("ref_op_slow_multiplier", Inf, 0),
      fast_decay = par("ref_op", "ref_op_fast_decay")) / count_scale
  }
  # Delphi's flu prod (method_parameters()): Fixed MultiQT on sqrt(counts per 100k) with a
  # constant rate; no S0 and no warmup. Unchanged forecasts keep the exact base.
  if (want("sqrt_constant")) {
    per_100k <- choose_count_scale(case, "population")$scale
    tracked <- multiqt(case$base / per_100k, case$feedback_observed / per_100k,
      case$feedback_by_week, case$alpha, par("sqrt_constant", "sqrt_constant_rate"), 0L, "sqrt")
    q <- tracked * per_100k
    unchanged <- tracked == case$base / per_100k
    q[unchanged] <- case$base[unchanged]
    predictions$sqrt_constant <- q / count_scale
  }
  # Slope rate: only the current base forecast's weekly slope in tracker units, never
  # future counts, with a floor and a cap relative to S0.
  if (want("absolute_slope")) {
    median_tracker <- to_tracker(case$base[, match(.5, case$alpha)], transform)
    slope <- c(0, diff(median_tracker) / (as.numeric(diff(case$reference_date)) / 7))
    predictions$absolute_slope <- adjust(case$base,
      pmin(par("absolute_slope", "slope_cap_multiplier") * s0$scale,
        par("absolute_slope", "slope_min_multiplier") * s0$scale + abs(slope)), s0$warmup_weeks)
  }
  # Studentized median MultiQT (median_multiqt()): Fixed MultiQT's coverage indicator
  # and update on the base median b_t, issued as b_t + sigma_t * offsets.
  # Studentized scales stay within sigma_bounds * S0 (missing scales stay missing).
  bound <- function(scale) if (!s0$initialized) scale else
    pmin(pmax(scale, params$sigma_bounds[1] * s0$scale), params$sigma_bounds[2] * s0$scale)
  bounds_key <- paste(params$sigma_bounds, collapse = ",")
  sigma <- memo(paste("sigma", params$median_window, params$median_min_rows,
    params$residual_estimator, bounds_key), {
    sigma <- median_sigma(case, "residuals", params$median_window, params$median_min_rows,
      settings$scale_floor, transform, estimator = params$residual_estimator)
    sigma$scale <- bound(sigma$scale)
    sigma
  })
  base_sigma <- memo(paste("base_sigma", paste(params$median_base_levels, collapse = ","),
    bounds_key), {
    base_sigma <- median_sigma(case, "base", floor = settings$scale_floor,
      transform = transform, levels = params$median_base_levels)
    base_sigma$scale <- bound(base_sigma$scale)
    base_sigma
  })
  # The original base is issued until sigma_t exists (and during the S0 warmup without replay).
  median_warmup_weeks <- max(s0$warmup_weeks, sigma$warmup_weeks)
  # Offsets before the first update: normal quantiles, or t quantiles of a residual divided
  # by the RMS of median_window other normal residuals (wider tails for small windows).
  start <- switch(params$median_start, normal = qnorm(case$alpha),
    t = qt(case$alpha, params$median_window), stop("Unknown median start ", params$median_start))
  if (want("studentized"))
    predictions$studentized <- median_multiqt(case, par("studentized", "median_rate"),
      sigma$scale, start, median_warmup_weeks, transform) / count_scale
  # sigma_t from the width of the base: no feedback lag, so no warmup beyond S0's. This
  # sigma_t is a normal sd, so it starts from normal quantiles: the base's own interval.
  if (want("studentized_base"))
    predictions$studentized_base <- median_multiqt(case, par("studentized_base", "median_rate"),
      base_sigma$scale, qnorm(case$alpha), max(s0$warmup_weeks, base_sigma$warmup_weeks),
      transform) / count_scale
  # Residual sigma frozen at its first issued value: additive MultiQT on the median.
  # It separates reducing the base to its median from studentizing.
  constant_sigma <- sigma$scale[median_warmup_weeks + 1L]
  if (want("median_constant"))
    predictions$median_constant <- median_multiqt(case, par("median_constant", "median_rate"),
      constant_sigma, start, median_warmup_weeks, transform) / count_scale
  # Studentized MultiQT on the full base: offsets start at 0 (the base) and are issued as
  # b_t^alpha + sigma_t * offsets, with one scale per week (the base width; isoreg) or one
  # per level (weighted isoreg): the normal sd of the base's central interval at that
  # level, or the RMS of that level's recent residuals y - b^alpha.
  full <- function(method, scales, warmup_weeks = s0$warmup_weeks) multiqt(case$base,
    case$feedback_observed, case$feedback_by_week, case$alpha, par(method, "median_rate"),
    warmup_weeks, transform, sigma = scales) / count_scale
  if (want("studentized_full"))
    predictions$studentized_full <- full("studentized_full", base_sigma$scale)
  level_scales <- if (want("studentized_levels")) memo(paste("level_scales", bounds_key),
    bound(level_sigma(case, settings$scale_floor, transform)))
  if (want("studentized_levels"))
    predictions$studentized_levels <- full("studentized_levels", level_scales)
  # The residual scales exist once median_min_rows residuals are released, as for the median.
  if (want("studentized_level_rms"))
    predictions$studentized_level_rms <- full("studentized_level_rms",
      memo(paste("level_rms", params$median_window, params$median_min_rows,
        params$residual_estimator, bounds_key),
        bound(level_residual_sigma(case, params$median_window, params$median_min_rows,
          settings$scale_floor, transform, params$residual_estimator))), median_warmup_weeks)

  # What if we knew the true peak? Hindsight benchmark: bypass correction for postpeak
  # targets, even if issued prepeak (each week's nearest peak). Setting eta=0 instead
  # would only freeze the accumulated offset.
  if (want("oracle_target")) {
    predictions$oracle_target <- predictions$fixed
    after_peak <- which(case$target_date > case$peak_date)
    predictions$oracle_target[after_peak, ] <- base[after_peak, ]
  }

  # Several trackers: weight their forecasts. Rate grid: one persistent expert per
  # multiplier. Restarts: a fresh copy per multiplier at every issue week after warmup;
  # no permanent base expert. Diverse experts: the base plus slow MultiQT on differently
  # corrected bases. Each weighting learns from released pinball losses of issued
  # forecasts, divided by the base's mean loss over the S0 rows, so the count scale does
  # not change the weights. Their deterministic average still uses LEARNED weights;
  # equal averaging is a control.
  loss_scale <- if (!is.null(case$frozen_loss_scale)) case$frozen_loss_scale else
    weighting_loss_scale(case, s0, count_scale, settings$scale_floor)
  weigh <- function(method, experts, ...) exp_weights(experts, base, feedback,
    case$feedback_by_week, case$alpha, H = par(method, "H"),
    lambda_start = par(method, "lambda_start"), warmup_weeks = s0$warmup_weeks,
    loss_scale = loss_scale, ...)
  # Key of a mixture: the method and every setting a mixture can depend on.
  mixture_key <- function(method) paste("mixture", method, paste(vapply(c("H", "lambda_start",
    "restart_lambda", "diverse_slow_multiplier", "switch_strength", "switch_width",
    "switch_threshold", "switch_weight_path", "switch_weight_nowcast",
    "switch_weight_revision"), function(name) format(par(method, name), digits = 15), ""),
    collapse = ","), paste(params$restart_multipliers, collapse = ","), params$restart_multiplier,
    params$diverse_leak, paste(params$diverse_damping, collapse = ","),
    paste(params$diverse_widths, collapse = ","))
  mix <- function(method, value) memo(mixture_key(method), value)
  mixtures <- list()
  rate_experts <- restart_experts <- fixed_restarts <- NULL
  if (want("uniform", "weighted", "weighted_common", "shared", "shared_common")) {
    rate_experts <- memo("rate_experts", make_rate_experts(case, settings$expert_multipliers,
      s0, count_scale, transform))
    if (want("uniform")) mixtures$uniform <- mix("uniform", weigh("uniform", rate_experts,
      lambda = 0))
    if (want("weighted")) mixtures$weighted <- mix("weighted", weigh("weighted", rate_experts,
      sharing = FALSE))
    if (want("weighted_common")) mixtures$weighted_common <- mix("weighted_common",
      weigh("weighted_common", rate_experts, sharing = FALSE, per_quantile = FALSE))
    if (want("shared")) mixtures$shared <- mix("shared", weigh("shared", rate_experts))
    if (want("shared_common")) mixtures$shared_common <- mix("shared_common",
      weigh("shared_common", rate_experts, per_quantile = FALSE))
  }
  if (want("restart_uniform", "restart_weighted", "restart_grid", "restart_persistent")) {
    restart_experts <- memo(paste("restart_experts", paste(params$restart_multipliers,
      collapse = ",")), make_restart_experts(case, params$restart_multipliers, s0, count_scale,
      transform))
    fixed_restarts <- select_experts(restart_experts,
      restart_experts$experts$multiplier == params$restart_multiplier)
    if (want("restart_uniform")) mixtures$restart_uniform <- mix("restart_uniform",
      weigh("restart_uniform", fixed_restarts, lambda = 0))
    if (want("restart_weighted")) mixtures$restart_weighted <- mix("restart_weighted",
      weigh("restart_weighted", fixed_restarts, lambda = par("restart_weighted", "restart_lambda")))
    if (want("restart_grid")) mixtures$restart_grid <- mix("restart_grid",
      weigh("restart_grid", restart_experts, lambda = par("restart_grid", "restart_lambda")))
  }
  # One diverse bank per slow rate, shared by the methods that use it.
  diverse_banks <- list()
  diverse <- function(method) {
    m <- par(method, "diverse_slow_multiplier")
    key <- format(m, digits = 15)
    if (is.null(diverse_banks[[key]])) diverse_banks[[key]] <<- memo(paste("diverse", key,
      params$diverse_leak, paste(params$diverse_damping, collapse = ","),
      paste(params$diverse_widths, collapse = ",")), make_diverse_experts(case, s0,
      count_scale, transform, m * s0$scale, params$diverse_leak, params$diverse_damping,
      params$diverse_widths))
    diverse_banks[[key]]
  }
  if (want("diverse_uniform")) mixtures$diverse_uniform <- mix("diverse_uniform",
    weigh("diverse_uniform", diverse("diverse_uniform"), lambda = 0))
  if (want("diverse_shared")) mixtures$diverse_shared <- mix("diverse_shared",
    weigh("diverse_shared", diverse("diverse_shared")))
  if (want("diverse_shared_common")) mixtures$diverse_shared_common <-
    mix("diverse_shared_common", weigh("diverse_shared_common", diverse("diverse_shared_common"),
      per_quantile = FALSE))
  soft_switch <- NULL
  if (want("diverse_soft_switch")) {
    soft_switch <- mix("diverse_soft_switch", {
      gate <- transition_gate(case, par("diverse_soft_switch", "switch_strength"),
        par("diverse_soft_switch", "switch_width"), par("diverse_soft_switch", "switch_threshold"),
        c(par("diverse_soft_switch", "switch_weight_path"),
          par("diverse_soft_switch", "switch_weight_nowcast"),
          par("diverse_soft_switch", "switch_weight_revision")))
      switched <- soft_switch_weights(diverse("diverse_soft_switch"), base, feedback,
        case$feedback_by_week, case$alpha, gate$gate, H = par("diverse_soft_switch", "H"),
        lambda_start = par("diverse_soft_switch", "lambda_start"),
        warmup_weeks = s0$warmup_weeks, loss_scale = loss_scale)
      switched$gate <- gate
      switched
    })
    mixtures$diverse_soft_switch <- soft_switch
  }
  for (method in names(mixtures)) predictions[[method]] <- mixtures[[method]]$q
  # The first-start control has exactly the same birth rule as the restart copies.
  if (want("restart_persistent")) {
    predictions$restart_persistent <- base
    if (length(fixed_restarts$birth_week)) {
      issued <- seq.int(fixed_restarts$birth_week[1], nrow(base))
      predictions$restart_persistent[issued, ] <- fixed_restarts$forecasts[issued, 1, ]
    }
  }
  # Each outer tracker has its own offsets and learns from its OWN final forecasts.
  # Its state stays continuous when new restart copies join the underlying mixture.
  for (method in intersect(names(outer_inner), methods))
    predictions[[method]] <- adjust(predictions[[outer_inner[[method]]]] * count_scale,
      par(method, "outer_multiplier") * s0$scale, s0$warmup_weeks)
  predictions <- predictions[intersect(method_names, names(predictions))]
  # expert_forecasts: the experts each mixture weighs [week x expert x level] in counts /
  # count_scale (04's weight pages; the diverse banks differ by slow rate).
  list(case = case, scaling = scaling, s0 = s0, sigma = sigma, base_sigma = base_sigma,
    level_scales = level_scales, median_warmup_weeks = median_warmup_weeks,
    constant_sigma = constant_sigma, predictions = predictions,
    weights = lapply(mixtures, `[[`, "weights"), lambda = lapply(mixtures, `[[`, "lambda"),
    expert_forecasts = lapply(mixtures, `[[`, "experts"),
    rate_experts = rate_experts, restart_experts = restart_experts,
    diverse_experts = if (length(diverse_banks)) diverse_banks[[1]], soft_switch = soft_switch)
}

# Re-run design, as in Delphi's calibration reports: at every issue week t of `weeks`
# (default: all), every method is re-run from the first week with the reports current at
# t's cutoff (case_as_of()): the same rows are consumed in the same rounds as in the
# stateful design with this delay, but each with its report as known at t. Only the week-t
# forecasts, weights and expert forecasts are issued; the earlier forecasts of a re-run are
# reconstructions and generally differ from those issued then. S0 and the weighting loss
# scale are those of the stateful design with this delay (from its stateful feedback),
# frozen for every re-run. Rows outside `weeks` are NA; experts not yet started have
# weight 0. The other results (scales, experts) are those of the run at the last week;
# rerun_omitted: the rows left out of the re-run at each week because their current report
# was invalid.
# cache: one environment per week (see apply_methods()).
rerun_methods <- function(case, settings, params, methods, cache = NULL, weeks = NULL) {
  n <- nrow(case$base)
  weeks <- sort(if (is.null(weeks)) seq_len(n) else weeks)
  stateful <- case
  stateful$feedback_revised <- FALSE
  s0 <- error_scale(stateful, settings$warmup_rounds, settings$scale_floor,
    settings$transform, settings$replay_warmup)
  loss_scale <- weighting_loss_scale(stateful, s0,
    choose_count_scale(case, settings$scale_method)$scale, settings$scale_floor)
  run_at <- function(t) {
    week_cache <- NULL
    if (!is.null(cache)) {
      key <- paste("re-run week", t)
      if (is.null(cache[[key]])) cache[[key]] <- new.env()
      week_cache <- cache[[key]]
    }
    apply_methods(case_as_of(case, t, s0, loss_scale), settings, params, methods, week_cache)
  }
  # n rows of NA with the other dimensions of x (one row per issue week).
  blank <- function(x) {
    if (is.data.frame(x)) {
      x <- x[rep(NA_integer_, n), , drop = FALSE]
      rownames(x) <- NULL
      return(x)
    }
    d <- if (is.null(dim(x))) c(length(x), 1L) else dim(x)
    array(NA_real_, c(n, d[-1]), c(list(NULL), dimnames(x)[-1]))
  }
  # Row t of x from the run at t (experts matched by name: weekly restarts join over time);
  # empty: the value of the experts missing at t (weight 0, forecast NA).
  fill <- function(target, x, t, empty = 0) {
    if (is.data.frame(target)) {
      target[t, ] <- x[t, ]
    } else if (length(dim(target)) == 3L) {
      target[t, , ] <- empty
      target[t, match(dimnames(x)[[2]], dimnames(target)[[2]]), ] <- x[t, , ]
    } else {
      target[t, ] <- if (is.null(dim(x))) x[t] else x[t, ]
    }
    target
  }
  last <- run_at(max(weeks))
  parts <- c("predictions", "weights", "lambda", "expert_forecasts")
  switch_parts <- c("q", "weights", "lambda", "regime_weights", "regime_lambda", "gate")
  issued <- lapply(last[parts], function(part) lapply(part, blank))
  switched <- if (!is.null(last$soft_switch)) lapply(last$soft_switch[switch_parts], blank)
  omitted <- list()
  for (t in weeks) {
    run <- if (t == max(weeks)) last else run_at(t)
    for (part in parts) for (method in names(issued[[part]]))
      issued[[part]][[method]] <- fill(issued[[part]][[method]], run[[part]][[method]], t,
        if (part == "expert_forecasts") NA else 0)
    for (part in names(switched)) switched[[part]] <- fill(switched[[part]],
      run$soft_switch[[part]], t)
    if (length(run$case$rerun_omitted))
      omitted[[length(omitted) + 1L]] <- data.frame(week = t, row = run$case$rerun_omitted)
  }
  last[parts] <- issued
  if (!is.null(switched)) last$soft_switch[switch_parts] <- switched
  last$case <- case
  last$rerun_omitted <- do.call(rbind, c(list(data.frame(week = integer(), row = integer())),
    omitted))
  last$frozen <- list(s0 = s0, loss_scale = loss_scale)
  last
}

# 8. SAVED RESULTS OF LONG RUNS --------------------------------------------
# R/04 and R/05 save every finished piece of work (one location's run, one tuning stage) in
# output/cache/<name>/, under a key of its arguments and of everything the function reads:
# the code of every function it calls (without comments) and the values of the global
# objects these use (codetools::findGlobals()). A restarted run, after a crash or a change,
# loads every piece whose key is unchanged and recomputes the others; changes to plots,
# colours or labels recompute nothing. Not in the key: the R and package versions and files
# read inside the function. Delete output/cache/ to recompute everything.
cache_dir <- "output/cache"
cache_ignore <- "cores"                   # Global objects that never change a result.
# The plotted locations (settings$location_ids) only choose pages: left out of keys.
key_value <- function(x) {
  if (is.list(x) && all(c("location_ids", "summary_exclude") %in% names(x)))
    x$location_ids <- NULL
  x
}
# The code and the global values a function reads, followed through the functions it calls
# (also functions inside lists).
code_key <- function(fun) {
  key <- list(deparse(fun))
  queue <- list(fun)
  while (length(queue)) {
    names_used <- setdiff(codetools::findGlobals(queue[[1]]), c(names(key), cache_ignore))
    queue <- queue[-1]
    for (name in names_used) if (exists(name, globalenv(), inherits = FALSE)) {
      value <- get(name, globalenv())
      functions <- if (is.function(value)) list(value) else if (is.list(value))
        Filter(is.function, value) else list()
      key[[name]] <- if (is.function(value)) deparse(value) else
        rapply(list(key_value(value)), deparse, classes = "function", how = "replace")
      queue <- c(queue, functions)
    }
  }
  key
}
# fun(...), loaded from output/cache/<name>/ if saved under the same key, else computed and
# saved (written to a temporary file first, so an interrupted write leaves no result).
cached <- function(name, fun, ...) {
  key <- tools::md5sum(bytes = serialize(list(code_key(fun), lapply(list(...), key_value)),
    NULL))
  file <- file.path(cache_dir, name, paste0(key, ".rds"))
  if (file.exists(file)) return(readRDS(file))
  value <- fun(...)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  temporary <- file.path(dirname(file), paste0("tmp_", Sys.getpid(), "_", basename(file)))
  saveRDS(value, temporary)
  file.rename(temporary, file)
  value
}
