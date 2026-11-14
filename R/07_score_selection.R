# Which score of the base forecast keeps its alpha-quantile stable over time? A tracker
# learns the alpha-quantile of its score; if that quantile stays put (quantile stability),
# a constant threshold stays calibrated and a small rate is enough. The candidates
# (score_candidates in R/03_functions.R; notes, "Scores and scales tested") are compared
# with three diagnostics, pooled over the 10 largest states without Florida
# (settings$pooled_locations) and the horizons, against the latest reports:
#   drift ratio: the hindsight alpha-quantile of each score is the best constant threshold;
#     the mean |coverage - alpha| of its hits in windows of 8 consecutive weeks, divided by
#     the same with the weeks permuted (1: no drift beyond chance);
#   peak-aligned coverage of that threshold, by week relative to the nearest peak;
#   running quantile: each week the empirical alpha-quantile of all scores reported before
#     the cutoff, inverted into a forecast (static recalibration, no rate), scored like R/04;
#     also with the scores of the same phase only, the phase taken from the peaks in
#     hindsight (what two scores, one for the rise and one for the decline, could give);
#   up to and after the peak: each score's alpha-quantile in the weeks up to the peak and
#     in the weeks after it (does the threshold change sign?);
#   MultiQT on each candidate's score at the adaptive rate (as the methods of R/04).
# Also draws each candidate's scores over time (multiqt_score_stationarity/).
# Run from the project root after R/01_data_overview.R: Rscript R/07_score_selection.R
source("R/03_functions.R")

# 1. SETTINGS AND PREPARED DATA -------------------------------------------------
settings <- multiqt_settings()
prepared <- read_prepared(settings)       # The base model's prepared data and output folder.
params <- method_parameters()
diseases <- c("covid", "flu")
disease_names <- c(covid = "COVID", flu = "Flu")
drift_window <- 8L                        # Weeks per window of the drift ratio.
permutations <- 30L                       # Permuted orders for its noise.
min_scores <- 5L                          # Running quantile: the base until 5 scores exist.
recent_scores <- 26L                      # Its window variant: the last 26 scores only.
weeks_relative <- -7:8                    # Peak-aligned coverage.
cores <- max(1L, parallel::detectCores() - 2L)
result_dir <- file.path(settings$output_dir, "results", "score_selection")
figure_dir <- file.path(settings$output_dir, "pdf")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figure_dir, "multiqt_score_stationarity"), recursive = TRUE,
  showWarnings = FALSE)
candidates <- score_candidates
origin_colours <- c("Trackers in use" = "black", "Level of the forecast" = "#d95f02",
  "Recent errors" = "#7570b3", "Spread of the forecast" = "#1b9e77",
  "Revisions of the forecast" = "#66a61e", "Uniejewski (2025)" = "#a6761d")
candidate_colours <- setNames(origin_colours[candidates$origin], candidates$id)
labels <- setNames(candidates$label, candidates$id)
# The candidates of the peak-aligned pages and of the scores over time.
shown <- c("identity", "sqrt", "log1p", "power_50", "power_75", "level_rms", "base_width",
  "spread", "pit", "revisions_target", "window_asinh", "tpit_25")
levels_shown <- c(.1, .5, .9)
phase_window <- list(before = (1L - settings$weeks_before_peak):0L,
  after = seq_len(settings$weeks_after_peak))
level_colours <- c("#0072B2", "black", "#D55E00")

# 2. THE THREE DIAGNOSTICS OF ONE CASE ------------------------------------------
# Drift ratio of a score matrix [week x level]: drift and noise are the mean over the windows
# of drift_window consecutive weeks and over the levels of |coverage - alpha| of the hindsight
# threshold, with the weeks in their order and permuted.
drift_ratio <- function(S, alpha) {
  n <- nrow(S)
  if (n < 2L * drift_window) return(NA_real_)
  u <- vapply(seq_along(alpha), function(j) quantile(S[, j], alpha[j], type = 1,
    names = FALSE), numeric(1))
  covered <- sweep(S, 2, u, "<=") * 1
  windowed <- function(x) {
    sums <- rbind(0, apply(x, 2, cumsum))
    rolling <- (sums[(drift_window + 1L):(n + 1L), , drop = FALSE] -
      sums[seq_len(n - drift_window + 1L), , drop = FALSE]) / drift_window
    mean(abs(sweep(rolling, 2, alpha)))
  }
  windowed(covered) / mean(replicate(permutations,
    windowed(covered[sample.int(n), , drop = FALSE])))
}
# The scores of the reports of a case under one candidate (rows in release order).
scores_of <- function(candidate, case) {
  released <- unlist(case$feedback_by_week)
  candidate_score(candidate, case, case$feedback_observed[released], released)
}
# Running-quantile forecasts [week x level] of one candidate: at each issue week the
# empirical alpha-quantile of the scores reported by then (each with the report learned and
# the scale of its own cutoff, as MultiQT compares a report with its historical forecast),
# inverted; the base until min_scores scores exist. window: Inf or the last `window` scores.
# phase (optional): a label per row; only the past scores of rows with week t's label count.
running_quantile <- function(candidate, case, window = Inf, phase = NULL) {
  n_released <- cumsum(lengths(case$feedback_by_week))
  released <- unlist(case$feedback_by_week)
  past <- scores_of(candidate, case)
  q <- case$base
  for (t in seq_len(nrow(q))) {
    k <- n_released[t]
    if (k < min_scores) next
    use <- if (is.finite(window)) max(1L, k - window + 1L):k else seq_len(k)
    if (!is.null(phase)) use <- use[phase[released[use]] == phase[t]]
    S <- past[use, , drop = FALSE]
    S <- S[apply(is.finite(S), 1, all), , drop = FALSE]
    if (nrow(S) < min_scores) next
    u <- vapply(seq_along(case$alpha), function(j) quantile(S[, j], case$alpha[j], type = 1,
      names = FALSE), numeric(1))
    forecast <- candidate_forecast(candidate, case, u, t)
    if (all(is.finite(forecast))) q[t, ] <- forecast
  }
  q
}
# MultiQT on a candidate's score (as the trackers of R/04): thresholds u, one per level,
# start at 0; the forecast of week t inverts S_t^alpha(y) <= u^alpha (candidate_forecast(),
# levels sorted; the base where the score is not defined yet), and a report on row s moves u
# by eta_t (alpha - covered_s). eta: Delphi's adaptive rate of the candidate's own scores
# (rolling_rate(), adaptive_multiplier and adaptive_rounds of method_parameters()).
candidate_multiqt <- function(candidate, case) {
  eta <- params$adaptive_multiplier * rolling_rate(case, function(rows)
    candidate_score(candidate, case, case$feedback_observed[rows], rows), 1,
    params$adaptive_rounds, settings$replay_warmup)
  u <- rep(0, length(case$alpha))
  q <- case$base
  for (week in seq_len(nrow(q))) {
    for (s in case$feedback_by_week[[week]])
      u <- u + eta[week] * (case$alpha - (case$feedback_observed[s] <= q[s, ]))
    if (any(u != 0)) {
      forecast <- candidate_forecast(candidate, case, u, week)
      if (all(is.finite(forecast))) q[week, ] <- forecast
    }
  }
  q
}
# One location and horizon of a disease: the drift ratio and the hindsight peak-aligned hits
# of every candidate (delay 0, all scored weeks), the running quantile under delay 0 and 1
# (evaluation weeks), and the cases' scores for the plots.
evaluate_case <- function(disease, location, horizon) {
  rows <- list(drift = list(), peak = list(), running = list(), phase = list())
  for (delay in c(0L, 1L)) {
    case <- get_case(prepared, disease, location, horizon, delay, FALSE)
    s0 <- error_scale(case, settings$warmup_rounds, settings$scale_floor, settings$transform,
      settings$replay_warmup)
    scores <- candidate_scores(case, s0, params)
    truth <- score_truth(case, settings$score_target, settings$score_filled)
    alpha <- case$alpha
    if (delay == 0L) {
      # Every candidate on the same weeks: a finite score for all of them (the own-target
      # revisions, undefined at the longest horizon, aside) and no filled round.
      y <- score_truth(case, "latest", FALSE)
      S <- lapply(scores, function(x) candidate_score(x, case, y, seq_len(nrow(case$base))))
      finite <- function(x) apply(is.finite(x), 1, all)
      keep <- is.finite(y) & Reduce(`&`, lapply(S[names(S) != "revisions_target"], finite))
      if (horizon < max(settings$horizons)) keep <- keep & finite(S$revisions_target)
      for (id in names(S)) {
        x <- S[[id]][keep, , drop = FALSE]
        ok <- all(is.finite(x))
        rows$drift[[id]] <- data.frame(disease, location, horizon, candidate = id,
          weeks = sum(keep), ratio = if (ok) drift_ratio(x, alpha) else NA_real_)
        if (!ok || sum(keep) < 2L * drift_window) next
        u <- vapply(seq_along(alpha), function(j) quantile(x[, j], alpha[j], type = 1,
          names = FALSE), numeric(1))
        # The score's alpha-quantile up to the peak and after it, in units of its median
        # absolute value (all weeks and levels): a change of sign is a bias that flips.
        spread <- median(abs(x))
        week <- case$relative_week[keep]
        phase_quantile <- function(w) if (sum(week %in% w) < 4L) rep(NA_real_, length(alpha))
          else vapply(seq_along(alpha), function(j) quantile(x[week %in% w, j], alpha[j],
            type = 1, names = FALSE), numeric(1)) / spread
        rows$phase[[id]] <- data.frame(disease, location, horizon, candidate = id, alpha,
          before = phase_quantile(phase_window$before), after = phase_quantile(phase_window$after))
        r <- which(case$relative_week[keep] %in% weeks_relative)
        rows$peak[[id]] <- data.frame(disease, location, horizon, candidate = id,
          threshold = "hindsight constant", week = rep(case$relative_week[keep][r],
            length(alpha)), alpha = rep(alpha, each = length(r)),
          covered = as.vector(sweep(x[r, , drop = FALSE], 2, u, "<=")))
      }
      plot_case <- list(case = case, scores = S)
    }
    count_scale <- choose_count_scale(case, settings$scale_method)$scale
    base <- score_quantiles(case$base, truth, alpha, count_scale)
    windows <- lapply(setNames(nm = c("overall", "before", "after")), function(phase)
      which(is.finite(truth) & phase_weeks(case, phase, settings$weeks_before_peak,
        settings$weeks_after_peak)))
    # The phase of each row in hindsight: the decline after the nearest peak, else the rise.
    row_phase <- ifelse(is.finite(case$relative_week) & case$relative_week > 0, "decline",
      "rise")
    variants <- c("all past scores", "last 26 scores", "same phase (hindsight)", "MultiQT")
    for (id in names(scores)) for (variant in variants) {
      # A candidate without any score in this case (the own-target revisions at the longest
      # horizon) would issue the base throughout: left out, not counted as no change.
      if (!any(is.finite(scores_of(scores[[id]], case)))) next
      q <- switch(variant, "all past scores" = running_quantile(scores[[id]], case),
        "last 26 scores" = running_quantile(scores[[id]], case, recent_scores),
        "same phase (hindsight)" = running_quantile(scores[[id]], case, phase = row_phase),
        MultiQT = candidate_multiqt(scores[[id]], case))
      s <- score_quantiles(q, truth, alpha, count_scale)
      for (phase in names(windows)) {
        r <- windows[[phase]]
        if (!length(r)) next
        rows$running[[length(rows$running) + 1L]] <- data.frame(disease, location, horizon,
          design = paste("delay", delay), candidate = id, variant, phase, weeks = length(r),
          loss = mean(s$loss[r, ]), base_loss = mean(base$loss[r, ]),
          ce = mean(abs(colMeans(s$covered[r, , drop = FALSE]) - alpha)),
          base_ce = mean(abs(colMeans(base$covered[r, , drop = FALSE]) - alpha)))
      }
      # Out of sample: the running quantile's hits by week relative to the peak.
      r <- which(is.finite(truth) & case$period == "evaluation" &
        case$relative_week %in% weeks_relative)
      if (delay == 0L && variant == "all past scores" && length(r))
        rows$peak[[paste(id, "running")]] <- data.frame(disease, location, horizon,
          candidate = id, threshold = "running quantile", week = rep(case$relative_week[r],
            length(alpha)), alpha = rep(alpha, each = length(r)),
          covered = as.vector(s$covered[r, , drop = FALSE]))
    }
  }
  list(drift = do.call(rbind, rows$drift), peak = do.call(rbind, rows$peak),
    running = do.call(rbind, rows$running), phase = do.call(rbind, rows$phase),
    plot = plot_case)
}

# 3. ALL POOLED CASES -----------------------------------------------------------
cases <- expand.grid(disease = diseases, location = union(settings$pooled_locations,
  settings$location_ids), horizon = settings$horizons, stringsAsFactors = FALSE)
results <- parallel::mclapply(seq_len(nrow(cases)), function(i) {
  set.seed(20261113L + i)                 # Reproducible permutations in every worker.
  with(cases[i, ], evaluate_case(disease, location, horizon))
}, mc.cores = cores, mc.preschedule = FALSE)
failed <- vapply(results, inherits, logical(1), "try-error")
if (any(failed)) stop(results[[which(failed)[1]]])
pooled_case <- cases$location %in% settings$pooled_locations
drift <- do.call(rbind, lapply(results[pooled_case], `[[`, "drift"))
peak <- do.call(rbind, lapply(results[pooled_case], `[[`, "peak"))
running <- do.call(rbind, lapply(results[pooled_case], `[[`, "running"))
write.csv(drift, file.path(result_dir, "drift_by_case.csv"), row.names = FALSE)
write.csv(running, file.path(result_dir, "running_quantile_by_case.csv"), row.names = FALSE)
phases <- do.call(rbind, lapply(results[pooled_case], `[[`, "phase"))
write.csv(phases, file.path(result_dir, "phase_quantiles_by_case.csv"), row.names = FALSE)

# Drift ratio per candidate: mean over the pooled cases, the share of cases with less drift
# than the identity score, and the mean over the horizons with own-target revisions.
identity_ratio <- drift[drift$candidate == "identity", c("disease", "location", "horizon",
  "ratio")]
drift <- merge(drift, identity_ratio, by = c("disease", "location", "horizon"),
  suffixes = c("", "_identity"))
drift_summary <- merge(aggregate(cbind(ratio, less_drift = ratio < ratio_identity) ~ disease +
  candidate, drift, mean), aggregate(cbind(cases = ratio) ~ disease + candidate, drift, length))
short <- aggregate(ratio ~ disease + candidate, drift[drift$horizon < max(settings$horizons), ],
  mean)
drift_summary$ratio_short_horizons <- short$ratio[match(paste(drift_summary$disease,
  drift_summary$candidate), paste(short$disease, short$candidate))]
drift_summary <- merge(candidates[c("id", "label", "origin", "methods")], drift_summary,
  by.x = "id", by.y = "candidate")
drift_summary <- drift_summary[order(drift_summary$disease, drift_summary$ratio), ]
write.csv(drift_summary, file.path(result_dir, "drift.csv"), row.names = FALSE)
drift_by_horizon <- aggregate(ratio ~ disease + horizon + candidate, drift, mean)
write.csv(drift_by_horizon, file.path(result_dir, "drift_by_horizon.csv"), row.names = FALSE)

# Peak-aligned coverage, mean over cases and peaks, and its drift: the mean over the weeks
# and levels of |coverage - alpha|.
peak_coverage <- aggregate(covered ~ disease + candidate + threshold + week + alpha, peak, mean)
write.csv(peak_coverage, file.path(result_dir, "peak_coverage.csv"), row.names = FALSE)
peak_drift <- aggregate(cbind(drift = abs(covered - alpha)) ~ disease + candidate + threshold,
  peak_coverage, mean)
write.csv(peak_drift, file.path(result_dir, "peak_drift.csv"), row.names = FALSE)

# Up to and after the peak: the mean over cases of each score's normalized alpha-quantile in
# both windows, and per candidate the share of cases and levels in which the two have
# opposite signs (the threshold flips at the peak) and their mean difference.
phase_quantiles <- aggregate(cbind(before, after) ~ disease + candidate + alpha, phases, mean)
write.csv(phase_quantiles, file.path(result_dir, "phase_quantiles.csv"), row.names = FALSE)
phase_flips <- aggregate(cbind(flip = sign(before) != sign(after) & before != 0 & after != 0,
  gap = before - after) ~ disease + candidate, phases, mean)
write.csv(phase_flips, file.path(result_dir, "phase_flips.csv"), row.names = FALSE)

# Running quantile: pooled change of the mean loss and CE versus the base (%), as R/04.
pooled_running <- aggregate(cbind(loss, base_loss, ce, base_ce) ~ disease + design + candidate +
  variant + phase, running, mean)
pooled_running$loss_change <- 100 * (pooled_running$loss / pooled_running$base_loss - 1)
pooled_running$ce_change <- 100 * (pooled_running$ce / pooled_running$base_ce - 1)
write.csv(pooled_running, file.path(result_dir, "running_quantile.csv"), row.names = FALSE)

sink(file.path(result_dir, "summary.txt"), split = TRUE)
options(width = 200)
cat("Drift ratio (1: no drift beyond chance), the", pooled_label(
  length(settings$pooled_locations), settings), "and horizons",
  paste(range(settings$horizons), collapse = " to "), "\n")
for (d in diseases) {
  cat("\n==", disease_names[[d]], "==\n")
  x <- drift_summary[drift_summary$disease == d, ]
  print(data.frame(candidate = x$label, ratio = round(x$ratio, 3),
    short_horizons = round(x$ratio_short_horizons, 3), less_drift = round(x$less_drift, 2),
    cases = x$cases), row.names = FALSE)
}
cat("\nPeak-aligned drift, mean |coverage - alpha| over weeks -7 to 8 and the levels:\n")
print(reshape(transform(peak_drift, drift = round(drift, 3)), idvar = c("disease", "candidate"),
  timevar = "threshold", direction = "wide"), row.names = FALSE)
cat("\nUp to and after the peak: share of cases and levels with a change of sign of the",
  "score's alpha-quantile, and the mean gap (in units of the median |score|):\n")
print(transform(phase_flips, flip = round(flip, 2), gap = round(gap, 2))[order(
  phase_flips$disease, -phase_flips$flip), ], row.names = FALSE)
for (variant in unique(pooled_running$variant)) {
  cat("\n", variant, ": pooled change versus the base (%):\n", sep = "")
  x <- pooled_running[pooled_running$variant == variant, ]
  x <- reshape(transform(x, loss = round(loss_change, 1), ce = round(ce_change, 1))[,
    c("disease", "design", "candidate", "phase", "loss", "ce")], idvar = c("disease",
    "design", "candidate"), timevar = "phase", direction = "wide")
  print(x[order(x$disease, x$design, x$loss.overall), ], row.names = FALSE)
}
sink()

# 4. PLOTS: multiqt_score_selection.pdf ----------------------------------------
cairo_pdf(file.path(figure_dir, "multiqt_score_selection.pdf"), width = 16, height = 11,
  onefile = TRUE, family = "DejaVu Sans")
# Page 1: the drift ratio of every candidate, both diseases; the bar spans the middle half of
# the cases.
order_ids <- with(aggregate(ratio ~ id, drift_summary, mean), id[order(ratio)])
par(mfrow = c(1, 2), mar = c(4.5, 16, 3, 1), oma = c(4, 0, 3, 0))
for (d in diseases) {
  x <- drift_summary[drift_summary$disease == d, ]
  x <- x[match(rev(order_ids), x$id), ]
  spread <- t(vapply(x$id, function(id) quantile(drift$ratio[drift$disease == d &
    drift$candidate == id], c(.25, .75), na.rm = TRUE), numeric(2)))
  plot(x$ratio, seq_len(nrow(x)), xlim = range(1, spread, na.rm = TRUE), yaxt = "n",
    pch = 16, col = candidate_colours[x$id], xlab = "Drift ratio (1: no drift beyond chance)",
    ylab = "", main = disease_names[[d]])
  segments(spread[, 1], seq_len(nrow(x)), spread[, 2], seq_len(nrow(x)),
    col = adjustcolor(candidate_colours[x$id], .5), lwd = 3)
  axis(2, at = seq_len(nrow(x)), labels = labels[x$id], las = 1, cex.axis = .8)
  abline(v = x$ratio[x$id == "identity"], col = "gray60", lty = 2)
}
legend_origins <- function() {
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("bottom", names(origin_colours), col = origin_colours, pch = 16, horiz = TRUE,
    bty = "n", cex = .9, inset = .01)
}
mtext(paste("Drift of the hindsight constant threshold of each score: the",
  pooled_label(length(settings$pooled_locations), settings), "and horizons",
  paste(range(settings$horizons), collapse = " to ")), outer = TRUE, font = 2, cex = 1.2)
mtext(paste("Dot: mean over the cases; bar: their middle half; dashed: the identity score.",
  "Every score of a case on the same weeks; own-target revisions only below the longest horizon."),
  side = 1, outer = TRUE, line = 1, cex = .8)
legend_origins()
# Page 2: the drift ratio against the power p of the forecast median, (1 + b^.5)^p.
powers <- c(power_m100 = -1, identity = 0, power_25 = .25, power_50 = .5, power_75 = .75,
  power_100 = 1, power_150 = 1.5, power_200 = 2)
references <- c("log1p", "sqrt", "level_rms", "base_width", "spread", "pit", "window_asinh",
  "tpit_25")
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 12), oma = c(4, 0, 3, 0), xpd = FALSE)
for (d in diseases) {
  x <- drift_summary[drift_summary$disease == d, ]
  ratio <- x$ratio[match(names(powers), x$id)]
  reference <- x$ratio[match(references, x$id)]
  plot(powers, ratio, type = "o", pch = 16, col = origin_colours[["Level of the forecast"]],
    lwd = 2, ylim = range(ratio, reference, na.rm = TRUE), xlab = "Power p of (1 + b^0.5)^p",
    ylab = "Drift ratio", main = disease_names[[d]])
  abline(h = reference, col = candidate_colours[references], lty = 3)
  axis(4, at = reference, labels = labels[references], las = 1, cex.axis = .7, tick = FALSE,
    line = -.5)
}
mtext("Drift ratio of the scale (1 + b^0.5)^p against p (p = 0: the identity); dotted: other scores",
  outer = TRUE, font = 2, cex = 1.2)
# Pages 3-4: coverage by week relative to the peak, hindsight threshold and running quantile.
for (d in diseases) {
  par(mfrow = c(3, 4), mar = c(3.8, 4, 3, 1), oma = c(4, 0, 3, 0), xpd = FALSE)
  for (id in shown) {
    x <- peak_coverage[peak_coverage$disease == d & peak_coverage$candidate == id, ]
    drifts <- peak_drift$drift[peak_drift$disease == d & peak_drift$candidate == id]
    plot(range(weeks_relative), c(0, 1), type = "n", xlab = "Week relative to the peak",
      ylab = "Coverage", main = sprintf("%s\n(drift %s)", labels[[id]],
        paste(sprintf("%.3f", drifts), collapse = " / ")), cex.main = .9)
    abline(v = .5, col = "gray70", lty = 3)
    for (j in seq_along(levels_shown)) {
      abline(h = levels_shown[j], col = level_colours[j], lty = 3)
      for (threshold in c("hindsight constant", "running quantile")) {
        z <- x[x$alpha == levels_shown[j] & x$threshold == threshold, ]
        z <- z[order(z$week), ]
        lines(z$week, z$covered, col = level_colours[j], lwd = 1.6,
          lty = if (threshold == "hindsight constant") 1 else 2)
      }
    }
  }
  mtext(paste0(disease_names[[d]], ": coverage by week relative to the peak, the ",
    pooled_label(length(settings$pooled_locations), settings), ", horizons ",
    paste(range(settings$horizons), collapse = " to ")), outer = TRUE, font = 2, cex = 1.2)
  mtext(paste("Levels 0.1 (blue), 0.5 (black), 0.9 (red); dotted: the target. Solid: the",
    "hindsight constant threshold (all weeks); dashed: the running quantile (evaluation weeks,",
    "delay 0).\nDrift: mean |coverage - alpha| over the weeks and levels, solid / dashed."),
    side = 1, outer = TRUE, line = 2, cex = .8)
}
# Pages 5-6: the running quantile, loss change against CE change, overall and after peaks.
for (phase in c("overall", "after")) {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), oma = c(4, 0, 3, 0), xpd = FALSE)
  for (d in diseases) for (design in c("delay 0", "delay 1")) {
    x <- pooled_running[pooled_running$variant == "all past scores" &
      pooled_running$disease == d & pooled_running$design == design &
      pooled_running$phase == phase, ]
    plot(x$loss_change, x$ce_change, pch = 16, col = candidate_colours[x$candidate],
      xlab = "Loss change versus the base (%)", ylab = "CE change versus the base (%)",
      main = paste0(disease_names[[d]], ", ", design), xlim = extendrange(x$loss_change, f = .3),
      ylim = extendrange(x$ce_change, f = .1))
    abline(h = 0, v = 0, col = "gray75", lty = 2)
    text(x$loss_change, x$ce_change, labels[x$candidate], pos = 4, cex = .6,
      col = candidate_colours[x$candidate])
  }
  mtext(paste0("Running quantile of each score's past values, ", c(overall = "overall",
    after = paste("after peaks (weeks 1 to", settings$weeks_after_peak, ")"))[[phase]],
    ": change versus the base, the ", pooled_label(length(settings$pooled_locations),
      settings)), outer = TRUE, font = 2, cex = 1.2)
  mtext(paste("Each week the threshold is the empirical alpha-quantile of all scores reported",
    "before the cutoff; no rate. Evaluation weeks. Lower left is better."), side = 1,
    outer = TRUE, line = 1, cex = .8)
  legend_origins()
}
# Pages 7-8: each score's alpha-quantile up to the peak (filled) and after it (open), in
# units of its median |score|, at three levels; a change of sign is a bias that flips.
for (d in diseases) {
  par(mfrow = c(1, 3), mar = c(4.5, if (TRUE) 15 else 4, 3, 1), oma = c(4, 0, 3, 0))
  x <- phase_quantiles[phase_quantiles$disease == d, ]
  ids <- rev(order_ids)
  for (j in seq_along(levels_shown)) {
    z <- x[abs(x$alpha - levels_shown[j]) < 1e-9, ]
    z <- z[match(ids, z$candidate), ]
    plot(c(z$before, z$after), rep(seq_along(ids), 2), type = "n", yaxt = "n",
      xlab = "alpha-quantile of the score / median |score|", ylab = "",
      main = paste("Level", levels_shown[j]))
    abline(v = 0, col = "gray50", lty = 2)
    segments(z$before, seq_along(ids), z$after, seq_along(ids),
      col = adjustcolor(candidate_colours[ids], .5))
    points(z$before, seq_along(ids), pch = 16, col = candidate_colours[ids])
    points(z$after, seq_along(ids), pch = 1, col = candidate_colours[ids])
    axis(2, at = seq_along(ids), labels = if (j == 1L) labels[ids] else FALSE, las = 1,
      cex.axis = .75)
  }
  mtext(paste0(disease_names[[d]], ": the score's quantile up to the peak (filled, weeks ",
    min(phase_window$before), " to 0) and after it (open, weeks 1 to ",
    max(phase_window$after), "), mean over the ",
    pooled_label(length(settings$pooled_locations), settings), " and horizons"), outer = TRUE,
    font = 2, cex = 1.1)
  mtext(paste("A threshold that is positive up to the peak and negative after it: the base",
    "forecasts too low on the rise and too high in the decline, whatever the scale."),
    side = 1, outer = TRUE, line = 1, cex = .8)
  legend_origins()
}
# Pages 9-12: MultiQT on every score at the adaptive rate, and the running quantile of the
# scores of the same phase only (phase from the peaks in hindsight), overall and after peaks.
for (variant in c("MultiQT", "same phase (hindsight)")) for (phase in c("overall", "after")) {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), oma = c(4, 0, 3, 0), xpd = FALSE)
  for (d in diseases) for (design in c("delay 0", "delay 1")) {
    x <- pooled_running[pooled_running$variant == variant & pooled_running$disease == d &
      pooled_running$design == design & pooled_running$phase == phase, ]
    plot(x$loss_change, x$ce_change, pch = 16, col = candidate_colours[x$candidate],
      xlab = "Loss change versus the base (%)", ylab = "CE change versus the base (%)",
      main = paste0(disease_names[[d]], ", ", design), xlim = extendrange(x$loss_change, f = .3),
      ylim = extendrange(x$ce_change, f = .1))
    abline(h = 0, v = 0, col = "gray75", lty = 2)
    text(x$loss_change, x$ce_change, labels[x$candidate], pos = 4, cex = .6,
      col = candidate_colours[x$candidate])
  }
  mtext(paste0(if (variant == "MultiQT") "MultiQT on each score at the adaptive rate" else
    "Running quantile of the past scores of the same phase (phase in hindsight)", ", ",
    c(overall = "overall", after = paste("after peaks (weeks 1 to",
      settings$weeks_after_peak, ")"))[[phase]], ": change versus the base, the ",
    pooled_label(length(settings$pooled_locations), settings)), outer = TRUE, font = 2,
    cex = 1.1)
  mtext(if (variant == "MultiQT") paste("Thresholds start at 0 and learn at .03 times the",
    "rolling 90% quantile of the score's absolute values (last 20 learned rows). Evaluation weeks.")
    else paste("Each week the alpha-quantile of the past scores of rows in the same phase as",
      "the target week (rise or decline, from the peaks in hindsight). Evaluation weeks."),
    side = 1, outer = TRUE, line = 1, cex = .8)
  legend_origins()
}
invisible(dev.off())

# 5. SCORES OVER TIME: multiqt_score_stationarity/ -------------------------------
# For the shown candidates and three levels: the scores of the latest counts divided by one
# constant per case (the median |score| over its weeks and levels), so the panels share an
# axis. A quantile-stable score keeps its level's quantile, the ideal threshold of a
# tracker, flat over time. One PDF per disease and horizon: first the pooled locations by
# week relative to the peak (the level's quantile over the pooled locations and all peaks),
# then each location of settings$location_ids over time, with the peaks dotted and the
# first season shaded. Scales from reports use the feedback of delay 0.
mark_time <- function(case) {
  tuning <- case$target_date[case$period == "tuning"]
  usr <- par("usr")
  if (length(tuning)) rect(min(tuning) - 3.5, usr[3], max(tuning) + 3.5, usr[4],
    col = "gray93", border = NA)
  abline(v = as.numeric(case$peaks$peak_date), col = "gray40", lty = 3)
}
normalized <- function(S, alpha) {
  x <- S[, match(levels_shown, alpha), drop = FALSE]
  x / median(abs(S[is.finite(S)]))
}
score_weeks <- -(settings$weeks_before_peak - 1L):settings$weeks_after_peak
plots <- lapply(results, `[[`, "plot")
for (d in diseases) for (horizon in settings$horizons) {
  index <- which(cases$disease == d & cases$horizon == horizon)
  path <- file.path(figure_dir, "multiqt_score_stationarity", paste0(d, "_horizon", horizon,
    ".pdf"))
  cairo_pdf(path, width = 16, height = 11, family = "DejaVu Sans", onefile = TRUE)
  panel_layout <- function() {
    layout(matrix(1:15, 3, 5, byrow = TRUE))
    par(oma = c(3.5, 0, 3.5, 0), mar = c(4, 4.2, 3, 1), mgp = c(2.4, .7, 0), cex = .8)
  }
  legend_panel <- function() {
    plot.new()
    legend("center", legend = paste("Level", levels_shown), col = level_colours, lwd = 2,
      bty = "n", cex = 1.1)
  }
  title_line <- function(where) mtext(paste(disease_names[[d]], "| horizon", horizon, "|",
    where, "| scores of the base forecast against the latest counts"), side = 3,
    outer = TRUE, line = 1.5, font = 2, cex = 1.1)
  # Pooled page: the level's quantile of the normalized scores by week relative to each peak.
  panel_layout()
  for (id in shown) {
    weekly <- do.call(rbind, lapply(index[cases$location[index] %in%
      settings$pooled_locations], function(i) {
      case <- plots[[i]]$case
      S <- normalized(plots[[i]]$scores[[id]], case$alpha)
      keep <- case$relative_week %in% score_weeks & apply(is.finite(S), 1, all)
      data.frame(week = case$relative_week[keep], S[keep, , drop = FALSE])
    }))
    paths <- sapply(seq_along(levels_shown), function(j) vapply(score_weeks, function(w)
      quantile(weekly[weekly$week == w, j + 1L], levels_shown[j], na.rm = TRUE,
        names = FALSE), numeric(1)))
    if (!any(is.finite(paths))) {                # Own-target revisions, longest horizon.
      plot.new()
      title(main = labels[[id]])
      text(.5, .5, "Not defined at this horizon", col = "gray40")
      next
    }
    plot(NA, xlim = range(score_weeks), ylim = range(-1.5, 1.5, paths, na.rm = TRUE),
      xlab = "Week relative to the peak", ylab = "Score / median |score|", main = labels[[id]])
    abline(h = 0, col = "gray65", lty = 2)
    abline(v = .5, col = "gray40", lty = 3)
    matlines(score_weeks, paths, col = level_colours, lty = 1, lwd = 2)
  }
  legend_panel()
  title_line(paste("the", pooled_label(length(settings$pooled_locations), settings),
    "pooled, all peaks"))
  mtext(paste("Each line: the level's quantile of the normalized scores over the pooled",
    "locations and peaks in that week, the threshold a tracker at that level would need;",
    "flat = quantile-stable.\nScores divided by their case's median |score| (all weeks and",
    "levels). Dotted: between the peak week and the week after."), side = 1, outer = TRUE,
    line = 2, cex = .75)
  # One page per plotted location: the normalized scores over time.
  for (i in index[match(settings$location_ids, cases$location[index])]) {
    case <- plots[[i]]$case
    dates <- case$target_date
    ticks <- seq(as.Date(format(min(dates), "%Y-%m-01")), max(dates), by = "6 months")
    panel_layout()
    for (id in shown) {
      values <- normalized(plots[[i]]$scores[[id]], case$alpha)
      if (!any(is.finite(values))) {             # Own-target revisions, longest horizon.
        plot.new()
        title(main = labels[[id]])
        text(.5, .5, "Not defined at this horizon", col = "gray40")
        next
      }
      limit <- max(1.5, quantile(abs(values), .98, na.rm = TRUE))
      broken <- break_gaps(dates, pmin(pmax(values, -limit), limit))
      plot(NA, xlim = range(dates), ylim = c(-limit, limit), xaxt = "n", xaxs = "i",
        xlab = "Target date", ylab = "Score / median |score|", main = labels[[id]])
      mark_time(case)
      abline(h = 0, col = "gray65", lty = 2)
      matlines(broken$dates, broken[[2]], col = level_colours, lty = 1, lwd = 1.2)
      axis.Date(1, at = ticks, format = "%B\n%Y", cex.axis = .75)
    }
    legend_panel()
    title_line(paste0(case$location_name, " (", case$location, ")"))
    mtext(paste("Score of the latest count relative to the base quantile of each level,",
      "divided by the case's median |score|; beyond the axis: at its edge.\nDotted lines:",
      "peaks (local maxima of the latest counts). Gray: first season. Scales from reports:",
      "feedback delay 0."), side = 1, outer = TRUE, line = 2, cex = .75)
  }
  invisible(dev.off())
}
