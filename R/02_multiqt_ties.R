# Run R/01_data_overview.R once to prepare output/data/peak_aligned.rds.

# 1. SETTINGS AND PREPARED DATA -------------------------------------
disease_name <- "covid"                    # "covid" or "flu"
location_id <- "US"                       # One location code; ignored if combined=TRUE.
horizon_id <- 0
combined <- FALSE                         # TRUE averages independently fitted locations.
make_phase_overview <- TRUE               # Also compare both diseases, h-1:3, locations and US.
use_log_scale <- TRUE                     # TRUE: log1p(counts); FALSE: raw counts for MultiQT.
# Script 03 inherits this switch when sourced afterward in the same R session.
scale_method <- "max"                     # Scoring/plots only: "median", "max", or "none".
expert_multipliers <- c(0, .01, .015, .02, .03, .05, .075, .1, .15, .2, .3)
expert_warmup_rounds <- 20L                 # Shared with script 03's ensemble.
expert_scale_floor <- 1e-6                 # Guard on the error scale, not each rate.
multiplier_min <- min(expert_multipliers[expert_multipliers > 0])
multiplier_max <- max(expert_multipliers)
n_positive_multipliers <- 81L
display_multipliers <- c(0, .01, .03, .075, .15, .3)
# No forecasts were issued 2025-10-25 to 2025-11-15, so a longer pre-peak window
# drops locations with December/early-January peaks. Only evaluation is windowed.
weeks_before_peak <- 3L                    # Before window also scores the peak week (0).
weeks_after_peak <- 8L
combined_weeks <- seq.int(-weeks_before_peak, weeks_after_peak)
positive_grid <- 10^seq(log10(multiplier_min), log10(multiplier_max), length.out = n_positive_multipliers)
rate_multipliers <- sort(unique(c(expert_multipliers, signif(positive_grid, 12), display_multipliers)))
tracker_scale <- if (use_log_scale) "log1p_counts" else "counts"
# The dense plot grid includes every script 03 multiplier. For each location, eta=m*S0:
# S0 is the frozen Q90 absolute base error in tracker units from the first 20 released rounds.
# Script 03 uses the same switch and error scale for its matching experts.
# The multiplier is shared across locations; actual eta is recorded separately.
peak_data <- readRDS("output/data/peak_aligned.rds")
if (make_phase_overview &&
    !all(-1:3 %in% vapply(peak_data$cases, function(x) x$horizon, integer(1))))
  stop("Phase overview needs prepared horizons -1:3. Rerun Rscript R/01_data_overview.R first.")
selected_cases <- Filter(function(x) x$disease == disease_name && x$horizon == horizon_id &&
  if (combined) x$location != "US" else x$location == location_id, peak_data$cases)
# Use the same complete calendar window for every location in a combined plot.
# Save all support flags, including excluded places, so the cohort is explicit.
get_location_support <- function(cases) do.call(rbind, lapply(cases, function(x) {
  i <- match(combined_weeks, x$relative_week)
  ok <- !is.na(i)
  ok[ok] <- x$within_winter[i[ok]] & is.finite(x$truth_latest[i[ok]]) &
    rowSums(is.finite(x$base[i[ok], , drop = FALSE])) == length(x$alpha)
  data.frame(location = x$location, location_name = x$location_name, peak_date = x$peak_date,
    boundary_peak = x$peak_boundary,
    n_observed = sum(ok), n_expected = length(combined_weeks), complete = all(ok))
}))
location_support <- get_location_support(selected_cases)
# A maximum on the winter edge is likely an earlier wave's tail, not a winter peak.
if (combined) selected_cases <- selected_cases[location_support$complete & !location_support$boundary_peak]

result_dir <- if (combined) file.path("output/results/combined", paste0(disease_name, "_h", horizon_id), "02") else "output/results"
figure_dir <- if (combined) file.path("output/pdf/combined", paste0(disease_name, "_h", horizon_id), "02") else "output/pdf"
overview_file <- "output/pdf/multiqt_phase_rate_tradeoffs.pdf"
overview_result_file <- "output/results/multiqt_phase_rate_tradeoffs.csv"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(location_support, file.path(result_dir, "multiqt_location_support.csv"), row.names = FALSE)

# 2. RESTORE EACH LOCATION'S INPUTS ----------------------------------
# 01 supplies full issue histories and first-use feedback in raw counts.
# Track log1p(counts) or raw counts; historical scaling is for scoring/plots only.
# Keep every location's offsets and past forecasts separate.
run_location <- function(case, case_name, combined) {
  scale_history <- case$scale_history[is.finite(case$scale_history)]
  count_scale <- 1
  scale_used <- scale_method
  if (scale_method != "none" && length(scale_history)) {
    count_scale <- if (scale_method == "median") median(scale_history) else max(scale_history)
    if (count_scale <= 0) {
      count_scale <- max(scale_history)
      scale_used <- "max"
    }
  }
  if (!is.finite(count_scale) || count_scale <= 0 || !length(scale_history)) {
    count_scale <- 1
    scale_used <- "none"
  }
  alpha <- case$alpha
  stopifnot(all(is.finite(case$base)), all(case$base >= 0),
            all(is.na(case$feedback_observed) |
              (is.finite(case$feedback_observed) & case$feedback_observed >= 0)))
  base_tracker <- if (use_log_scale) log1p(case$base) else case$base
  feedback_tracker <- if (use_log_scale) log1p(case$feedback_observed) else case$feedback_observed
  # These count-scale arrays are for scoring, saved forecasts and display only.
  base <- case$base / count_scale
  feedback_observed <- case$feedback_observed / count_scale
  truth_latest <- feedback_observed#case$truth_latest / count_scale # Scoring target; updates use feedback only.
  feedback_by_week <- case$feedback_by_week
  stopifnot(all(is.finite(feedback_tracker[unlist(feedback_by_week)])))
  feedback_date <- case$feedback_date
  example <- data.frame(reference_date = case$reference_date, target_end_date = case$target_date,
    truth_finalized = case$truth_finalized / count_scale, relative_week = case$relative_week)
  n_weeks <- nrow(base)
  n_levels <- length(alpha)
  # Same warmup as script 03: first 20 released rows, all quantile levels, no replay.
  # This freezes the notebook's error-scaled rate rather than updating it each week.
  # https://delphi-forecasting-reports.netlify.app/calibration_experiments/e01_eta_settings
  initialization_week <- match(TRUE, cumsum(lengths(feedback_by_week)) >= expert_warmup_rounds)
  initialized <- !is.na(initialization_week)
  warmup_weeks <- if (initialized) initialization_week else n_weeks
  # Consume the whole freeze-date batch, but use only the first 20 rows for its scale.
  released <- Reduce(c, feedback_by_week[seq_len(warmup_weeks)], init = integer())
  warmup_rows <- if (initialized) released[seq_len(expert_warmup_rounds)] else released
  rate_scale <- raw_rate_scale <- NA_real_
  if (initialized) {
    residuals <- sweep(base_tracker[warmup_rows, , drop = FALSE], 1,
      feedback_tracker[warmup_rows], "-")
    raw_rate_scale <- as.numeric(quantile(abs(residuals), .9, names = FALSE))
    rate_scale <- max(expert_scale_floor, raw_rate_scale)
  }
  # If fewer than 20 reports arrive, there are no post-initialization weeks to score.
  learning_rates <- rate_multipliers * rate_scale
  learning_rates[rate_multipliers == 0] <- 0
  rate_initialization <- list(initialized = initialized,
    week = if (initialized) warmup_weeks else NA_integer_,
    date = if (initialized) case$reference_date[warmup_weeks] else as.Date(NA),
    rows = warmup_rows, consumed_rows = released,
    raw_scale = raw_rate_scale, scale = rate_scale,
    required_rounds = expert_warmup_rounds, scale_floor = expert_scale_floor)
  rate_grid_row <- data.frame(case_id = case_name, disease = case$disease,
    location = case$location, horizon = case$horizon, tracker_scale, multiplier = rate_multipliers,
    eta = learning_rates, ensemble_expert = rate_multipliers %in% expert_multipliers,
    initialized, initialization_date = rate_initialization$date,
    initialization_rounds = length(warmup_rows), consumed_rounds = length(released),
    raw_error_scale = raw_rate_scale, error_scale = rate_scale, scale_floor = expert_scale_floor)
  # Evaluate only forecasts issued strictly after initialization, including m=0.
  # If initialization never completes, there are no eligible scoring weeks.
  post_warmup <- seq_len(n_weeks) > warmup_weeks
  score_weeks <- which(post_warmup & is.finite(truth_latest) &
    if (combined) case$relative_week %in% combined_weeks & case$within_winter else TRUE)
  if (!length(score_weeks)) stop("No post-initialization weeks to score for ", case_name)
  feedback_log <- data.frame(location = case$location, reference_date = example$reference_date,
    target_end_date = example$target_end_date, relative_week = case$relative_week,
    feedback_date, feedback_observed, truth_finalized = example$truth_finalized, truth_latest)
  scale_row <- data.frame(disease = case$disease, location = case$location,
    horizon = case$horizon, tracker_scale, scale_method, scale_used,
    count_scale, scale_as_of = case$scale_as_of,
    reference_weeks = length(scale_history), reference_start = case$scale_reference_start,
    reference_end = case$scale_reference_end)
  # 3. RUN MULTIQT AND SCORE EACH ALPHA --------------------------------
  forecasts_by_multiplier <- covered_by_multiplier <- loss_by_multiplier <- list()
  coverage <- quantile_loss_by_level <- matrix(NA_real_, n_levels, length(rate_multipliers))
  for (i in seq_along(rate_multipliers)) {
    eta <- learning_rates[i]
    offsets <- numeric(n_levels)               # Restart from zero for every rate.
    adjusted_tracker <- base_tracker
    if (initialized && rate_multipliers[i] != 0) for (week in which(post_warmup)) {
      for (past_week in feedback_by_week[[week]]) {
        covered <- feedback_tracker[past_week] <= adjusted_tracker[past_week, ]
        offsets <- offsets + eta * (alpha - covered)
        # Compare with the HISTORICAL ISSUED forecast in tracker units.
        # At alpha=.9, eta=.01, a miss adds .009; coverage subtracts .001.
      }
      # Project only the issued forecasts, NOT the offsets ("lazy" MultiQT).
      # isoreg restores quantile order: e.g., (12,10) becomes (11,11).
      adjusted_tracker[week, ] <- isoreg(base_tracker[week, ] + offsets)$yf
    }
    # Back-transform log mode before scoring. Pinball loss in these units is count loss
    # divided by count_scale exactly once (historical maximum by default).
    adjusted_counts <- if (use_log_scale) expm1(adjusted_tracker) else adjusted_tracker
    unchanged <- adjusted_tracker == base_tracker
    adjusted_counts[unchanged] <- case$base[unchanged] # Preserve exact base ties.
    adjusted <- adjusted_counts / count_scale
    key <- as.character(rate_multipliers[i])
    forecasts_by_multiplier[[key]] <- adjusted
    # Score against the latest vintage; the offsets above saw only first-use feedback.
    covered <- sweep(adjusted_counts, 1, case$truth_latest, ">=")
    error <- sweep(-adjusted, 1, truth_latest, "+")
    # Pinball loss penalizes underprediction by alpha, overprediction by 1-alpha.
    loss <- sweep(pmax(error, 0), 2, alpha, "*") +
      sweep(pmax(-error, 0), 2, 1 - alpha, "*")
    # Weekly score plots must also omit warmup and the initialization week.
    covered[!post_warmup, ] <- NA
    loss[!post_warmup, ] <- NA_real_
    coverage[, i] <- colMeans(covered[score_weeks, , drop = FALSE], na.rm = TRUE)
    quantile_loss_by_level[, i] <- colMeans(loss[score_weeks, , drop = FALSE], na.rm = TRUE)
    # Reuse these weekly diagnostics in the timeline and peak comparisons.
    covered_by_multiplier[[key]] <- covered
    loss_by_multiplier[[key]] <- loss
  }
  # Calibration is checked per alpha. Positive signed gap = coverage above alpha.
  signed_coverage_error <- sweep(coverage, 1, alpha, "-") # coverage - alpha
  coverage_error_by_level <- abs(signed_coverage_error)
  results <- data.frame(multiplier = rate_multipliers,
                        coverage_error = colMeans(coverage_error_by_level),
                        quantile_loss = colMeans(quantile_loss_by_level))
  # expand.grid varies alpha fastest, matching the matrices' column order.
  scores_by_level <- data.frame(
    expand.grid(alpha = alpha, multiplier = rate_multipliers)[c("multiplier", "alpha")],
    coverage = as.vector(coverage), signed_coverage_error = as.vector(signed_coverage_error),
    coverage_error = as.vector(coverage_error_by_level),
    quantile_loss = as.vector(quantile_loss_by_level), n_weeks = length(score_weeks))
  # Retain the full issued forecasts; only scoring/plotting uses the peak window.
  list(case = case, count_scale = count_scale, scale_used = scale_used,
    tracker_scale = tracker_scale, prediction_scale = "counts / count_scale",
    rate_multipliers = rate_multipliers, learning_rates = learning_rates,
    rate_initialization = rate_initialization, score_weeks = score_weeks,
    alpha = alpha, base = base, example = example, feedback_observed = feedback_observed,
    truth_latest = truth_latest, forecasts_by_multiplier = forecasts_by_multiplier, covered_by_multiplier = covered_by_multiplier,
    loss_by_multiplier = loss_by_multiplier, coverage = coverage,
    quantile_loss_by_level = quantile_loss_by_level, coverage_error_by_level = coverage_error_by_level,
    results = results, scores_by_level = scores_by_level, feedback_log = feedback_log,
    scale_row = scale_row, rate_grid_row = rate_grid_row)
}
location_runs <- setNames(lapply(names(selected_cases), function(case_name)
  run_location(selected_cases[[case_name]], case_name, combined)), names(selected_cases))
scale_rows <- lapply(location_runs, `[[`, "scale_row")
rate_grid_rows <- lapply(location_runs, `[[`, "rate_grid_row")
# Restore the selected run's plotting inputs; overview runs below stay in their own scope.
last_run <- location_runs[[length(location_runs)]]
case <- last_run$case
alpha <- last_run$alpha
n_levels <- length(alpha)
example <- last_run$example
count_scale <- last_run$count_scale
scale_used <- last_run$scale_used
coverage <- last_run$coverage
coverage_error_by_level <- last_run$coverage_error_by_level
quantile_loss_by_level <- last_run$quantile_loss_by_level
results <- last_run$results
scores_by_level <- last_run$scores_by_level
n_locations <- length(location_runs)
mean_location_field <- function(runs, field)
  Reduce(`+`, lapply(runs, `[[`, field)) / length(runs)
# Absolute gaps are formed INSIDE each location above, before this extra average.
if (combined) {
  coverage <- mean_location_field(location_runs, "coverage")
  coverage_error_by_level <- mean_location_field(location_runs, "coverage_error_by_level")
  quantile_loss_by_level <- mean_location_field(location_runs, "quantile_loss_by_level")
  results$coverage_error <- colMeans(coverage_error_by_level)
  results$quantile_loss <- colMeans(quantile_loss_by_level)
  scores_by_level$coverage <- as.vector(coverage)
  scores_by_level$signed_coverage_error <- as.vector(sweep(coverage, 1, alpha, "-"))
  scores_by_level$coverage_error <- as.vector(coverage_error_by_level)
  scores_by_level$quantile_loss <- as.vector(quantile_loss_by_level)
}
# Pooled CE: average signed gaps over locations first, then |gap| per level, then
# over levels. Location-specific errors cancel; it is always <= coverage_error.
results$pooled_coverage_error <- colMeans(abs(sweep(coverage, 1, alpha, "-")))
results$n_locations <- n_locations
scores_by_level$n_locations <- n_locations
feedback_log <- do.call(rbind, lapply(location_runs, `[[`, "feedback_log"))
scale_table <- do.call(rbind, scale_rows)
write.csv(results, file.path(result_dir, "multiqt_learning_rates.csv"), row.names = FALSE)
write.csv(scores_by_level, file.path(result_dir, "multiqt_scores_by_level.csv"), row.names = FALSE)
write.csv(feedback_log, file.path(result_dir, "multiqt_feedback.csv"), row.names = FALSE)
write.csv(scale_table, file.path(result_dir, "multiqt_scale.csv"), row.names = FALSE)
write.csv(do.call(rbind, rate_grid_rows), file.path(result_dir, "multiqt_rate_grid.csv"), row.names = FALSE)
local_scores <- do.call(rbind, lapply(location_runs, function(x)
  data.frame(location = x$case$location,
    eta = rep(x$learning_rates, each = length(x$alpha)), x$scores_by_level)))
write.csv(local_scores, file.path(result_dir, "multiqt_scores_by_location.csv"), row.names = FALSE)
scale_label <- paste(if (use_log_scale) "Frozen log1p tracker; eta = m * S0 |" else
  "Frozen count tracker; eta = m * S0 |", if (combined)
  paste(n_locations, "locations | scores/plots:", scale_method,
    "scale | weeks", min(combined_weeks), "to", max(combined_weeks)) else
  if (scale_used == "none") "Scores/plots: raw counts" else
    paste0("Scores/plots: counts / historical ", scale_used, " (", format(count_scale, trim = TRUE), ")"))
score_units <- if (all(scale_table$scale_used == "none")) "counts" else
  if (all(scale_table$scale_used == "max")) "counts / max" else
    if (all(scale_table$scale_used == "median")) "counts / median" else "counts / local scale"
count_axis_label <- paste0("Hospitalizations (", score_units, ")")
loss_axis_label <- paste0("Quantile loss (", score_units, ")")
mean_loss_axis_label <- paste("Mean", tolower(loss_axis_label))
time_axis_label <- if (combined) "Weeks from winter peak (peak = 0)" else "Target date"
multiplier_axis_label <- "Rate multiplier m (log axis)"
location_label <- if (combined) paste(n_locations, "locations combined (US excluded)") else location_id
display_indices <- match(display_multipliers, rate_multipliers)
print(results[display_indices, ])

# 4. COMPARE ERRORS: --------------------------------
# This is exploratory, not held-out rate selection. Full-period coverage can
# hide seasonal errors. Averaging levels also hides which quantiles struggle;
# every supplied level has equal weight, despite their unequal alpha spacing.
positive <- which(rate_multipliers > 0)
# log(0) is undefined: zero appears as a horizontal baseline or tradeoff cross.
baseline <- match(0, rate_multipliers)
log_multipliers <- log10(rate_multipliers[positive])
# Label both endpoints; omit nearby decade labels so the maximum stays readable.
multiplier_ticks <- 10^seq(floor(min(log_multipliers)), ceiling(max(log_multipliers)))
multiplier_ticks <- multiplier_ticks[log10(multiplier_ticks) > min(log_multipliers) + 0.5 &
                         log10(multiplier_ticks) < max(log_multipliers) - 0.5]
multiplier_ticks <- sort(unique(c(range(rate_multipliers[positive]), multiplier_ticks)))
multiplier_tick_labels <- trimws(formatC(multiplier_ticks, format = "fg", digits = 6))
colour_index <- 1 + round(255 * (log_multipliers - min(log_multipliers)) / diff(range(log_multipliers)))
multiplier_colours <- rep("orange", length(rate_multipliers))
multiplier_colours[positive] <- hcl.colors(320, "Viridis")[colour_index]

# One small helper avoids duplicating the same three panels for every alpha.
# It reads the rate grid/colours above; loss and gap are vectors over that grid.
plot_error_row <- function(loss, gap, label,
                           loss_range = range(loss), gap_range = range(gap), pooled_gap = NULL) {
  loss_label <- paste(strwrap(mean_loss_axis_label, width = 20), collapse = "\n")
  plot(rate_multipliers[positive], loss[positive], log = "x", type = "b", pch = 16,
       col = "steelblue4", cex = 0.45, ylim = loss_range, xaxt = "n",
       xlab = multiplier_axis_label, ylab = loss_label,
       main = paste(label, "| Quantile loss"))
  axis(1, at = multiplier_ticks, labels = multiplier_tick_labels)
  abline(h = loss[baseline], lty = 2, col = "orange")
  plot(rate_multipliers[positive], gap[positive], log = "x", type = "b", pch = 16,
       col = "steelblue4", cex = 0.45, xaxt = "n",  # Headroom for the pooled legend.
       ylim = range(0, gap_range) * if (is.null(pooled_gap)) 1 else c(1, 1.3),
       xlab = multiplier_axis_label, ylab = "Calibration error",
       main = paste(label, "| Calibration error"))
  axis(1, at = multiplier_ticks, labels = multiplier_tick_labels)
  abline(h = gap[baseline], lty = 2, col = "gray40")
  if (!is.null(pooled_gap)) {
    lines(rate_multipliers[positive], pooled_gap[positive], type = "b", pch = 1,
          col = "firebrick", cex = 0.45)
    abline(h = pooled_gap[baseline], lty = 2, col = "firebrick")
    legend("topright", c("Per-location CE", "Pooled CE"), col = c("steelblue4", "firebrick"),
           pch = c(16, 1), lty = 1, cex = 0.7, bty = "n")
  }
  plot(gap[positive], loss[positive], type = "l", col = "gray65",
       xlim = extendrange(gap_range), ylim = extendrange(loss_range),
       xlab = "Calibration error", ylab = loss_label,
       main = paste(label, "| Tradeoff"))
  points(gap[positive], loss[positive], pch = 16, cex = 0.6, col = multiplier_colours[positive])
  points(gap[baseline], loss[baseline], pch = 4, lwd = 2, col = "orange")
}
# pdf() starts writing vector graphics; dev.off() closes the file.
# Six alphas per page; PDF automatically starts a new page when panels fill.
cairo_pdf(file.path(figure_dir, "multiqt_by_level.pdf"), width = 14, height = 15, family = "DejaVu Sans")
par(mfrow = c(6, 3), mar = c(3.5, 5, 2.5, 1), oma = c(4, 0, 0, 0),
    mgp = c(2.2, 0.7, 0), cex = 0.85)
for (j in seq_along(alpha)) {
  plot_error_row(quantile_loss_by_level[j, ], coverage_error_by_level[j, ],
                 paste("alpha =", alpha[j]))
  if (j %% 6 == 0 || j == n_levels) {
    mtext(paste("m =", display_multipliers), side = 1, outer = TRUE, line = 1,
          at = seq(0.12, 0.88, length.out = length(display_multipliers)),
          col = multiplier_colours[display_indices], cex = 0.8)
    mtext(paste(scale_label, "| Dashed line / orange cross: m = 0. Lower is better."),
          side = 1, outer = TRUE, line = 2.5, cex = 0.8)
  }
}
dev.off()

cairo_pdf(file.path(figure_dir, "multiqt_error_tradeoff.pdf"), width = 12, height = 4.5, family = "DejaVu Sans")
par(mfrow = c(1, 3), mar = c(4, 4.5, 3, 1), oma = c(3, 0, 0, 0), mgp = c(2.3, 0.7, 0))
plot_error_row(results$quantile_loss, results$coverage_error, "Mean over quantile levels",
               pooled_gap = if (combined) results$pooled_coverage_error)
mtext(paste("m =", display_multipliers), side = 1, outer = TRUE, line = 0.5,
      at = seq(0.12, 0.88, length.out = length(display_multipliers)),
      col = multiplier_colours[display_indices], cex = 0.8)
mtext(paste(scale_label, "| Dashed line / orange cross: m = 0. Lower is better."),
      side = 1, outer = TRUE, line = 1.8, cex = 0.8)
dev.off()

cairo_pdf(file.path(figure_dir, "multiqt_learning_rates.pdf"), width = 7, height = 6, family = "DejaVu Sans")
matplot(alpha, coverage[, display_indices, drop = FALSE], type = "b", pch = 16, lty = 1,
        col = multiplier_colours[display_indices], xlim = c(0, 1), ylim = c(0, 1),
        xlab = "Quantile level (alpha)", ylab = "Coverage",
        main = paste(strwrap(scale_label, width = 48), collapse = "\n"))
abline(0, 1, lty = 2, col = "gray40")   # Diagonal = agreement with nominal coverage.
legend("topleft", paste("m =", display_multipliers), col = multiplier_colours[display_indices],
       pch = 16, lty = 1, bty = "n")
dev.off()

# 5. INSPECT THE TIME SERIES AROUND THE WINTER PEAK -------------------
plot_columns <- match(c(0.15, 0.5, 0.75), alpha)
plot_alpha <- alpha[plot_columns]
alpha_colours <- c("#0072B2", "#009E73", "#D55E00")
season_dates <- if (combined) combined_weeks else example$target_end_date
peak_dates <- if (combined) 0 else case$peak_date
date_ticks <- if (combined) pretty(season_dates) else
  seq(as.Date(format(min(season_dates), "%Y-%m-01")), max(season_dates), by = "3 months")
# Use the same display code for one location and an aligned location average.
season_mean <- function(extract) Reduce(`+`, lapply(location_runs, function(x) {
  i <- if (combined) match(combined_weeks, x$case$relative_week) else seq_len(nrow(x$base))
  extract(x, i)
})) / n_locations
season_observed <- season_mean(function(x, i) x$truth_latest[i])
base_season <- season_mean(function(x, i) x$base[i, plot_columns, drop = FALSE])
season_forecasts_by_multiplier <- season_coverage_by_multiplier <- season_loss_by_multiplier <- list()
for (multiplier in display_multipliers) {
  key <- as.character(multiplier)
  season_forecasts_by_multiplier[[key]] <- season_mean(function(x, i) x$forecasts_by_multiplier[[key]][i, plot_columns, drop = FALSE])
  covered <- season_mean(function(x, i) x$covered_by_multiplier[[key]][i, plot_columns, drop = FALSE])
  season_coverage_by_multiplier[[key]] <- sweep(covered, 2, plot_alpha, "-")
  season_loss_by_multiplier[[key]] <- season_mean(function(x, i) x$loss_by_multiplier[[key]][i, plot_columns, drop = FALSE])
}
# Averaged forecasts/outcomes are display profiles only. Scores were computed locally.
plot_limits <- range(0, season_observed, base_season, unlist(season_forecasts_by_multiplier), na.rm = TRUE)
loss_limits <- range(0, unlist(season_loss_by_multiplier), na.rm = TRUE)
season_xlim <- range(season_dates)
# Both weekly diagnostics use the same quantile styles and peak marker.
plot_season_diagnostic <- function(values, limits, y_label, title,
                                   show_time = FALSE, zero_line = FALSE) {
  par(mar = c(if (show_time) 3.5 else 1, 4.5, 2.5, 1))
  matplot(season_dates, values, type = "b",
          col = alpha_colours, lty = 1, pch = c(1, 16, 17), cex = 0.5,
          ylim = limits, xlim = season_xlim, xaxt = "n", xaxs = "i",
          xlab = if (show_time) time_axis_label else "", ylab = y_label, main = title)
  if (zero_line) abline(h = 0, col = "gray65", lty = 2)
  if (show_time) {
    if (combined) axis(1, at = date_ticks) else axis.Date(1, at = date_ticks, format = "%b\n%Y", cex.axis = 0.8)
  }
  abline(v = as.numeric(peak_dates), col = "gray40", lty = 3)
}
cairo_pdf(file.path(figure_dir, "multiqt_season_forecasts.pdf"), width = 12, height = 11, family = "DejaVu Sans")
layout(matrix(c(1, 4, 2, 5, 3, 6), nrow = 3, byrow = TRUE), heights = c(2, 1, 1))
par(oma = c(5, 0, 3, 0), mgp = c(2.5, 0.7, 0), cex = 0.9)
for (page_start in seq(1, length(display_indices), by = 2)) {
  page_indices <- display_indices[page_start:min(page_start + 1, length(display_indices))]
  for (i in page_indices) {
    multiplier <- rate_multipliers[i]
    key <- as.character(multiplier)
    q <- season_forecasts_by_multiplier[[key]]
    par(mar = c(1, 4.5, 2.5, 1))
    plot(season_dates, season_observed, type = "n", ylim = plot_limits,
         xlim = season_xlim, xaxt = "n", xaxs = "i", xlab = "",
         ylab = count_axis_label, main = paste("Rate multiplier m =", multiplier))
    polygon(c(season_dates, rev(season_dates)), c(q[, 1], rev(q[, 3])),
            col = adjustcolor(multiplier_colours[i], alpha.f = 0.2), border = NA)
    # Thick median lines; shading spans the selected lower/upper quantiles.
    matlines(season_dates, base_season, col = "gray50", lty = 2, lwd = c(1, 2, 1))
    matlines(season_dates, q, col = multiplier_colours[i], lty = 1, lwd = c(1, 2, 1))
    lines(season_dates, season_observed, type = "o", pch = 16, cex = 0.5)
    abline(v = as.numeric(peak_dates), lty = 3, col = "gray40")
    legend("topright", c("Latest reported counts", "Original", "Fixed MultiQT", "Winter peak (retrospective)"),
           col = c("black", "gray50", multiplier_colours[i], "gray40"),
           lty = c(1, 2, 1, 3), pch = c(16, NA, NA, NA), cex = 0.7, bty = "n")

    plot_season_diagnostic(season_coverage_by_multiplier[[key]], c(-1, 1),
      "Coverage - alpha", "Weekly coverage residual", zero_line = TRUE)
    plot_season_diagnostic(season_loss_by_multiplier[[key]], loss_limits,
      loss_axis_label, "Weekly quantile loss", show_time = TRUE)
  }
  if (length(page_indices) == 1L) {
    for (panel in 1:3) plot.new()
  }
  mtext(paste(disease_name, "| location", location_label, "| horizon", horizon_id,
              "|", min(season_dates), "to", max(season_dates)),
        side = 3, outer = TRUE, line = 1, cex = 1)
  mtext(paste(scale_label, "| Vertical line: largest count in the selected winter."),
        side = 3, outer = TRUE, line = -0.3, cex = 0.8)
  mtext(paste("alpha =", plot_alpha), side = 1, outer = TRUE,
        line = 0.5, at = c(0.2, 0.5, 0.8), col = alpha_colours, cex = 0.85)
  mtext(if (combined) "Profiles and coverage residuals are location averages; calibration error is calculated locally first." else "Coverage residual = 1(latest reported count <= quantile) - alpha",
        side = 1, outer = TRUE, line = 2, cex = 0.75)
}
dev.off()
# Compare loss/coverage with these curves: do larger multipliers leave forecasts too high,
# too low, or unstable after the peak? This motivates future adaptive rates.

# 6. COMPARE ERRORS BEFORE AND AFTER THE WINTER PEAK ------------------
# Each location has one winter peak. Score its two windows directly.
summarize_peak_windows <- function(location_runs, combined) {
  alpha <- location_runs[[1L]]$alpha
  case <- location_runs[[1L]]$case
  peak_rows <- list()
  for (window in c("before", "after")) {
    local <- lapply(location_runs, function(x) {
      relative <- x$case$relative_week
      weeks <- which(x$case$within_winter & is.finite(x$truth_latest) &
        if (window == "before") relative >= -weeks_before_peak & relative <= 0 else
          relative > 0 & relative <= weeks_after_peak)
      weeks <- intersect(weeks, x$score_weeks)
      gaps <- losses <- rep(NA_real_, length(rate_multipliers))
      signed <- matrix(NA_real_, length(rate_multipliers), length(alpha))  # multiplier x level
      if (length(weeks)) {
        for (i in seq_along(rate_multipliers)) {
          key <- as.character(rate_multipliers[i])
          signed[i, ] <- colMeans(x$covered_by_multiplier[[key]][weeks, , drop = FALSE]) - alpha
          gaps[i] <- mean(abs(signed[i, ]))
          losses[i] <- mean(x$loss_by_multiplier[[key]][weeks, , drop = FALSE])
        }
      }
      list(gap = gaps, loss = losses, signed = signed, weeks = weeks, dates = x$case$target_date[weeks])
    })
    first <- local[[1L]]
    bounds <- if (combined) {
      if (window == "before") c(-weeks_before_peak, 0) else c(1, weeks_after_peak)
    } else if (length(first$weeks)) range(first$dates) else rep(as.Date(NA), 2L)
    peak_rows[[window]] <- data.frame(peak_number = 1L, window,
      peak_date = if (combined) 0 else case$peak_date,
      target_start = bounds[1], target_end = bounds[2], n_weeks = length(first$weeks),
      n_locations = sum(vapply(local, function(x) length(x$weeks) > 0L, logical(1))), multiplier = rate_multipliers,
      coverage_error = mean_location_field(local, "gap"),
      pooled_coverage_error = rowMeans(abs(mean_location_field(local, "signed"))),
      quantile_loss = mean_location_field(local, "loss"))
  }
  do.call(rbind, peak_rows)
}
peak_results <- summarize_peak_windows(location_runs, combined)
write.csv(peak_results, file.path(result_dir, "multiqt_error_tradeoff_all_peaks.csv"), row.names = FALSE)

finite_loss <- peak_results$quantile_loss[is.finite(peak_results$quantile_loss)]
finite_gap <- peak_results$coverage_error[is.finite(peak_results$coverage_error)]
shared_loss_range <- if (length(finite_loss)) range(finite_loss) else c(0, 1)
shared_gap_range <- if (length(finite_gap)) range(finite_gap) else c(0, 1)
window_labels <- c(before = "Up to winter peak", after = "After winter peak")
cairo_pdf(file.path(figure_dir, "multiqt_error_tradeoff_all_peaks.pdf"), width = 12, height = 8.5, family = "DejaVu Sans")
layout(rbind(c(1, 1, 1), c(2, 3, 4), c(5, 5, 5), c(6, 7, 8)), heights = c(.14, 1, .14, 1))
par(oma = c(5, 0, 3, 0), mgp = c(2.3, .7, 0), cex = .9)
for (window in names(window_labels)) {
  row <- peak_results[peak_results$window == window, ]
  label <- window_labels[[window]]
  par(mar = c(0, 0, 0, 0))
  plot.new()
  support_note <- if (row$n_weeks[1] > 0L) paste0(
    if (combined) "Relative week: " else "Target date: ", row$target_start[1], " to ", row$target_end[1],
    " | ", row$n_weeks[1], " observed weeks/location | ", row$n_locations[1], " location(s)") else
    "No observations in this winter window"
  text(.5, .5, support_note, cex = 1)
  par(mar = c(4, 4.5, 2.5, 1), cex.main = .95)
  if (row$n_weeks[1] > 0L) {
    plot_error_row(row$quantile_loss, row$coverage_error, label,
      loss_range = shared_loss_range, gap_range = shared_gap_range,
      pooled_gap = if (combined) row$pooled_coverage_error)
  } else {
    for (panel in 1:3) {
      plot.new()
      text(.5, .5, paste(label, "- no observed weeks"))
    }
  }
}
mtext(paste(disease_name, "| location", location_label, "| horizon", horizon_id,
            "| Winter peak:", peak_results$peak_date[1]),
      side = 3, outer = TRUE, line = 1, cex = 1)
mtext(paste("m =", display_multipliers), side = 1, outer = TRUE, line = .5,
      at = seq(.12, .88, length.out = length(display_multipliers)),
      col = multiplier_colours[display_indices], cex = .8)
mtext(paste(scale_label, "| All quantile levels; shared axes. Dashed lines / orange cross: m = 0."),
      side = 1, outer = TRUE, line = 1.8, cex = .8)
mtext(paste(weeks_before_peak, "calendar weeks before plus the peak week;", weeks_after_peak,
            "weeks after. Short-window calibration is descriptive."),
      side = 1, outer = TRUE, line = 3.1, cex = .8)
dev.off()

# 7. OVERVIEW: DISEASE, HORIZON AND AGGREGATION -----------------------
# Refit the remaining cases with the SAME tracker, scaling and multiplier settings.
# This uses prepared inputs directly, so it never depends on earlier exported results.
# Loss compares each phase with its own m=0 baseline, averaging locations before
# forming percentages. Calibration is the mean per-location absolute coverage gap,
# averaged over quantile levels. Both metrics use the same weeks, including for US.
build_phase_overview <- function() {
  overview_rows <- list()
  for (disease in c("covid", "flu")) for (horizon in -1:3) {
    cases <- Filter(function(x) x$disease == disease && x$horizon == horizon, peak_data$cases)
    if (!length(cases)) next
    support <- get_location_support(cases)
    for (mode in c("combined", "US")) {
      keep <- if (mode == "combined") support$location != "US" & support$complete &
        !support$boundary_peak else support$location == "US"
      cases_in_mode <- cases[keep]
      if (!length(cases_in_mode)) next
      message("Phase overview: ", disease, " h", horizon, " | ", mode)
      runs <- setNames(lapply(names(cases_in_mode), function(id) {
        # The selected analysis above has already fitted some of these histories.
        if (id %in% names(location_runs)) location_runs[[id]] else
          run_location(cases_in_mode[[id]], id, combined = mode == "combined")
      }), names(cases_in_mode))
      scores <- summarize_peak_windows(runs, combined = mode == "combined")
      scores$relative_loss <- NA_real_
      for (window in c("before", "after")) {
        i <- which(scores$window == window)
        base_loss <- scores$quantile_loss[i[scores$multiplier[i] == 0]]
        if (is.finite(base_loss) && base_loss > 0)
          scores$relative_loss[i] <- 100 * (scores$quantile_loss[i] / base_loss - 1)
        else warning("No positive baseline loss for ", disease, " h", horizon,
                     " ", mode, " ", window, "; curve omitted.")
      }
      overview_rows[[paste(disease, horizon, mode)]] <- data.frame(
        disease, horizon, mode, tracker_scale, scale_method,
        scores[c("window", "multiplier", "n_weeks", "n_locations", "quantile_loss",
                 "relative_loss", "coverage_error")])
    }
  }
  overview <- do.call(rbind, overview_rows)
  dir.create(dirname(overview_result_file), recursive = TRUE, showWarnings = FALSE)
  write.csv(overview, overview_result_file, row.names = FALSE)
  dir.create(dirname(overview_file), recursive = TRUE, showWarnings = FALSE)
  cairo_pdf(overview_file, width = 20.3125, height = 15.5, family = "DejaVu Sans")
  on.exit(dev.off())
  par(mfrow = c(4, 5), mar = c(4, 4, 3, 1), oma = c(4.5, 0, 4, 0),
      mgp = c(2.5, .7, 0), cex = .95)
  # Keep the original overview ticks when they lie within the configured grid.
  ticks <- sort(unique(c(range(rate_multipliers[positive]), .01, .03, .1, .3)))
  ticks <- ticks[ticks >= min(rate_multipliers[positive]) &
                 ticks <= max(rate_multipliers[positive])]
  # Each disease's loss row is followed immediately by its calibration row.
  # Shared calibration limits within a disease make the five horizons comparable.
  for (disease in c("covid", "flu")) for (metric in c("relative_loss", "coverage_error"))
      for (horizon in -1:3) {
    z <- overview[overview$disease == disease & overview$horizon == horizon &
                    overview$multiplier > 0, ]
    loss_panel <- metric == "relative_loss"
    values <- if (loss_panel) z[[metric]] else
      overview$coverage_error[overview$disease == disease & overview$multiplier > 0]
    finite_values <- values[is.finite(values)]
    limits <- if (length(finite_values)) range(0, finite_values) else
      if (loss_panel) c(-1, 1) else c(0, 1)
    if (!loss_panel) limits <- c(0, max(.01, limits[2]) * 1.05)
    if (diff(limits) == 0) limits <- limits + c(-1, 1)
    plot(range(rate_multipliers[positive]), limits, type = "n", log = "x", xaxt = "n",
         xlab = "Rate multiplier m",
         ylab = if (loss_panel) "Loss change versus base (%)" else "Calibration error",
         main = paste(toupper(disease), "| horizon", horizon, "|",
                      if (loss_panel) "loss" else "calibration"))
    axis(1, ticks, trimws(formatC(ticks, format = "fg", digits = 4)))
    abline(h = 0, col = "gray40", lty = 3)
    for (mode in c("combined", "US")) for (window in c("before", "after")) {
      curve <- z[z$mode == mode & z$window == window, ]
      lines(curve$multiplier, curve[[metric]],
            col = if (window == "before") "#D55E00" else "#0072B2",
            lty = if (mode == "combined") 1 else 2, lwd = 1.8)
    }
    if (!any(is.finite(z[[metric]])))
      text(sqrt(prod(range(rate_multipliers[positive]))), mean(limits), "No usable phase scores")
  }
  mtext("Fixed-rate MultiQT: loss and calibration before and after the winter peak",
        side = 3, outer = TRUE, line = 2, cex = 1.3)
  mtext(paste(if (use_log_scale) "Log1p-count tracking" else "Raw-count tracking",
              "| before = weeks", -weeks_before_peak, "to 0; after = weeks 1 to",
              weeks_after_peak, "| eta = m * S0 | Loss scale:", scale_method),
        side = 3, outer = TRUE, line = .5, cex = .9)
  par(xpd = NA)
  legend(grconvertX(.5, "ndc", "user"), grconvertY(.035, "ndc", "user"),
         c("Locations: before", "Locations: after", "US: before", "US: after"),
         col = c("#D55E00", "#0072B2", "#D55E00", "#0072B2"),
         lty = c(1, 1, 2, 2), lwd = 2, bty = "n", horiz = TRUE, xjust = .5, cex = .85)
  mtext("Loss: % change from base (negative is better). Calibration: mean |coverage - alpha| over locations and levels (lower is better).",
        side = 1, outer = TRUE, line = 3, cex = .85)
  overview
}
if (make_phase_overview) phase_overview <- build_phase_overview()
