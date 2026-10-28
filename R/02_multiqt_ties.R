# Run R/01_data_overview.R once to prepare output/data/peak_aligned.rds.
# Run from the project root: Rscript R/02_multiqt_ties.R
# The method and plotting code stays in this file; no helper scripts are sourced.
.libPaths(c("agent/tmp/multiqt-lib", .libPaths()))

# 1. SETTINGS AND PREPARED DATA -------------------------------------
disease_name <- "covid"                    # "covid" or "flu"; run separately.
location_id <- "06"                        # One location code; ignored if combined=TRUE.
horizon_id <- 0L
combined <- TRUE                         # TRUE averages independently fitted locations.
scale_method <- "max"                # "median", "max", or "none".
rate_min <- 0.000001
rate_max <- 0.1
n_positive_rates <- 81L
display_rates <- c(0, 0.005, 0.01, 0.025, 0.05, 0.1)
# No forecasts were issued 2025-10-25 to 2025-11-15, so a longer pre-peak window
# drops locations with December/early-January peaks. Only evaluation is windowed.
weeks_before_peak <- 3L                    # Before window also scores the peak week (0).
weeks_after_peak <- 8L
combined_weeks <- seq.int(-weeks_before_peak, weeks_after_peak)
positive_grid <- 10^seq(log10(rate_min), log10(rate_max), length.out = n_positive_rates)
learning_rates <- sort(unique(signif(c(0, positive_grid, display_rates), 12)))
# Log spacing is BETWEEN rates; eta stays fixed WITHIN each run.
peak_data <- readRDS("output/data/peak_aligned.rds")
start_date <- peak_data$settings$start_date
winter_start <- peak_data$settings$winter_start
winter_end <- peak_data$settings$winter_end
selected_cases <- Filter(function(x) x$disease == disease_name && x$horizon == horizon_id &&
  if (combined) x$location != "US" else x$location == location_id, peak_data$cases)
# Use the same complete calendar window for every location in a combined plot.
# Save all support flags, including excluded places, so the cohort is explicit.
location_support <- do.call(rbind, lapply(selected_cases, function(x) {
  i <- match(combined_weeks, x$relative_week)
  ok <- !is.na(i)
  ok[ok] <- x$within_winter[i[ok]] & is.finite(x$truth_latest[i[ok]]) &
    rowSums(is.finite(x$base[i[ok], , drop = FALSE])) == length(x$alpha)
  data.frame(location = x$location, location_name = x$location_name, peak_date = x$peak_date,
    boundary_peak = x$peak_boundary,
    n_observed = sum(ok), n_expected = length(combined_weeks), complete = all(ok))
}))
# A maximum on the winter edge is likely an earlier wave's tail, not a winter peak.
if (combined) selected_cases <- selected_cases[location_support$complete & !location_support$boundary_peak]
if (!length(selected_cases)) stop("No complete combined windows. Reduce weeks_before_peak/weeks_after_peak or inspect an individual location.")
result_dir <- if (combined) file.path("output/results/combined", paste0(disease_name, "_h", horizon_id), "02") else "output/results"
figure_dir <- if (combined) file.path("output/pdf/combined", paste0(disease_name, "_h", horizon_id), "02") else "output/pdf"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(location_support, file.path(result_dir, "multiqt_location_support.csv"), row.names = FALSE)

# 2. RESTORE EACH LOCATION'S INPUTS ----------------------------------
# 01 supplies full issue histories and first-use feedback in raw counts.
# Scale once here, then keep every location's offsets and past forecasts separate.
location_runs <- scale_rows <- list()
for (case_name in names(selected_cases)) {
  case <- selected_cases[[case_name]]
  scale_history <- case$scale_history
  count_scale <- 1
  scale_used <- scale_method
  if (scale_method != "none" && length(scale_history)) {
    count_scale <- if (scale_method == "median") median(scale_history) else max(scale_history)
    if (count_scale <= 0) { count_scale <- max(scale_history); scale_used <- "max" }
  }
  if (!is.finite(count_scale) || count_scale <= 0 || !length(scale_history)) {
    count_scale <- 1; scale_used <- "none"
  }
  alpha <- case$alpha
  base <- case$base / count_scale
  feedback_observed <- case$feedback_observed / count_scale
  truth_latest <- case$truth_latest / count_scale  # Scoring target; updates use feedback only.
  feedback_by_week <- case$feedback_by_week
  feedback_date <- case$feedback_date
  example <- data.frame(reference_date = case$reference_date, target_end_date = case$target_date,
    truth_finalized = case$truth_finalized / count_scale, relative_week = case$relative_week)
  n_weeks <- nrow(base)
  n_levels <- length(alpha)
  score_weeks <- which(is.finite(truth_latest) &
    if (combined) case$relative_week %in% combined_weeks & case$within_winter else TRUE)
  feedback_log <- data.frame(location = case$location, reference_date = example$reference_date,
    target_end_date = example$target_end_date, relative_week = case$relative_week,
    feedback_date, feedback_observed, truth_finalized = example$truth_finalized, truth_latest)
  scale_rows[[case_name]] <- data.frame(disease = disease_name, location = case$location,
    horizon = horizon_id, scale_method, scale_used, count_scale, scale_as_of = case$scale_as_of,
    reference_weeks = length(scale_history), reference_start = case$scale_reference_start,
    reference_end = case$scale_reference_end)
  # 3. RUN MULTIQT AND SCORE EACH ALPHA --------------------------------
  forecasts_by_rate <- list()
  covered_by_rate <- list()
  loss_by_rate <- list()
  coverage <- matrix(NA_real_, n_levels, length(learning_rates))
  quantile_loss_by_level <- matrix(NA_real_, n_levels, length(learning_rates))
  for (i in seq_along(learning_rates)) {
    eta <- learning_rates[i]
    offsets <- numeric(n_levels)               # Restart from zero for every rate.
    adjusted <- base
    for (week in seq_len(n_weeks)) {
      for (past_week in feedback_by_week[[week]]) {
        covered <- feedback_observed[past_week] <= adjusted[past_week, ]
        offsets <- offsets + eta * (alpha - covered)
        # Compare with the HISTORICAL ISSUED forecast. At alpha=.9 and eta=10,
        # outcome above that forecast adds 9; outcome at/below subtracts 1.
      }
      # Project only the issued forecasts, NOT the offsets ("lazy" MultiQT).
      # isoreg restores quantile order: e.g., (12,10) becomes (11,11).
      adjusted[week, ] <- isoreg(base[week, ] + offsets)$yf
    }
    key <- as.character(eta)
    forecasts_by_rate[[key]] <- adjusted
    covered <- matrix(NA, n_weeks, n_levels, dimnames = dimnames(base))
    loss <- adjusted
    # Score against the latest vintage; the offsets above saw only first-use feedback.
    for (j in seq_along(alpha)) {
      covered[, j] <- truth_latest <= adjusted[, j]
      coverage[j, i] <- mean(covered[score_weeks, j], na.rm = TRUE) # average over weeks
      error <- truth_latest - adjusted[, j]
      # Pinball loss penalizes underprediction by alpha, overprediction by 1-alpha.
      loss[, j] <- ifelse(error >= 0, alpha[j] * error, (1 - alpha[j]) * abs(error))
      quantile_loss_by_level[j, i] <- mean(loss[score_weeks, j], na.rm = TRUE)
    }
    # Reuse these weekly diagnostics in the timeline and peak comparisons.
    covered_by_rate[[key]] <- covered
    loss_by_rate[[key]] <- loss
  }
  # Calibration is checked per alpha. Positive signed gap = coverage above alpha.
  signed_coverage_error <- sweep(coverage, 1, alpha, "-") # coverage - alpha
  coverage_error_by_level <- abs(signed_coverage_error)
  results <- data.frame(rate = learning_rates,
                        coverage_error = colMeans(coverage_error_by_level),
                        quantile_loss = colMeans(quantile_loss_by_level))
  # expand.grid varies alpha fastest, matching the matrices' column order.
  scores_by_level <- expand.grid(alpha = alpha, rate = learning_rates)
  scores_by_level$coverage <- as.vector(coverage)
  scores_by_level$signed_coverage_error <- as.vector(signed_coverage_error)
  scores_by_level$coverage_error <- as.vector(coverage_error_by_level)
  scores_by_level$quantile_loss <- as.vector(quantile_loss_by_level)
  scores_by_level$n_weeks <- length(score_weeks)
  scores_by_level <- scores_by_level[c("rate", "alpha", "coverage",
                                      "signed_coverage_error", "coverage_error",
                                      "quantile_loss", "n_weeks")]
  # Retain the full issued forecasts; only scoring/plotting uses the peak window.
  location_runs[[case_name]] <- list(case = case, count_scale = count_scale, scale_used = scale_used,
    alpha = alpha, base = base, example = example, feedback_observed = feedback_observed,
    truth_latest = truth_latest, forecasts_by_rate = forecasts_by_rate, covered_by_rate = covered_by_rate,
    loss_by_rate = loss_by_rate, coverage = coverage,
    quantile_loss_by_level = quantile_loss_by_level, coverage_error_by_level = coverage_error_by_level,
    results = results, scores_by_level = scores_by_level, feedback_log = feedback_log)
}
n_locations <- length(location_runs)
# Absolute gaps are formed INSIDE each location above, before this extra average.
if (combined) {
  average_matrix <- function(name) Reduce(`+`, lapply(location_runs, `[[`, name)) / n_locations
  coverage <- average_matrix("coverage")                 # Descriptive mean coverage.
  coverage_error_by_level <- average_matrix("coverage_error_by_level")
  quantile_loss_by_level <- average_matrix("quantile_loss_by_level")
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
local_scores <- do.call(rbind, lapply(location_runs, function(x)
  data.frame(location = x$case$location, x$scores_by_level)))
write.csv(local_scores, file.path(result_dir, "multiqt_scores_by_location.csv"), row.names = FALSE)
scale_label <- if (combined) paste(n_locations, "locations |", scale_method,
  "scale within location | weeks", min(combined_weeks), "to", max(combined_weeks)) else
  if (scale_used == "none") "Unscaled counts" else paste0("Counts / historical ", scale_used, " (", format(count_scale, trim = TRUE), ")")
count_axis_label <- if (scale_method == "none") "Hospitalizations" else "Hospitalizations / historical scale"
loss_axis_label <- if (scale_method == "none") "Mean quantile loss" else "Mean quantile loss (scaled)"
rate_axis_label <- if (scale_method == "none") "Fixed eta (log scale)" else "Fixed eta (scaled units; log)"
location_label <- if (combined) paste(n_locations, "locations combined (US excluded)") else location_id
display_indices <- match(display_rates, learning_rates)
print(results[display_indices, ])

# 4. COMPARE ERRORS: --------------------------------
# This is exploratory, not held-out rate selection. Full-period coverage can
# hide seasonal errors. Averaging levels also hides which quantiles struggle;
# every supplied level has equal weight, despite their unequal alpha spacing.
positive <- which(learning_rates > 0)
# log(0) is undefined: zero appears as a horizontal baseline or tradeoff cross.
baseline <- match(0, learning_rates)
log_rates <- log10(learning_rates[positive])
# Label both endpoints; omit nearby decade labels so the maximum stays readable.
rate_ticks <- 10^seq(floor(min(log_rates)), ceiling(max(log_rates)))
rate_ticks <- rate_ticks[log10(rate_ticks) > min(log_rates) + 0.5 &
                         log10(rate_ticks) < max(log_rates) - 0.5]
rate_ticks <- sort(unique(c(range(learning_rates[positive]), rate_ticks)))
rate_tick_labels <- trimws(formatC(rate_ticks, format = "fg", digits = 6))
colour_index <- 1 + round(255 * (log_rates - min(log_rates)) / diff(range(log_rates)))
rate_colours <- rep("orange", length(learning_rates))
rate_colours[positive] <- hcl.colors(320, "Viridis")[colour_index]

# One small helper avoids duplicating the same three panels for every alpha.
# It reads the rate grid/colours above; loss and gap are vectors over that grid.
plot_error_row <- function(loss, gap, label,
                           loss_range = range(loss), gap_range = range(gap),
                           gap_axis_label = "Absolute coverage gap", pooled_gap = NULL) {
  plot(learning_rates[positive], loss[positive], log = "x", type = "b", pch = 16,
       col = "steelblue4", cex = 0.45, ylim = loss_range, xaxt = "n",
       xlab = rate_axis_label, ylab = loss_axis_label,
       main = paste(label, "| Quantile error"))
  axis(1, at = rate_ticks, labels = rate_tick_labels)
  abline(h = loss[baseline], lty = 2, col = "orange")
  plot(learning_rates[positive], gap[positive], log = "x", type = "b", pch = 16,
       col = "steelblue4", cex = 0.45, xaxt = "n",  # Headroom for the pooled legend.
       ylim = range(0, gap_range) * if (is.null(pooled_gap)) 1 else c(1, 1.3),
       xlab = rate_axis_label, ylab = gap_axis_label,
       main = paste(label, "| Calibration error"))
  axis(1, at = rate_ticks, labels = rate_tick_labels)
  abline(h = gap[baseline], lty = 2, col = "gray40")
  if (!is.null(pooled_gap)) {
    lines(learning_rates[positive], pooled_gap[positive], type = "b", pch = 1,
          col = "firebrick", cex = 0.45)
    abline(h = pooled_gap[baseline], lty = 2, col = "firebrick")
    legend("topright", c("Per-location CE", "Pooled CE"), col = c("steelblue4", "firebrick"),
           pch = c(16, 1), lty = 1, cex = 0.7, bty = "n")
  }
  plot(gap[positive], loss[positive], type = "l", col = "gray65",
       xlim = extendrange(gap_range), ylim = extendrange(loss_range),
       xlab = gap_axis_label, ylab = loss_axis_label,
       main = paste(label, "| Tradeoff"))
  points(gap[positive], loss[positive], pch = 16, cex = 0.6, col = rate_colours[positive])
  points(gap[baseline], loss[baseline], pch = 4, lwd = 2, col = "orange")
}
# pdf() starts writing vector graphics; dev.off() closes the file.
# Six alphas per page; PDF automatically starts a new page when panels fill.
cairo_pdf(file.path(figure_dir, "multiqt_by_level.pdf"), width = 14, height = 15, family = "DejaVu Sans")
par(mfrow = c(6, 3), mar = c(3.5, 4, 2.5, 1), oma = c(4, 0, 0, 0),
    mgp = c(2.2, 0.7, 0), cex = 0.85)
for (j in seq_along(alpha)) {
  plot_error_row(quantile_loss_by_level[j, ], coverage_error_by_level[j, ],
                 paste("alpha =", alpha[j]))
  if (j %% 6 == 0 || j == n_levels) {
    mtext(paste("eta =", display_rates), side = 1, outer = TRUE, line = 1,
          at = seq(0.12, 0.88, length.out = length(display_rates)),
          col = rate_colours[display_indices], cex = 0.8)
    mtext(paste(scale_label, "| Dashed line / orange cross: eta = 0. Lower is better."),
          side = 1, outer = TRUE, line = 2.5, cex = 0.8)
  }
}
dev.off()

cairo_pdf(file.path(figure_dir, "multiqt_error_tradeoff.pdf"), width = 12, height = 4.5, family = "DejaVu Sans")
par(mfrow = c(1, 3), mar = c(4, 4.5, 3, 1), oma = c(3, 0, 0, 0), mgp = c(2.3, 0.7, 0))
plot_error_row(results$quantile_loss, results$coverage_error, "Mean over levels",
               gap_axis_label = "Mean absolute coverage gap",
               pooled_gap = if (combined) results$pooled_coverage_error)
mtext(paste("eta =", display_rates), side = 1, outer = TRUE, line = 0.5,
      at = seq(0.12, 0.88, length.out = length(display_rates)),
      col = rate_colours[display_indices], cex = 0.8)
mtext(paste(scale_label, "| Dashed line / orange cross: eta = 0. Lower is better."),
      side = 1, outer = TRUE, line = 1.8, cex = 0.8)
dev.off()

cairo_pdf(file.path(figure_dir, "multiqt_learning_rates.pdf"), width = 7, height = 6, family = "DejaVu Sans")
matplot(alpha, coverage[, display_indices, drop = FALSE], type = "b", pch = 16, lty = 1,
        col = rate_colours[display_indices], xlim = c(0, 1), ylim = c(0, 1),
        xlab = "Nominal quantile level", ylab = "Coverage of latest vintage",
        main = paste(strwrap(scale_label, width = 48), collapse = "\n"))
abline(0, 1, lty = 2, col = "gray40")   # Diagonal = agreement with nominal coverage.
legend("topleft", paste("eta =", display_rates), col = rate_colours[display_indices],
       pch = 16, lty = 1, bty = "n")
dev.off()

# 5. INSPECT THE FULL TIME SERIES AND ITS PEAKS -----------------------
plot_columns <- match(c(0.15, 0.5, 0.75), alpha)
plot_alpha <- alpha[plot_columns]
alpha_colours <- c("#0072B2", "#009E73", "#D55E00")
season_forecasts_by_rate <- season_coverage_by_rate <- season_loss_by_rate <- list()
if (combined) {
  season_dates <- combined_weeks
  peak_dates <- 0
  aligned_mean <- function(extract) Reduce(`+`, lapply(location_runs, function(x) {
    i <- match(combined_weeks, x$case$relative_week)
    extract(x, i)
  })) / n_locations
  season_observed <- aligned_mean(function(x, i) x$truth_latest[i])
  base_season <- aligned_mean(function(x, i) x$base[i, plot_columns, drop = FALSE])
  for (eta in display_rates) {
    key <- as.character(eta)
    season_forecasts_by_rate[[key]] <- aligned_mean(function(x, i) x$forecasts_by_rate[[key]][i, plot_columns, drop = FALSE])
    covered <- aligned_mean(function(x, i) x$covered_by_rate[[key]][i, plot_columns, drop = FALSE])
    season_coverage_by_rate[[key]] <- sweep(covered, 2, plot_alpha, "-")
    season_loss_by_rate[[key]] <- aligned_mean(function(x, i) x$loss_by_rate[[key]][i, plot_columns, drop = FALSE])
  }
  date_ticks <- pretty(season_dates)
} else {
  season_weeks <- seq_len(n_weeks)
  season_dates <- example$target_end_date
  season_observed <- truth_latest
  peak_dates <- case$peak_date
  base_season <- base[, plot_columns, drop = FALSE]
  for (eta in display_rates) {
    key <- as.character(eta)
    season_forecasts_by_rate[[key]] <- forecasts_by_rate[[key]][, plot_columns, drop = FALSE]
    season_coverage_by_rate[[key]] <- sweep(covered_by_rate[[key]][, plot_columns, drop = FALSE], 2, plot_alpha, "-")
    season_loss_by_rate[[key]] <- loss_by_rate[[key]][, plot_columns, drop = FALSE]
  }
  date_ticks <- seq(as.Date(format(min(season_dates), "%Y-%m-01")), max(season_dates), by = "3 months")
}
# Averaged forecasts/outcomes are display profiles only. Scores were computed locally.
plot_limits <- range(0, season_observed, base_season, unlist(season_forecasts_by_rate), na.rm = TRUE)
loss_limits <- range(0, unlist(season_loss_by_rate), na.rm = TRUE)
season_xlim <- range(season_dates)
cairo_pdf(file.path(figure_dir, "multiqt_season_forecasts.pdf"), width = 12, height = 11, family = "DejaVu Sans")
layout(matrix(c(1, 4, 2, 5, 3, 6), nrow = 3, byrow = TRUE), heights = c(2, 1, 1))
par(oma = c(5, 0, 3, 0), mgp = c(2.5, 0.7, 0), cex = 0.9)
for (page_start in seq(1, length(display_indices), by = 2)) {
  page_indices <- display_indices[page_start:min(page_start + 1, length(display_indices))]
  for (i in page_indices) {
    eta <- learning_rates[i]
    key <- as.character(eta)
    q <- season_forecasts_by_rate[[key]]
    par(mar = c(1, 4.5, 2.5, 1))
    plot(season_dates, season_observed, type = "n", ylim = plot_limits,
         xlim = season_xlim, xaxt = "n", xaxs = "i", xlab = "",
         ylab = count_axis_label, main = paste("Fixed eta =", eta))
    polygon(c(season_dates, rev(season_dates)), c(q[, 1], rev(q[, 3])),
            col = adjustcolor(rate_colours[i], alpha.f = 0.2), border = NA)
    # Thick median lines; shading spans the selected lower/upper quantiles.
    matlines(season_dates, base_season, col = "gray50", lty = 2, lwd = c(1, 2, 1))
    matlines(season_dates, q, col = rate_colours[i], lty = 1, lwd = c(1, 2, 1))
    lines(season_dates, season_observed, type = "o", pch = 16, cex = 0.5)
    abline(v = as.numeric(peak_dates), lty = 3, col = "gray40")
    legend("topright", c("Latest vintage", "Original", "Adjusted", "Peaks (retrospective)"),
           col = c("black", "gray50", rate_colours[i], "gray40"),
           lty = c(1, 2, 1, 3), pch = c(16, NA, NA, NA), cex = 0.7, bty = "n")

    par(mar = c(1, 4.5, 2.5, 1))
    matplot(season_dates, season_coverage_by_rate[[key]], type = "b",
            col = alpha_colours, lty = 1, pch = c(1, 16, 17), cex = 0.5,
            ylim = c(-1, 1), xlim = season_xlim, xaxt = "n", xaxs = "i",
            xlab = "", ylab = "Covered - alpha", main = "Weekly coverage residual")
    abline(h = 0, col = "gray65", lty = 2)
    abline(v = as.numeric(peak_dates), col = "gray40", lty = 3)

    par(mar = c(3.5, 4.5, 2.5, 1))
    matplot(season_dates, season_loss_by_rate[[key]], type = "b",
            col = alpha_colours, lty = 1, pch = c(1, 16, 17), cex = 0.5,
            ylim = loss_limits, xlim = season_xlim, xaxt = "n", xaxs = "i",
            xlab = if (combined) "Weeks from local winter maximum" else "target_end_date", ylab = if (scale_used == "none") "Quantile loss" else "Quantile loss (scaled)",
            main = "Weekly quantile loss")
    if (combined) axis(1, at = date_ticks) else axis.Date(1, at = date_ticks, format = "%b\n%Y", cex.axis = 0.8)
    abline(v = as.numeric(peak_dates), col = "gray40", lty = 3)
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
  mtext(if (combined) "Profiles and signed residuals are descriptive location averages; CE is calculated locally first." else "Coverage residual = 1(latest vintage <= quantile) - alpha",
        side = 1, outer = TRUE, line = 2, cex = 0.75)
}
dev.off()
# Compare loss/coverage with these curves: do larger eta leave forecasts too high,
# too low, or unstable after the peak? This motivates future adaptive rates.

# 6. COMPARE ERRORS BEFORE AND AFTER THE WINTER MAXIMUM ----------------------
empty_peak_results <- data.frame(peak_number = integer(), window = character(),
  peak_date = as.Date(character()), target_start = as.Date(character()),
  target_end = as.Date(character()), n_weeks = integer(), rate = numeric(),
  coverage_error = numeric(), pooled_coverage_error = numeric(), quantile_loss = numeric())
score_peak_windows <- function(peak_number, peak_windows = NULL) {
  rows <- lapply(c("before", "after"), function(window) {
    local <- lapply(location_runs, function(x) {
      relative <- x$case$relative_week
      weeks <- which(x$case$within_winter & is.finite(x$truth_latest) &
        if (window == "before") relative >= -weeks_before_peak & relative <= 0 else
          relative > 0 & relative <= weeks_after_peak)
      if (!is.null(peak_windows) && !combined) weeks <- peak_windows[[window]]
      gaps <- losses <- rep(NA_real_, length(learning_rates))
      signed <- matrix(NA_real_, length(learning_rates), length(alpha))  # rate x level
      if (length(weeks)) for (i in seq_along(learning_rates)) {
        key <- as.character(learning_rates[i])
        signed[i, ] <- colMeans(x$covered_by_rate[[key]][weeks, , drop = FALSE]) - alpha
        gaps[i] <- mean(abs(signed[i, ]))
        losses[i] <- mean(x$loss_by_rate[[key]][weeks, , drop = FALSE])
      }
      list(gap = gaps, loss = losses, signed = signed, weeks = weeks, dates = x$case$target_date[weeks])
    })
    first <- local[[1L]]
    bounds <- if (combined) if (window == "before") c(-weeks_before_peak, 0) else c(1, weeks_after_peak) else
      if (length(first$weeks)) range(first$dates) else rep(as.Date(NA), 2L)
    data.frame(peak_number, window, peak_date = if (combined) 0 else case$peak_date,
      target_start = bounds[1], target_end = bounds[2], n_weeks = length(first$weeks),
      n_locations = sum(vapply(local, function(x) length(x$weeks) > 0L, logical(1))), rate = learning_rates,
      coverage_error = Reduce(`+`, lapply(local, `[[`, "gap")) / length(local),
      pooled_coverage_error = rowMeans(abs(Reduce(`+`, lapply(local, `[[`, "signed")) / length(local))),
      quantile_loss = Reduce(`+`, lapply(local, `[[`, "loss")) / length(local))
  })
  do.call(rbind, rows)
}

plot_peak_tradeoff <- function(peak_results, peak_label,
  window_labels = c(before = "Up to peak", after = "After peak"),
  window_note = paste(weeks_before_peak, "calendar weeks before plus the peak week;", weeks_after_peak,
    "weeks after. Short-window calibration is descriptive.")) {
  finite_loss <- peak_results$quantile_loss[is.finite(peak_results$quantile_loss)]
  finite_gap <- peak_results$coverage_error[is.finite(peak_results$coverage_error)]
  shared_loss_range <- if (length(finite_loss)) range(finite_loss) else c(0, 1)
  shared_gap_range <- if (length(finite_gap)) range(finite_gap) else c(0, 1)
  layout(rbind(c(1, 1, 1), c(2, 3, 4), c(5, 5, 5), c(6, 7, 8)),
         heights = c(.14, 1, .14, 1))
  par(oma = c(5, 0, 3, 0), mgp = c(2.3, .7, 0), cex = .9)
  for (window in c("before", "after")) {
    row <- peak_results[peak_results$window == window, ]
    label <- window_labels[[window]]
    par(mar = c(0, 0, 0, 0))
    plot.new()
    support_note <- if (row$n_weeks[1] > 0L) paste0(
      if (combined) "Relative week: " else "target_end_date: ", row$target_start[1], " to ", row$target_end[1],
      " | ", row$n_weeks[1], " observed weeks/location | ", row$n_locations[1], " location(s)") else
      "No observations in this winter window"
    text(.5, .5, support_note, cex = 1)
    par(mar = c(4, 4.5, 2.5, 1))
    if (row$n_weeks[1] > 0L) {
      plot_error_row(row$quantile_loss, row$coverage_error, label,
                     loss_range = shared_loss_range, gap_range = shared_gap_range,
                     gap_axis_label = "Mean absolute coverage gap",
                     pooled_gap = if (combined) row$pooled_coverage_error)
    } else {
      for (panel in 1:3) {
        plot.new()
        text(.5, .5, paste(label, "- no observed weeks"))
      }
    }
  }
  mtext(paste(disease_name, "| location", location_label, "| horizon", horizon_id,
              "|", peak_label, peak_results$peak_date[1]),
        side = 3, outer = TRUE, line = 1, cex = 1)
  mtext(paste("eta =", display_rates), side = 1, outer = TRUE, line = .5,
        at = seq(.12, .88, length.out = length(display_rates)),
        col = rate_colours[display_indices], cex = .8)
  mtext(paste(scale_label, "| All levels; shared axes. Dashed lines / orange cross: eta = 0."),
        side = 1, outer = TRUE, line = 1.8, cex = .8)
  mtext(window_note,
        side = 1, outer = TRUE, line = 3.1, cex = .8)
}

all_peak_tradeoff_results <- if (length(peak_dates)) {
  do.call(rbind, lapply(seq_along(peak_dates), score_peak_windows))
} else empty_peak_results

# One chronological peak per page, using the same before/after panel layout.
# Axes are shared between the two windows on each page, as in the original.
cairo_pdf(file.path(figure_dir, "multiqt_error_tradeoff_all_peaks.pdf"), width = 12, height = 8.5, family = "DejaVu Sans")
if (!length(peak_dates)) {
  plot.new()
  text(.5, .55, "Peak comparison unavailable", font = 2)
  text(.5, .45, "No peaks were detected in feedback_observed.")
} else {
  for (peak_number in seq_along(peak_dates)) {
    plot_peak_tradeoff(
      all_peak_tradeoff_results[all_peak_tradeoff_results$peak_number == peak_number, ],
      paste0("peak ", peak_number, " of ", length(peak_dates), ":"))
  }
}
dev.off()
write.csv(all_peak_tradeoff_results, file.path(result_dir, "multiqt_error_tradeoff_all_peaks.csv"),
          row.names = FALSE)
