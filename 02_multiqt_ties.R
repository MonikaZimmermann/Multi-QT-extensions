# install.packages(c("dplyr", "tidyr", "nanoparquet", "pak"))
# pak::pkg_install("cmu-delphi/epiprocess@main")

# 1. SETTINGS AND DATA ----------------------------------------------
.libPaths(c("tmp/multiqt-lib", .libPaths()))
library(epiprocess)
library(dplyr)
library(tidyr)

disease_name <- "covid"                       # "covid" or "flu"
location_id <- "18"                         
horizon_id <- 0                           
start_date <- as.Date("2024-11-23") # this is where the covid forecasts actually start, but might be needed later

# Divide this location's forecasts and weekly counts by one fixed historical scale.
# "median" or "max" uses reports available before the first forecast;
scale_method <- "median"

# With scaling, eta is a fraction of the location's reference count, not a count.
rate_min <- 0.000001
rate_max <- 0.1
n_positive_rates <- 81L
display_rates <- c(0, 0.005, 0.01,0.025,0.05, 0.1)

peak_window_weeks <- 6L                    # Each side of the second peak; exclude its week.

positive_grid <- 10^seq(log10(rate_min), log10(rate_max), length.out = n_positive_rates)
learning_rates <- sort(unique(signif(c(0, positive_grid, display_rates), 12)))
# Include displayed rates outside the log-spaced grid, plus the zero baseline.
# Log spacing means equal ratios BETWEEN rates; eta stays fixed WITHIN each run.

forecasts <- nanoparquet::read_parquet("cmu_timeseries_submitted_forecasts.parquet")
reports <- nanoparquet::read_parquet("nhsn_vintages.parquet")
selected_forecasts <- forecasts |>
  filter(disease == disease_name, location == location_id,
         horizon == horizon_id, reference_date >= start_date)
selected_reports <- reports |>
  filter(disease == disease_name, location == location_id)
# Organize the same reports for historical lookups
report_archive <- as_epi_archive(selected_reports, time_value = "target_end_date")

# Freeze the divisor before this run starts. Use one latest vintage per target
# week, so frequently revised weeks do not receive extra weight in the median.
# These reference counts are used only for scaling, not to initialize offsets.
scale_as_of <- min(selected_forecasts$reference_date) - 1L
scale_history <- epix_as_of(report_archive, version = scale_as_of)
scale_history <- scale_history[scale_history$time_value <= scale_as_of &
  is.finite(scale_history$value), ]
count_scale <- 1
scale_used <- scale_method
if (scale_method != "none" && nrow(scale_history)) {
  count_scale <- if (scale_method == "median") median(scale_history$value) else max(scale_history$value)
  # Small locations can have a zero median; fall back to their maximum, then 1.
  if (count_scale <= 0) {
    count_scale <- max(scale_history$value)
    scale_used <- "max"
  }
}
if (!is.finite(count_scale) || count_scale <= 0 || !nrow(scale_history)) {
  count_scale <- 1
  scale_used <- "none"
}
scale_label <- if (scale_used == "none") "Unscaled counts" else paste0("Counts / historical ", scale_used, " (", format(count_scale, trim = TRUE), ")")
count_axis_label <- if (scale_used == "none") "Hospitalizations" else paste("Hospitalizations /", scale_used)
loss_axis_label <- if (scale_used == "none") "Mean quantile loss" else "Mean quantile loss (scaled)"
rate_axis_label <- if (scale_used == "none") "Fixed eta (log scale)" else "Fixed eta (scaled units; log)"

# Keep the archive in raw counts; divide its feedback by this same fixed scale below.
selected_forecasts$value <- selected_forecasts$value / count_scale
selected_forecasts$truth_finalized <- selected_forecasts$truth_finalized / count_scale
alpha <- sort(unique(selected_forecasts$quantile_level))
example <- selected_forecasts |>
  filter(is.finite(truth_finalized)) |>
  select(reference_date, target_end_date, truth_finalized, quantile_level, value) |>
  pivot_wider(names_from = quantile_level, values_from = value) |>
  arrange(reference_date)
base <- as.matrix(example[, as.character(alpha)])  # Rows = weeks; columns = levels.
finalized <- example$truth_finalized               # Later truth, retained in feedback_log.
n_weeks <- nrow(base)
n_levels <- length(alpha)                          # 23 thresholds for ONE outcome.
cat(n_weeks, "weeks;", n_levels, "levels;", length(learning_rates), "rates\n")

# 2. WHAT FEEDBACK WAS AVAILABLE? ------------------------------------
# The numerical examples in the original comments below use raw count units.
# We treat example$reference_date[week] as the current forecast's issue date.
# Use only reports with reports$version < example$reference_date[week].
# For example, reference_date = Nov 30 allows reports$version <= Nov 29.
#
# Example timeline: Alabama COVID, horizon = 0; all dates below are in 2024.
# Horizon 0 means target_end_date = reference_date, not that truth is known yet.
# Nov 19: Archive begins with older counts, for target weeks through Nov 9.
# Nov 23: First forecast predicts the week ending Nov 23. No earlier forecast
#         in this experiment can provide feedback, so the adjustment is zero.
# Nov 27: The Nov 23 week's count is first reported as 70.
# Nov 29: That count is revised to 71. No forecast/update occurs on Nov 27 or 29.
# Nov 30: Look up reports through Nov 29. Compare 71 with the quantiles ISSUED
#         ON NOV 23, update the offsets, THEN issue the forecast for Nov 30.
# Dec  4: The Nov 30 week's count is first reported as 91.
# Dec  7: Look up reports through Dec 6. Compare 91 with the quantiles ISSUED
#         ON NOV 30, update the offsets, THEN issue the forecast for Dec 7.
# Thus reports for an earlier week adjust a later forecast; they do not rewrite
# the earlier issued forecast. feedback_observed is also used for scoring.
# Older reports already available on Nov 23 do not initialize the tracker:
# feedback requires an earlier issued forecast in this experiment to evaluate.
# If reports or forecast weeks are missing, feedback can take longer to be used.
feedback_by_week <- vector("list", n_weeks) # feedback_by_week[[week]] contains the row indices of earlier forecasts whose reported counts can now be used for their one update
feedback_observed <- rep(NA_real_, n_weeks) # stores the reported count used for each earlier forecast update
feedback_date <- rep(as.Date(NA), n_weeks) # stores the reference_date of the later forecast receiving that update—not the report’s version date. It is when a report is USED, not its version/release date.
# Targets without feedback_observed are omitted from scoring.

report_used <- rep(FALSE, n_weeks) # prevents repeated updates from the same earlier forecast
# It starts as FALSE. When that forecast is added to feedback_by_week, its entry becomes TRUE, 
# so subsequent weeks skip it—even if its reported count is revised.
for (week in seq_len(n_weeks)) {
  available <- epix_as_of(report_archive, version = example$reference_date[week] - 1L)
  counts <- available$value[match(example$target_end_date, available$time_value)] / count_scale
  ready <- which(seq_len(n_weeks) < week & !report_used & is.finite(counts))
  feedback_by_week[[week]] <- ready
  feedback_observed[ready] <- counts[ready]
  feedback_date[ready] <- example$reference_date[week]
  report_used[ready] <- TRUE
}
feedback_log <- data.frame(reference_date = example$reference_date,
                          target_end_date = example$target_end_date,
                          feedback_date, feedback_observed, truth_finalized = finalized)

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
  for (j in seq_along(alpha)) {
    covered[, j] <- feedback_observed <= adjusted[, j]
    coverage[j, i] <- mean(covered[, j], na.rm = TRUE) # average over weeks
    error <- feedback_observed - adjusted[, j]
    # Pinball loss penalizes underprediction by alpha, overprediction by 1-alpha.
    loss[, j] <- ifelse(error >= 0, alpha[j] * error, (1 - alpha[j]) * abs(error))
    quantile_loss_by_level[j, i] <- mean(loss[, j], na.rm = TRUE)
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
scores_by_level$n_weeks <- sum(is.finite(feedback_observed))
scores_by_level <- scores_by_level[c("rate", "alpha", "coverage",
                                    "signed_coverage_error", "coverage_error",
                                    "quantile_loss", "n_weeks")]
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
dir.create("output/pdf", recursive = TRUE, showWarnings = FALSE)
write.csv(results, "output/results/multiqt_learning_rates.csv", row.names = FALSE)
write.csv(scores_by_level, "output/results/multiqt_scores_by_level.csv", row.names = FALSE)
write.csv(feedback_log, "output/results/multiqt_feedback.csv", row.names = FALSE)
# Save the units alongside the existing CSVs without changing their column schemas.
write.csv(data.frame(disease = disease_name, location = location_id, horizon = horizon_id,
  scale_method, scale_used, count_scale, scale_as_of, reference_weeks = nrow(scale_history),
  reference_start = if (nrow(scale_history)) min(scale_history$time_value) else as.Date(NA),
  reference_end = if (nrow(scale_history)) max(scale_history$time_value) else as.Date(NA)),
  "output/results/multiqt_scale.csv", row.names = FALSE)
display_indices <- match(display_rates, learning_rates)
results[display_indices, ]
# Explore interactively: head(forecasts_by_rate[["10"]])

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
                           gap_axis_label = "Absolute coverage gap") {
  plot(learning_rates[positive], loss[positive], log = "x", type = "b", pch = 16,
       col = "steelblue4", cex = 0.45, ylim = loss_range, xaxt = "n",
       xlab = rate_axis_label, ylab = loss_axis_label,
       main = paste(label, "| Quantile error"))
  axis(1, at = rate_ticks, labels = rate_tick_labels)
  abline(h = loss[baseline], lty = 2, col = "orange")
  plot(learning_rates[positive], gap[positive], log = "x", type = "b", pch = 16,
       col = "steelblue4", cex = 0.45, ylim = range(0, gap_range), xaxt = "n",
       xlab = rate_axis_label, ylab = gap_axis_label,
       main = paste(label, "| Calibration error"))
  axis(1, at = rate_ticks, labels = rate_tick_labels)
  abline(h = gap[baseline], lty = 2, col = "gray40")
  plot(gap[positive], loss[positive], type = "l", col = "gray65",
       xlim = extendrange(gap_range), ylim = extendrange(loss_range),
       xlab = gap_axis_label, ylab = loss_axis_label,
       main = paste(label, "| Tradeoff"))
  points(gap[positive], loss[positive], pch = 16, cex = 0.6, col = rate_colours[positive])
  points(gap[baseline], loss[baseline], pch = 4, lwd = 2, col = "orange")
}
# pdf() starts writing vector graphics; dev.off() closes the file.
# Six alphas per page; PDF automatically starts a new page when panels fill.
pdf("output/pdf/multiqt_by_level.pdf", width = 14, height = 15)
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

pdf("output/pdf/multiqt_error_tradeoff.pdf", width = 12, height = 4.5)
par(mfrow = c(1, 3), mar = c(4, 4.5, 3, 1), oma = c(3, 0, 0, 0), mgp = c(2.3, 0.7, 0))
plot_error_row(results$quantile_loss, results$coverage_error, "Mean over levels",
               gap_axis_label = "Mean absolute coverage gap")
mtext(paste("eta =", display_rates), side = 1, outer = TRUE, line = 0.5,
      at = seq(0.12, 0.88, length.out = length(display_rates)),
      col = rate_colours[display_indices], cex = 0.8)
mtext(paste(scale_label, "| Dashed line / orange cross: eta = 0. Lower is better."),
      side = 1, outer = TRUE, line = 1.8, cex = 0.8)
dev.off()

pdf("output/pdf/multiqt_learning_rates.pdf", width = 7, height = 6)
matplot(alpha, coverage[, display_indices, drop = FALSE], type = "b", pch = 16, lty = 1,
        col = rate_colours[display_indices], xlim = c(0, 1), ylim = c(0, 1),
        xlab = "Nominal quantile level", ylab = "Coverage of feedback_observed",
        main = scale_label)
abline(0, 1, lty = 2, col = "gray40")   # Diagonal = agreement with nominal coverage.
legend("topleft", paste("eta =", display_rates), col = rate_colours[display_indices],
       pch = 16, lty = 1, bty = "n")
dev.off()

# 5. INSPECT THE FULL TIME SERIES AND ITS PEAKS -----------------------
season_weeks <- seq_len(n_weeks)
season_dates <- example$target_end_date[season_weeks]
season_observed <- feedback_observed[season_weeks]
# Peaks of feedback_observed within eight calendar weeks on either side.
peak_weeks <- integer()
for (week in season_weeks) {
  if (!is.finite(season_observed[week])) next
  days_apart <- as.numeric(season_dates - season_dates[week])
  nearby <- which(abs(days_apart) <= 8L * 7L & is.finite(season_observed))
  is_maximum <- week == nearby[which.max(season_observed[nearby])]
  lower_before <- any(days_apart[nearby] < 0 & season_observed[nearby] < season_observed[week])
  lower_after <- any(days_apart[nearby] > 0 & season_observed[nearby] < season_observed[week])
  if (is_maximum && lower_before && lower_after) {
    peak_weeks <- c(peak_weeks, week)
  }
}
peak_dates <- season_dates[peak_weeks]
plot_columns <- match(c(0.15, 0.5, 0.75), alpha)
plot_alpha <- alpha[plot_columns]
alpha_colours <- c("#0072B2", "#009E73", "#D55E00")
base_season <- base[season_weeks, plot_columns]
# Shared axes include negative forecasts: do not hide overcorrection by clipping.
plot_limits <- range(0, feedback_observed[season_weeks], base_season, na.rm = TRUE)
loss_limits <- c(0, 0)
season_coverage_by_rate <- list()
season_loss_by_rate <- list()
for (eta in display_rates) {
  key <- as.character(eta)
  q <- forecasts_by_rate[[key]][season_weeks, plot_columns, drop = FALSE]
  plot_limits <- range(plot_limits, q)
  covered <- covered_by_rate[[key]][season_weeks, plot_columns, drop = FALSE]
  season_coverage_by_rate[[key]] <- sweep(covered, 2, plot_alpha, "-")
  season_loss_by_rate[[key]] <- loss_by_rate[[key]][season_weeks, plot_columns, drop = FALSE]
  loss_limits <- range(loss_limits, season_loss_by_rate[[key]], na.rm = TRUE)
}
season_xlim <- range(season_dates)
date_ticks <- seq(as.Date(format(min(season_dates), "%Y-%m-01")),
                  max(season_dates), by = "3 months")
pdf("output/pdf/multiqt_season_forecasts.pdf", width = 12, height = 11)
layout(matrix(c(1, 4, 2, 5, 3, 6), nrow = 3, byrow = TRUE), heights = c(2, 1, 1))
par(oma = c(5, 0, 3, 0), mgp = c(2.5, 0.7, 0), cex = 0.9)
for (page_start in seq(1, length(display_indices), by = 2)) {
  page_indices <- display_indices[page_start:min(page_start + 1, length(display_indices))]
  for (i in page_indices) {
    eta <- learning_rates[i]
    key <- as.character(eta)
    q <- forecasts_by_rate[[key]][season_weeks, plot_columns, drop = FALSE]
    par(mar = c(1, 4.5, 2.5, 1))
    plot(season_dates, season_observed, type = "n", ylim = plot_limits,
         xlim = season_xlim, xaxt = "n", xaxs = "i", xlab = "",
         ylab = count_axis_label, main = paste("Fixed eta =", eta))
    polygon(c(season_dates, rev(season_dates)), c(q[, 1], rev(q[, 3])),
            col = adjustcolor(rate_colours[i], alpha.f = 0.2), border = NA)
    # Thick median lines; shaded 20--80% band = nominal central 60% interval.
    matlines(season_dates, base_season, col = "gray50", lty = 2, lwd = c(1, 2, 1))
    matlines(season_dates, q, col = rate_colours[i], lty = 1, lwd = c(1, 2, 1))
    lines(season_dates, season_observed, type = "o", pch = 16, cex = 0.5)
    abline(v = as.numeric(peak_dates), lty = 3, col = "gray40")
    legend("topright", c("feedback_observed", "Original", "Adjusted", "Peaks (retrospective)"),
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
            xlab = "target_end_date", ylab = if (scale_used == "none") "Quantile loss" else "Quantile loss (scaled)",
            main = "Weekly quantile loss")
    axis.Date(1, at = date_ticks, format = "%b\n%Y", cex.axis = 0.8)
    abline(v = as.numeric(peak_dates), col = "gray40", lty = 3)
  }
  if (length(page_indices) == 1L) {
    for (panel in 1:3) plot.new()
  }
  mtext(paste(disease_name, "| location", location_id, "| horizon", horizon_id,
              "|", min(season_dates), "to", max(season_dates)),
        side = 3, outer = TRUE, line = 1, cex = 1)
  mtext(paste(scale_label, "peaks within +/-8 weeks."),
        side = 3, outer = TRUE, line = -0.3, cex = 0.8)
  mtext(paste("alpha =", plot_alpha), side = 1, outer = TRUE,
        line = 0.5, at = c(0.2, 0.5, 0.8), col = alpha_colours, cex = 0.85)
  mtext("Coverage residual = 1(feedback_observed <= quantile) - alpha",
        side = 1, outer = TRUE, line = 2, cex = 0.75)
}
dev.off()
# Compare loss/coverage with these curves: do larger eta leave forecasts too high,
# too low, or unstable after the peak? This motivates future adaptive rates.

# 6. COMPARE ERRORS BEFORE AND AFTER EVERY PEAK ----------------------
empty_peak_results <- data.frame(peak_number = integer(), window = character(),
  peak_date = as.Date(character()), target_start = as.Date(character()),
  target_end = as.Date(character()), n_weeks = integer(), rate = numeric(),
  coverage_error = numeric(), quantile_loss = numeric())
score_peak_windows <- function(peak_number, peak_windows = NULL) {
  peak_date <- peak_dates[peak_number]
  window_results <- list()
  if (is.null(peak_windows)) peak_windows <- list(
    before = which(example$target_end_date >= peak_date - 7L * peak_window_weeks &
      example$target_end_date < peak_date & is.finite(feedback_observed)),
    after = which(example$target_end_date > peak_date &
      example$target_end_date <= peak_date + 7L * peak_window_weeks & is.finite(feedback_observed))
  )
  for (window in names(peak_windows)) {
    weeks <- peak_windows[[window]]
    scores <- data.frame(
      peak_number = peak_number, window = window, peak_date = peak_date,
      target_start = if (length(weeks)) min(example$target_end_date[weeks]) else as.Date(NA),
      target_end = if (length(weeks)) max(example$target_end_date[weeks]) else as.Date(NA),
      n_weeks = length(weeks), rate = learning_rates,
      coverage_error = NA_real_, quantile_loss = NA_real_)

    # Score the historical issued forecasts, retaining their full update history.
    # The tracker is never restarted at a peak or at a window boundary.
    if (length(weeks)) {
      for (i in seq_along(learning_rates)) {
        key <- as.character(learning_rates[i])
        covered <- covered_by_rate[[key]][weeks, , drop = FALSE]
        loss <- loss_by_rate[[key]][weeks, , drop = FALSE]
        scores$coverage_error[i] <- mean(abs(colMeans(covered) - alpha))
        scores$quantile_loss[i] <- mean(loss)
      }
    }
    window_results[[window]] <- scores
  }
  do.call(rbind, unname(window_results))
}

plot_peak_tradeoff <- function(peak_results, peak_label,
  window_labels = c(before = "Before peak", after = "After peak"),
  window_note = paste(peak_window_weeks, "calendar weeks on each side; peak week excluded. Short-window calibration is descriptive.")) {
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
    text(.5, .5, paste0("target_end_date: ", row$target_start[1], " to ", row$target_end[1],
                        " | ", row$n_weeks[1], " observed weeks"), cex = 1)
    par(mar = c(4, 4.5, 2.5, 1))
    if (row$n_weeks[1] > 0L) {
      plot_error_row(row$quantile_loss, row$coverage_error, label,
                     loss_range = shared_loss_range, gap_range = shared_gap_range,
                     gap_axis_label = "Mean absolute coverage gap")
    } else {
      for (panel in 1:3) {
        plot.new()
        text(.5, .5, paste(label, "- no observed weeks"))
      }
    }
  }
  mtext(paste(disease_name, "| location", location_id, "| horizon", horizon_id,
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
pdf("output/pdf/multiqt_error_tradeoff_all_peaks.pdf", width = 12, height = 8.5)
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
write.csv(all_peak_tradeoff_results, "output/results/multiqt_error_tradeoff_all_peaks.csv",
          row.names = FALSE)


