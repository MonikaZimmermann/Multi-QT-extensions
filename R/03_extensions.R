# Run R/01_data_overview.R once to prepare the data, then run this script.
# Run from the project root: Rscript R/03_extensions.R
# Reuse script 02's selectors when run in the same R session, or set them here.
# Each location keeps its own complete forecasting and feedback history.

# 1. SETTINGS AND PREPARED DATA -------------------------------------
combined_weeks <- -weeks_before_peak:weeks_after_peak
reference_rate <- .05
outer_rate <- .005
slope_rate_min <- .005
expert_rates <- c(0, .1 * 2^seq(-17, 0, by = .5))  # Base plus 35 fixed-rate experts.
H <- 4L                                  # Fixed lookback for expert-weight updates.
plot_levels <- (1:9) / 10
lambda_start <- 2.72
stopifnot(length(disease_name) == 1L, length(location_id) == 1L,
          length(horizon_id) == 1L, scale_method %in% c("median", "max", "none"))
if (!file.exists("output/data/peak_aligned.rds"))
  stop("Run R/01_data_overview.R first to create output/data/peak_aligned.rds.")
prepared <- readRDS("output/data/peak_aligned.rds")
selected_cases <- Filter(function(case) case$disease == disease_name &&
  case$horizon == horizon_id && if (combined) case$location != "US" else
    case$location == location_id, prepared$cases)
if (!length(selected_cases)) stop("No prepared cases match disease, location and horizon.")

# A common set of calendar weeks gives every retained location the same weight.
# Individual mode retains the full history, even for locations missing a peak-window week.
location_support <- do.call(rbind, lapply(names(selected_cases), function(case_id) {
  case <- selected_cases[[case_id]]
  indices <- match(combined_weeks, case$relative_week)
  observed <- !is.na(indices)
  observed[observed] <- is.finite(case$truth_latest[indices[observed]]) &
    rowSums(is.finite(case$base[indices[observed], , drop = FALSE])) == length(case$alpha) &
    case$within_winter[indices[observed]]
  data.frame(case_id, location = case$location, location_name = case$location_name,
    peak_date = case$peak_date, boundary_peak = case$peak_boundary,
    expected_weeks = length(combined_weeks),
    observed_weeks = sum(observed), complete = all(observed),
    missing_relative_weeks = paste(combined_weeks[!observed], collapse = ","))
}))
if (combined) {
  # A maximum on the winter edge is likely an earlier wave's tail, not a winter peak.
  selected_cases <- selected_cases[location_support$case_id[location_support$complete &
    !location_support$boundary_peak]]
  if (!length(selected_cases))
    stop("No location has complete feedback in the selected peak window; reduce weeks_before_peak/weeks_after_peak or choose individual mode.")
}
n_locations <- length(selected_cases)
output_dir <- if (combined) file.path("output/results/combined",
  paste0(disease_name, "_h", horizon_id), "03") else "output/results/peak_reasoning"
figure_dir <- if (combined) file.path("output/pdf/combined",
  paste0(disease_name, "_h", horizon_id), "03") else "output/pdf"
figure_file <- file.path(figure_dir, "multiqt_peak_reasoning.pdf")
overall_figure_file <- file.path(figure_dir, "multiqt_peak_reasoning_overall.pdf")
tradeoff_figure_file <- file.path(figure_dir, "multiqt_peak_reasoning_tradeoff.pdf")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(location_support, file.path(output_dir, "location_support.csv"), row.names = FALSE)

# 2. RUN EACH LOCATION INDEPENDENTLY --------------------------------
# Peaks are the largest latest-reported count in the prepared winter window.
# Only the target oracle uses that hindsight date to choose forecasts.
# Alignment changes evaluation coordinates, never when feedback becomes available.
extension_runs <- list()
error_rows <- list()
scale_rows <- list()
for (case_id in names(selected_cases)) {
  case <- selected_cases[[case_id]]
  alpha <- case$alpha
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
  base <- case$base / count_scale
  feedback_observed <- case$feedback_observed / count_scale
  truth_latest <- case$truth_latest / count_scale  # Scoring target; updates use feedback only.
  feedback_by_week <- case$feedback_by_week
  example <- data.frame(reference_date = case$reference_date,
                        target_end_date = case$target_date)
  dates <- example$target_end_date
  peak_date <- case$peak_date
  stopifnot(!is.na(peak_date), all(is.finite(base)))
  n_weeks <- nrow(base)
  n_levels <- length(alpha)
  observed_weeks <- which(is.finite(truth_latest))
  # 3. ONE TRACKER: FIXED RATE OR A RATE BASED ON ABSOLUTE SLOPE --------
  run_multiqt <- function(rates, base_forecasts = base) {
    if (length(rates) == 1) rates <- rep(rates, n_weeks)
    offsets <- numeric(n_levels)
    q <- base_forecasts
    history <- matrix(0, n_weeks, n_levels)
    for (week in seq_len(n_weeks)) {
      for (past_week in feedback_by_week[[week]]) {
        covered <- feedback_observed[past_week] <= q[past_week, ]
        offsets <- offsets + rates[week] * (alpha - covered)
      }
      q[week, ] <- isoreg(base_forecasts[week, ] + offsets)$yf
      history[week, ] <- offsets
    }
    list(q = q, offsets = history)
  }
  fixed <- run_multiqt(reference_rate)
  median_column <- match(0.5, alpha)
  slope <- c(0, diff(base[, median_column]) / (as.numeric(diff(example$reference_date)) / 7))
  # A deliberately simple heuristic: more movement means a larger current rate.
  # This is a slope of ORIGINAL forecasts available now, not a slope of future counts.
  slope_rates <- pmin(rate_max, slope_rate_min + abs(slope))
  adaptive <- run_multiqt(slope_rates)
  # A larger rate can unwind an existing offset faster once reversing feedback arrives.
  # But it can also build a larger offset earlier or amplify stale feedback now.
  # Large forecast slope does not itself imply a large forecast error.

  # 4. WHAT IF WE KNEW THE TRUE PEAK? ----------------------------------
  predictions <- list(original = base, fixed = fixed$q, absolute_slope = adaptive$q)
  # Bypass the correction completely; eta=0 would only freeze the accumulated offset.
  # Target oracle: replace every postpeak TARGET, including forecasts issued prepeak.
  # This is a hindsight benchmark that knows the peak before reports arrive.
  predictions$oracle_target <- fixed$q
  predictions$oracle_target[dates > peak_date, ] <- base[dates > peak_date, ]

  # 5. SEVERAL TRACKERS: AVERAGE THEIR FORECASTS ------------------------
  # Each expert keeps its OWN state and learns from its OWN historical forecast.
  # Averaging their rates and updating one state would be a different algorithm.
  # Pinball loss is convex in the issued quantile. This motivates combining
  # expert forecasts using past losses, but does not guarantee better calibration.
  experts <- list()
  for (eta in expert_rates) {
    experts[[as.character(eta)]] <- if (eta == reference_rate) fixed else run_multiqt(eta)
  }
  n_experts <- length(experts)
  pinball <- function(error, level) pmax(error, 0) * level + pmax(-error, 0) * (1 - level)
  # Use ordinary pinball loss in the count units already scaled for this location.
  # Each quantile has its own weights. The fixed H is the same for all quantiles.
  run_weights <- function(sharing = TRUE) {
    p <- matrix(1 / n_experts, n_experts, n_levels) # initial equal weights for each quantile
    log_weights <- log(p)
    lambda <- rep(lambda_start, n_levels)
    share <- if (sharing) 1 / (2 * H) else 0
    moments <- matrix(0, n_weeks, n_levels) # store loss summaries used to adapt lambda
    n_feedback <- 0L
    issued_weights <- array(NA_real_, c(n_weeks, n_experts, n_levels)) # store weights for each week, expert, and quantile used for each forecast
    q <- base
    for (week in seq_len(n_weeks)) {
      for (past_week in feedback_by_week[[week]]) {
        # Use the PREVIOUS H loss moments, excluding the feedback being processed now.
        # week = current forecast issue week; past_week = feedback week for that forecast.
        if (n_feedback >= H) {
          recent <- seq.int(n_feedback - H + 1L, n_feedback)
          denominator <- colSums(moments[recent, , drop = FALSE])
          positive <- denominator > 0
          lambda[positive] <- sqrt((log(2 * n_experts * H) + 1) / denominator[positive])
          # If recent losses are all zero, keep the previous finite learning rate.
        }
        loss <- matrix(0, n_experts, n_levels)
        for (i in seq_along(experts)) {
          error <- feedback_observed[past_week] - experts[[i]]$q[past_week, ]
          loss[i, ] <- pinball(error, alpha)
        }
        # Paper's second moment: average squared expert losses using pre-update weights.
        n_feedback <- n_feedback + 1L
        moments[n_feedback, ] <- colSums(p * loss^2)
        for (j in seq_along(alpha)) {
          # Exponential weights: an expert with a smaller released loss gains weight.
          log_weights[, j] <- log_weights[, j] - lambda[j] * loss[, j]
          log_weights[, j] <- log_weights[, j] - max(log_weights[, j])  # Numerical stability.
          p[, j] <- exp(log_weights[, j])
          p[, j] <- p[, j] / sum(p[, j])
          if (sharing) {
            # Fixed sharing gives every expert some weight again, including eta=0.
            # It does not reset offsets or simply discard the oldest losses.
            p[, j] <- (1 - share) * p[, j] + share / n_experts
            log_weights[, j] <- log(p[, j])
          }
          # Without sharing, retain log_weights even if exp() underflows.
        }
      }
      issued_weights[week, , ] <- p
      for (j in seq_along(alpha)) {
        values <- vapply(experts, function(expert) expert$q[week, j], numeric(1))
        q[week, j] <- sum(p[, j] * values)
      }
      # Quantile-specific weights can cause crossings; repair the issued forecast.
      q[week, ] <- isoreg(q[week, ])$yf
    }
    list(q = q, weights = issued_weights)
  }
  # Gibbs and Candes (2024), Algorithms 1-2, motivate exponential weighting and sharing.
  # Their deterministic average still uses LEARNED weights; equal averaging is a control.
  # H sets both the lookback and share=1/(2H); lambda starts at 2.72 for H updates.
  # https://jmlr.org/papers/v25/22-1218.html
  # https://github.com/isgibbs/DtACI/blob/main/DtACI.R

  shared <- run_weights()
  weighted <- run_weights(sharing = FALSE)  # Same H, with sharing removed.
  # Final layer: use the ordered sharing ensemble as a new base forecast.
  # This tracker has its own hidden offsets and learns from its OWN final forecasts.
  shared_outer <- run_multiqt(outer_rate, base_forecasts = shared$q)
  predictions$uniform <- Reduce(`+`, lapply(experts, function(expert) expert$q)) / n_experts
  predictions$weighted <- weighted$q
  predictions$shared <- shared$q
  predictions$shared_outer <- shared_outer$q



  # 6. KEEP THE ISSUED FORECASTS AND SIGNED WEEKLY DIAGNOSTICS --------
  # Score against the latest vintage; the online updates above saw only first-use feedback.
  # Counts are already scaled: do not divide quantile loss a second time.
  extension_runs[[case_id]] <- list(location = case$location, base = base, alpha = alpha,
    feedback_observed = feedback_observed, truth_latest = truth_latest,
    feedback_by_week = feedback_by_week,
    reference_date = case$reference_date, target_date = dates, peak_date = peak_date,
    relative_week = case$relative_week, predictions = predictions,
    weights = list(weighted = weighted$weights, shared = shared$weights),
    slope_rates = slope_rates, count_scale = count_scale, scale_used = scale_used)
  scale_rows[[case_id]] <- data.frame(location = case$location, scale_method, scale_used,
    count_scale, scale_as_of = case$scale_as_of, reference_weeks = length(scale_history))
  for (j in seq_along(alpha)) {
    for (method in names(predictions)) {
      q <- predictions[[method]][observed_weeks, j]
      covered <- truth_latest[observed_weeks] <= q
      error_rows[[paste(case_id, j, method)]] <- data.frame(location = case$location,
        target_date = dates[observed_weeks], relative_week = case$relative_week[observed_weeks],
        within_winter = case$within_winter[observed_weeks],
        alpha = alpha[j], method, calibration_residual = covered - alpha[j],
        quantile_error = pinball(truth_latest[observed_weeks] - q, alpha[j]))
    }
  }
}
location_errors <- do.call(rbind, error_rows)
rownames(location_errors) <- NULL
write.csv(location_errors, file.path(output_dir, "errors_by_location.csv"), row.names = FALSE)
write.csv(do.call(rbind, scale_rows), file.path(output_dir, "scales.csv"), row.names = FALSE)
saveRDS(extension_runs, file.path(output_dir, "issued_forecasts.rds"))

# 7. AVERAGE LOCATIONS ONLY AFTER CALCULATING THEIR OWN ERRORS -------
# Per-quantile weekly curves remain signed; the overall curve averages their absolute values.
# These weekly diagnostics are not time-window calibration-error scores.
# CE = mean over locations and quantiles of |mean residual over that location's weeks|.
# This prevents overcoverage in one location cancelling undercoverage elsewhere.
evaluation_errors <- if (combined)
  location_errors[location_errors$relative_week %in% combined_weeks, ] else location_errors
local_time_averages <- aggregate(cbind(calibration_residual, quantile_error) ~
  location + alpha + method, evaluation_errors, mean)
local_time_averages$calibration_error <- abs(local_time_averages$calibration_residual)
time_averages <- aggregate(cbind(calibration_residual, calibration_error, quantile_error) ~
  alpha + method, local_time_averages, mean)
# Pooled CE: |location-averaged signed residual| per level. Location-specific errors
# cancel, so it measures miscalibration shared by locations (always <= CE above).
time_averages$pooled_calibration_error <- abs(time_averages$calibration_residual)
overall_time_average <- aggregate(cbind(calibration_residual, calibration_error,
  pooled_calibration_error, quantile_error) ~ method, time_averages, mean)
write.csv(local_time_averages, file.path(output_dir, "time_averages_by_location.csv"), row.names = FALSE)
write.csv(rbind(time_averages, data.frame(alpha = NA_real_, overall_time_average)),
  file.path(output_dir, "time_averages.csv"), row.names = FALSE)

# The complete combined cohort contributes once at every relative week.
# In individual mode retain calendar dates and the original full-history view.
if (combined) {
  errors <- aggregate(cbind(calibration_residual, quantile_error) ~
    relative_week + alpha + method, evaluation_errors, mean)
  errors$plot_time <- errors$relative_week
} else {
  errors <- evaluation_errors
  errors$plot_time <- as.numeric(errors$target_date)
}
# Signed mean over locations and quantiles: the overall direction of miscoverage.
# Opposite gaps (e.g. too-narrow intervals) can cancel; see the per-quantile panels.
overall <- aggregate(cbind(mean_calibration_residual = calibration_residual,
  quantile_error) ~ plot_time + method, errors, mean)
write.csv(errors, file.path(output_dir, "errors_over_time.csv"), row.names = FALSE)
write.csv(overall, file.path(output_dir, "errors_average_quantiles.csv"), row.names = FALSE)

# Pre ends with the peak week (0); post starts one week later. Whole means full
# history for one location, or the common aligned window for the combined cohort.
summary_rows <- list()
summary_by_location <- list()
for (window in c("whole", "pre", "post")) {
  rows <- if (window == "whole") evaluation_errors else location_errors[
    location_errors$within_winter & location_errors$relative_week %in%
      if (window == "pre") -weeks_before_peak:0 else seq_len(weeks_after_peak), ]
  if (!nrow(rows)) {
    summary_rows[[window]] <- data.frame(method = names(predictions),
      calibration_residual = NA_real_, calibration_error = NA_real_, quantile_error = NA_real_,
      pooled_calibration_error = NA_real_, window, n_locations = 0L, min_weeks = 0L, max_weeks = 0L)
    next
  }
  local <- aggregate(cbind(calibration_residual, quantile_error) ~ location + alpha + method, rows, mean)
  local$calibration_error <- abs(local$calibration_residual)
  local$window <- window
  summary_by_location[[window]] <- local
  score <- aggregate(cbind(calibration_residual, calibration_error, quantile_error) ~ method, local, mean)
  # Pooled CE: average signed residuals over locations, then |.| per level, then levels.
  pooled <- aggregate(calibration_residual ~ alpha + method, local, mean)
  pooled <- aggregate(cbind(pooled_calibration_error = abs(calibration_residual)) ~ method, pooled, mean)
  score <- merge(score, pooled, by = "method")
  support <- aggregate(target_date ~ location, unique(rows[c("location", "target_date")]), length)
  score$window <- window
  score$n_locations <- nrow(support)
  score$min_weeks <- min(support$target_date)
  score$max_weeks <- max(support$target_date)
  summary_rows[[window]] <- score
}
summary_scores <- do.call(rbind, summary_rows)
write.csv(summary_scores, file.path(output_dir, "summary_scores.csv"), row.names = FALSE)
write.csv(do.call(rbind, summary_by_location),
  file.path(output_dir, "summary_scores_by_location_quantile.csv"), row.names = FALSE)

# 8. ORIGINAL RESIDUAL/LOSS PLOTS, WITH ONE WINTER PEAK ---------------
labels <- c(original = "Original", fixed = "Fixed MultiQT",
  absolute_slope = "Absolute-slope rate", oracle_target = "Target oracle (hindsight)",
  uniform = "Equal average", weighted = "Per-quantile, no sharing",
  shared = paste0("Per-quantile + sharing (H=", H, ")"),
  shared_outer = paste0("Sharing + outer MultiQT (rate=", outer_rate, ")"))
colours <- c(original = "gray40", fixed = "#0072B2", absolute_slope = "#009E73",
  oracle_target = "#CC79A7", uniform = "#E69F00", weighted = "#56B4E9", shared = "#882255",
  shared_outer = "#D55E00")
method_names <- names(labels)
axis_labels <- c(calibration_residual = "Coverage -\nquantile level",
  mean_calibration_residual = "Coverage - level,\nall quantiles", quantile_error = "Quantile loss")
mean_axis_labels <- c(calibration_residual = "Mean coverage -\nquantile level",
  mean_calibration_residual = "Time mean of\nsigned gap", quantile_error = "Mean quantile loss")
plot_levels <- plot_levels[plot_levels %in% alpha]
peak_marker <- if (combined) 0 else as.numeric(peak_date)
ticks <- if (combined) combined_weeks else
  pretty(as.Date(range(errors$plot_time), origin = "1970-01-01"), n = 3)
mean_limit <- 1.15 * max(.01, abs(time_averages$calibration_residual))
mean_limits <- list(calibration_residual = c(-mean_limit, mean_limit),
                    quantile_error = c(0, 1.15 * max(1e-8, time_averages$quantile_error)))
cohort_label <- if (combined) paste(n_locations, "locations with all", length(combined_weeks),
  "aligned weeks; US excluded") else paste(case$location_name, "| full history")
plot_title <- paste(disease_name, "| horizon", horizon_id, "|", cohort_label)
time_label <- if (combined) "Weeks from each location's winter maximum (peak = 0)" else "Target dates"
loss_axis_label <- if (scale_method == "none") "Quantile loss in counts" else "Quantile loss in locally scaled counts"
plot_errors <- function(rows, title, calibration_metric = "calibration_residual") {
  for (metric in c(calibration_metric, "quantile_error")) {
    residual <- metric != "quantile_error"
    # Axes follow the plotted values, like the loss panels; residual axes keep 0 visible.
    limits <- if (residual) 1.05 * range(-.01, .01, rows[[metric]]) else
      c(0, 1.05 * max(1e-8, rows[[metric]]))
    par(mar = if (residual) c(.5, 5.2, 2, 9) else c(3, 5.2, .5, 9))
    plot(rows$plot_time, rows[[metric]], type = "n", xaxt = "n", ylim = limits,
      xlab = "", ylab = axis_labels[metric], main = if (residual) title else "")
    if (residual) abline(h = 0, col = "gray85") else {
      if (combined) axis(1, at = ticks) else axis.Date(1, at = ticks, format = "%b %Y")
    }
    abline(v = peak_marker, col = "gray60", lty = 3)
    for (method in method_names) {
      line <- rows[rows$method == method, ]
      line <- line[order(line$plot_time), ]
      lines(line$plot_time, line[[metric]], col = colours[method], lwd = 1.5)
    }
    means <- vapply(method_names, function(method) mean(rows[[metric]][rows$method == method]), numeric(1))
    # Right-hand scale: shared across quantile panels; the overall panel uses its own means.
    mean_range <- if (metric == "mean_calibration_residual")
      c(-1, 1) * 1.15 * max(.01, abs(means)) else mean_limits[[metric]]
    mean_ticks <- pretty(mean_range, n = 3)
    mean_ticks <- mean_ticks[mean_ticks >= mean_range[1] & mean_ticks <= mean_range[2]]
    scale <- diff(limits) / diff(mean_range)
    offset <- limits[1] - scale * mean_range[1]
    edge <- par("usr")[2]
    inch <- diff(par("usr")[1:2]) / par("pin")[1]
    rect(edge + .04 * inch, limits[1], edge + .68 * inch, limits[2], col = "gray97", border = NA, xpd = NA)
    segments(edge + .04 * inch, offset, edge + .68 * inch, offset, col = "gray80", xpd = NA)
    points(edge + seq(.12, .62, length.out = length(means)) * inch,
      offset + scale * means, col = colours[names(means)], pch = 16, cex = .85, xpd = NA)
    axis(4, at = offset + scale * mean_ticks, labels = mean_ticks,
      pos = edge + .75 * inch, las = 1, cex.axis = .8, xpd = NA)
    mtext(mean_axis_labels[metric], side = 4, line = 7, cex = .8)
  }
}
panels <- matrix(0L, nrow = 2 * ceiling(length(plot_levels) / 3), ncol = 3)
for (j in seq_along(plot_levels)) {
  row <- 2 * ((j - 1L) %/% 3) + 1L
  column <- (j - 1L) %% 3 + 1L
  panels[row:(row + 1L), column] <- c(2 * j - 1L, 2 * j)
}
last_panel <- 2 * length(plot_levels)
panels <- rbind(panels, rep(last_panel + 1L, 3), rep(last_panel + 2L, 3), rep(last_panel + 3L, 3))
cairo_pdf(figure_file, width = 16, height = 16, family = "DejaVu Sans")
layout(panels, heights = c(rep(1, nrow(panels) - 1L), .7))
par(oma = c(2, 0, 4, 0), mgp = c(2.3, .6, 0), cex = 1)
for (level in plot_levels) plot_errors(errors[errors$alpha == level, ], paste("Quantile", level))
plot_errors(overall, paste("Signed average across all", n_levels, "quantiles"),
  calibration_metric = "mean_calibration_residual")
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", labels, col = colours, lty = 1, lwd = 2, ncol = 4, bty = "n")
mtext(plot_title, side = 3, outer = TRUE, line = 2, font = 2, cex = 1.3)
mtext("Gap = coverage - quantile level | All gaps signed | Overall: mean gap across quantiles (opposite gaps can cancel)",
  side = 3, outer = TRUE, line = .5)
mtext(paste(time_label, "|", loss_axis_label, "| Dotted line: winter peak"),
  side = 1, outer = TRUE, line = .5)
dev.off()

cairo_pdf(overall_figure_file, width = 12, height = 8, family = "DejaVu Sans")
layout(matrix(1:3, ncol = 1), heights = c(1, 1, .25))
par(oma = c(2, 0, 3, 0), mgp = c(2.5, .6, 0), cex = .9)
plot_errors(overall, paste("Signed average across all", n_levels, "quantiles"),
  calibration_metric = "mean_calibration_residual")
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", labels, col = colours, lty = 1, lwd = 2, ncol = 3, bty = "n", cex = .85)
mtext(plot_title, side = 3, outer = TRUE, line = 1, font = 2)
mtext(paste(time_label, "| Gap = coverage - quantile level; signed, averaged over locations and quantiles"),
  side = 1, outer = TRUE, line = .5, cex = .85)
dev.off()

# 9. TRUE CALIBRATION ERROR AND QUANTILE LOSS ON MATCHED WINDOWS -----
# Do not average data or signed gaps before the absolute value. Every location
# and quantile contributes one calibration gap; scale is already applied above.
cairo_pdf(tradeoff_figure_file, width = 14, height = 5.5, family = "DejaVu Sans")
layout(rbind(1:3, c(4, 4, 4)), heights = c(1, .22))
par(oma = c(2, 0, 3, 0), mgp = c(2.6, .7, 0))
window_titles <- c(whole = if (combined) "Whole aligned window" else "Full observed history",
                   pre = "Up to and including the winter peak", post = "After the winter peak")
xlim <- c(0, 1.1 * max(1e-8, summary_scores$quantile_error, na.rm = TRUE))
ylim <- c(0, 1.1 * max(.01, summary_scores$calibration_error, na.rm = TRUE))
for (window in c("whole", "pre", "post")) {
  rows <- summary_scores[summary_scores$window == window, ]
  par(mar = c(4, 4, 4, 1))
  plot(NA, xlim = xlim, ylim = ylim, xlab = "Mean quantile loss", ylab = "Calibration error",
    main = "")
  mtext(window_titles[window], side = 3, line = 1.8, font = 2, cex = .95)
  grid(col = "gray90")
  if (nrow(rows)) {
    support <- if (rows$min_weeks[1] == rows$max_weeks[1]) rows$min_weeks[1] else
      paste0(rows$min_weeks[1], "-", rows$max_weeks[1])
    mtext(paste(rows$n_locations[1], "locations;", support, "weeks each"), side = 3, line = .25, cex = .8)
    points(rows$quantile_error, rows$calibration_error, col = colours[rows$method], pch = 16, cex = 1.2)
    if (combined) {
      # Open circles: pooled CE of the same method, directly below its local CE.
      segments(rows$quantile_error, rows$pooled_calibration_error, rows$quantile_error,
        rows$calibration_error, col = adjustcolor(colours[rows$method], .4))
      points(rows$quantile_error, rows$pooled_calibration_error, col = colours[rows$method], pch = 1, cex = 1.2)
    }
    if (!any(is.finite(rows$calibration_error)))
      text(mean(xlim), mean(ylim), "No observations inside winter", cex = .8)
  }
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", labels, col = colours, pch = 16, ncol = 4, bty = "n", cex = .8)
mtext(plot_title, side = 3, outer = TRUE, line = 1, font = 2)
mtext(paste("CE: mean absolute local/quantile coverage gap.",
  if (combined) "Open circles: pooled CE (signed gaps averaged over locations first).",
  "Compare methods within matched windows; window length affects CE."),
  side = 1, outer = TRUE, line = .6, cex = .8)
dev.off()
cat("Script 03:", disease_name, "horizon", horizon_id, "|", cohort_label, "\n")
