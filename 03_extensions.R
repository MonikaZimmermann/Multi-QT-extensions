# Run sections 1 and 2 of 02_multiqt_ties.R first

# 1. EXTENSION SETTINGS ---------------------------------------------
peak_number <- 2L
reference_rate <- 0.05                     # One of the rates displayed in script 02.
outer_rate <- 0.0050
# Fixed small rate for the final MultiQT layer.
slope_rate_min <- 0.005
expert_rates <- c(0, .1 * 2^seq(-17, 0, by = .5))  # Base plus 35 fixed-rate experts.
H <- 4L                                  # Fixed lookback for the expert-weight updates.
plot_levels <- (1:9) / 10                 # Show 0.1, 0.2, ..., 0.9 and the overall average.
lambda_start <- 2.72                      # Initial strength of the weight updates.
output_dir <- "output/results/peak_reasoning"
figure_file <- "output/pdf/multiqt_peak_reasoning.pdf"
overall_figure_file <- "output/pdf/multiqt_peak_reasoning_overall.pdf"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)

# 2. THE PEAK FOR THE TARGET ORACLE ---------------------------------
# These outcome-based labels are retrospective. Only the target oracle
# uses the peak when choosing forecasts; the other methods never see it.
dates <- example$target_end_date
peak_weeks <- integer()
for (week in which(is.finite(feedback_observed))) {
  distance <- as.numeric(dates - dates[week]) # time diff in days
  nearby <- which(abs(distance) <= 56 & is.finite(feedback_observed))
  if (week == nearby[which.max(feedback_observed[nearby])] &&
      any(distance[nearby] < 0 & feedback_observed[nearby] < feedback_observed[week]) &&
      any(distance[nearby] > 0 & feedback_observed[nearby] < feedback_observed[week])) {
    peak_weeks <- c(peak_weeks, week)
  }
}
stopifnot(length(peak_weeks) >= peak_number)
peak_week <- peak_weeks[peak_number]
peak_date <- dates[peak_week]
peak_dates <- dates[peak_weeks]            # Mark all peaks; the oracle still uses peak_number.
observed_weeks <- which(is.finite(feedback_observed))

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
# Use ordinary pinball loss in the count units already scaled by script 02.
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


# 6. WEEKLY RESIDUAL AND LOSS, FOR EACH QUANTILE ---------------------
# There is no calibration split: every method is evaluated on the same full history.
# Calibration residual = 1{y <= q} - alpha, retaining its sign at each target date.
# Its time average is the signed coverage gap; one residual is not a calibration estimate.
# Quantile error = pinball loss, in scaled counts, for that SAME target date.
# These are retrospective evaluation curves; they are not inputs to the weight updates.
errors <- list()
for (j in seq_along(alpha)) {
  for (method in names(predictions)) {
    q <- predictions[[method]][observed_weeks, j]
    covered <- feedback_observed[observed_weeks] <= q
    loss <- pinball(feedback_observed[observed_weeks] - q, alpha[j])
    errors[[paste(j, method)]] <- data.frame(target_date = dates[observed_weeks],
      alpha = alpha[j], method, calibration_residual = covered - alpha[j], quantile_error = loss)
  }
}
errors <- do.call(rbind, errors)
write.csv(errors, file.path(output_dir, "errors_over_time.csv"), row.names = FALSE)

# 7. AVERAGE ACROSS QUANTILES AND ACROSS TIME ------------------------
# Give every quantile equal weight, with no averaging over time in these curves.
# The residual remains signed, so opposite residuals can cancel across quantiles.
overall <- aggregate(errors[c("calibration_residual", "quantile_error")],
  by = errors[c("target_date", "method")], FUN = mean)
write.csv(overall, file.path(output_dir, "errors_average_quantiles.csv"), row.names = FALSE)
# The right-hand markers summarize the entire observed period, for each method.
time_averages <- aggregate(cbind(calibration_residual, quantile_error) ~ alpha + method, errors, mean)
overall_time_average <- aggregate(cbind(calibration_residual, quantile_error) ~ method, overall, mean)
write.csv(rbind(time_averages, data.frame(alpha = NA_real_, overall_time_average)),
  file.path(output_dir, "time_averages.csv"), row.names = FALSE)

# 8. SELECTED QUANTILES AND THE OVERALL AVERAGE ON ONE PAGE ----------
labels <- c(original = "Original", fixed = "Fixed MultiQT",
  absolute_slope = "Absolute-slope rate", oracle_target = "Target oracle (hindsight)",
  uniform = "Equal average", weighted = "Per-quantile, no sharing",
  shared = paste0("Per-quantile + sharing (H=", H, ")"),
  shared_outer = paste0("Sharing + outer MultiQT (rate=", outer_rate, ")"))
colours <- c(original = "gray40", fixed = "#0072B2", absolute_slope = "#009E73",
  oracle_target = "#CC79A7", uniform = "#E69F00", weighted = "#56B4E9", shared = "#882255",
  shared_outer = "#D55E00")
axis_labels <- c(calibration_residual = "Residual", quantile_error = "Quantile loss")
ticks <- pretty(dates[observed_weeks], n = 3)
# Each metric has the same zoomed right-hand scale across all quantiles.
mean_limit <- 1.15 * max(.01, abs(time_averages$calibration_residual))
mean_limits <- list(calibration_residual = c(-mean_limit, mean_limit),
                    quantile_error = c(0, 1.15 * max(1e-8, time_averages$quantile_error)))
plot_errors <- function(rows, title) {
  for (metric in names(axis_labels)) {
    residual <- metric == "calibration_residual"
    limits <- if (residual) c(-1, 1) else c(0, 1.05 * max(1e-8, rows[[metric]]))
    par(mar = if (residual) c(.5, 4, 2, 7) else c(3, 4, .5, 7))
    plot(rows$target_date, rows[[metric]], type = "n", xaxt = "n",
      ylim = limits,
      xlab = "", ylab = axis_labels[metric], main = if (residual) title else "")
    if (residual) abline(h = 0, col = "gray85") else axis.Date(1, at = ticks, format = "%b %Y")
    abline(v = as.numeric(peak_dates), col = "gray60", lty = 3)
    for (method in names(predictions)) {
      line <- rows[rows$method == method, ]
      lines(line$target_date, line[[metric]], col = colours[method], lwd = 1.5)
    }
    means <- vapply(names(predictions), function(method)
      mean(rows[[metric]][rows$method == method]), numeric(1))
    # Place the means outside the timeline, with one horizontal position per method.
    # Map the separate right-hand scale to the plotting coordinates of this panel.
    mean_range <- mean_limits[[metric]]
    mean_ticks <- pretty(mean_range, n = 3)
    mean_ticks <- mean_ticks[mean_ticks >= mean_range[1] & mean_ticks <= mean_range[2]]
    scale <- diff(limits) / diff(mean_range)
    offset <- limits[1] - scale * mean_range[1]
    edge <- par("usr")[2]
    inch <- diff(par("usr")[1:2]) / par("pin")[1]
    rect(edge + .04 * inch, limits[1], edge + .68 * inch, limits[2],
      col = "gray97", border = NA, xpd = NA)
    segments(edge + .04 * inch, offset, edge + .68 * inch, offset, col = "gray80", xpd = NA)
    points(edge + seq(.12, .62, length.out = length(means)) * inch,
      offset + scale * means, col = colours[names(means)], pch = 16, cex = .85, xpd = NA)
    axis(4, at = offset + scale * mean_ticks, labels = mean_ticks,
      pos = edge + .75 * inch, las = 1, cex.axis = .8, xpd = NA)
    mtext("Time mean", side = 4, line = 5.8, cex = .8)
  }
}
# Three quantiles per row; the overall pair spans the full width below them.
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
plot_errors(overall, paste("Overall average across all", n_levels, "quantiles"))
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", labels, col = colours, lty = 1, lwd = 2, ncol = 4, bty = "n")
mtext(paste(disease_name, location_id, "| horizon", horizon_id, "| H =", H),
  side = 3, outer = TRUE, line = 2, font = 2, cex = 1.4)
mtext("Weekly residual (top) and quantile loss (bottom) | Right-hand dots: full-period means",
  side = 3, outer = TRUE, line = .5)
mtext(paste("Target dates |", sub("^Mean ", "", loss_axis_label), "| Dotted lines: all peak dates"),
  side = 1, outer = TRUE, line = .5)
dev.off()

# Also keep the overall comparison in a larger, separate figure.
cairo_pdf(overall_figure_file, width = 12, height = 8, family = "DejaVu Sans")
layout(matrix(1:3, ncol = 1), heights = c(1, 1, .25))
par(oma = c(2, 0, 3, 0), mgp = c(2.5, .7, 0), cex = 1)
plot_errors(overall, paste("Overall average across all", n_levels, "quantiles"))
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", labels, col = colours, lty = 1, lwd = 2, ncol = 3, bty = "n", cex = .85)
mtext(paste(disease_name, location_id, "| horizon", horizon_id, "| H =", H),
  side = 3, outer = TRUE, line = 1, font = 2)
mtext("Target dates | Dotted lines: all peak dates | Right-hand dots: full-period means",
  side = 1, outer = TRUE, line = .5, cex = .85)
dev.off()
