# Apply the recalibration methods of R/03_functions.R to each location and plot them.
# Run R/01_data_overview.R once to prepare output/data/peak_aligned.rds.
# Run from the project root: Rscript R/04_methods.R
# Offsets and ensemble weights are fitted independently for each location.
source("R/03_functions.R")

# 1. SETTINGS AND PREPARED DATA -------------------------------------
disease_name <- "covid"                    # "covid" or "flu"
horizon_id <- 0                            # -1, 0, 1, 2 or 3
settings <- multiqt_settings()             # Locations, tracker, experts, H and windows (R/03_functions.R).
transform <- settings$transform
H <- settings$H
# Other rates and caps below are in tracker units (log1p counts or raw counts).
outer_rate <- .005
slope_rate_min <- .005
slope_rate_max <- .1                       # Tracker-unit cap; independent of the multiplier grid.
# Weekly restart comparison: each copy votes from birth, initially forecasting the base.
restart_multiplier <- .05
restart_multipliers <- c(.01, .05, .2)     # Smaller grid at each new start.
restart_lambda <- 2.72                     # Fixed meta-rate on scaled-count pinball loss.
plot_levels <- (1:9) / 10
phases <- c("overall", "before", "after")
prepared <- readRDS("output/data/peak_aligned.rds")
result_dir <- "output/results/peak_reasoning"
figure_dir <- "output/pdf"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write_result <- function(rows, filename)
  write.csv(rows, file.path(result_dir, filename), row.names = FALSE)

# 2. RUN EACH LOCATION INDEPENDENTLY --------------------------------
# Peaks are the largest latest-reported count in the prepared winter window.
# Only the target oracle uses that hindsight date to choose forecasts.
# Alignment changes evaluation coordinates, never when feedback becomes available.
# Every location starts fresh: offsets, expert weights, lambda and loss history.
run_methods <- function(case) {
  scaling <- choose_count_scale(case, settings$scale_method)
  count_scale <- scaling$scale
  s0 <- error_scale(case, settings$warmup_rounds, settings$scale_floor, transform)
  # Keep stored predictions and the scoring/weighting inputs in counts / count_scale.
  base <- case$base / count_scale
  feedback <- case$feedback_observed / count_scale
  adjust <- function(counts, eta, warmup_weeks) multiqt(counts, case$feedback_observed,
    case$feedback_by_week, case$alpha, eta, warmup_weeks, transform) / count_scale

  # 3. ONE TRACKER: FIXED RATE OR A RATE BASED ON ABSOLUTE SLOPE --------
  # Baseline eta = reference_multiplier * S0, using the same frozen scale as the experts.
  predictions <- list(original = base)
  predictions$fixed <- adjust(case$base, settings$reference_multiplier * s0$scale, s0$warmup_weeks)
  # Use only the current base forecast's weekly slope in tracker units, never future counts.
  # This rate needs no S0, so it learns from the first released report on.
  median_tracker <- to_tracker(case$base[, match(.5, case$alpha)], transform)
  slope <- c(0, diff(median_tracker) / (as.numeric(diff(case$reference_date)) / 7))
  slope_rates <- pmin(slope_rate_max, slope_rate_min + abs(slope))
  predictions$absolute_slope <- adjust(case$base, slope_rates, 0L)
  # Larger rates can unwind offsets faster, but can also amplify stale feedback.

  # 4. WHAT IF WE KNEW THE TRUE PEAK? ----------------------------------
  # Hindsight benchmark: bypass correction for postpeak targets, even if issued prepeak.
  # Setting eta=0 instead would only freeze the accumulated offset.
  predictions$oracle_target <- predictions$fixed
  after_peak <- case$target_date > case$peak_date
  predictions$oracle_target[after_peak, ] <- base[after_peak, ]

  # 5. SEVERAL TRACKERS: WEIGHT THEIR FORECASTS ------------------------
  # Rate grid: one persistent expert per multiplier. Restarts: a fresh copy per
  # multiplier at every issue week after warmup; no permanent base expert.
  rate_experts <- make_rate_experts(case, settings$expert_multipliers, s0, count_scale, transform)
  restart_experts <- make_restart_experts(case, restart_multipliers, s0, count_scale, transform)
  fixed_restarts <- select_experts(restart_experts,
    restart_experts$experts$multiplier == restart_multiplier)
  # Each weighting learns from released count-scale pinball losses of issued forecasts.
  # Their deterministic average still uses LEARNED weights; equal averaging is a control.
  weigh <- function(experts, ...) exp_weights(experts, base, feedback, case$feedback_by_week,
    case$alpha, H = H, lambda_start = settings$lambda_start, warmup_weeks = s0$warmup_weeks, ...)
  mixtures <- list(
    uniform = weigh(rate_experts, lambda = 0),
    weighted = weigh(rate_experts, sharing = FALSE),
    weighted_common = weigh(rate_experts, sharing = FALSE, per_quantile = FALSE),
    shared = weigh(rate_experts),
    shared_common = weigh(rate_experts, per_quantile = FALSE),
    restart_uniform = weigh(fixed_restarts, lambda = 0),
    restart_weighted = weigh(fixed_restarts, lambda = restart_lambda),
    restart_grid = weigh(restart_experts, lambda = restart_lambda))
  for (method in names(mixtures)) predictions[[method]] <- mixtures[[method]]$q
  # The first-start control has exactly the same birth rule as the restart copies.
  predictions$restart_persistent <- base
  if (length(fixed_restarts$birth_week)) {
    issued <- seq.int(fixed_restarts$birth_week[1], nrow(base))
    predictions$restart_persistent[issued, ] <- fixed_restarts$forecasts[issued, 1, ]
  }
  # Each outer tracker has its own offsets and learns from its OWN final forecasts.
  # Its state stays continuous when new restart copies join the underlying mixture.
  for (method in c("shared", "shared_common", "restart_uniform", "restart_weighted", "restart_grid"))
    predictions[[paste0(method, "_outer")]] <- adjust(predictions[[method]] * count_scale,
      outer_rate, s0$warmup_weeks)

  # 6. SCORE AGAINST THE LATEST VINTAGE ---------------------------------
  # The online updates above saw only first-use feedback.
  # Counts are already scaled: do not divide quantile loss a second time.
  truth <- case$truth_latest / count_scale
  windows <- lapply(setNames(phases, phases), function(phase) which(is.finite(truth) &
    phase_weeks(case, phase, settings$weeks_before_peak, settings$weeks_after_peak)))
  weekly <- lapply(predictions, score_quantiles, truth = truth, alpha = case$alpha)
  scores <- scores_by_level <- list()
  for (method in names(predictions)) for (phase in phases) {
    s <- summarize_scores(weekly[[method]], case$alpha, windows[[phase]])
    scores[[paste(method, phase)]] <- data.frame(method, phase, n_weeks = s$n_weeks,
      quantile_loss = s$quantile_loss, coverage_error = s$coverage_error)
    scores_by_level[[paste(method, phase)]] <- data.frame(method, phase, n_weeks = s$n_weeks,
      s$by_level)
  }
  id <- data.frame(disease = case$disease, location = case$location,
    location_name = case$location_name, horizon = case$horizon)
  list(case = case, scaling = scaling, s0 = s0, truth = truth, predictions = predictions,
    weights = lapply(mixtures, `[[`, "weights"), lambda = lapply(mixtures, `[[`, "lambda"),
    rate_experts = rate_experts, restart_experts = restart_experts, weekly = weekly,
    scores = data.frame(id, do.call(rbind, scores), row.names = NULL),
    scores_by_level = data.frame(id, do.call(rbind, scores_by_level), row.names = NULL))
}
runs <- list()
for (location in settings$location_ids) {
  case <- get_case(prepared, disease_name, location, horizon_id)
  check_peak_window(case, settings$weeks_before_peak, settings$weeks_after_peak)
  message("Methods: ", disease_name, " | ", location, " | horizon ", horizon_id)
  runs[[location]] <- run_methods(case)
}
saveRDS(runs, file.path(result_dir, "issued_forecasts.rds"))
summary_scores <- do.call(rbind, lapply(runs, `[[`, "scores"))
rownames(summary_scores) <- NULL
write_result(summary_scores, "summary_scores.csv")
write_result(do.call(rbind, lapply(runs, `[[`, "scores_by_level")), "summary_scores_by_level.csv")
write_result(do.call(rbind, lapply(runs, function(run) data.frame(location = run$case$location,
  scale_method = run$scaling$method, scale_used = run$scaling$used, count_scale = run$scaling$scale,
  transform, scale_as_of = run$case$scale_as_of, reference_weeks = run$scaling$reference_weeks,
  initialization_date = run$s0$date, warmup_weeks = run$s0$warmup_weeks,
  raw_error_scale = run$s0$raw_scale, error_scale = run$s0$scale))), "scales.csv")
write_result(do.call(rbind, lapply(runs, function(run)
  data.frame(location = run$case$location, run$rate_experts$experts))), "expert_grid.csv")
write_result(do.call(rbind, lapply(runs, function(run)
  data.frame(location = rep(run$case$location, nrow(run$restart_experts$experts)),
    run$restart_experts$experts, birth_week = run$restart_experts$birth_week))),
  "restart_experts.csv")

# 7. PLOT SETTINGS ---------------------------------------------------
labels <- c(original = "Original",
  fixed = paste0("Fixed MultiQT (eta=", settings$reference_multiplier, "*S0)"),
  absolute_slope = if (transform == "log1p") "Log-slope rate" else "Count-slope rate",
  oracle_target = "Target oracle (hindsight)",
  uniform = "Rate grid: equal",
  weighted = "Rate grid: learned, no sharing",
  weighted_common = "Rate grid: learned, no sharing, common weights",
  shared = paste0("Rate grid: learned + sharing (H=", H, ")"),
  shared_common = "Rate grid: learned + sharing, common weights",
  shared_outer = paste0("Rate grid: sharing + outer (eta=", outer_rate, ")"),
  shared_common_outer = "Rate grid: sharing, common weights + outer",
  restart_persistent = paste0("First start only (m=", restart_multiplier, ")"),
  restart_uniform = paste0("Weekly starts: equal (m=", restart_multiplier, ")"),
  restart_weighted = paste0("Weekly starts: learned (m=", restart_multiplier, ")"),
  restart_grid = "Starts x rate grid: learned",
  restart_uniform_outer = "Starts equal + outer",
  restart_weighted_outer = "Starts learned + outer",
  restart_grid_outer = "Starts x rate grid + outer")
colours <- c(original = "gray40", fixed = "#0072B2", absolute_slope = "#009E73",
  oracle_target = "#CC79A7", uniform = "#E69F00", weighted = "#56B4E9",
  weighted_common = "#56B4E9", shared = "#882255", shared_common = "#882255",
  shared_outer = "#D55E00", shared_common_outer = "#D55E00",
  restart_persistent = "#0072B2", restart_uniform = "#E69F00", restart_weighted = "#009E73",
  restart_grid = "#882255", restart_uniform_outer = "#56B4E9",
  restart_weighted_outer = "#D55E00", restart_grid_outer = "#CC79A7")
# Common weights across quantile levels: dashed lines and triangles.
common <- setNames(grepl("_common", names(labels)), names(labels))
line_types <- ifelse(common, 2, 1)
point_shapes <- ifelse(common, 17, 16)
method_groups <- list(
  "Single trackers" = c("original", "fixed", "absolute_slope", "oracle_target"),
  "Rate-grid experts" = c("original", "fixed", "uniform", "weighted", "weighted_common",
    "shared", "shared_common", "shared_outer", "shared_common_outer"),
  "Weekly starts: immediate admission" = c("original", "restart_persistent", "restart_uniform",
    "restart_weighted", "restart_grid", "restart_uniform_outer", "restart_weighted_outer",
    "restart_grid_outer"))
score_units <- if (settings$scale_method == "none") "counts" else
  paste0("counts / ", settings$scale_method)
loss_axis_label <- paste0("Quantile loss (", score_units, ")")
mean_loss_axis_label <- paste("Mean", tolower(loss_axis_label))
window_titles <- c(overall = "Full observed history",
  before = paste0("Up to winter peak (weeks -", settings$weeks_before_peak, " to 0)"),
  after = paste0("After winter peak (weeks 1 to ", settings$weeks_after_peak, ")"))
plot_title <- function(run) paste(disease_name, "| horizon", horizon_id, "|",
  run$case$location_name, paste0("(", run$case$location, ")"), "| full history |",
  if (transform == "log1p") "log1p(counts) updates" else "raw-count updates")
date_ticks <- function(dates) pretty(as.Date(range(dates), origin = "1970-01-01"), n = 4)

# 8. FORECASTS AND WEEKLY QUANTILE LOSS PER LEVEL ----------------------
# Each level panel shows the latest reported counts and every method's issued
# alpha-quantile; the panel below shows that level's weekly pinball loss.
blocks <- matrix(0L, 6, 3)
for (b in seq_along(plot_levels)) {
  r <- (b - 1L) %/% 3L
  blocks[2L * r + 1:2, (b - 1L) %% 3L + 1L] <- c(2L * b - 1L, 2L * b)
}
cairo_pdf(file.path(figure_dir, "multiqt_method_forecasts.pdf"), width = 18, height = 16, family = "DejaVu Sans")
for (run in runs) for (group_name in names(method_groups)) {
  method_names <- method_groups[[group_name]]
  layout(rbind(blocks, 2L * length(plot_levels) + 1L), heights = c(rep(c(2, 1), 3), .45))
  par(oma = c(1.5, 0, 4, 0), mgp = c(2.4, .6, 0), cex = .85)
  dates <- run$case$target_date
  ticks <- date_ticks(dates)
  for (level in plot_levels) {
    j <- match(level, run$case$alpha)
    values <- vapply(method_names, function(method) run$predictions[[method]][, j], numeric(length(dates)))
    par(mar = c(1, 4.6, 2, 1))
    plot(dates, run$truth, type = "n", ylim = range(0, run$truth, values, na.rm = TRUE),
      xaxt = "n", xlab = "", ylab = paste0("Hospitalizations\n(", score_units, ")"),
      main = paste("alpha =", level))
    abline(v = as.numeric(run$case$peak_date), col = "gray60", lty = 3)
    for (method in method_names)
      lines(dates, values[, method], col = colours[method], lty = line_types[method], lwd = 1.3)
    lines(dates, run$truth, col = "black", lwd = .7)
    points(dates, run$truth, pch = 16, cex = .45)
    losses <- vapply(method_names, function(method) run$weekly[[method]]$loss[, j], numeric(length(dates)))
    par(mar = c(2.6, 4.6, .7, 1))
    plot(dates, losses[, 1], type = "n", ylim = c(0, max(1e-8, losses, na.rm = TRUE)),
      xaxt = "n", xlab = "", ylab = "Quantile\nloss")
    axis.Date(1, at = ticks, format = "%b %Y", cex.axis = .85)
    abline(v = as.numeric(run$case$peak_date), col = "gray60", lty = 3)
    for (method in method_names)
      lines(dates, losses[, method], col = colours[method], lty = line_types[method], lwd = 1.1)
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", c("Latest reported counts", labels[method_names]),
    col = c("black", colours[method_names]), lty = c(NA, line_types[method_names]),
    pch = c(16, rep(NA, length(method_names))), lwd = 2, ncol = 4, bty = "n")
  mtext(paste(plot_title(run), "|", group_name), side = 3, outer = TRUE, line = 2, font = 2, cex = 1.15)
  mtext(paste("Issued alpha-quantile forecasts and their weekly pinball loss in", score_units,
    "against the latest reported counts | Dotted line: winter peak"), side = 3, outer = TRUE, line = .5)
}
dev.off()

# 9. CALIBRATION ERROR AND QUANTILE LOSS ON MATCHED WINDOWS ----------
# CE = mean over levels of |coverage - alpha| over the window's weeks. Axes are
# shared by all pages of a location; up to winter peak includes the peak week.
cairo_pdf(file.path(figure_dir, "multiqt_peak_reasoning_tradeoff.pdf"), width = 14, height = 5.5, family = "DejaVu Sans")
for (run in runs) {
  xlim <- c(0, 1.1 * max(1e-8, run$scores$quantile_loss, na.rm = TRUE))
  ylim <- c(0, 1.1 * max(.01, run$scores$coverage_error, na.rm = TRUE))
  for (group_name in names(method_groups)) {
    method_names <- method_groups[[group_name]]
    layout(rbind(1:3, c(4, 4, 4)), heights = c(1, .22))
    par(oma = c(3, 0, 3, 0), mgp = c(2.6, .7, 0))
    for (phase in phases) {
      rows <- run$scores[run$scores$phase == phase & run$scores$method %in% method_names, ]
      par(mar = c(4, 4, 4, 1))
      plot(NA, xlim = xlim, ylim = ylim, xlab = mean_loss_axis_label, ylab = "Calibration error",
        main = "")
      mtext(window_titles[phase], side = 3, line = 1.8, font = 2, cex = .95)
      grid(col = "gray90")
      mtext(paste(rows$n_weeks[1], "weeks"), side = 3, line = .25, cex = .8)
      points(rows$quantile_loss, rows$coverage_error, col = colours[rows$method],
        pch = point_shapes[rows$method], cex = 1.3)
      if (!any(is.finite(rows$coverage_error)))
        text(mean(xlim), mean(ylim), "No observations inside winter", cex = .8)
    }
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", labels[method_names], col = colours[method_names],
      pch = point_shapes[method_names], ncol = 3, bty = "n", cex = .8)
    mtext(paste(plot_title(run), "|", group_name), side = 3, outer = TRUE, line = 1, font = 2, cex = .9)
    mtext("Calibration error (CE): mean absolute coverage residual after averaging weeks within each quantile level.",
      side = 1, outer = TRUE, line = .5, cex = .8)
    mtext("Triangles: one weight vector shared by all quantile levels. Up to winter peak includes the peak week.",
      side = 1, outer = TRUE, line = 1.8, cex = .8)
  }
}
dev.off()

# 10. ISSUED RATE-GRID WEIGHTS OVER TIME -----------------------------
# These are the weights used to mix forecasts before crossing repair or an outer
# correction. The outer methods keep their underlying ensemble's weights.
# Rows 0.1-0.9: one weight vector per level; last row: common weights for all levels.
weight_columns <- data.frame(method = c("uniform", "weighted", "shared"),
  common_method = c("uniform", "weighted_common", "shared_common"),
  title = c("Equal weights", "Learned, no sharing", "Learned + sharing"))
rate_colours <- setNames(c("gray45", hcl.colors(length(settings$expert_multipliers) - 1L, "Viridis")),
  as.character(settings$expert_multipliers))
plot_weight_distribution <- function(run, weights, title, row_label, first_column, last_row) {
  times <- as.numeric(run$case$target_date)
  par(mar = c(if (last_row) 2.4 else .3, if (first_column) 3.7 else 2.2, 1.5, .6))
  plot(NA, xlim = range(times) + c(-.5, .5) * 7, ylim = c(0, 1),
    xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "",
    ylab = if (first_column) row_label else "", main = title)
  bottom <- numeric(length(times))
  for (j in seq_len(ncol(weights))) {
    top <- bottom + weights[, j]
    rect(times - .45 * 7, bottom, times + .45 * 7, top, col = rate_colours[j], border = NA)
    bottom <- top
  }
  abline(h = .5, col = adjustcolor("white", .6), lty = 3)
  abline(v = as.numeric(run$case$peak_date), col = "gray30", lty = 3)
  axis(2, at = c(0, .5, 1), las = 1, cex.axis = .75)
  if (last_row) axis.Date(1, at = date_ticks(run$case$target_date), format = "%b %Y", cex.axis = .8)
  box(col = "gray70")
}
cairo_pdf(file.path(figure_dir, "multiqt_weight_timelines.pdf"), width = 18,
  height = 1.45 * (length(plot_levels) + 1) + 3.7, family = "DejaVu Sans")
for (run in runs) {
  stopifnot(identical(dimnames(run$weights$shared)[[2]], paste0("m=", names(rate_colours))))
  n_rows <- length(plot_levels) + 1L
  panel_ids <- matrix(seq_len(n_rows * nrow(weight_columns)), ncol = nrow(weight_columns), byrow = TRUE)
  layout(rbind(panel_ids, max(panel_ids) + 1L), heights = c(rep(1, n_rows), .8))
  par(oma = c(4.8, .8, 5.2, .8), mgp = c(2, .4, 0), cex = .9)
  for (row in seq_len(n_rows)) for (j in seq_len(nrow(weight_columns))) {
    level_row <- row <= length(plot_levels)
    weights <- if (level_row) run$weights[[weight_columns$method[j]]][, , match(plot_levels[row], run$case$alpha)] else
      run$weights[[weight_columns$common_method[j]]][, , 1]
    plot_weight_distribution(run, weights, if (row == 1L) weight_columns$title[j] else "",
      if (level_row) paste("alpha =", plot_levels[row]) else "Common\n(all levels)",
      first_column = j == 1L, last_row = row == n_rows)
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", paste("m =", names(rate_colours)), fill = rate_colours,
    border = NA, ncol = 6, bty = "n", cex = .95)
  mtext("Frozen rate-grid weights", side = 3, outer = TRUE, line = 2.4, font = 2, cex = 1.3)
  mtext(plot_title(run), side = 3, outer = TRUE, line = .9, cex = .9)
  mtext("Target date | Weight at issuance for each target | Dotted line: winter peak",
    side = 1, outer = TRUE, line = .4, cex = .85)
  mtext("Rows alpha = 0.1-0.9: one weight vector per quantile level. Last row: one weight vector for all levels, learned from the mean pinball loss.",
    side = 1, outer = TRUE, line = 1.6, cex = .8)
}
dev.off()
