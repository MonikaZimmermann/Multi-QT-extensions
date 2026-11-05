# First exploratory analysis of the fixed-rate MultiQT rate grid.
# Run R/01_data_overview.R once to prepare output/data/peak_aligned.rds.
# Run from the project root: Rscript R/02_multiqt_ties.R
source("R/03_functions.R")

# 1. SETTINGS AND PREPARED DATA -------------------------------------
settings <- multiqt_settings()            # Locations, tracker, S0 warmup and windows (R/03_functions.R).
season_disease <- "covid"                 # Season forecast plot only: "covid" or "flu".
season_horizon <- 2                       # Season forecast plot only: -1, 0, 1, 2 or 3.
diseases <- c("covid", "flu")
expert_multipliers <- settings$expert_multipliers
multiplier_min <- min(expert_multipliers[expert_multipliers > 0])
multiplier_max <- max(expert_multipliers)
n_positive_multipliers <- 81L
display_multipliers <- c(0, .01, .03, .075, .15, .3)
positive_grid <- 10^seq(log10(multiplier_min), log10(multiplier_max), length.out = n_positive_multipliers)
rate_multipliers <- sort(unique(c(expert_multipliers, signif(positive_grid, 12), display_multipliers)))
# The dense plot grid includes every script 04 multiplier. For each location, eta=m*S0:
# S0 is the frozen Q90 absolute base error in tracker units from the first 20 released rounds.
# The multiplier is shared across locations; actual eta is recorded separately.
phases <- c("before", "after", "overall")
weeks_before_peak <- settings$weeks_before_peak   # Before window also scores the peak week (0).
weeks_after_peak <- settings$weeks_after_peak
prepared <- readRDS("output/data/peak_aligned.rds")
result_dir <- "output/results"
figure_dir <- "output/pdf"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# 2. RUN MULTIQT ON THE RATE GRID, ONE LOCATION AT A TIME ----------------
# 01 supplies full issue histories and the online feedback in raw counts.
# MultiQT tracks log1p(counts) or counts and learns only from the feedback released
# so far: the report available feedback_delay weeks after the first report.
# Scores use score_truth() (latest vintage by default), in counts / count_scale
# (historical maximum by default).
# Every location keeps its own offsets, error scale S0 and count scale.
# multipliers: the grid to run (default: the dense plot grid; it must contain 0).
run_rate_grid <- function(case, keep_weekly = FALSE, multipliers = rate_multipliers) {
  scaling <- choose_count_scale(case, settings$scale_method)
  s0 <- error_scale(case, settings$warmup_rounds, settings$scale_floor, settings$transform,
    settings$replay_warmup)
  # If fewer than 20 reports arrive, every multiplier keeps the base.
  eta <- multipliers * s0$scale
  eta[multipliers == 0] <- 0
  # Score every issued week, including the S0 warmup, where every m issues the base
  # unless the warmup is replayed.
  truth <- score_truth(case, settings$score_target, settings$score_filled)
  scored <- is.finite(truth)
  phase_rows <- lapply(setNames(phases, phases), function(phase) which(scored &
    phase_weeks(case, phase, weeks_before_peak, weeks_after_peak)))
  scores <- scores_by_level <- weekly <- list()
  for (i in seq_along(multipliers)) {
    q <- multiqt(case$base, case$feedback_observed, case$feedback_by_week, case$alpha,
      eta[i], s0$warmup_weeks, settings$transform)
    # Coverage compares counts, so exact ties with integer reports stay ties.
    # The pinball loss is count loss divided by count_scale exactly once.
    week_scores <- score_quantiles(q, truth, case$alpha, scaling$scale)
    for (phase in phases) {
      s <- summarize_scores(week_scores, case$alpha, phase_rows[[phase]])
      key <- paste(i, phase)
      scores[[key]] <- data.frame(phase, multiplier = multipliers[i], eta = eta[i],
        n_weeks = s$n_weeks, quantile_loss = s$quantile_loss, coverage_error = s$coverage_error)
      scores_by_level[[key]] <- data.frame(phase, multiplier = multipliers[i],
        eta = eta[i], n_weeks = s$n_weeks, s$by_level)
    }
    if (keep_weekly && multipliers[i] %in% display_multipliers) {
      weekly[[as.character(multipliers[i])]] <- list(q = q / scaling$scale,
        covered = week_scores$covered, loss = week_scores$loss)
    }
  }
  id <- data.frame(disease = case$disease, location = case$location,
    location_name = case$location_name, horizon = case$horizon)
  scores <- do.call(rbind, scores)
  # Loss change versus the same phase's m=0 baseline (the original forecasts).
  zero <- scores[scores$multiplier == 0, ]
  base_loss <- zero$quantile_loss[match(scores$phase, zero$phase)]
  scores$relative_loss <- ifelse(is.finite(base_loss) & base_loss > 0,
    100 * (scores$quantile_loss / base_loss - 1), NA_real_)
  list(case = case, scaling = scaling, s0 = s0, weekly = weekly,
    scores = data.frame(id, scores, row.names = NULL),
    scores_by_level = data.frame(id, do.call(rbind, scores_by_level), row.names = NULL),
    rate_grid = data.frame(id, transform = settings$transform,
      feedback_delay = case$feedback_delay, score_target = settings$score_target,
      multiplier = multipliers,
      eta, ensemble_expert = multipliers %in% expert_multipliers,
      initialized = s0$initialized, initialization_date = s0$date,
      warmup_weeks = s0$warmup_weeks, initialization_rounds = length(s0$rows),
      consumed_rounds = length(s0$consumed_rows), raw_error_scale = s0$raw_scale,
      error_scale = s0$scale, scale_floor = s0$floor, scale_method = scaling$method,
      scale_used = scaling$used, count_scale = scaling$scale, scale_as_of = case$scale_as_of,
      reference_weeks = scaling$reference_weeks))
}
runs <- list()
for (disease in diseases) for (location in settings$location_ids) for (horizon in settings$horizons) {
  case <- get_case(prepared, disease, location, horizon, settings$feedback_delay)
  check_peak_window(case, weeks_before_peak, weeks_after_peak,
    score_truth(case, settings$score_target, settings$score_filled))
  message("Rate grid: ", disease, " | ", location, " | horizon ", horizon)
  runs[[paste(disease, location, horizon, sep = "_")]] <- run_rate_grid(case,
    keep_weekly = disease == season_disease && horizon == season_horizon)
}
scores <- do.call(rbind, lapply(runs, `[[`, "scores"))
scores_by_level <- do.call(rbind, lapply(runs, `[[`, "scores_by_level"))
rownames(scores) <- rownames(scores_by_level) <- NULL
write.csv(do.call(rbind, lapply(runs, `[[`, "rate_grid")),
  file.path(result_dir, "multiqt_rate_grid.csv"), row.names = FALSE)
write.csv(scores, file.path(result_dir, "multiqt_phase_rate_tradeoffs.csv"), row.names = FALSE)
# About 300,000 rows (cases x multipliers x phases x levels): stored as RDS.
saveRDS(scores_by_level, file.path(result_dir, "multiqt_scores_by_level.rds"))
# Every prepared location (not only location_ids), on the multipliers drawn as dots, for
# the first row of multiqt_error_tradeoff.pdf. Only the scores are kept.
all_scores <- list()
for (prepared_case in prepared$cases) {
  key <- with(prepared_case, paste(disease, location, horizon, sep = "_"))
  if (!prepared_case$disease %in% diseases || !prepared_case$horizon %in% settings$horizons) next
  all_scores[[key]] <- if (key %in% names(runs))
    runs[[key]]$scores[runs[[key]]$scores$multiplier %in% expert_multipliers, ] else
    run_rate_grid(get_case(prepared, prepared_case$disease, prepared_case$location,
      prepared_case$horizon, settings$feedback_delay), multipliers = expert_multipliers)$scores
}
all_scores <- do.call(rbind, all_scores)
rownames(all_scores) <- NULL
write.csv(all_scores, file.path(result_dir, "multiqt_phase_rate_tradeoffs_all_locations.csv"),
  row.names = FALSE)
print(scores[scores$disease == season_disease & scores$horizon == season_horizon &
  scores$phase == "overall" & scores$multiplier %in% display_multipliers,
  c("location", "multiplier", "n_weeks", "quantile_loss", "coverage_error")], row.names = FALSE)

# 3. PLOT SETTINGS ---------------------------------------------------
score_units <- if (settings$scale_method == "none") "counts" else
  paste0("counts / historical ", settings$scale_method)
tracker_label <- if (settings$transform == "log1p") "Log1p-count tracking" else "Count tracking"
settings_label <- paste(tracker_label, "| eta = m * S0, S0 frozen after",
  settings$warmup_rounds, if (settings$replay_warmup) "released rows (warmup replayed)" else
  "released rows", "| Feedback delay", settings$feedback_delay, "week(s) | Scored against",
  score_truth_label(settings$score_target, settings$feedback_delay), "in",
  score_units, "per location")
phase_labels <- c(before = paste0("Before peak (weeks -", weeks_before_peak - 1, " to 0)"),
  after = paste0("After peak (weeks 1 to ", weeks_after_peak, ")"),
  overall = "Overall (full history)")
phase_colours <- c(before = "#D55E00", after = "#0072B2", overall = "gray25")
phase_pch <- c(before = 17, after = 16, overall = 15)   # m = 0: open version of the symbol.
phase_open_pch <- c(before = 2, after = 1, overall = 0)
location_names <- setNames(vapply(settings$location_ids, function(location)
  get_case(prepared, diseases[1], location, settings$horizons[1])$location_name, ""),
  settings$location_ids)
positive <- rate_multipliers > 0
# log(0) is undefined: m = 0 appears as an open symbol or a horizontal baseline.
dot_multipliers <- expert_multipliers[expert_multipliers > 0]
dot_size <- function(m) .45 + 1.35 * (log10(m) - log10(multiplier_min)) /
  (log10(multiplier_max) - log10(multiplier_min))
multiplier_ticks <- sort(unique(c(range(rate_multipliers[positive]), .01, .03, .1, .3)))
multiplier_ticks <- multiplier_ticks[multiplier_ticks >= min(rate_multipliers[positive]) &
  multiplier_ticks <= max(rate_multipliers[positive])]
multiplier_tick_labels <- trimws(formatC(multiplier_ticks, format = "fg", digits = 4))

# 4. CALIBRATION ERROR VS MEAN QUANTILE LOSS ON THE MULTIPLIER GRID ----------
# Each panel is one location and horizon. Dots mark the script 04 expert multipliers
# and grow with m; the axes span only the plotted points.
# This is exploratory, not held-out rate selection. Short-window calibration is descriptive.
plot_tradeoff_panel <- function(rows, title) {
  rows <- rows[rows$multiplier %in% c(0, dot_multipliers), ]
  finite <- is.finite(rows$quantile_loss) & is.finite(rows$coverage_error)
  if (!any(finite)) {
    plot.new()
    title(main = title)
    text(.5, .5, "No scored weeks")
    return(invisible())
  }
  plot(NA, xlim = extendrange(rows$quantile_loss[finite]),
    ylim = extendrange(rows$coverage_error[finite]),
    xlab = "Mean quantile loss", ylab = "Calibration error", main = title)
  support <- vapply(phases, function(phase) as.character(rows$n_weeks[rows$phase == phase][1]), "")
  draw_tradeoff_points(rows)
  mtext(paste("Weeks before / after / overall:", paste(support, collapse = " / ")),
    side = 3, line = .1, cex = .6)
}
# Shades of a phase colour from light (smallest m) to dark (largest m), on the log scale.
m_shade <- function(colour, m) {
  ramp <- colorRampPalette(c(colorRampPalette(c("white", colour))(5)[2], colour,
    colorRampPalette(c(colour, "black"))(5)[3]))(100)
  ramp[1 + round(99 * (log10(m) - log10(multiplier_min)) /
    (log10(multiplier_max) - log10(multiplier_min)))]
}
# The symbols of a panel: dots for the expert multipliers, growing with m, and the open
# m = 0 symbol, in each phase's colour; alpha and size shrink them for single locations.
# shade = TRUE darkens the dots with m (m = 0 in the lightest shade, thin); halo = TRUE
# draws a white rim around every symbol so it stands out from a cloud.
draw_tradeoff_points <- function(rows, alpha = 1, size = 1, shade = FALSE, halo = FALSE) {
  for (phase in phases) {
    z <- rows[rows$phase == phase & is.finite(rows$quantile_loss), ]
    colour <- adjustcolor(if (shade) m_shade(phase_colours[phase], multiplier_min) else
      phase_colours[phase], alpha.f = alpha)
    dots <- z[z$multiplier %in% dot_multipliers, ]
    dot_colour <- if (shade) adjustcolor(m_shade(phase_colours[phase], dots$multiplier),
      alpha.f = alpha) else colour
    zero <- z[z$multiplier == 0, ]
    if (halo) {
      points(dots$quantile_loss, dots$coverage_error, pch = phase_pch[phase], col = "white",
        cex = 1.4 * size * dot_size(dots$multiplier))
      points(zero$quantile_loss, zero$coverage_error, pch = phase_open_pch[phase],
        col = "white", cex = 2 * size, lwd = 5 * size)
    }
    points(dots$quantile_loss, dots$coverage_error, pch = phase_pch[phase], col = dot_colour,
      cex = size * dot_size(dots$multiplier))
    points(zero$quantile_loss, zero$coverage_error, pch = phase_open_pch[phase],
      col = colour, cex = 2 * size, lwd = if (shade) size else 2 * size)
  }
}
# All prepared locations in one panel: small, light symbols for single locations and the
# usual symbols for the unweighted means over locations. Axes: 1st to 99th percentile.
plot_all_locations_panel <- function(rows, title) {
  rows <- rows[rows$multiplier %in% c(0, dot_multipliers) & is.finite(rows$quantile_loss) &
    is.finite(rows$coverage_error), ]
  means <- aggregate(cbind(quantile_loss, coverage_error) ~ phase + multiplier, rows, mean)
  plot(NA, xlim = extendrange(c(quantile(rows$quantile_loss, c(.01, .99)), means$quantile_loss)),
    ylim = extendrange(c(quantile(rows$coverage_error, c(.01, .99)), means$coverage_error)),
    xlab = "Mean quantile loss", ylab = "Calibration error", main = title)
  draw_tradeoff_points(rows, alpha = .5, size = .7, shade = TRUE)
  draw_tradeoff_points(means, halo = TRUE)
  mtext("Small: locations, darker = larger m; large: means", side = 3, line = .1, cex = .6)
}
# One row per location and one column per horizon, legend below.
n_rows <- length(settings$location_ids)
n_columns <- length(settings$horizons)
page_height <- 3.1 * n_rows + 1.5   # 17 inches for 5 locations.
# all_rows: scores of every prepared location; if given, they fill a first row.
plot_tradeoff_page <- function(rows, disease, title, all_rows = NULL) {
  page_rows <- n_rows + !is.null(all_rows)
  panels <- page_rows * n_columns
  layout(rbind(matrix(seq_len(panels), page_rows, n_columns, byrow = TRUE), panels + 1L),
    heights = c(rep(1, page_rows), .4))
  par(oma = c(1, 0, 4, 0), mar = c(3.6, 3.6, 2.6, .8), mgp = c(2.2, .6, 0), cex = .8)
  if (!is.null(all_rows)) for (horizon in settings$horizons)
    plot_all_locations_panel(all_rows[all_rows$horizon == horizon, ],
      paste("All", length(unique(all_rows$location)), "locations | horizon", horizon))
  for (location in settings$location_ids) for (horizon in settings$horizons)
    plot_tradeoff_panel(rows[rows$location == location & rows$horizon == horizon, ],
      paste(location_names[location], "| horizon", horizon))
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", c(phase_labels, "Open symbol: m = 0 (original forecasts)",
    paste("m =", c(.01, .03, .1, .3))),
    pch = c(phase_pch, 1, rep(16, 4)), col = c(phase_colours, "gray25", rep("gray50", 4)),
    pt.cex = c(rep(1.2, 3), 2, dot_size(c(.01, .03, .1, .3))), pt.lwd = c(1, 1, 1, 2, rep(1, 4)),
    ncol = 4, bty = "n", cex = 1)
  mtext(paste(toupper(disease), "|", title), side = 3, outer = TRUE, line = 2.2, font = 2, cex = 1.2)
  mtext(paste(settings_label, "| Lower left is better."), side = 3, outer = TRUE, line = .7, cex = .85)
}
# Calibration error = mean over levels of |coverage - alpha|; loss = mean over levels.
# First row: every prepared location (all_scores).
cairo_pdf(file.path(figure_dir, "multiqt_error_tradeoff.pdf"), width = 17, height = page_height + 3.1, family = "DejaVu Sans")
for (disease in diseases)
  plot_tradeoff_page(scores[scores$disease == disease, ], disease,
    "Fixed-rate MultiQT: calibration error vs mean quantile loss, mean over quantile levels",
    all_scores[all_scores$disease == disease, ])
dev.off()

# The same comparison for each quantile level: one page per disease and level.
# Averaging levels hides which quantiles struggle; CE here is |coverage - alpha| itself.
cairo_pdf(file.path(figure_dir, "multiqt_by_level.pdf"), width = 17, height = page_height, family = "DejaVu Sans")
for (disease in diseases) for (level in sort(unique(scores_by_level$alpha))) {
  rows <- scores_by_level[scores_by_level$disease == disease & scores_by_level$alpha == level, ]
  plot_tradeoff_page(rows, disease, paste0("Fixed-rate MultiQT at quantile level alpha = ", level,
    ": |coverage - alpha| vs quantile loss"))
}
dev.off()

# 5. LOSS AND CALIBRATION ALONG THE GRID: ONE PAGE PER LOCATION ----------
# Loss compares each phase with its own m=0 baseline. Calibration is the mean
# absolute coverage gap over levels; dashed lines show its m=0 value.
cairo_pdf(file.path(figure_dir, "multiqt_phase_rate_tradeoffs.pdf"), width = 20.3125, height = 15.5, family = "DejaVu Sans")
par(mfrow = c(4, 5), mar = c(4, 4, 3, 1), oma = c(4.5, 0, 4, 0), mgp = c(2.5, .7, 0), cex = .95)
for (location in settings$location_ids) {
  # Each disease's loss row is followed immediately by its calibration row.
  # Shared calibration limits within a disease make the five horizons comparable.
  for (disease in diseases) for (metric in c("relative_loss", "coverage_error"))
    for (horizon in settings$horizons) {
      z <- scores[scores$disease == disease & scores$location == location &
        scores$horizon == horizon, ]
      loss_panel <- metric == "relative_loss"
      values <- if (loss_panel) z[[metric]][z$multiplier > 0] else
        scores$coverage_error[scores$disease == disease & scores$location == location]
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
      axis(1, multiplier_ticks, multiplier_tick_labels)
      abline(h = 0, col = "gray40", lty = 3)
      for (phase in phases) {
        curve <- z[z$phase == phase, ]
        lines(curve$multiplier[curve$multiplier > 0], curve[[metric]][curve$multiplier > 0],
          col = phase_colours[phase], lwd = 1.8)
        if (!loss_panel) abline(h = curve$coverage_error[curve$multiplier == 0],
          col = phase_colours[phase], lty = 2)
      }
      if (!any(is.finite(z[[metric]][z$multiplier > 0])))
        text(sqrt(prod(range(rate_multipliers[positive]))), mean(limits), "No usable phase scores")
    }
  mtext(paste(location_names[location], "| Fixed-rate MultiQT: loss and calibration before and after the peak"),
    side = 3, outer = TRUE, line = 2, cex = 1.3)
  mtext(settings_label, side = 3, outer = TRUE, line = .5, cex = .9)
  legend(grconvertX(.5, "ndc", "user"), grconvertY(.035, "ndc", "user"), phase_labels,
    col = phase_colours, lwd = 2, bty = "n", horiz = TRUE, xjust = .5, cex = .85, xpd = NA)
  mtext("Loss: % change from m = 0 in the same phase (negative is better). Calibration: mean |coverage - alpha| over levels (lower is better); dashed: m = 0.",
    side = 1, outer = TRUE, line = 3, cex = .85)
}
dev.off()

# 6. INSPECT THE TIME SERIES AROUND THE PEAK -------------------------
plot_alpha <- c(0.15, 0.5, 0.75)
alpha_colours <- c("#0072B2", "#009E73", "#D55E00")
colour_index <- 1 + round(255 * (log10(display_multipliers[-1]) - log10(min(display_multipliers[-1]))) /
  diff(range(log10(display_multipliers[-1]))))
multiplier_colours <- c("orange", hcl.colors(320, "Viridis")[colour_index])
# Both weekly diagnostics use the same quantile styles and peak marker.
plot_season_diagnostic <- function(dates, values, limits, y_label, title, ticks,
                                   peak_date, show_time = FALSE, zero_line = FALSE) {
  par(mar = c(if (show_time) 3.5 else 1, 4.5, 2.5, 1))
  matplot(dates, values, type = "b", col = alpha_colours, lty = 1, pch = c(1, 16, 17),
    cex = 0.5, ylim = limits, xlim = range(dates), xaxt = "n", xaxs = "i",
    xlab = if (show_time) "Target date" else "", ylab = y_label, main = title)
  if (zero_line) abline(h = 0, col = "gray65", lty = 2)
  if (show_time) axis.Date(1, at = ticks, format = "%b\n%Y", cex.axis = 0.8)
  abline(v = as.numeric(peak_date), col = "gray40", lty = 3)
}
truth_label <- score_truth_label(settings$score_target, settings$feedback_delay)
truth_legend <- paste0(toupper(substr(truth_label, 1, 1)), substring(truth_label, 2))
cairo_pdf(file.path(figure_dir, "multiqt_season_forecasts.pdf"), width = 18, height = 15, family = "DejaVu Sans")
# Six multipliers per page; each column block holds forecast, coverage and loss panels.
layout(rbind(cbind(1:3, 4:6, 7:9), cbind(10:12, 13:15, 16:18)), heights = c(2, 1, 1, 2, 1, 1))
for (location in settings$location_ids) {
  run <- runs[[paste(season_disease, location, season_horizon, sep = "_")]]
  case <- run$case
  columns <- match(plot_alpha, case$alpha)
  dates <- case$target_date
  ticks <- seq(as.Date(format(min(dates), "%Y-%m-01")), max(dates), by = "3 months")
  observed <- score_truth(case, settings$score_target, settings$score_filled) / run$scaling$scale
  base <- case$base[, columns] / run$scaling$scale
  forecasts <- lapply(run$weekly, function(x) x$q[, columns, drop = FALSE])
  plot_limits <- range(0, observed, base, unlist(forecasts), na.rm = TRUE)
  loss_limits <- range(0, unlist(lapply(run$weekly, function(x) x$loss[, columns])), na.rm = TRUE)
  par(oma = c(5, 0, 3, 0), mgp = c(2.5, 0.7, 0), cex = 0.9)
  for (i in seq_along(display_multipliers)) {
    key <- as.character(display_multipliers[i])
    q <- forecasts[[key]]
    par(mar = c(1, 4.5, 2.5, 1))
    plot(dates, observed, type = "n", ylim = plot_limits, xlim = range(dates), xaxt = "n",
      xaxs = "i", xlab = "", ylab = paste0("Hospitalizations (", score_units, ")"),
      main = paste("Rate multiplier m =", key))
    polygon(c(dates, rev(dates)), c(q[, 1], rev(q[, 3])),
      col = adjustcolor(multiplier_colours[i], alpha.f = 0.2), border = NA)
    # Thick median lines; shading spans the selected lower/upper quantiles.
    matlines(dates, base, col = "gray50", lty = 2, lwd = c(1, 2, 1))
    matlines(dates, q, col = multiplier_colours[i], lty = 1, lwd = c(1, 2, 1))
    lines(dates, observed, type = "o", pch = 16, cex = 0.5)
    abline(v = as.numeric(case$peak_date), lty = 3, col = "gray40")
    legend("topright", c(truth_legend, "Original", "Fixed MultiQT", "Peak (retrospective)"),
      col = c("black", "gray50", multiplier_colours[i], "gray40"),
      lty = c(1, 2, 1, 3), pch = c(16, NA, NA, NA), cex = 0.7, bty = "n")
    residual <- sweep(run$weekly[[key]]$covered[, columns], 2, plot_alpha, "-")
    plot_season_diagnostic(dates, residual, c(-1, 1), "Coverage - alpha",
      "Weekly coverage residual", ticks, case$peak_date, zero_line = TRUE)
    plot_season_diagnostic(dates, run$weekly[[key]]$loss[, columns], loss_limits,
      "Quantile loss", "Weekly quantile loss", ticks, case$peak_date, show_time = TRUE)
  }
  mtext(paste(season_disease, "|", case$location_name, paste0("(", location, ")"),
    "| horizon", season_horizon, "|", min(dates), "to", max(dates)),
    side = 3, outer = TRUE, line = 1, cex = 1)
  mtext(paste0(tracker_label, " | eta = m * S0, S0 = ", signif(run$s0$scale, 3),
    " | Scores/plots: counts / historical ", run$scaling$used, " (",
    format(run$scaling$scale, trim = TRUE), ") | Vertical line: largest count in the peak window."),
    side = 3, outer = TRUE, line = -0.3, cex = 0.8)
  mtext(paste("alpha =", plot_alpha), side = 1, outer = TRUE,
    line = 0.5, at = c(0.2, 0.5, 0.8), col = alpha_colours, cex = 0.85)
  mtext(paste0("Coverage residual = 1(", truth_label, " <= quantile) - alpha."),
    side = 1, outer = TRUE, line = 2, cex = 0.75)
}
dev.off()
# Compare loss/coverage with these curves: do larger multipliers leave forecasts too high,
# too low, or unstable after the peak? This motivates future adaptive rates.
