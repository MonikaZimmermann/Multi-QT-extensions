# First exploratory analysis of the fixed-rate MultiQT rate grid.
# Run R/01_data_overview.R once to prepare output/data/peak_aligned.rds.
# Run from the project root: Rscript R/02_multiqt_ties.R
# Rebuild only the base timeline and indicators: add --forecast-paths-only.
.libPaths(c("agent/tmp/multiqt-lib", .libPaths()))  # nanoparquet, as in script 01.
source("R/03_functions.R")

# 1. SETTINGS AND PREPARED DATA -------------------------------------
settings <- multiqt_settings()            # Locations, tracker, S0 warmup and windows (R/03_functions.R).
# The plots go through every disease, horizon and feedback design (designs), the locations
# of multiqt_settings() and, where a plot can pool them, all locations except the US.
diseases <- c("covid", "flu")
# Feedback designs, as in R/04 and R/05: the delay in weeks after the first report, after
# which each forecast is learned once from its report current then, or the weekly re-run on
# the reports known then.
designs <- data.frame(delay = c(0L, 1L, 0L), revisions = c(FALSE, FALSE, TRUE))
design_names <- with(designs, ifelse(revisions, "re-run", paste("delay", delay)))
design_files <- gsub("[ -]", "", design_names)          # File names: "delay0", "rerun".
design_of <- function(rows) ifelse(rows$feedback_revisions, "re-run",
  paste("delay", rows$feedback_delay))
cores <- max(1L, parallel::detectCores() - 2L)
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

# 2. BASE FORECAST PATHS AND THE EARLIER THREE INDICATORS ------------------
# One page per disease; each location has a full count timeline and an indicator
# panel on the same calendar axis. Peaks come directly from script 01. Curves are
# selected by reference weeks -7..8, then trimmed to those target-date windows.
# Use the selected prepared base for both curves and indicators, without filled
# forecasts. This keeps the diagnostic consistent when the seasons design uses
# our own base. The indicators are descriptive; no recalibration method changes.
plot_base_forecast_paths <- function(prepared, settings, figure_dir, diseases) {
  path_weeks <- seq.int(1L - settings$weeks_before_peak, settings$weeks_after_peak)
  path_colours <- colorRampPalette(c("#B7D9EB", "#448FBE", "#073B70"))(length(path_weeks))
  before_fill <- "#EFF5FA"
  after_fill <- "#F5F5F5"
  indicator_colours <- c(active = "#246C9D", inactive = "white", missing = "#D5D9DC")
  # Full latest-count timelines include weeks with no base forecast (e.g. flu summer).
  reports <- as.data.frame(nanoparquet::read_parquet("nhsn_vintages_epidata.parquet"))
  reports <- reports[reports$disease %in% diseases &
    reports$location %in% settings$location_ids, ]
  reports <- reports[order(reports$disease, reports$location,
    reports$target_end_date, reports$version), ]
  reports <- reports[!duplicated(reports[c("disease", "location", "target_end_date")],
    fromLast = TRUE), ]
  pdf_file <- file.path(figure_dir, "multiqt_base_forecast_paths.pdf")
  cairo_pdf(pdf_file, width = 15, height = 20, family = "DejaVu Sans")
  on.exit(dev.off())
  for (disease in diseases) {
    cases <- lapply(settings$location_ids, function(location)
      get_case(prepared, disease, location, 0, settings$feedback_delay))
    range_dates <- range(do.call(c, lapply(cases, function(case)
      case$context_forecasts$target_date)))
    dates <- seq(range_dates[1], range_dates[2], by = "week")
    times <- as.numeric(dates)
    ticks <- seq(as.Date(format(range_dates[1], "%Y-%m-01")), range_dates[2], by = "3 months")
    ticks <- ticks[ticks >= range_dates[1]]
    n_panels <- 2L * length(cases)
    layout(matrix(seq_len(n_panels + 1L), ncol = 1),
      heights = c(rep(c(2.2, .8), length(cases)), .9))
    par(oma = c(4.2, .5, 4.2, .5), family = "DejaVu Sans", mgp = c(3.4, .7, 0))
    for (case in cases) {
      peaks <- sort(unique(case$peaks$peak_date))
      peaks <- peaks[peaks >= range_dates[1] & peaks <= range_dates[2]]
      # Nearest peak, earlier at a tie, as in the evaluation windows of script 01.
      week_from_peak <- function(dates) {
        if (!length(peaks)) return(rep(NA_real_, length(dates)))
        vapply(seq_along(dates), function(i) {
          distance <- as.numeric(dates[i] - peaks) / 7
          distance[which.min(abs(distance))]
        }, numeric(1))
      }
      shade_windows <- function() {
        usr <- par("usr")
        for (p in as.numeric(peaks)) {
          rect(p + 7 * min(path_weeks) - 3.5, usr[3], p + 3.5, usr[4],
            col = before_fill, border = NA)
          rect(p + 3.5, usr[3], p + 7 * max(path_weeks) + 3.5, usr[4],
            col = after_fill, border = NA)
        }
        abline(v = as.numeric(peaks), col = "gray50", lty = 3)
      }
      scale <- choose_count_scale(case, settings$scale_method)$scale
      truth <- reports[reports$disease == disease & reports$location == case$location, ]
      observed <- truth$value[match(dates, truth$target_end_date)] / scale
      observed[!is.finite(observed) | observed < 0] <- NA_real_
      # Retained reference dates and their own cutoffs; lookup rejects filled,
      # late, invalid or absent points. Blue curves and signal inputs agree.
      references <- sort(unique(case$context_forecasts$reference_date))
      reference_weeks <- week_from_peak(references)
      references <- references[reference_weeks %in% path_weeks]
      reference_weeks <- reference_weeks[reference_weeks %in% path_weeks]
      lookup <- forecast_lookup(case$context_forecasts)
      curves <- lapply(seq_along(references), function(i) {
        at <- match(references[i], case$reference_date)
        cutoff <- if (is.na(at)) as.Date(NA) else case$as_of[at]
        target <- references[i] + 7 * (-1:3)
        median <- expm1(lookup(references[i], -1:3, cutoff)) / scale
        median[!week_from_peak(target) %in% path_weeks] <- NA_real_
        list(target = target, median = median, week = reference_weeks[i])
      })
      limits <- range(0, observed, unlist(lapply(curves, `[[`, "median")), finite = TRUE)
      if (diff(limits) == 0) limits <- c(0, 1)
      par(mar = c(.25, 6.3, 1.8, 1), cex = .9)
      plot(times, observed, type = "n", ylim = limits, xlim = range(times),
        xaxt = "n", xaxs = "i", xlab = "", ylab = "Hospitalizations", bty = "n")
      shade_windows()
      abline(h = axTicks(2), col = "gray92")
      for (curve in curves) lines(as.numeric(curve$target), curve$median, type = "o",
        pch = c(21, rep(16, 4)), bg = "white", cex = c(.65, rep(.35, 4)),
        lwd = 1.15, col = path_colours[match(curve$week, path_weeks)])
      lines(times, observed, col = "gray20", lwd = 1.8, lty = 2)
      mtext(paste0(case$location_name, " | ", length(peaks), " detected peaks"),
        side = 3, line = .4, adj = 0, font = 2, cex = 1)

      # Indicators are dated when the information is available (reference week),
      # not backdated to the nowcast target. Unavailable evidence stays missing.
      indicator_case <- case
      indicator_case$reference_date <- dates
      indicator_case$as_of <- case$as_of[match(dates, case$reference_date)]
      indicators <- forecast_path_indicators(indicator_case, threshold = .05)
      par(mar = c(2.1, 6.3, .1, 1), cex = .9)
      plot(range(times), c(.5, 3.5), type = "n", axes = FALSE, xaxs = "i", yaxs = "i",
        xlab = "", ylab = "")
      shade_windows()
      for (j in seq_along(c("path", "nowcast", "revision"))) {
        flag <- indicators[[paste0(c("path", "nowcast", "revision")[j], "_indicator")]]
        fill <- ifelse(is.na(flag), indicator_colours["missing"],
          ifelse(flag, indicator_colours["active"], indicator_colours["inactive"]))
        rect(times - 2.7, 4 - j - .28, times + 2.7, 4 - j + .28,
          col = fill, border = "#B5BEC5", lwd = .35)
      }
      axis(2, at = 3:1, labels = c("Path", "Nowcast", "Revision"), las = 1,
        tick = FALSE, cex.axis = .9)
      axis.Date(1, at = ticks, format = "%b %Y", cex.axis = .85)
    }
    par(mar = c(0, 6.3, 0, 1))
    plot.new()
    legend("topleft", legend = c("Latest reported counts", "Retrospective peak",
      "Weeks -7 to 0", "Weeks 1 to 8", "Indicator = 1", "Indicator = 0", "Unavailable"),
      col = c("gray20", "gray50", "#CCD6DF", "gray80", "#246C9D", "#B5BEC5", "#B5BEC5"),
      lty = c(2, 3, rep(NA, 5)), pch = c(NA, NA, rep(22, 5)),
      pt.bg = c(NA, NA, before_fill, after_fill, indicator_colours),
      lwd = 1.5, pt.cex = 1.3, ncol = 4, bty = "n", cex = .88)
    legend("bottomleft", legend = paste0("Reference week ", c(-7, 0, 8)),
      col = path_colours[match(c(-7, 0, 8), path_weeks)], lty = 1, lwd = 2,
      horiz = TRUE, bty = "n", cex = .88)
    base_names <- unique(vapply(cases, `[[`, "", "base_name"))
    mtext(paste(toupper(disease), "| Base forecast paths and three indicators"),
      side = 3, outer = TRUE, line = 2.6, font = 2, cex = 1.25)
    mtext(paste(paste(base_names, collapse = ", "), "|", count_scale_label(settings$scale_method),
      "| All peaks from script 01 | Indicators: score > log(1.05)"),
      side = 3, outer = TRUE, line = 1.1, cex = .9)
    mtext("Blue curves: one reference week, horizons -1 to 3; darker = later within each peak window. Open point: h = -1.",
      side = 3, outer = TRUE, line = 0, cex = .8)
    mtext("Timeline x-axis: target date. Indicator x-axis: reference date, using its cutoff. Same calendar scale; gray cells mean unavailable, not zero.",
      side = 1, outer = TRUE, line = .5, cex = .78)
    mtext("Blue curves use reference weeks -7 to 8 and are clipped to the same target-date windows. Counts and peak dates are retrospective.",
      side = 1, outer = TRUE, line = 1.6, cex = .78)
    mtext("Path: fall across horizons. Nowcast: fall or slowing growth. Revision: downward same-target revisions that become more negative.",
      side = 1, outer = TRUE, line = 2.7, cex = .78)
  }
  invisible(pdf_file)
}
plot_base_forecast_paths(prepared, settings, figure_dir, diseases)
if ("--forecast-paths-only" %in% commandArgs(trailingOnly = TRUE)) quit(save = "no")


# 3. RUN MULTIQT ON THE RATE GRID, ONE LOCATION AT A TIME ----------------
# 01 supplies full issue histories and the online feedback in raw counts.
# MultiQT tracks log1p(counts) or counts and learns only from the feedback released
# so far: each forecast once, at the first later round feedback_delay weeks after its
# first valid report, from its report current then (feedback_protocol() in
# R/03_functions.R). Re-run design (get_case(revisions = TRUE)): at every issue week t,
# MultiQT is re-run from the first week on the reports current at t's cutoff (case_as_of())
# and issues its week-t forecast; S0 stays that of the stateful design, as in R/04.
# Every run is saved (cached() in R/03_functions.R): a restarted run recomputes only what
# changed.
# Scores use score_truth() (latest vintage by default), in counts / count_scale
# (per 100k residents by default).
# Every location keeps its own offsets, error scale S0 and count scale.
# multipliers: the grid to run (default: the dense plot grid; it must contain 0).
run_rate_grid <- function(case, keep_weekly = FALSE, multipliers = rate_multipliers) {
  scaling <- choose_count_scale(case, settings$scale_method)
  s0 <- error_scale(case, settings$warmup_rounds, settings$scale_floor, settings$transform,
    settings$replay_warmup)
  # If fewer than 20 reports arrive, every multiplier keeps the base.
  eta <- multipliers * s0$scale
  eta[multipliers == 0] <- 0
  # The issued forecasts of one rate: stateful, or re-run at every issue week.
  replays <- if (isTRUE(case$feedback_revised))
    lapply(seq_len(nrow(case$base)), function(t) case_as_of(case, t, s0))
  issue <- function(eta) {
    if (is.null(replays)) return(multiqt(case$base, case$feedback_observed,
      case$feedback_by_week, case$alpha, eta, s0$warmup_weeks, settings$transform))
    t(vapply(seq_along(replays), function(t) {
      replay <- replays[[t]]
      multiqt(replay$base, replay$feedback_observed, replay$feedback_by_week, replay$alpha,
        eta, s0$warmup_weeks, settings$transform)[t, ]
    }, numeric(length(case$alpha))))
  }
  # Score every issued week, including the S0 warmup, where every m issues the base
  # unless the warmup is replayed.
  truth <- score_truth(case, settings$score_target, settings$score_filled)
  scored <- is.finite(truth)
  phase_rows <- lapply(setNames(phases, phases), function(phase) which(scored &
    phase_weeks(case, phase, weeks_before_peak, weeks_after_peak)))
  scores <- scores_by_level <- weekly <- list()
  for (i in seq_along(multipliers)) {
    q <- issue(eta[i])
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
    location_name = case$location_name, horizon = case$horizon,
    feedback_delay = case$feedback_delay, feedback_revisions = case$feedback_revised,
    feedback_protocol = case$feedback_protocol)
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
      score_target = settings$score_target,
      multiplier = multipliers,
      eta, ensemble_expert = multipliers %in% expert_multipliers,
      initialized = s0$initialized, initialization_date = s0$date,
      warmup_weeks = s0$warmup_weeks, initialization_rounds = length(s0$rows),
      consumed_rounds = length(s0$consumed_rows), raw_error_scale = s0$raw_scale,
      error_scale = s0$scale, scale_floor = s0$floor, scale_method = scaling$method,
      scale_used = scaling$used, count_scale = scaling$scale, scale_as_of = case$scale_as_of,
      reference_weeks = scaling$reference_weeks))
}
# The plotted locations on the dense grid, with the weekly forecasts of display_multipliers
# (season plots), for every disease, horizon and design, in parallel.
run_key <- function(disease, location, horizon, k) paste(disease, location, horizon,
  design_files[k], sep = "_")
plotted <- expand.grid(horizon = settings$horizons, location = settings$location_ids,
  disease = diseases, design = seq_len(nrow(designs)), stringsAsFactors = FALSE)
plotted <- plotted[order(!designs$revisions[plotted$design]), ]     # Re-runs first (longest).
runs <- parallel::mclapply(seq_len(nrow(plotted)), function(i) {
  k <- plotted$design[i]
  case <- get_case(prepared, plotted$disease[i], plotted$location[i], plotted$horizon[i],
    designs$delay[k], designs$revisions[k])
  check_peak_window(case, weeks_before_peak, weeks_after_peak,
    score_truth(case, settings$score_target, settings$score_filled))
  message("Rate grid: ", plotted$disease[i], " | ", plotted$location[i], " | horizon ",
    plotted$horizon[i], " | ", design_names[k])
  cached("02_rate_grid", run_rate_grid, case, TRUE)
}, mc.cores = cores, mc.preschedule = FALSE)
failed <- vapply(runs, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
if (any(failed)) stop("Rate grid failed: ", runs[[which(failed)[1]]])
names(runs) <- with(plotted, run_key(disease, location, horizon, design))
scores <- do.call(rbind, lapply(runs, `[[`, "scores"))
scores_by_level <- do.call(rbind, lapply(runs, `[[`, "scores_by_level"))
rownames(scores) <- rownames(scores_by_level) <- NULL
scores$design <- design_of(scores)
scores_by_level$design <- design_of(scores_by_level)
write.csv(do.call(rbind, lapply(runs, `[[`, "rate_grid")),
  file.path(result_dir, "multiqt_rate_grid.csv"), row.names = FALSE)
write.csv(scores, file.path(result_dir, "multiqt_phase_rate_tradeoffs.csv"), row.names = FALSE)
# About a million rows (cases x multipliers x phases x levels): stored as RDS.
saveRDS(scores_by_level, file.path(result_dir, "multiqt_scores_by_level.rds"))
# Every prepared location (not only location_ids), on the multipliers drawn as dots, for
# the pooled row and page of the plots. Only the scores are kept.
dots <- do.call(rbind, lapply(seq_len(nrow(designs)), function(k) do.call(rbind,
  lapply(prepared$cases, function(case) if (case$disease %in% diseases &&
    case$horizon %in% settings$horizons) data.frame(disease = case$disease,
      location = case$location, horizon = case$horizon, design = k)))))
dot_scores <- parallel::mclapply(seq_len(nrow(dots)), function(i) {
  k <- dots$design[i]
  key <- run_key(dots$disease[i], dots$location[i], dots$horizon[i], k)
  if (key %in% names(runs)) return(lapply(runs[[key]][c("scores", "scores_by_level")],
    function(x) x[x$multiplier %in% expert_multipliers, ]))
  cached("02_rate_grid_dots", function(case) run_rate_grid(case,
    multipliers = expert_multipliers)[c("scores", "scores_by_level")], get_case(prepared,
      dots$disease[i], dots$location[i], dots$horizon[i], designs$delay[k], designs$revisions[k]))
}, mc.cores = cores)
failed <- vapply(dot_scores, inherits, logical(1), "try-error")
if (any(failed)) stop(dot_scores[[which(failed)[1]]])
all_scores <- do.call(rbind, lapply(dot_scores, `[[`, "scores"))
all_scores_by_level <- do.call(rbind, lapply(dot_scores, `[[`, "scores_by_level"))
rownames(all_scores) <- rownames(all_scores_by_level) <- NULL
all_scores$design <- design_of(all_scores)
all_scores_by_level$design <- design_of(all_scores_by_level)
write.csv(all_scores, file.path(result_dir, "multiqt_phase_rate_tradeoffs_all_locations.csv"),
  row.names = FALSE)

# 4. PLOT SETTINGS ---------------------------------------------------
score_units <- count_scale_label(settings$scale_method)
tracker_label <- if (settings$transform == "log1p") "Log1p-count tracking" else "Count tracking"
# The settings line of a page for feedback design k.
settings_label <- function(k) paste(tracker_label, "| eta = m * S0, S0 frozen after",
  settings$warmup_rounds, if (settings$replay_warmup)
  "released rows (in hindsight: warmup replayed)" else "released rows", "| Feedback",
  design_names[k], "| Scored against",
  score_truth_label(settings$score_target, designs$delay[k]), "in", score_units)
phase_labels <- c(before = paste0("Before peak (weeks -", weeks_before_peak - 1, " to 0)"),
  after = paste0("After peak (weeks 1 to ", weeks_after_peak, ")"),
  overall = paste0("Overall (", prepared$settings$evaluation_label, ")"))
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

# 5. CALIBRATION ERROR VS MEAN QUANTILE LOSS ON THE MULTIPLIER GRID ----------
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
# All pooled locations (without the US) in one panel: small, light symbols for single locations and the
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
# all_rows: scores of every pooled location (pooled_rows()); if given, they fill a first row.
# k: the feedback design.
plot_tradeoff_page <- function(rows, disease, title, k, all_rows = NULL) {
  page_rows <- n_rows + !is.null(all_rows)
  panels <- page_rows * n_columns
  layout(rbind(matrix(seq_len(panels), page_rows, n_columns, byrow = TRUE), panels + 1L),
    heights = c(rep(1, page_rows), .4))
  par(oma = c(1, 0, 4, 0), mar = c(3.6, 3.6, 2.6, .8), mgp = c(2.2, .6, 0), cex = .8)
  if (!is.null(all_rows)) for (horizon in settings$horizons)
    plot_all_locations_panel(all_rows[all_rows$horizon == horizon, ],
      paste("All", pooled_label(length(unique(all_rows$location)), settings), "| horizon",
        horizon))
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
  mtext(paste(toupper(disease), "|", design_names[k], "|", title), side = 3, outer = TRUE,
    line = 2.2, font = 2, cex = 1.2)
  mtext(paste(settings_label(k), "| Lower left is better."), side = 3, outer = TRUE, line = .7,
    cex = .85)
}
# Calibration error = mean over levels of |coverage - alpha|; loss = mean over levels.
# One page per disease and feedback design. First row: every prepared location except the
# US (all_scores, pooled_rows()).
in_page <- function(rows, disease, k) rows[rows$disease == disease &
  rows$design == design_names[k], ]
cairo_pdf(file.path(figure_dir, "multiqt_error_tradeoff.pdf"), width = 17, height = page_height + 3.1, family = "DejaVu Sans")
for (disease in diseases) for (k in seq_len(nrow(designs)))
  plot_tradeoff_page(in_page(scores, disease, k), disease,
    "Fixed-rate MultiQT: calibration error vs mean quantile loss, mean over quantile levels", k,
    pooled_rows(in_page(all_scores, disease, k), settings))
dev.off()

# The same comparison for each quantile level: one PDF per feedback design in
# multiqt_by_level/, one page per disease and level, first row every location except the US.
# Averaging levels hides which quantiles struggle; CE here is |coverage - alpha| itself.
dir.create(file.path(figure_dir, "multiqt_by_level"), showWarnings = FALSE)
for (k in seq_len(nrow(designs))) {
  cairo_pdf(file.path(figure_dir, "multiqt_by_level", paste0(design_files[k], ".pdf")),
    width = 17, height = page_height + 3.1, family = "DejaVu Sans")
  for (disease in diseases) for (level in sort(unique(scores_by_level$alpha))) {
    at_level <- function(rows) in_page(rows[rows$alpha == level, ], disease, k)
    plot_tradeoff_page(at_level(scores_by_level), disease, paste0("Fixed-rate MultiQT at ",
      "quantile level alpha = ", level, ": |coverage - alpha| vs quantile loss"), k,
      pooled_rows(at_level(all_scores_by_level), settings))
  }
  dev.off()
}

# 6. LOSS AND CALIBRATION ALONG THE GRID: ONE PAGE PER LOCATION ----------
# Per feedback design: first all locations except the US pooled (on the multipliers drawn
# as dots: each phase's loss and calibration error pooled over the locations, each location
# weighted by its weeks), then one page per location of location_ids (the dense grid).
# Loss compares each phase with its own m=0 baseline. Calibration is the mean
# absolute coverage gap over levels; dashed lines show its m=0 value.
pooled_curves <- function(rows) {
  rows <- pooled_rows(rows, settings)
  rows <- rows[is.finite(rows$quantile_loss) & rows$n_weeks > 0, ]
  pooled <- do.call(rbind, lapply(split(rows, rows[c("horizon", "phase", "multiplier")],
    drop = TRUE), function(x) data.frame(x[1, c("horizon", "phase", "multiplier")],
      quantile_loss = weighted.mean(x$quantile_loss, x$n_weeks),
      coverage_error = weighted.mean(x$coverage_error, x$n_weeks))))
  zero <- pooled[pooled$multiplier == 0, ]
  base_loss <- zero$quantile_loss[match(paste(pooled$horizon, pooled$phase),
    paste(zero$horizon, zero$phase))]
  pooled$relative_loss <- 100 * (pooled$quantile_loss / base_loss - 1)
  pooled[order(pooled$multiplier), ]
}
cairo_pdf(file.path(figure_dir, "multiqt_phase_rate_tradeoffs.pdf"), width = 20.3125, height = 15.5, family = "DejaVu Sans")
par(mfrow = c(4, 5), mar = c(4, 4, 3, 1), oma = c(4.5, 0, 4, 0), mgp = c(2.5, .7, 0), cex = .95)
for (k in seq_len(nrow(designs))) for (location in c(NA, settings$location_ids)) {
  pooled <- is.na(location)
  page <- if (pooled) do.call(rbind, lapply(diseases, function(disease)
    data.frame(disease, pooled_curves(in_page(all_scores, disease, k))))) else
    scores[scores$location == location & scores$design == design_names[k], ]
  # Each disease's loss row is followed immediately by its calibration row.
  # Shared calibration limits within a disease make the five horizons comparable.
  for (disease in diseases) for (metric in c("relative_loss", "coverage_error"))
    for (horizon in settings$horizons) {
      z <- page[page$disease == disease & page$horizon == horizon, ]
      loss_panel <- metric == "relative_loss"
      values <- if (loss_panel) z[[metric]][z$multiplier > 0] else
        page$coverage_error[page$disease == disease]
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
          col = phase_colours[phase], lwd = 1.8, type = if (pooled) "o" else "l", pch = 16,
          cex = .6)
        if (!loss_panel) abline(h = curve$coverage_error[curve$multiplier == 0],
          col = phase_colours[phase], lty = 2)
      }
      if (!any(is.finite(z[[metric]][z$multiplier > 0])))
        text(sqrt(prod(range(rate_multipliers[positive]))), mean(limits), "No usable phase scores")
    }
  where <- if (pooled) paste("All", pooled_label(length(unique(pooled_rows(in_page(all_scores,
    diseases[1], k), settings)$location)), settings), "pooled (multipliers of the dots)") else
    location_names[location]
  mtext(paste(where, "|", design_names[k], "| Fixed-rate MultiQT: loss and calibration before",
    "and after the peak"), side = 3, outer = TRUE, line = 2, cex = 1.3)
  mtext(settings_label(k), side = 3, outer = TRUE, line = .5, cex = .9)
  legend(grconvertX(.5, "ndc", "user"), grconvertY(.035, "ndc", "user"), phase_labels,
    col = phase_colours, lwd = 2, bty = "n", horiz = TRUE, xjust = .5, cex = .85, xpd = NA)
  mtext(paste("Loss: % change from m = 0 in the same phase (negative is better). Calibration:",
    "mean |coverage - alpha| over levels (lower is better); dashed: m = 0.", if (pooled)
      "Pooled: each location weighted by its weeks."), side = 1, outer = TRUE, line = 3,
    cex = .85)
}
dev.off()

# 7. INSPECT THE TIME SERIES AROUND THE PEAK -------------------------
plot_alpha <- c(0.15, 0.5, 0.75)
alpha_colours <- c("#0072B2", "#009E73", "#D55E00")
colour_index <- 1 + round(255 * (log10(display_multipliers[-1]) - log10(min(display_multipliers[-1]))) /
  diff(range(log10(display_multipliers[-1]))))
multiplier_colours <- c("orange", hcl.colors(320, "Viridis")[colour_index])
# Both weekly diagnostics use the same quantile styles and peak marker.
# Background of the time plots: the tuning period shaded, dotted lines at the peaks.
mark_time <- function(case) {
  tuning <- case$target_date[case$period == "tuning"]
  usr <- par("usr")
  if (length(tuning)) rect(min(tuning) - 3.5, usr[3], max(tuning) + 3.5, usr[4], col = "gray93",
    border = NA)
  abline(v = as.numeric(case$peaks$peak_date), col = "gray40", lty = 3)
}
plot_season_diagnostic <- function(case, values, limits, y_label, title, ticks,
                                   show_time = FALSE, zero_line = FALSE) {
  par(mar = c(if (show_time) 3.5 else 1, 4.5, 2.5, 1))
  broken <- break_gaps(case$target_date, values)   # Lines break at the summer 2024 gap.
  matplot(broken$dates, broken[[2]], type = "n", ylim = limits, xlim = range(case$target_date),
    xaxt = "n", xaxs = "i", xlab = if (show_time) "Target date" else "", ylab = y_label,
    main = title)
  mark_time(case)
  matplot(broken$dates, broken[[2]], type = "b", col = alpha_colours, lty = 1, pch = c(1, 16, 17),
    cex = 0.5, add = TRUE)
  if (zero_line) abline(h = 0, col = "gray65", lty = 2)
  if (show_time) axis.Date(1, at = ticks, format = "%b\n%Y", cex.axis = 0.8)
}
truth_label <- score_truth_label(settings$score_target, settings$feedback_delay)
truth_legend <- paste0(toupper(substr(truth_label, 1, 1)), substring(truth_label, 2))
# One PDF per disease, horizon and feedback design in multiqt_season_forecasts/, one page
# per location.
draw_season_pages <- function(disease, horizon, k) {
  cairo_pdf(file.path(figure_dir, "multiqt_season_forecasts", paste0(disease, "_horizon",
    horizon, "_", design_files[k], ".pdf")), width = 18, height = 15, family = "DejaVu Sans")
  # Six multipliers per page; each column block holds forecast, coverage and loss panels.
  layout(rbind(cbind(1:3, 4:6, 7:9), cbind(10:12, 13:15, 16:18)), heights = c(2, 1, 1, 2, 1, 1))
  for (location in settings$location_ids) {
    run <- runs[[run_key(disease, location, horizon, k)]]
    case <- run$case
    columns <- match(plot_alpha, case$alpha)
    dates <- case$target_date
    ticks <- seq(as.Date(format(min(dates), "%Y-%m-01")), max(dates), by = "6 months")
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
      mark_time(case)
      # Thick median lines; shading spans the selected lower/upper quantiles. The polygon
      # and lines break at the summer 2024 gap.
      broken <- break_gaps(dates, q, base, observed)
      for (part in split(seq_along(broken$dates), cumsum(is.na(broken$dates)))) {
        part <- part[!is.na(broken$dates[part])]
        polygon(c(broken$dates[part], rev(broken$dates[part])),
          c(broken[[2]][part, 1], rev(broken[[2]][part, 3])),
          col = adjustcolor(multiplier_colours[i], alpha.f = 0.2), border = NA)
      }
      matlines(broken$dates, broken[[3]], col = "gray50", lty = 2, lwd = c(1, 2, 1))
      matlines(broken$dates, broken[[2]], col = multiplier_colours[i], lty = 1, lwd = c(1, 2, 1))
      lines(broken$dates, broken[[4]], type = "o", pch = 16, cex = 0.5)
      legend("topright", c(truth_legend, "Original", "Fixed MultiQT", "Peaks (retrospective)"),
        col = c("black", "gray50", multiplier_colours[i], "gray40"),
        lty = c(1, 2, 1, 3), pch = c(16, NA, NA, NA), cex = 0.7, bty = "n")
      residual <- sweep(run$weekly[[key]]$covered[, columns], 2, plot_alpha, "-")
      plot_season_diagnostic(case, residual, c(-1, 1), "Coverage - alpha",
        "Weekly coverage residual", ticks, zero_line = TRUE)
      plot_season_diagnostic(case, run$weekly[[key]]$loss[, columns], loss_limits,
        "Quantile loss", "Weekly quantile loss", ticks, show_time = TRUE)
    }
    mtext(paste(disease, "|", case$location_name, paste0("(", location, ")"),
      "| horizon", horizon, "|", design_names[k], "|", min(dates), "to", max(dates)),
      side = 3, outer = TRUE, line = 1, cex = 1)
    mtext(paste0(tracker_label, " | eta = m * S0, S0 = ", signif(run$s0$scale, 3),
      " | Scores/plots: ", count_scale_label(run$scaling$used), " | Dotted lines: peaks (local ",
      "maxima of the latest counts) | Gray: tuning period (not scored)"),
      side = 3, outer = TRUE, line = -0.3, cex = 0.8)
    mtext(paste("alpha =", plot_alpha), side = 1, outer = TRUE,
      line = 0.5, at = c(0.2, 0.5, 0.8), col = alpha_colours, cex = 0.85)
    mtext(paste0("Coverage residual = 1(", truth_label, " <= quantile) - alpha."),
      side = 1, outer = TRUE, line = 2, cex = 0.75)
  }
  invisible(dev.off())
}
dir.create(file.path(figure_dir, "multiqt_season_forecasts"), showWarnings = FALSE)
for (k in seq_len(nrow(designs))) for (disease in diseases) for (horizon in settings$horizons)
  draw_season_pages(disease, horizon, k)
# Compare loss/coverage with these curves: do larger multipliers leave forecasts too high,
# too low, or unstable after the peak? This motivates future adaptive rates.
