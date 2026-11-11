# Apply the recalibration methods of R/03_functions.R to each location and plot them.
# Run R/01_data_overview.R once to prepare output/data/peak_aligned.rds.
# Run from the project root: Rscript R/04_methods.R
# Offsets and ensemble weights are fitted independently for each location. The forecast
# and weight plots go through every disease, horizon and feedback design (one
# PDF per combination in a folder per plot), their pages through the locations of
# multiqt_settings(); every location, horizon and feedback design of both diseases is scored
# for the trade-off pages and R/06_evaluation_tree.R (about 1.5 hours on 18 cores).
# Every location's run is saved (cached() in R/03_functions.R): a restarted run, after a
# crash or a change, recomputes only the runs whose settings, data or code changed; changes
# to the plots recompute nothing.
source("R/03_functions.R")

# 1. SETTINGS AND PREPARED DATA -------------------------------------
# The combination whose runs are saved in issued_forecasts.rds (every one is plotted).
disease_name <- "covid"                    # "covid" or "flu"
horizon_id <- 2                            # -1, 0, 1, 2 or 3
settings <- multiqt_settings()             # Locations, tracker, experts and windows (R/03_functions.R).
transform <- settings$transform
# Feedback designs of the plots and trade-off pages: the delay in weeks after the first
# report, and the weekly re-run on the reports known then (get_case() in R/03_functions.R).
tradeoff_feedback <- data.frame(delay = c(0L, 1L, 0L), revisions = c(FALSE, FALSE, TRUE))
# Diseases plotted and scored for every location, horizon and feedback design.
score_diseases <- c("covid", "flu")
cores <- max(1L, parallel::detectCores() - 2L)
plot_levels <- (1:9) / 10
phases <- c("overall", "before", "after")
prepared <- readRDS("output/data/peak_aligned.rds")
# The methods' settings: method_parameters() in R/03_functions.R, with the values tuned
# on the tuning period by R/05_tuning.R for this disease, tracker, design and feedback.
tuned_for <- function(delay, revisions, disease = disease_name) tuned_parameters(disease,
  transform, design = prepared$settings$design, feedback_delay = delay, revisions = revisions)
seasons <- c("all", prepared$settings$evaluation_seasons)   # Scored: each season and all.
result_dir <- "output/results/peak_reasoning"
figure_dir <- "output/pdf"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write_result <- function(rows, filename)
  write.csv(rows, file.path(result_dir, filename), row.names = FALSE)

# 2. RUN EACH LOCATION INDEPENDENTLY --------------------------------
# apply_methods() (R/03_functions.R) issues every method's forecasts from the online
# feedback only. Peaks are local maxima of the latest-reported counts (R/01);
# only the target oracle uses that hindsight date to choose forecasts. Alignment
# changes evaluation coordinates, never when feedback becomes available. Every
# location starts fresh: offsets, expert weights, lambda and loss history.
# run_params: tuned_for() of the case's disease and feedback design; issue_weeks: the weeks a
# re-run case re-runs (apply_methods(); default all).
run_methods <- function(case, run_params, issue_weeks = NULL) {
  run <- apply_methods(case, settings, run_params, weeks = issue_weeks)
  predictions <- run$predictions

  # SCORE AGAINST THE CHOSEN TARGET -------------------------------------
  # The online updates saw only the feedback released so far: each forecast's report
  # current at the first later round feedback_delay weeks after its first valid report
  # (re-run: the same rows, each as known at the issue week). Scores use score_truth():
  # the latest vintage by default, the same truth for every feedback design.
  # Counts are already scaled: do not divide quantile loss a second time. The log-scale
  # loss (score_quantiles(score_scale = "log")) needs counts. Windows lie in the evaluation
  # seasons: each season, and "all" of them together (the default of the plots).
  count_scale <- run$scaling$scale
  truth <- score_truth(case, settings$score_target, settings$score_filled) / count_scale
  weekly <- lapply(predictions, score_quantiles, truth = truth, alpha = case$alpha)
  weekly_log <- lapply(predictions, function(q) score_quantiles(q * count_scale,
    truth * count_scale, case$alpha, score_scale = "log"))
  # A method whose issued forecasts equal an earlier method's is that method's duplicate
  # here, e.g. the first restart copy equals Fixed MultiQT when the warmup is replayed
  # and restart_multiplier = fixed_multiplier.
  methods <- names(predictions)
  duplicate_of <- character()
  for (i in seq_along(methods)[-1]) {
    twin <- Find(function(j) identical(unname(predictions[[methods[i]]]),
      unname(predictions[[methods[j]]])), seq_len(i - 1L))
    if (!is.null(twin)) duplicate_of[methods[i]] <- methods[twin]
  }
  scores <- scores_by_level <- list()
  for (season in seasons) for (phase in phases) {
    weeks <- which(is.finite(truth) & phase_weeks(case, phase, settings$weeks_before_peak,
      settings$weeks_after_peak) & (season == "all" | case$season == season))
    for (method in methods) {
      s <- summarize_scores(weekly[[method]], case$alpha, weeks)
      key <- paste(method, season, phase)
      scores[[key]] <- data.frame(method, season, phase, n_weeks = s$n_weeks,
        quantile_loss = s$quantile_loss,
        log_quantile_loss = summarize_scores(weekly_log[[method]], case$alpha, weeks)$quantile_loss,
        coverage_error = s$coverage_error, duplicate_of = unname(duplicate_of[method]))
      if (season == "all") scores_by_level[[key]] <- data.frame(method, phase,
        n_weeks = s$n_weeks, s$by_level)
    }
  }
  id <- data.frame(disease = case$disease, location = case$location,
    location_name = case$location_name, horizon = case$horizon,
    feedback_delay = case$feedback_delay, feedback_revisions = case$feedback_revised,
    feedback_protocol = case$feedback_protocol)
  c(run, list(truth = truth, weekly = weekly, duplicate_of = duplicate_of,
    scores = data.frame(id, do.call(rbind, scores), row.names = NULL),
    scores_by_level = data.frame(id, do.call(rbind, scores_by_level), row.names = NULL)))
}
locations_of <- function(disease) unique(unlist(lapply(prepared$cases, function(case)
  if (case$disease == disease) case$location)))
# The settings tuned with each disease and feedback design (rows of tradeoff_feedback).
tradeoff_params <- lapply(setNames(score_diseases, score_diseases), function(disease)
  Map(tuned_for, tradeoff_feedback$delay, tradeoff_feedback$revisions, disease))
# Short names for the panel headings: "delay d" or "re-run".
feedback_names <- with(tradeoff_feedback, ifelse(revisions,
  ifelse(delay > 0, paste("re-run, delay", delay), "re-run"), paste("delay", delay)))
feedback_files <- gsub("[ ,-]", "", feedback_names)    # File names: "delay0", "rerun".

# 3. PLOT SETTINGS ---------------------------------------------------
# Labels show each method's own value of its tuned setting (tuned_setting in R/03) in the
# settings p of one disease and feedback design (tuned_for()).
method_labels <- function(p) {
  par_of <- function(method, name) method_parameter(p, method, name)
  residual_label <- paste("RMS of", p$median_window)
  base_width <- paste0(round(100 * diff(p$median_base_levels)), "%")
  value <- function(method) par_of(method, tuned_setting[[method]])
  eta <- function(method) paste0("eta=", value(method))
  outer <- function(method) paste0("outer m=", value(method))
  H_of <- function(method) paste0("H=", value(method))
  lambda_of <- function(method) paste0("lambda=", value(method))
  learned <- paste0("H=", p$score_H, ", sharing 1/", 2 * p$score_H, ", lambda on the variance")
  c(original = "Original",
    fixed = paste0("Fixed MultiQT (eta=", value("fixed"), "*S0)"),
    leaky = paste0("Leaky MultiQT (eta=", value("leaky"), "*S0, leak ", par_of("leaky", "leak"),
      ")"),
    ref_op = paste0("Delphi REF-op (sqrt, eta=", par_of("ref_op", "ref_op_multiplier"), "*Q90 of ",
      par_of("ref_op", "ref_op_window"), ", slow ", par_of("ref_op", "ref_op_slow_multiplier"),
      ", leak ", par_of("ref_op", "ref_op_fast_decay"), ")"),
    sqrt_constant = paste0("Delphi flu prod (sqrt per 100k, eta=",
      par_of("sqrt_constant", "sqrt_constant_rate"), ")"),
    absolute_slope = paste0(if (transform == "log1p") "Log-slope rate" else "Count-slope rate",
      " (min m=", value("absolute_slope"), ")"),
    oracle_target = "Target oracle (hindsight)",
    uniform = "Rate grid: equal",
    weighted = paste0("Rate grid: learned, no sharing (", H_of("weighted"), ")"),
    weighted_common = paste0("Rate grid: learned, no sharing, common weights (",
      H_of("weighted_common"), ")"),
    shared = paste0("Rate grid: learned + sharing (", H_of("shared"), ")"),
    shared_common = paste0("Rate grid: learned + sharing, common weights (", H_of("shared_common"),
      ")"),
    shared_outer = paste0("Rate grid: sharing + outer (", outer("shared_outer"), ")"),
    shared_common_outer = paste0("Rate grid: sharing, common weights + outer (",
      outer("shared_common_outer"), ")"),
    restart_persistent = paste0("First start only (m=", p$restart_multiplier, ")"),
    restart_uniform = paste0("Weekly starts: equal (m=", p$restart_multiplier, ")"),
    restart_weighted = paste0("Weekly starts: learned (m=", p$restart_multiplier, ", ",
      lambda_of("restart_weighted"), ")"),
    restart_grid = paste0("Starts x rate grid: learned (", lambda_of("restart_grid"), ")"),
    restart_uniform_outer = paste0("Starts equal + outer (", outer("restart_uniform_outer"), ")"),
    restart_weighted_outer = paste0("Starts learned + outer (", outer("restart_weighted_outer"), ")"),
    restart_grid_outer = paste0("Starts x rate grid + outer (", outer("restart_grid_outer"), ")"),
    restart_recent = paste0("Weekly starts: most recent ", value("restart_recent"),
      ", equal (m=", p$restart_multiplier, ")"),
    decline = paste0("Decline MultiQT (eta=", par_of("decline", "decline_multiplier"),
      "*S0, declining context, delta=", p$context_delta, ")"),
    score_equal = "Score experts: equal weights",
    score_learned = paste0("Score experts: learned weights per level (", learned, ")"),
    score_gated = paste0("Score experts: learned, (1 - G) x mixture + G x original (delta=",
      p$context_delta, ")"),
    studentized = paste0("Studentized median (", residual_label, ", ", p$median_start,
      " start, ", eta("studentized"), ")"),
    studentized_base = paste0("Studentized median (base ", base_width, " width, normal start, ",
      eta("studentized_base"), ")"),
    studentized_full = paste0("Studentized full base (base ", base_width, " width, ",
      eta("studentized_full"), ")"),
    studentized_level_rms = paste0("Studentized full base (level RMS of ", p$median_window,
      ", ", eta("studentized_level_rms"), ")"),
    median_constant = paste0("Median, constant sigma (", residual_label, ", ", p$median_start,
      " start, ", eta("median_constant"), ")"),
    # The median-start variants use the rate of their normal-start method (parameter_source).
    studentized_zero = paste0("Studentized median (", residual_label, ", median start, eta=",
      par_of("studentized_zero", "median_rate"), ")"),
    median_constant_zero = paste0("Median, constant sigma (", residual_label,
      ", median start, eta=", par_of("median_constant_zero", "median_rate"), ")"))
}
# Colours as in the trade-off plots and 06's tree (tradeoff_colours in R/03); the original
# in gray, as the scoring target is drawn in black.
colours <- replace(tradeoff_colours, "original", "gray40")
# Lines as the symbols of the trade-off plots: common weights across quantile levels dashed
# (triangles), + outer MultiQT dotted (squares), common weights + outer dot-dashed (diamonds).
common <- setNames(grepl("_common", names(method_short_labels)), names(method_short_labels))
common["score_equal"] <- TRUE
with_outer <- setNames(grepl("_outer$", names(method_short_labels)), names(method_short_labels))
line_types <- ifelse(with_outer, ifelse(common, 4, 3), ifelse(common, 2, 1))
point_shapes <- ifelse(common, 17, 16)
# The method families of the plots: method_groups in R/03_functions.R (shared with 06).
score_units <- count_scale_label(settings$scale_method)
loss_axis_label <- paste0("Quantile loss (", score_units, ")")
mean_loss_axis_label <- paste("Mean", tolower(loss_axis_label))
evaluation_label <- prepared$settings$evaluation_label
window_titles <- c(overall = paste0(toupper(substr(evaluation_label, 1, 1)),
  substring(evaluation_label, 2)),
  before = paste0("Up to peak (weeks -", settings$weeks_before_peak - 1, " to 0)"),
  after = paste0("After peak (weeks 1 to ", settings$weeks_after_peak, ")"))
plot_title <- function(run) paste(run$case$disease, "| horizon", run$case$horizon, "|",
  run$case$location_name, paste0("(", run$case$location, ")"), "| evaluation:",
  evaluation_label, "|",
  if (transform == "log1p") "log1p(counts) updates" else "raw-count updates",
  "| feedback", feedback_label(run$case$feedback_delay, run$case$feedback_revised),
  paste0("(", feedback_protocol()$id, ")"))
# Plots leave out a method whose identical twin is drawn in the same plot.
shown_methods <- function(run, methods) {
  twins <- run$duplicate_of[names(run$duplicate_of) %in% methods]
  setdiff(methods, names(twins)[twins %in% methods])
}
date_ticks <- function(dates) pretty(as.Date(range(dates), origin = "1970-01-01"), n = 6)
# Background of the time plots: the tuning period shaded, dotted lines at the peaks.
mark_time <- function(case, col = "gray60") {
  tuning <- case$target_date[case$period == "tuning"]
  usr <- par("usr")
  if (length(tuning)) rect(min(tuning) - 3.5, usr[3], max(tuning) + 3.5, usr[4], col = "gray93",
    border = NA)
  abline(v = as.numeric(case$peaks$peak_date), col = col, lty = 3)
}

# 4. FORECASTS AND WEEKLY QUANTILE LOSS PER LEVEL ----------------------
# Each level panel shows the scoring target (score_truth()) and every method's issued
# alpha-quantile; the panel below shows that level's weekly quantile loss. One page per run
# and method group; labels: method_labels() of the runs' settings.
blocks <- matrix(0L, 6, 3)
for (b in seq_along(plot_levels)) {
  r <- (b - 1L) %/% 3L
  blocks[2L * r + 1:2, (b - 1L) %% 3L + 1L] <- c(2L * b - 1L, 2L * b)
}
draw_forecast_pages <- function(runs, labels) for (run in runs) for (group_name in
  names(method_groups)) {
  truth_label <- score_truth_label(settings$score_target, run$case$feedback_delay)
  truth_legend <- paste0(toupper(substr(truth_label, 1, 1)), substring(truth_label, 2))
  method_names <- shown_methods(run, method_groups[[group_name]])
  layout(rbind(blocks, 2L * length(plot_levels) + 1L), heights = c(rep(c(2, 1), 3), .45))
  par(oma = c(1.5, 0, 4, 0), mgp = c(2.4, .6, 0), cex = .85)
  dates <- run$case$target_date
  ticks <- date_ticks(dates)
  for (level in plot_levels) {
    j <- match(level, run$case$alpha)
    values <- vapply(method_names, function(method) run$predictions[[method]][, j], numeric(length(dates)))
    losses <- vapply(method_names, function(method) run$weekly[[method]]$loss[, j], numeric(length(dates)))
    broken <- break_gaps(dates, values, run$truth, losses)   # Lines break at the summer 2024 gap.
    par(mar = c(1, 4.6, 2, 1))
    plot(dates, run$truth, type = "n", ylim = range(0, run$truth, values, na.rm = TRUE),
      xaxt = "n", xlab = "", ylab = paste0("Hospitalizations\n(", score_units, ")"),
      main = paste("alpha =", level))
    mark_time(run$case)
    for (method in method_names) lines(broken$dates, broken[[2]][, method], col = colours[method],
      lty = line_types[method], lwd = 1.3)
    lines(broken$dates, broken[[3]], col = "black", lwd = .7)
    points(dates, run$truth, pch = 16, cex = .45)
    par(mar = c(2.6, 4.6, .7, 1))
    plot(dates, losses[, 1], type = "n", ylim = c(0, max(1e-8, losses, na.rm = TRUE)),
      xaxt = "n", xlab = "", ylab = "Quantile\nloss")
    axis.Date(1, at = ticks, format = "%b %Y", cex.axis = .85)
    mark_time(run$case)
    for (method in method_names) lines(broken$dates, broken[[4]][, method], col = colours[method],
      lty = line_types[method], lwd = 1.1)
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  # Column-wise widths, so long labels fit; trailing spaces keep a gap between columns.
  legend("center", paste0(c(truth_legend, labels[method_names]), "    "),
    col = c("black", colours[method_names]), lty = c(NA, line_types[method_names]),
    pch = c(16, rep(NA, length(method_names))), lwd = 2, ncol = 4, bty = "n", text.width = NA)
  mtext(paste(plot_title(run), "|", group_name), side = 3, outer = TRUE, line = 2, font = 2, cex = 1.15)
  mtext(paste("Issued alpha-quantile forecasts and their weekly quantile loss in", score_units,
    "against the", truth_label, "| Vertical dotted lines: peaks | Gray: tuning period (not scored)"),
    side = 3, outer = TRUE, line = .5)
}

# 5. ISSUED WEIGHTS OVER TIME ---------------------------------------
# One page per location and expert bank, one column per mixture. Rows alpha = 0.1-0.9: one
# weight vector per quantile level; below alpha = 0.2, 0.5 and 0.8, that level's forecasts:
# the latest reported counts (black), every expert in its weight colour and the issued
# mixture (thick, in its colour of the other plots; dashed: the common-weight variant of the
# last row). Last row: one weight vector for all levels. The weights mix the forecasts before
# crossing repair or an outer correction; the outer methods keep their underlying ensemble's
# weights. Expert colours: the base gray, the rate grid by multiplier, the weekly starts in
# one palette per multiplier from dark (first start) to light (last start), the score
# experts in one qualitative palette. Score experts: equal weights for all levels (row 'All
# levels'), learned weights per level, and the learned weights backing off to the original
# with G (the context row: G and c).
weight_forecast_levels <- c(.2, .5, .8)
weight_pages <- list(
  "Frozen rate-grid weights" = data.frame(method = c("uniform", "weighted", "shared"),
    common_method = c("uniform", "weighted_common", "shared_common"),
    title = c("Equal weights", "Learned, no sharing", "Learned + sharing")),
  "Score-based expert weights" = data.frame(
    method = c("score_equal", "score_learned", "score_gated"), common_method = NA,
    title = c("Equal weights", "Learned weights", "Learned, original in decline (G)")),
  "Weekly-start weights" = data.frame(
    method = c("restart_uniform", "restart_weighted", "restart_grid", "restart_recent"),
    common_method = NA, title = c("Equal weights", "Learned", "Starts x rate grid: learned",
      "Most recent starts: equal")))
# Experts named by their methods where they are methods (the score experts).
expert_names <- function(experts) ifelse(experts %in% names(method_short_labels),
  method_short_labels[experts], experts)
restart_palettes <- c("Purples 3", "Blues 3", "Reds 3", "Oranges", "Greens 3")  # By multiplier.
# The colour of every expert of a page, named by expert, in stacking order (bottom first).
expert_fills <- function(run, page) {
  if (page != "Weekly-start weights") {
    experts <- dimnames(run$expert_forecasts[[weight_pages[[page]]$method[1]]])[[2]]
    return(setNames(c("gray45", hcl.colors(length(experts) - 1L,
      if (page == "Score-based expert weights") "Dark 3" else "Viridis")), experts))
  }
  experts <- run$restart_experts$experts
  experts <- experts[order(experts$multiplier, experts$start_date), ]
  multipliers <- unique(experts$multiplier)
  stopifnot(length(multipliers) <= length(restart_palettes))
  # Sequential palettes run dark to light; drop the lightest quarter.
  setNames(unlist(lapply(seq_along(multipliers), function(i) {
    n <- sum(experts$multiplier == multipliers[i])
    hcl.colors(ceiling(4 / 3 * n), restart_palettes[i])[seq_len(n)]
  })), experts$expert)
}
plot_weight_distribution <- function(run, weights, fills, title, row_label, first_column, last_row,
                                     dates = run$case$target_date) {
  times <- as.numeric(dates)
  par(mar = c(if (last_row) 2.4 else .3, if (first_column) 3.7 else 2.2, 1.5, .6))
  plot(NA, xlim = range(times) + c(-.5, .5) * 7, ylim = c(0, 1),
    xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "",
    ylab = if (first_column) row_label else "", main = title)
  # Full-week bars, widened by 2% so neighbours overlap and no seam shows between them.
  # Filled rounds are drawn lighter.
  if (ncol(weights) <= 40L) {
    bottom <- numeric(length(times))
    for (j in seq_len(ncol(weights))) {
      top <- bottom + weights[, j]
      rect(times - .52 * 7, bottom, times + .52 * 7, top,
        col = ifelse(run$case$filled, adjustcolor(fills[j], .45), fills[j]), border = NA)
      bottom <- top
    }
  } else {
    # Many experts (weekly starts): one image with a column of 500 pixels per target week
    # (blank without a forecast) instead of a rectangle per expert and week; an expert below
    # 1/500 of the weight drops out.
    pixels <- (500:1 - .5) / 500                         # Top pixel first.
    week <- round((times - min(times)) / 7) + 1
    stopifnot(all(abs(times - min(times) - 7 * (week - 1)) < 1e-6))   # Weekly targets.
    image <- matrix("transparent", length(pixels), max(week))
    for (t in which(rowSums(is.na(weights)) == 0)) {
      shades <- if (run$case$filled[t]) adjustcolor(fills, .45) else fills
      column <- shades[findInterval(pixels, c(0, cumsum(weights[t, ])), left.open = TRUE)]
      column[is.na(column)] <- "transparent"             # Above the total weight.
      image[, week[t]] <- column
    }
    rasterImage(as.raster(image), min(times) - .5 * 7, 0, max(times) + .5 * 7, 1,
      interpolate = FALSE)
  }
  abline(h = .5, col = adjustcolor("white", .6), lty = 3)
  abline(v = as.numeric(run$case$peaks$peak_date), col = "gray30", lty = 3)
  axis(2, at = c(0, .5, 1), las = 1, cex.axis = .75)
  if (last_row) axis.Date(1, at = date_ticks(dates), format = "%b %Y", cex.axis = .8)
  box(col = "gray70")
}
# One level's forecasts in a mixture's column: the latest reported counts, every expert in
# its weight colour and the issued mixture (thick; dashed: its common-weight variant). key:
# name the mixture's lines in the panel; title: the column title (first row of a page).
plot_mixture_forecasts <- function(run, method, common_method, level, fills, ylim,
                                   first_column, key, title = "") {
  dates <- run$case$target_date
  times <- as.numeric(dates)
  a <- match(level, run$case$alpha)
  experts <- run$expert_forecasts[[method]]
  experts <- matrix(experts[, , a], nrow(experts), dimnames = dimnames(experts)[1:2])
  mixtures <- c(method, if (!is.na(common_method) && common_method != method) common_method)
  issued <- sapply(mixtures, function(m) run$predictions[[m]][, a])
  truth <- run$case$truth_latest / run$scaling$scale
  par(mar = c(2.2, if (first_column) 3.7 else 2.2, if (nzchar(title)) 1.5 else .5, .6))
  plot(NA, xlim = range(times) + c(-.5, .5) * 7, ylim = ylim, xaxs = "i", xaxt = "n",
    xlab = "", ylab = if (first_column) paste0("Forecast, alpha = ", level, "\n(",
      score_units, ")") else "", las = 1, cex.axis = .75, main = title)
  abline(v = as.numeric(run$case$peaks$peak_date), col = "gray30", lty = 3)
  broken <- break_gaps(dates, experts, issued, truth)    # Lines break at gaps.
  many <- ncol(experts) > 40L
  for (k in seq_len(ncol(experts))) lines(broken$dates, broken[[2]][, k],
    col = adjustcolor(fills[colnames(experts)[k]], if (many) .6 else .9), lwd = if (many) .35 else .8)
  lines(broken$dates, broken[[4]], col = "black", lwd = .8)
  points(dates, truth, pch = 16, cex = .4)
  for (k in seq_along(mixtures)) {
    lines(broken$dates, broken[[3]][, k], col = "white", lwd = 4.2)
    lines(broken$dates, broken[[3]][, k], col = colours[mixtures[k]], lwd = 2, lty = k)
  }
  if (key) legend("topleft", method_short_labels[mixtures], col = colours[mixtures],
    lwd = 2, lty = seq_along(mixtures), bty = "o", box.col = NA,
    bg = adjustcolor("white", .85), cex = .75, inset = c(.005, .02))
  axis.Date(1, at = date_ticks(dates), format = "%b %Y", cex.axis = .75)
  box(col = "gray70")
}
# Key of the weekly starts: one ramp per multiplier, first start (dark) to last start (light).
draw_start_ramps <- function(run, fills, single) {
  experts <- run$restart_experts$experts
  multipliers <- sort(unique(experts$multiplier))
  for (i in seq_along(multipliers)) {
    starts <- experts[experts$multiplier == multipliers[i], ]
    names <- intersect(names(fills), starts$expert)          # First start first.
    y <- .8 - (i - 1) * .3
    x <- seq(.17, .37, length.out = length(names) + 1L)
    rect(x[-length(x)], y - .1, x[-1], y + .1, col = fills[names], border = NA)
    text(.165, y, paste0("Weekly starts, m=", multipliers[i],
      if (multipliers[i] %in% single) " (all columns)" else " (grid only)"), adj = 1, cex = .95)
    text(.375, y, paste("first start", format(min(starts$start_date), "%b %Y"), "(dark) to last",
      format(max(starts$start_date), "%b %Y"), "(light)"), adj = 0, cex = .85)
  }
}
# Rows of a page: each level's weights (only if a mixture of the page has one weight vector
# per level; empty for the others), below 0.2, 0.5 and 0.8 the forecasts at that level, then
# the weights that are not quantile specific: the vector for all levels of each mixture (or
# of its common-weight variant), and for score_gated the target signal (G: line; c: points).
row_inches <- c(weights = 1.25, forecasts = 2.3, common = 1.25, context = .9)
row_titles <- c(common = "All levels:\nweights", context = "Context\n(G; c)")
weight_page_height <- sum(row_inches[c(rep("weights", length(plot_levels)),
  rep("forecasts", length(weight_forecast_levels)), "common", "context")]) + 1.2 + 3.7
per_level <- function(run, method) dim(run$weights[[method]])[3] > 1L
page_rows <- function(run, columns) {
  common <- !vapply(columns$method, function(m) per_level(run, m), logical(1))
  kind <- level <- NULL
  for (a in plot_levels) {
    shown <- c(if (!all(common)) "weights", if (a %in% weight_forecast_levels) "forecasts")
    kind <- c(kind, shown)
    level <- c(level, rep(a, length(shown)))
  }
  extra <- c(if (any(common) || any(!is.na(columns$common_method))) "common",
    if ("score_gated" %in% columns$method) "context")
  data.frame(kind = c(kind, extra), level = c(level, rep(NA, length(extra))))
}
# The context of score_gated at each target: G (line) and c (points), 0 without a signal.
plot_context <- function(run, method, first_column, last_row) {
  signal <- if (method == "score_gated") run$score_gate
  if (is.null(signal)) return(plot.new())
  dates <- run$case$target_date
  times <- as.numeric(dates)
  par(mar = c(if (last_row) 2.4 else .3, if (first_column) 3.7 else 2.2, .5, .6))
  plot(NA, xlim = range(times) + c(-.5, .5) * 7, ylim = c(0, 1.05), xaxs = "i", yaxs = "i",
    xaxt = "n", yaxt = "n", xlab = "", ylab = if (first_column) row_titles[["context"]] else "")
  abline(h = .5, col = "gray80", lty = 3)
  abline(v = as.numeric(run$case$peaks$peak_date), col = "gray30", lty = 3)
  points(times, signal$context, pch = 15, cex = .35, col = "gray35")
  broken <- break_gaps(dates, signal$gate)
  lines(broken$dates, broken[[2]], col = colours[method], lwd = 1.5)
  axis(2, at = c(0, .5, 1), las = 1, cex.axis = .75)
  if (last_row) axis.Date(1, at = date_ticks(dates), format = "%b %Y", cex.axis = .8)
  box(col = "gray70")
}
draw_weight_pages <- function(runs) for (run in runs) for (page in names(weight_pages)) {
  columns <- weight_pages[[page]]
  columns <- columns[columns$method %in% names(run$weights), ]
  if (!nrow(columns)) next
  fills <- expert_fills(run, page)
  rows <- page_rows(run, columns)
  has_common <- any(!is.na(columns$common_method) & columns$common_method != columns$method)
  last_row <- nrow(rows)
  # One axis per forecast row, shared by its columns.
  limits <- lapply(setNames(weight_forecast_levels, weight_forecast_levels), function(level) {
    a <- match(level, run$case$alpha)
    range(run$case$truth_latest / run$scaling$scale, lapply(seq_len(nrow(columns)), function(j)
      c(run$expert_forecasts[[columns$method[j]]][, , a], sapply(na.omit(c(columns$method[j],
        columns$common_method[j])), function(m) run$predictions[[m]][, a]))), finite = TRUE)
  })
  panel_ids <- matrix(seq_len(nrow(rows) * nrow(columns)), ncol = nrow(columns), byrow = TRUE)
  # Rows in inches; a shorter page keeps the free space below its key.
  heights <- c(row_inches[rows$kind], 1.2)
  layout(rbind(panel_ids, max(panel_ids) + 1L, max(panel_ids) + 2L),
    heights = c(lcm(2.54 * heights), 1))
  par(oma = c(4.8, .8, 5.2, .8), mgp = c(2, .4, 0), cex = .9)
  for (row in seq_len(nrow(rows))) for (j in seq_len(nrow(columns))) {
    kind <- rows$kind[row]
    level <- rows$level[row]
    method <- columns$method[j]
    if (kind == "forecasts") {
      plot_mixture_forecasts(run, method, columns$common_method[j], level, fills,
        limits[[as.character(level)]], j == 1L, level == weight_forecast_levels[1],
        if (row == 1L) columns$title[j] else "")
      next
    }
    if (kind == "context") {
      plot_context(run, method, j == 1L, row == last_row)
      next
    }
    weights <- switch(kind,
      weights = if (per_level(run, method)) run$weights[[method]][, , match(level, run$case$alpha)],
      common = if (!per_level(run, method)) run$weights[[method]][, , 1] else
        if (!is.na(columns$common_method[j])) run$weights[[columns$common_method[j]]][, , 1])
    if (is.null(weights)) {                              # Not this mixture's kind of weights.
      plot.new()
      if (row == 1L) title(columns$title[j])
      next
    }
    weights <- matrix(weights, nrow(run$weights[[method]]),
      dimnames = dimnames(run$weights[[method]])[1:2])
    stack <- intersect(names(fills), colnames(weights))
    stopifnot(length(stack) == ncol(weights))
    plot_weight_distribution(run, weights[, stack, drop = FALSE], fills[stack],
      if (row == 1L) columns$title[j] else "",
      if (kind == "weights") paste("alpha =", level) else row_titles[[kind]],
      first_column = j == 1L, last_row = row == last_row)
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  if (page == "Weekly-start weights") {
    single <- unique(run$restart_experts$experts$multiplier[run$restart_experts$experts$expert %in%
      colnames(run$weights[[columns$method[1]]])])
    draw_start_ramps(run, fills, single)
  } else legend(0, .5, expert_names(names(fills)), fill = fills, border = NA, yjust = .5,
    ncol = ceiling(length(fills) / 2), bty = "n", cex = .95,
    title = "Experts: weights and forecasts", title.adj = 0, title.font = 2)
  legend(1, .5, c("Latest reported counts", "Each expert: its weight colour",
    "Mixture: thick, in the colour of its method", if (has_common)
      "Common-weight variant: dashed"), xjust = 1, yjust = .5,
    col = c("black", "gray55", "gray25", if (has_common) "gray25"),
    lty = c(1, 1, 1, if (has_common) 2), lwd = c(.8, .8, 2, if (has_common) 2),
    pch = c(16, NA, NA, if (has_common) NA), pt.cex = .6, bty = "n", cex = .95)
  plot.new()                                             # The free space of a shorter page.
  mtext(page, side = 3, outer = TRUE, line = 2.4, font = 2, cex = 1.3)
  mtext(plot_title(run), side = 3, outer = TRUE, line = .9, cex = .9)
  mtext("Target date | Weight at issuance for each target | Vertical dotted lines: peaks",
    side = 1, outer = TRUE, line = .4, cex = .85)
  mtext(paste("Rows alpha = 0.1-0.9 (mixtures with one weight vector per level): the weights of",
    "each level; below alpha = 0.2, 0.5 and 0.8 the forecasts at that level (per-level mixtures",
    "after crossing repair, before any outer correction)."), side = 1, outer = TRUE, line = 1.6,
    cex = .8)
  mtext(paste("Rows 'All levels': one weight vector for all levels (the mixture, or its",
    "common-weight variant), learned from the mean quantile loss over the levels (equal",
    "weights: not learned); context: G (line) and c (points) of the target signal."),
    side = 1, outer = TRUE, line = 2.7, cex = .8)
  if (any(run$case$filled)) mtext(paste("Lighter bars: rounds without a CMU forecast, filled with",
    "the latest earlier CMU forecast of the same target week."), side = 1, outer = TRUE,
    line = 3.8, cex = .8)
}

# 6. EVERY DISEASE, HORIZON AND FEEDBACK DESIGN ---------------------
# One PDF per combination and plot (one folder per plot), its pages going through the
# locations of multiqt_settings(); the combinations run in parallel, re-runs first (the
# longest). A re-run combination re-runs every week, so the plots show the tuning period too;
# its scores (evaluation weeks only) equal those of section 7. Each combination returns its
# scores and result tables; the runs of disease_name, horizon_id and the feedback design of
# multiqt_settings() also go to issued_forecasts.rds.
plot_folders <- c(forecasts = "multiqt_method_forecasts", weights = "multiqt_weight_timelines")
for (folder in plot_folders) dir.create(file.path(figure_dir, folder), showWarnings = FALSE)
combinations <- expand.grid(feedback = seq_len(nrow(tradeoff_feedback)),
  horizon = settings$horizons, disease = score_diseases, stringsAsFactors = FALSE)
combinations <- combinations[order(!tradeoff_feedback$revisions[combinations$feedback]), ]
saved <- which(combinations$disease == disease_name & combinations$horizon == horizon_id &
  tradeoff_feedback$delay[combinations$feedback] == settings$feedback_delay &
  tradeoff_feedback$revisions[combinations$feedback] == settings$feedback_revisions)
stopifnot(length(saved) == 1L)
combination_file <- function(i) paste0(combinations$disease[i], "_horizon",
  combinations$horizon[i], "_", feedback_files[combinations$feedback[i]], ".pdf")
# The case of each row of a result table.
with_case <- function(run, rows) data.frame(disease = run$case$disease,
  location = run$case$location, location_name = run$case$location_name,
  horizon = run$case$horizon, feedback_delay = run$case$feedback_delay,
  feedback_revisions = run$case$feedback_revised, rows, row.names = NULL)
# The result tables of one combination's runs.
run_results <- function(runs, run_params) list(
  scores = do.call(rbind, lapply(runs, `[[`, "scores")),
  scores_by_level = do.call(rbind, lapply(runs, `[[`, "scores_by_level")),
  scales = do.call(rbind, lapply(runs, function(run) with_case(run, data.frame(
    feedback_protocol = run$case$feedback_protocol,
    scale_method = run$scaling$method, scale_used = run$scaling$used, count_scale = run$scaling$scale,
    transform, scale_as_of = run$case$scale_as_of, reference_weeks = run$scaling$reference_weeks,
    initialization_date = run$s0$date, warmup_weeks = run$s0$warmup_weeks,
    raw_error_scale = run$s0$raw_scale, error_scale = run$s0$scale,
    sigma_window = run$sigma$window, sigma_warmup_weeks = run$median_warmup_weeks,
    constant_sigma = run$constant_sigma,
    base_sigma_levels = paste(run_params$median_base_levels, collapse = "-"),
    median_start = run_params$median_start)))),
  expert_grid = do.call(rbind, lapply(runs, function(run)
    with_case(run, run$rate_experts$experts))),
  restart_experts = do.call(rbind, lapply(runs, function(run) with_case(run,
    data.frame(run$restart_experts$experts, birth_week = run$restart_experts$birth_week)))),
  # The target signal of the score experts (target_signal() in R/03: D, c and G, 0 without
  # a signal), by which Decline MultiQT and score_gated (G) act.
  score_gate = do.call(rbind, lapply(runs, function(run) if (!is.null(run$score_gate))
    with_case(run, data.frame(reference_date = run$case$reference_date,
      target_date = run$case$target_date, signal = run$score_gate$signal,
      context = run$score_gate$context, gate = run$score_gate$gate)))))
plotted <- parallel::mclapply(seq_len(nrow(combinations)), function(i) {
  disease <- combinations$disease[i]
  horizon <- combinations$horizon[i]
  k <- combinations$feedback[i]
  run_params <- tradeoff_params[[disease]][[k]]
  labels <- method_labels(run_params)
  runs <- list()
  for (location in intersect(settings$location_ids, locations_of(disease))) {
    case <- get_case(prepared, disease, location, horizon, tradeoff_feedback$delay[k],
      tradeoff_feedback$revisions[k])
    check_peak_window(case, settings$weeks_before_peak, settings$weeks_after_peak,
      score_truth(case, settings$score_target, settings$score_filled))
    runs[[location]] <- cached("04_runs", run_methods, case, run_params)
  }
  message("Plots: ", disease, " | horizon ", horizon, " | ", feedback_names[k])
  file <- combination_file(i)
  cairo_pdf(file.path(figure_dir, plot_folders[["forecasts"]], file), width = 18, height = 16,
    family = "DejaVu Sans")
  draw_forecast_pages(runs, labels)
  dev.off()
  cairo_pdf(file.path(figure_dir, plot_folders[["weights"]], file), width = 18,
    height = weight_page_height, family = "DejaVu Sans")
  draw_weight_pages(runs)
  dev.off()
  list(results = run_results(runs, run_params), runs = if (i == saved) runs)
}, mc.cores = cores, mc.preschedule = FALSE)
failed <- vapply(plotted, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
if (any(failed)) stop("Plots of ", paste(unlist(combinations[which(failed)[1], ]),
  collapse = " "), " failed: ", if (is.null(plotted[[which(failed)[1]]])) "no result (memory?)"
  else plotted[[which(failed)[1]]])
saveRDS(plotted[[saved]]$runs, file.path(result_dir, "issued_forecasts.rds"))
# The saved combination also keeps the single files (multiqt_method_forecasts.pdf, ...),
# whose pages notes/ loads.
for (folder in plot_folders) if (file.exists(file.path(figure_dir, folder, combination_file(saved))))
  file.copy(file.path(figure_dir, folder, combination_file(saved)),
    file.path(figure_dir, paste0(folder, ".pdf")), overwrite = TRUE)
results <- lapply(setNames(names(plotted[[1]]$results), names(plotted[[1]]$results)),
  function(table) do.call(rbind, lapply(plotted, function(p) p$results[[table]])))

# 7. SCORES OF EVERY LOCATION ----------------------------------------
# Scores of every prepared location, horizon and feedback design of score_diseases for the
# trade-off pages (section 8) and R/06_evaluation_tree.R, in parallel; only the scores are
# kept. Each disease and design uses the settings tuned with it; a re-run re-runs only the
# evaluation weeks. The plotted locations (section 6) are reused.
tradeoff_cases <- do.call(rbind, lapply(score_diseases, function(disease)
  expand.grid(disease = disease, location = setdiff(locations_of(disease),
    settings$location_ids), horizon = settings$horizons,
    feedback = seq_len(nrow(tradeoff_feedback)), stringsAsFactors = FALSE)))
for (disease in score_diseases) message("Trade-off scores: ", disease, " | ",
  length(locations_of(disease)), " locations x ", length(settings$horizons), " horizons x ",
  nrow(tradeoff_feedback), " feedback designs")
# The scores of one case.
evaluation_scores <- function(case, run_params)
  run_methods(case, run_params, which(case$period == "evaluation"))$scores
other_scores <- parallel::mclapply(seq_len(nrow(tradeoff_cases)), function(i) {
  disease <- tradeoff_cases$disease[i]
  k <- tradeoff_cases$feedback[i]
  case <- get_case(prepared, disease, tradeoff_cases$location[i], tradeoff_cases$horizon[i],
    tradeoff_feedback$delay[k], tradeoff_feedback$revisions[k])
  cached("04_scores", evaluation_scores, case, tradeoff_params[[disease]][[k]])
}, mc.cores = cores)
failed <- vapply(other_scores, inherits, logical(1), "try-error")
if (any(failed)) stop(other_scores[[which(failed)[1]]])
all_scores <- rbind(results$scores, do.call(rbind, other_scores))
rownames(all_scores) <- NULL
write_result(results$scores, "summary_scores.csv")                # The plotted locations.
write_result(all_scores, "summary_scores_all_locations.csv")    # Every disease of score_diseases.
write_result(results$scores_by_level, "summary_scores_by_level.csv")
write_result(results$scales, "scales.csv")
write_result(results$expert_grid, "expert_grid.csv")
write_result(results$restart_experts, "restart_experts.csv")
write_result(results$score_gate, "score_gate.csv")

# 8. CALIBRATION ERROR AND QUANTILE LOSS ON MATCHED WINDOWS ----------
# Per disease of score_diseases: one page per location, after a first page with all prepared
# locations except the US. Columns:
# horizons; rows: the feedback designs of tradeoff_feedback, each with its three windows.
# Every method in each panel. The loss is scored on the tracker's scale: log(counts + 1)
# for the log1p tracker, counts per 100k residents for raw counts. Each feedback design
# uses the settings tuned (R/05_tuning.R) with it.
# Colours: tradeoff_colours in R/03_functions.R (one colour family per method group, one
# shade per method; common-weight and outer variants share their method's shade; also
# 06's colours). CE = mean over levels of |coverage - alpha| over the window's weeks.
# Every panel has its own axis ranges (panel_limits()); up to peak includes the peak week.
# Original: open circle. Filled circle: weights per level; triangle: common weights;
# square: + outer MultiQT; diamond: common weights + outer.
tradeoff_methods <- names(tradeoff_colours)
outer_variant <- grepl("_outer$", tradeoff_methods)
tradeoff_shapes <- setNames(ifelse(outer_variant, ifelse(common[tradeoff_methods], 18, 15),
  ifelse(common[tradeoff_methods], 17, 16)), tradeoff_methods)
tradeoff_shapes["original"] <- 1
tradeoff_sizes <- ifelse(tradeoff_methods == "original", 2.2, 1.6)
names(tradeoff_sizes) <- tradeoff_methods
# Axis range of one panel: the values, with the loss axis capped at `cap` times the
# original's loss (reference), so one diverging method does not squeeze the others.
panel_limits <- function(x, reference = NA, cap = 4) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0, 1))
  upper <- if (is.finite(reference) && reference > 0) min(max(x), cap * reference) else max(x)
  extendrange(c(min(x), upper))
}
tradeoff_loss <- if (transform == "log1p") "log_quantile_loss" else "quantile_loss"
tradeoff_loss_label <- if (transform == "log1p") "Mean quantile loss of log(counts + 1)" else
  mean_loss_axis_label
panel_rows <- expand.grid(phase = phases, feedback = seq_len(nrow(tradeoff_feedback)),
  stringsAsFactors = FALSE)                 # Rows: the windows of each feedback design.
panel_height <- 2.6                         # Inches per row of panels.
legend_height <- 2.6
# One disease's pages: all prepared locations except the US, then each location of
# settings$location_ids. Legend values: the settings of that disease and the feedback design
# of multiqt_settings().
draw_tradeoff_pages <- function(disease) {
  page_params <- tuned_for(settings$feedback_delay, settings$feedback_revisions, disease)
  labels <- method_labels(page_params)
  base_width <- paste0(round(100 * diff(page_params$median_base_levels)), "%")
  all_locations <- setdiff(locations_of(disease), settings$summary_exclude)
  tradeoff_scores <- all_scores[all_scores$disease == disease & all_scores$season == "all" &
    all_scores$method %in% tradeoff_methods, ]
  # Duplicates (identical forecasts) of a drawn method are not drawn; those hidden in
  # every case also leave the legend, and a footnote names each pair of the location pages.
  tradeoff_scores$hidden <- tradeoff_scores$duplicate_of %in% tradeoff_methods
  hidden_everywhere <- names(which(tapply(tradeoff_scores$hidden, tradeoff_scores$method, all)))
  hidden_rows <- tradeoff_scores[tradeoff_scores$hidden &
    tradeoff_scores$location %in% settings$location_ids, ]
  hidden_pairs <- unique(paste(labels[hidden_rows$method], "=", labels[hidden_rows$duplicate_of],
    recycle0 = TRUE))
  # Legend: one column per group, headed by the group name; the two studentized groups
  # share a column, so the legend fits the page width.
  legend_groups <- lapply(setNames(names(method_groups), names(method_groups)), function(group_name)
    setdiff(c(if (group_name == names(method_groups)[1]) "original",
      names(home_group)[home_group == group_name]), hidden_everywhere))
  legend_columns <- c(as.list(setdiff(names(legend_groups), c("Studentized median",
    "Studentized full base"))), list(c("Studentized median", "Studentized full base")))
  legend_columns <- lapply(legend_columns, function(groups)
    unlist(lapply(groups, function(group_name) c(group_name, legend_groups[[group_name]]))))
  legend_rows <- max(lengths(legend_columns))
  legend_entries <- unlist(lapply(legend_columns, function(column)
    c(column, rep("", legend_rows - length(column)))))
  legend_method <- legend_entries %in% tradeoff_methods
  legend_header <- legend_entries %in% names(legend_groups)
  tradeoff_title <- function(location = NULL) paste(disease, "|",
    if (is.null(location)) paste("all", pooled_label(length(all_locations), settings)) else
      paste0(tradeoff_scores$location_name[match(location, tradeoff_scores$location)],
        " (", location, ")"), "| evaluation:", evaluation_label, "|",
    if (transform == "log1p") "log1p(counts) updates, loss of log(counts + 1)" else
      paste("raw-count updates, loss in", score_units))
  # One page; location NULL: all prepared locations except the US (pooled_rows()). Every
  # panel has its own axis ranges; a loss beyond the axis is drawn at its edge as ">".
  draw_tradeoff_page <- function(location = NULL) {
    all_page <- is.null(location)
    page <- if (all_page) pooled_rows(tradeoff_scores[!tradeoff_scores$method %in%
      hidden_everywhere, ], settings) else
      tradeoff_scores[tradeoff_scores$location == location & !tradeoff_scores$hidden, ]
    page <- page[order(page$method == "original"), ]   # Original last, on top of the others.
    # Notes below the legend, one per line; the bottom margin fits them.
    notes <- c("Calibration error (CE): mean absolute coverage residual after averaging weeks within each quantile level.",
      paste("Open circle: original forecasts. Filled circle: weights per quantile level; triangle: one",
        "weight vector shared by all levels; square / diamond: the same + outer MultiQT."),
      # A long list of identical pairs on several lines (labels only: no relabelled words).
      if (length(hidden_pairs)) strwrap(paste("Not drawn (identical forecasts):",
        paste(hidden_pairs, collapse = "; ")), width = 200),
      paste0("Each feedback design uses the settings tuned with it; the values in the legend are ",
        "those of feedback ", feedback_label(settings$feedback_delay, settings$feedback_revisions),
        " (all: output/results/tuning)."),
      paste0("Feedback protocol ", feedback_protocol()$id, ": Wednesday issue cutoffs (actual ",
        "generation dates where known); S0 and the weighting loss scale from the first ",
        settings$warmup_rounds, " released rows in hindsight, frozen for the re-runs."),
      # Studentized median: what sigma_t is, and the first weeks without a residual sigma_t,
      # which are the original base (e.g. most of the flu up-to-peak window).
      paste0("Studentized median methods: sigma_t = RMS of the last ", page_params$median_window,
        " released median residuals, or the base's ", base_width, " interval width / ",
        round(diff(qnorm(page_params$median_base_levels)), 3), " (a normal sd)."),
      paste("With residuals, the original base is issued until sigma_t exists; these weeks are",
        "scored as the base. Median start: all quantiles start at the median (offsets 0)."),
      paste0("Studentized full base: offsets start at the base; level RMS = the RMS of each ",
        "level's last ", page_params$median_window,
        " released residuals (weighted isotonic projection)."))
    n_panels <- nrow(panel_rows) * length(settings$horizons)
    layout(rbind(matrix(seq_len(n_panels), nrow(panel_rows), byrow = TRUE), n_panels + 1L),
      heights = c(rep(panel_height, nrow(panel_rows)), legend_height))
    par(oma = c(1.8 + 1.3 * (length(notes) - 1), 3.4, 7.8, 0), mgp = c(2.2, .6, 0))
    for (r in seq_len(nrow(panel_rows))) for (horizon in settings$horizons) {
      phase <- panel_rows$phase[r]
      k <- panel_rows$feedback[r]
      rows <- page[page$horizon == horizon & page$feedback_delay == tradeoff_feedback$delay[k] &
        page$feedback_revisions == tradeoff_feedback$revisions[k] & page$phase == phase, ]
      # All locations: every location as a small, light symbol per method (legend shapes)
      # and the mean over locations as the usual symbol; the axes span the 1st to 99th
      # percentile of the single locations and the means.
      drawn <- if (all_page) aggregate(rows[c(tradeoff_loss, "coverage_error")], rows["method"],
        mean, na.rm = TRUE) else rows
      drawn <- drawn[order(drawn$method == "original"), ]
      spread <- function(x) if (all_page) quantile(x, c(.01, .99), na.rm = TRUE)
      xlim <- panel_limits(c(spread(rows[[tradeoff_loss]]), drawn[[tradeoff_loss]]),
        drawn[[tradeoff_loss]][drawn$method == "original"])
      ylim <- panel_limits(c(spread(rows$coverage_error), drawn$coverage_error))
      beyond <- is.finite(drawn[[tradeoff_loss]]) & drawn[[tradeoff_loss]] > xlim[2]
      par(mar = c(3.4, 3.6, 2.9, .8))
      plot(NA, xlim = xlim, ylim = ylim, xlab = tradeoff_loss_label, ylab = "Calibration error")
      mtext(window_titles[phase], side = 3, line = 1.4, font = 2, cex = .75)
      mtext(paste0("Horizon ", horizon, " | ", feedback_names[k], " | ",
        if (all_page) paste(length(unique(rows$location)), "locations") else
          paste(rows$n_weeks[1], "weeks"), if (any(beyond)) " | >: beyond 4x the original's"),
        side = 3, line = .3, cex = .65)
      grid(col = "gray90")
      if (all_page) points(rows[[tradeoff_loss]], rows$coverage_error,
        pch = tradeoff_shapes[rows$method], cex = .5 * tradeoff_sizes[rows$method],
        col = adjustcolor(tradeoff_colours[rows$method], alpha.f = .35))
      points(pmin(drawn[[tradeoff_loss]], xlim[2]), drawn$coverage_error,
        col = tradeoff_colours[drawn$method], pch = ifelse(beyond, 62, tradeoff_shapes[drawn$method]),
        cex = tradeoff_sizes[drawn$method], lwd = ifelse(drawn$method == "original", 2, 1))
      if (!any(is.finite(drawn$coverage_error)))
        text(mean(xlim), mean(ylim), "No observations inside the peak window", cex = .8)
    }
    par(mar = c(0, 0, 0, 0))
    plot.new()
    key <- function(cex, plot = TRUE) legend("center", ifelse(legend_method,
      labels[legend_entries], legend_entries),
      col = ifelse(legend_method, tradeoff_colours[legend_entries], NA),
      pch = ifelse(legend_method, tradeoff_shapes[legend_entries], NA),
      pt.cex = ifelse(legend_method, tradeoff_sizes[legend_entries], 1),
      pt.lwd = ifelse(legend_entries == "original", 2, 1),
      text.font = ifelse(legend_header, 2, 1), ncol = length(legend_columns), bty = "n", cex = cex,
      text.width = NA, plot = plot)            # Column-wise widths, so every group column fits.
    # The largest text size up to 1.1 at which the key fits the page width.
    legend_cex <- 1.1
    while (legend_cex > .6 && key(legend_cex, FALSE)$rect$w > diff(par("usr")[1:2]))
      legend_cex <- legend_cex - .05
    key(legend_cex)
    # Feedback blocks: a label in the left margin and a line between blocks; horizons on top.
    inner_height <- nrow(panel_rows) * panel_height + legend_height
    block <- length(phases) * panel_height / inner_height
    for (k in seq_len(nrow(tradeoff_feedback))) mtext(paste("Feedback", feedback_label(
      tradeoff_feedback$delay[k], tradeoff_feedback$revisions[k])), side = 2, outer = TRUE,
      line = 1.6, at = 1 - (k - .5) * block, font = 2, cex = 1.2)
    par(xpd = NA)
    for (k in seq_len(nrow(tradeoff_feedback))[-1]) segments(grconvertX(0, "nic"),
      grconvertY(1 - (k - 1) * block, "nic"), grconvertX(1, "nic"),
      grconvertY(1 - (k - 1) * block, "nic"), col = "gray45", lwd = 1.5)
    par(xpd = FALSE)
    for (j in seq_along(settings$horizons)) mtext(paste("Horizon", settings$horizons[j]),
      side = 3, outer = TRUE, line = .5, at = (j - .5) / length(settings$horizons), font = 2,
      cex = 1.2)
    mtext(tradeoff_title(location), side = 3, outer = TRUE, line = 5.9, font = 2, cex = 1.45)
    mtext(paste0("Rows: feedback design (protocol ", feedback_protocol()$id, "; delay d: learn ",
      "once, d weeks after the first report, from the report current then; re-run: re-run ",
      "every week on the reports known then) and window; columns: horizon."), side = 3,
      outer = TRUE, line = 3.9, cex = 1.05)
    if (all_page) mtext(paste("Small: single locations; large: means over locations;",
      "axes: 1st to 99th percentile of both."), side = 3, outer = TRUE, line = 2.5, cex = 1.05)
    mtext(notes, side = 1, outer = TRUE, line = .6 + 1.3 * (seq_along(notes) - 1), cex = .95)
  }
  draw_tradeoff_page()
  for (location in intersect(settings$location_ids, tradeoff_scores$location))
    draw_tradeoff_page(location)
}
cairo_pdf(file.path(figure_dir, "multiqt_peak_reasoning_tradeoff.pdf"), width = 19,
  height = nrow(panel_rows) * panel_height + legend_height + 2.8, family = "DejaVu Sans")
for (disease in score_diseases) draw_tradeoff_pages(disease)
dev.off()
