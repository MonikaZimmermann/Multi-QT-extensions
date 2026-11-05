# Apply the recalibration methods of R/03_functions.R to each location and plot them.
# Run R/01_data_overview.R once to prepare output/data/peak_aligned.rds.
# Run from the project root: Rscript R/04_methods.R
# Offsets and ensemble weights are fitted independently for each location.
source("R/03_functions.R")

# 1. SETTINGS AND PREPARED DATA -------------------------------------
disease_name <- "covid"                    # "covid" or "flu"
horizon_id <- 2                            # -1, 0, 1, 2 or 3
settings <- multiqt_settings()             # Locations, tracker, experts, H and windows (R/03_functions.R).
transform <- settings$transform
H <- settings$H
# Other rates and caps below are in tracker units (log1p counts or raw counts), except median_rate.
outer_rate <- .005
slope_rate_min <- .005
slope_rate_max <- .1                       # Tracker-unit cap; independent of the multiplier grid.
# Weekly restart comparison: each copy votes from birth, initially forecasting the base.
restart_multiplier <- .05
restart_multipliers <- c(.01, .05, .2)     # Smaller grid at each new start.
restart_lambda <- 2.72                     # Fixed meta-rate on scaled-count pinball loss.
# Diverse experts: slow MultiQT on differently corrected starting forecasts.
diverse_slow_rate <- .01                   # Tracker units, the same for every expert.
diverse_leak <- .2                         # Leaky shift per report; halves every week, cap +-.2.
diverse_damping <- c(-.5, .5, 1)           # Start from the base minus k * implied growth.
diverse_widths <- c(.75, 1.5)              # Start from the base intervals scaled by w.
# Studentized MultiQT: offsets in units of sigma_t. On the base median, sigma_t from released
# median residuals (standard) or from the width of the base; on the full base, the width of
# the base or quantile-specific scales from the base (level_sigma()) or from recent
# residuals of each level (level_residual_sigma(), same window and warmup as the median).
median_window <- 5L                        # Residuals: RMS of this many latest released median residuals.
median_min_rows <- 5L                      # Residuals: the original base until this many are released.
median_base_levels <- c(.1, .9)            # Base width: the 80% interval / 2.563, a normal sd.
median_rate <- .1                          # eta in units of sigma_t: every studentized method and location.
median_start <- "normal"                   # Residual sigma_t: "normal" or "t" start (df = median_window).
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
# Peaks are the largest latest-reported count in the disease's peak window.
# Only the target oracle uses that hindsight date to choose forecasts.
# Alignment changes evaluation coordinates, never when feedback becomes available.
# Every location starts fresh: offsets, expert weights, lambda and loss history.
run_methods <- function(case) {
  scaling <- choose_count_scale(case, settings$scale_method)
  count_scale <- scaling$scale
  s0 <- error_scale(case, settings$warmup_rounds, settings$scale_floor, transform,
    settings$replay_warmup)
  # Keep stored predictions and the scoring/weighting inputs in counts / count_scale.
  base <- case$base / count_scale
  feedback <- case$feedback_observed / count_scale
  adjust <- function(counts, eta, warmup_weeks) multiqt(counts, case$feedback_observed,
    case$feedback_by_week, case$alpha, eta, warmup_weeks, transform) / count_scale

  # 3. ONE TRACKER: FIXED RATE, SLOPE RATE OR STUDENTIZED MEDIAN --------
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

  # Studentized median MultiQT (median_multiqt()): Fixed MultiQT's coverage indicator
  # and update on the base median b_t, issued as b_t + sigma_t * offsets.
  sigma <- median_sigma(case, "residuals", median_window, median_min_rows,
    settings$scale_floor, transform)
  base_sigma <- median_sigma(case, "base", floor = settings$scale_floor, transform = transform,
    levels = median_base_levels)
  # The original base is issued until sigma_t exists (and during the S0 warmup without replay).
  median_warmup_weeks <- max(s0$warmup_weeks, sigma$warmup_weeks)
  # Offsets before the first update: normal quantiles, or t quantiles of a residual divided
  # by the RMS of median_window other normal residuals (wider tails for small windows).
  start <- switch(median_start, normal = qnorm(case$alpha), t = qt(case$alpha, median_window),
    stop("Unknown median start ", median_start))
  predictions$studentized <- median_multiqt(case, median_rate, sigma$scale, start,
    median_warmup_weeks, transform) / count_scale
  # sigma_t from the width of the base: no feedback lag, so no warmup beyond S0's. This
  # sigma_t is a normal sd, so it starts from normal quantiles: the base's own interval.
  predictions$studentized_base <- median_multiqt(case, median_rate, base_sigma$scale,
    qnorm(case$alpha), max(s0$warmup_weeks, base_sigma$warmup_weeks), transform) / count_scale
  # Residual sigma frozen at its first issued value: additive MultiQT on the median.
  # It separates reducing the base to its median from studentizing.
  constant_sigma <- sigma$scale[median_warmup_weeks + 1L]
  predictions$median_constant <- median_multiqt(case, median_rate, constant_sigma, start,
    median_warmup_weeks, transform) / count_scale
  # Studentized MultiQT on the full base: offsets start at 0 (the base) and are issued as
  # b_t^alpha + sigma_t * offsets, with one scale per week (the base width; isoreg) or one
  # per level (weighted isoreg): the normal sd of the base's central interval at that
  # level, or the RMS of that level's recent residuals y - b^alpha.
  level_scales <- level_sigma(case, settings$scale_floor, transform)
  level_residual_scales <- level_residual_sigma(case, median_window, median_min_rows,
    settings$scale_floor, transform)
  full <- function(scales, warmup_weeks = s0$warmup_weeks) multiqt(case$base,
    case$feedback_observed, case$feedback_by_week, case$alpha, median_rate, warmup_weeks,
    transform, sigma = scales) / count_scale
  predictions$studentized_full <- full(base_sigma$scale)
  predictions$studentized_levels <- full(level_scales)
  # The residual scales exist once median_min_rows residuals are released, as for the median.
  predictions$studentized_level_rms <- full(level_residual_scales, median_warmup_weeks)

  # 4. WHAT IF WE KNEW THE TRUE PEAK? ----------------------------------
  # Hindsight benchmark: bypass correction for postpeak targets, even if issued prepeak.
  # Setting eta=0 instead would only freeze the accumulated offset.
  predictions$oracle_target <- predictions$fixed
  after_peak <- case$target_date > case$peak_date
  predictions$oracle_target[after_peak, ] <- base[after_peak, ]

  # 5. SEVERAL TRACKERS: WEIGHT THEIR FORECASTS ------------------------
  # Rate grid: one persistent expert per multiplier. Restarts: a fresh copy per
  # multiplier at every issue week after warmup; no permanent base expert.
  # Diverse experts: the base plus slow MultiQT on differently corrected bases.
  rate_experts <- make_rate_experts(case, settings$expert_multipliers, s0, count_scale, transform)
  diverse_experts <- make_diverse_experts(case, s0, count_scale, transform, diverse_slow_rate,
    diverse_leak, diverse_damping, diverse_widths)
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
    restart_grid = weigh(restart_experts, lambda = restart_lambda),
    diverse_uniform = weigh(diverse_experts, lambda = 0),
    diverse_shared = weigh(diverse_experts),
    diverse_shared_common = weigh(diverse_experts, per_quantile = FALSE))
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

  # 6. SCORE AGAINST THE CHOSEN TARGET ----------------------------------
  # The online updates above saw only the feedback released so far: the report
  # available feedback_delay weeks after the first report. Scores use
  # score_truth(): the latest vintage by default.
  # Counts are already scaled: do not divide quantile loss a second time.
  truth <- score_truth(case, settings$score_target, settings$score_filled) / count_scale
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
  list(case = case, scaling = scaling, s0 = s0, sigma = sigma, base_sigma = base_sigma,
    level_scales = level_scales,
    median_warmup_weeks = median_warmup_weeks, constant_sigma = constant_sigma,
    truth = truth, predictions = predictions,
    weights = lapply(mixtures, `[[`, "weights"), lambda = lapply(mixtures, `[[`, "lambda"),
    rate_experts = rate_experts, restart_experts = restart_experts,
    diverse_experts = diverse_experts, weekly = weekly,
    scores = data.frame(id, do.call(rbind, scores), row.names = NULL),
    scores_by_level = data.frame(id, do.call(rbind, scores_by_level), row.names = NULL))
}
runs <- list()
for (location in settings$location_ids) {
  case <- get_case(prepared, disease_name, location, horizon_id, settings$feedback_delay)
  check_peak_window(case, settings$weeks_before_peak, settings$weeks_after_peak,
    score_truth(case, settings$score_target, settings$score_filled))
  message("Methods: ", disease_name, " | ", location, " | horizon ", horizon_id)
  runs[[location]] <- run_methods(case)
}
# A method whose issued forecasts equal an earlier method's at a location is that
# method's duplicate there, e.g. the first restart copy equals Fixed MultiQT when
# the warmup is replayed and restart_multiplier = reference_multiplier.
for (location in names(runs)) {
  p <- runs[[location]]$predictions
  methods <- names(p)
  runs[[location]]$duplicate_of <- unlist(lapply(seq_along(methods)[-1], function(i) {
    twin <- Find(function(j) identical(unname(p[[methods[i]]]), unname(p[[methods[j]]])),
      seq_len(i - 1L))
    if (is.null(twin)) NULL else setNames(methods[twin], methods[i])
  }))
}
saveRDS(runs, file.path(result_dir, "issued_forecasts.rds"))
summary_scores <- do.call(rbind, lapply(runs, `[[`, "scores"))
rownames(summary_scores) <- NULL
write_result(summary_scores, "summary_scores.csv")
# Scores of every prepared location of the disease and horizon, for the first row of the
# trade-off page; only the scores are kept for locations without their own plots.
all_locations <- unique(unlist(lapply(prepared$cases, function(case)
  if (case$disease == disease_name && case$horizon == horizon_id) case$location)))
all_scores <- do.call(rbind, lapply(all_locations, function(location) {
  if (location %in% names(runs)) return(runs[[location]]$scores)
  message("Scores: ", disease_name, " | ", location, " | horizon ", horizon_id)
  run_methods(get_case(prepared, disease_name, location, horizon_id, settings$feedback_delay))$scores
}))
rownames(all_scores) <- NULL
write_result(all_scores, "summary_scores_all_locations.csv")
write_result(do.call(rbind, lapply(runs, `[[`, "scores_by_level")), "summary_scores_by_level.csv")
write_result(do.call(rbind, lapply(runs, function(run) data.frame(location = run$case$location,
  scale_method = run$scaling$method, scale_used = run$scaling$used, count_scale = run$scaling$scale,
  transform, scale_as_of = run$case$scale_as_of, reference_weeks = run$scaling$reference_weeks,
  initialization_date = run$s0$date, warmup_weeks = run$s0$warmup_weeks,
  raw_error_scale = run$s0$raw_scale, error_scale = run$s0$scale,
  sigma_window = run$sigma$window, sigma_warmup_weeks = run$median_warmup_weeks,
  constant_sigma = run$constant_sigma, base_sigma_levels = paste(median_base_levels, collapse = "-"),
  median_start, median_rate))), "scales.csv")
write_result(do.call(rbind, lapply(runs, function(run)
  data.frame(location = run$case$location, run$rate_experts$experts))), "expert_grid.csv")
write_result(do.call(rbind, lapply(runs, function(run)
  data.frame(location = rep(run$case$location, nrow(run$restart_experts$experts)),
    run$restart_experts$experts, birth_week = run$restart_experts$birth_week))),
  "restart_experts.csv")

# 7. PLOT SETTINGS ---------------------------------------------------
residual_label <- paste("RMS of", median_window)
base_width <- paste0(round(100 * diff(median_base_levels)), "%")
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
  restart_grid_outer = "Starts x rate grid + outer",
  diverse_uniform = "Diverse experts: equal",
  diverse_shared = paste0("Diverse experts: learned + sharing (H=", H, ")"),
  diverse_shared_common = "Diverse experts: learned + sharing, common weights",
  studentized = paste0("Studentized median (", residual_label, ", ", median_start, " start, eta=",
    median_rate, ")"),
  studentized_base = paste0("Studentized median (base ", base_width, " width, normal start, eta=",
    median_rate, ")"),
  studentized_full = paste0("Studentized full base (base ", base_width, " width, eta=",
    median_rate, ")"),
  studentized_levels = paste0("Studentized full base (level scales, eta=", median_rate, ")"),
  studentized_level_rms = paste0("Studentized full base (level RMS of ", median_window, ", eta=",
    median_rate, ")"),
  median_constant = paste0("Median, constant sigma (", residual_label, ", ", median_start,
    " start, eta=", median_rate, ")"))
colours <- c(original = "gray40", fixed = "#0072B2", absolute_slope = "#009E73",
  oracle_target = "#CC79A7", uniform = "#E69F00", weighted = "#56B4E9",
  weighted_common = "#56B4E9", shared = "#882255", shared_common = "#882255",
  shared_outer = "#D55E00", shared_common_outer = "#D55E00",
  restart_persistent = "#0072B2", restart_uniform = "#E69F00", restart_weighted = "#009E73",
  restart_grid = "#882255", restart_uniform_outer = "#56B4E9",
  restart_weighted_outer = "#D55E00", restart_grid_outer = "#CC79A7",
  diverse_uniform = "#E69F00", diverse_shared = "#009E73", diverse_shared_common = "#009E73",
  studentized = "#D55E00", studentized_base = "#009E73", median_constant = "#CC79A7",
  studentized_full = "#D55E00", studentized_levels = "#009E73", studentized_level_rms = "#CC79A7")
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
    "restart_grid_outer"),
  "Diverse experts" = c("original", "shared", "shared_common", "diverse_uniform",
    "diverse_shared", "diverse_shared_common"),
  "Studentized median" = c("original", "fixed", "studentized", "studentized_base",
    "median_constant"),
  "Studentized full base" = c("original", "fixed", "studentized_full", "studentized_levels",
    "studentized_level_rms"))
score_units <- if (settings$scale_method == "none") "counts" else
  paste0("counts / ", settings$scale_method)
loss_axis_label <- paste0("Quantile loss (", score_units, ")")
mean_loss_axis_label <- paste("Mean", tolower(loss_axis_label))
window_titles <- c(overall = "Full observed history",
  before = paste0("Up to peak (weeks -", settings$weeks_before_peak - 1, " to 0)"),
  after = paste0("After peak (weeks 1 to ", settings$weeks_after_peak, ")"))
plot_title <- function(run = NULL) paste(disease_name, "| horizon", horizon_id, "|",
  if (is.null(run)) "all locations" else
    paste(run$case$location_name, paste0("(", run$case$location, ")")), "| full history |",
  if (transform == "log1p") "log1p(counts) updates" else "raw-count updates",
  "| feedback delay", settings$feedback_delay, "week(s)")
# Plots leave out a method whose identical twin is drawn in the same plot.
shown_methods <- function(run, methods) {
  twins <- run$duplicate_of[names(run$duplicate_of) %in% methods]
  setdiff(methods, names(twins)[twins %in% methods])
}
truth_label <- score_truth_label(settings$score_target, settings$feedback_delay)
truth_legend <- paste0(toupper(substr(truth_label, 1, 1)), substring(truth_label, 2))
date_ticks <- function(dates) pretty(as.Date(range(dates), origin = "1970-01-01"), n = 4)

# 8. FORECASTS AND WEEKLY QUANTILE LOSS PER LEVEL ----------------------
# Each level panel shows the scoring target (score_truth()) and every method's issued
# alpha-quantile; the panel below shows that level's weekly pinball loss.
blocks <- matrix(0L, 6, 3)
for (b in seq_along(plot_levels)) {
  r <- (b - 1L) %/% 3L
  blocks[2L * r + 1:2, (b - 1L) %% 3L + 1L] <- c(2L * b - 1L, 2L * b)
}
cairo_pdf(file.path(figure_dir, "multiqt_method_forecasts.pdf"), width = 18, height = 16, family = "DejaVu Sans")
for (run in runs) for (group_name in names(method_groups)) {
  method_names <- shown_methods(run, method_groups[[group_name]])
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
  # Column-wise widths, so long labels fit; trailing spaces keep a gap between columns.
  legend("center", paste0(c(truth_legend, labels[method_names]), "    "),
    col = c("black", colours[method_names]), lty = c(NA, line_types[method_names]),
    pch = c(16, rep(NA, length(method_names))), lwd = 2, ncol = 4, bty = "n", text.width = NA)
  mtext(paste(plot_title(run), "|", group_name), side = 3, outer = TRUE, line = 2, font = 2, cex = 1.15)
  mtext(paste("Issued alpha-quantile forecasts and their weekly pinball loss in", score_units,
    "against the", truth_label, "| Dotted line: peak"), side = 3, outer = TRUE, line = .5)
}
dev.off()

# 9. CALIBRATION ERROR AND QUANTILE LOSS ON MATCHED WINDOWS ----------
# One page: one row per location, one column per window, every method in each panel.
# A method is coloured by the first group that lists it: one colour family per group,
# one shade per method; common-weight and outer variants share their method's shade.
# CE = mean over levels of |coverage - alpha| over the window's weeks. Axes are
# shared by the three windows of a location; up to peak includes the peak week.
# Every method group needs a palette: methods of a group without one are not drawn.
group_palettes <- c("Single trackers" = "Blues 3", "Rate-grid experts" = "Oranges",
  "Weekly starts: immediate admission" = "Greens 3", "Diverse experts" = "Purples 3",
  "Studentized median" = "Burg", "Studentized full base" = "BluGrn")
home_group <- character()
for (group_name in names(method_groups)) for (method in method_groups[[group_name]])
  if (!method %in% c("original", names(home_group))) home_group[method] <- group_name
base_method <- sub("_outer$", "", sub("_common", "", names(home_group)))
tradeoff_colours <- c(original = "black")
for (group_name in names(group_palettes)) {
  members <- home_group == group_name
  bases <- unique(base_method[members])
  # Sequential palettes run dark to light; drop the two lightest shades.
  shades <- hcl.colors(length(bases) + 2L, group_palettes[[group_name]])[seq_along(bases)]
  tradeoff_colours[names(home_group)[members]] <- shades[match(base_method[members], bases)]
}
# Original: open circle. Filled circle: weights per level; triangle: common weights;
# square: + outer MultiQT; diamond: common weights + outer.
tradeoff_methods <- names(tradeoff_colours)
outer_variant <- grepl("_outer$", tradeoff_methods)
tradeoff_shapes <- setNames(ifelse(outer_variant, ifelse(common[tradeoff_methods], 18, 15),
  ifelse(common[tradeoff_methods], 17, 16)), tradeoff_methods)
tradeoff_shapes["original"] <- 1
tradeoff_sizes <- ifelse(tradeoff_methods == "original", 2.2, 1.6)
names(tradeoff_sizes) <- tradeoff_methods
# Duplicates (identical forecasts) are not drawn; those hidden at every location
# also leave the legend, and a footnote names each pair.
hidden_everywhere <- Reduce(intersect, lapply(runs, function(run)
  setdiff(tradeoff_methods, shown_methods(run, tradeoff_methods))))
hidden_pairs <- unique(unlist(lapply(runs, function(run) {
  hidden <- setdiff(tradeoff_methods, shown_methods(run, tradeoff_methods))
  if (length(hidden)) paste(labels[hidden], "=", labels[run$duplicate_of[hidden]])
})))
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

tradeoff_rows <- length(runs) + 1L          # First row: all prepared locations.
cairo_pdf(file.path(figure_dir, "multiqt_peak_reasoning_tradeoff.pdf"), width = 15,
  height = 3.3 * tradeoff_rows + 3, family = "DejaVu Sans")
layout(rbind(matrix(seq_len(3 * tradeoff_rows), tradeoff_rows, 3, byrow = TRUE),
  3 * tradeoff_rows + 1L), heights = c(rep(1, tradeoff_rows), .9))
par(oma = c(7.6, 0, 3, 0), mgp = c(2.4, .7, 0))
# First row: every prepared location as a small, light symbol per method (legend shapes),
# and the mean over locations as the usual symbol. The axes end at the 99th percentile of
# the single locations.
cloud <- all_scores[all_scores$method %in% setdiff(tradeoff_methods, hidden_everywhere), ]
cloud <- cloud[order(cloud$method == "original"), ]
xlim <- c(0, 1.05 * max(1e-8, quantile(cloud$quantile_loss, .99, na.rm = TRUE)))
ylim <- c(0, 1.05 * max(.01, quantile(cloud$coverage_error, .99, na.rm = TRUE)))
for (phase in phases) {
  rows <- cloud[cloud$phase == phase, ]
  means <- aggregate(cbind(quantile_loss, coverage_error) ~ method, rows, mean)
  means <- means[order(means$method == "original"), ]
  par(mar = c(3.8, 4, 3.4, 1))
  plot(NA, xlim = xlim, ylim = ylim, xlab = mean_loss_axis_label, ylab = "Calibration error",
    main = "")
  mtext(paste("All", length(all_locations), "locations |", window_titles[phase]), side = 3,
    line = 1.6, font = 2, cex = .85)
  mtext("Small: single locations; large: means; axes to the 99th percentile", side = 3,
    line = .3, cex = .75)
  grid(col = "gray90")
  points(rows$quantile_loss, rows$coverage_error, pch = tradeoff_shapes[rows$method],
    cex = .5 * tradeoff_sizes[rows$method],
    col = adjustcolor(tradeoff_colours[rows$method], alpha.f = .35))
  points(means$quantile_loss, means$coverage_error, col = tradeoff_colours[means$method],
    pch = tradeoff_shapes[means$method], cex = tradeoff_sizes[means$method],
    lwd = ifelse(means$method == "original", 2, 1))
}
for (run in runs) {
  scores <- run$scores[run$scores$method %in% shown_methods(run, tradeoff_methods), ]
  # Original last, so it stays visible on top of the other methods.
  scores <- scores[order(scores$method == "original"), ]
  xlim <- c(0, 1.1 * max(1e-8, scores$quantile_loss, na.rm = TRUE))
  ylim <- c(0, 1.1 * max(.01, scores$coverage_error, na.rm = TRUE))
  for (phase in phases) {
    rows <- scores[scores$phase == phase, ]
    par(mar = c(3.8, 4, 3.4, 1))
    plot(NA, xlim = xlim, ylim = ylim, xlab = mean_loss_axis_label, ylab = "Calibration error",
      main = "")
    mtext(paste(run$case$location_name, "|", window_titles[phase]), side = 3, line = 1.6,
      font = 2, cex = .85)
    mtext(paste(rows$n_weeks[1], "weeks"), side = 3, line = .3, cex = .75)
    grid(col = "gray90")
    points(rows$quantile_loss, rows$coverage_error, col = tradeoff_colours[rows$method],
      pch = tradeoff_shapes[rows$method], cex = tradeoff_sizes[rows$method],
      lwd = ifelse(rows$method == "original", 2, 1))
    if (!any(is.finite(rows$coverage_error)))
      text(mean(xlim), mean(ylim), "No observations inside the peak window", cex = .8)
  }
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", ifelse(legend_method, labels[legend_entries], legend_entries),
  col = ifelse(legend_method, tradeoff_colours[legend_entries], NA),
  pch = ifelse(legend_method, tradeoff_shapes[legend_entries], NA),
  pt.cex = ifelse(legend_method, tradeoff_sizes[legend_entries], 1),
  pt.lwd = ifelse(legend_entries == "original", 2, 1),
  text.font = ifelse(legend_header, 2, 1), ncol = length(legend_columns), bty = "n", cex = 1.05,
  text.width = NA)                         # Column-wise widths, so every group column fits.
mtext(plot_title(), side = 3, outer = TRUE, line = 1, font = 2, cex = 1.05)
mtext("Calibration error (CE): mean absolute coverage residual after averaging weeks within each quantile level.",
  side = 1, outer = TRUE, line = .5, cex = .8)
mtext(paste("Open circle: original forecasts. Filled circle: weights per quantile level; triangle: one",
  "weight vector shared by all levels; square / diamond: the same + outer MultiQT."),
  side = 1, outer = TRUE, line = 1.7, cex = .8)
if (length(hidden_pairs)) mtext(paste("Not drawn (identical forecasts):",
  paste(hidden_pairs, collapse = "; ")), side = 1, outer = TRUE, line = 2.9, cex = .8)
# Studentized median: what sigma_t is, and the first weeks without a residual sigma_t,
# which are the original base (e.g. most of the flu up-to-peak window).
mtext(paste0("Studentized median methods: sigma_t = RMS of the last ", median_window,
  " released median residuals, or the base's ", base_width, " interval width / ",
  round(diff(qnorm(median_base_levels)), 3), " (a normal sd)."), side = 1, outer = TRUE,
  line = 4.1, cex = .8)
median_warmup <- range(sapply(runs, function(run) run$sigma$warmup_weeks))
if (median_warmup[2] > 0) mtext(paste0("With residuals, the original base is issued in the first ",
  paste(unique(median_warmup), collapse = "-"), " weeks, until sigma_t exists; these weeks are ",
  "scored as the base."), side = 1, outer = TRUE, line = 5.3, cex = .8)
mtext(paste0("Studentized full base: offsets start at the base; level scales = normal sd of the ",
  "base's central interval at each level, or the RMS of each level's last ", median_window,
  " released residuals (weighted isotonic projection)."), side = 1, outer = TRUE, line = 6.5,
  cex = .8)
dev.off()

# 10. ISSUED WEIGHTS OVER TIME ---------------------------------------
# These are the weights used to mix forecasts before crossing repair or an outer
# correction. The outer methods keep their underlying ensemble's weights.
# Rows 0.1-0.9: one weight vector per level; last row: common weights for all levels.
# One page per location and expert bank; the base expert (m = 0) is gray.
weight_pages <- list(
  "Frozen rate-grid weights" = list(palette = "Viridis", columns = data.frame(
    method = c("uniform", "weighted", "shared"),
    common_method = c("uniform", "weighted_common", "shared_common"),
    title = c("Equal weights", "Learned, no sharing", "Learned + sharing"))),
  "Diverse expert weights" = list(palette = "Dark 3", columns = data.frame(
    method = c("diverse_uniform", "diverse_shared"),
    common_method = c("diverse_uniform", "diverse_shared_common"),
    title = c("Equal weights", "Learned + sharing"))))
plot_weight_distribution <- function(run, weights, fills, title, row_label, first_column, last_row) {
  times <- as.numeric(run$case$target_date)
  par(mar = c(if (last_row) 2.4 else .3, if (first_column) 3.7 else 2.2, 1.5, .6))
  plot(NA, xlim = range(times) + c(-.5, .5) * 7, ylim = c(0, 1),
    xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "",
    ylab = if (first_column) row_label else "", main = title)
  bottom <- numeric(length(times))
  for (j in seq_len(ncol(weights))) {
    top <- bottom + weights[, j]
    # Full-week bars, widened by 2% so neighbours overlap and no seam shows between them.
    # Filled rounds are drawn lighter.
    rect(times - .52 * 7, bottom, times + .52 * 7, top,
      col = ifelse(run$case$filled, adjustcolor(fills[j], .45), fills[j]), border = NA)
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
for (run in runs) for (page in names(weight_pages)) {
  columns <- weight_pages[[page]]$columns
  experts <- dimnames(run$weights[[columns$method[1]]])[[2]]
  fills <- c("gray45", hcl.colors(length(experts) - 1L, weight_pages[[page]]$palette))
  n_rows <- length(plot_levels) + 1L
  panel_ids <- matrix(seq_len(n_rows * nrow(columns)), ncol = nrow(columns), byrow = TRUE)
  layout(rbind(panel_ids, max(panel_ids) + 1L), heights = c(rep(1, n_rows), .8))
  par(oma = c(4.8, .8, 5.2, .8), mgp = c(2, .4, 0), cex = .9)
  for (row in seq_len(n_rows)) for (j in seq_len(nrow(columns))) {
    level_row <- row <= length(plot_levels)
    weights <- if (level_row) run$weights[[columns$method[j]]][, , match(plot_levels[row], run$case$alpha)] else
      run$weights[[columns$common_method[j]]][, , 1]
    stopifnot(identical(colnames(weights), experts))
    plot_weight_distribution(run, weights, fills, if (row == 1L) columns$title[j] else "",
      if (level_row) paste("alpha =", plot_levels[row]) else "Common\n(all levels)",
      first_column = j == 1L, last_row = row == n_rows)
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", experts, fill = fills, border = NA, ncol = ceiling(length(experts) / 2),
    bty = "n", cex = .95)
  mtext(page, side = 3, outer = TRUE, line = 2.4, font = 2, cex = 1.3)
  mtext(plot_title(run), side = 3, outer = TRUE, line = .9, cex = .9)
  mtext("Target date | Weight at issuance for each target | Dotted line: peak",
    side = 1, outer = TRUE, line = .4, cex = .85)
  mtext("Rows alpha = 0.1-0.9: one weight vector per quantile level. Last row: one weight vector for all levels, learned from the mean pinball loss.",
    side = 1, outer = TRUE, line = 1.6, cex = .8)
  if (any(run$case$filled)) mtext(paste("Lighter bars: rounds without a CMU forecast, filled with",
    "the latest earlier CMU forecast of the same target week."), side = 1, outer = TRUE,
    line = 2.8, cex = .8)
}
dev.off()
