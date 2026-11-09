# Apply the recalibration methods of R/03_functions.R to each location and plot them.
# Run R/01_data_overview.R once to prepare output/data/peak_aligned.rds.
# Run from the project root: Rscript R/04_methods.R
# Offsets and ensemble weights are fitted independently for each location.
source("R/03_functions.R")

# 1. SETTINGS AND PREPARED DATA -------------------------------------
disease_name <- "covid"                    # "covid" or "flu"
horizon_id <- 2                            # -1, 0, 1, 2 or 3: forecast and weight plots.
settings <- multiqt_settings()             # Locations, tracker, experts and windows (R/03_functions.R).
transform <- settings$transform
# Feedback designs of the trade-off pages (all horizons): the delay in weeks after the first
# report, and the weekly re-run on the reports known then (get_case() in R/03_functions.R).
tradeoff_feedback <- data.frame(delay = c(0L, 1L, 0L), revisions = c(FALSE, FALSE, TRUE))
cores <- max(1L, parallel::detectCores() - 2L)
plot_levels <- (1:9) / 10
phases <- c("overall", "before", "after")
prepared <- readRDS("output/data/peak_aligned.rds")
# The methods' settings: method_parameters() in R/03_functions.R, with the values tuned
# on the tuning period by R/05_tuning.R for this disease, tracker, design and feedback.
tuned_for <- function(delay, revisions) tuned_parameters(disease_name, transform,
  design = prepared$settings$design, feedback_delay = delay, revisions = revisions)
params <- tuned_for(settings$feedback_delay, settings$feedback_revisions)
seasons <- c("all", prepared$settings$evaluation_seasons)   # Scored: each season and all.
result_dir <- "output/results/peak_reasoning"
figure_dir <- "output/pdf"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write_result <- function(rows, filename)
  write.csv(rows, file.path(result_dir, filename), row.names = FALSE)
par_of <- function(method, name) method_parameter(params, method, name)

# 2. RUN EACH LOCATION INDEPENDENTLY --------------------------------
# apply_methods() (R/03_functions.R) issues every method's forecasts from the online
# feedback only. Peaks are local maxima of the latest-reported counts (R/01);
# only the target oracle uses that hindsight date to choose forecasts. Alignment
# changes evaluation coordinates, never when feedback becomes available. Every
# location starts fresh: offsets, expert weights, lambda and loss history.
# issue_weeks: the weeks a re-run case re-runs (apply_methods(); default all).
run_methods <- function(case, run_params = params, issue_weeks = NULL) {
  run <- apply_methods(case, settings, run_params, weeks = issue_weeks)
  predictions <- run$predictions

  # SCORE AGAINST THE CHOSEN TARGET -------------------------------------
  # The online updates saw only the feedback released so far: the report
  # available feedback_delay weeks after the first report (re-run: as known at each
  # issue week). Scores use score_truth(): the latest vintage by default.
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
    feedback_delay = case$feedback_delay, feedback_revisions = case$feedback_revised)
  c(run, list(truth = truth, weekly = weekly, duplicate_of = duplicate_of,
    scores = data.frame(id, do.call(rbind, scores), row.names = NULL),
    scores_by_level = data.frame(id, do.call(rbind, scores_by_level), row.names = NULL)))
}
runs <- list()
for (location in settings$location_ids) {
  case <- get_case(prepared, disease_name, location, horizon_id, settings$feedback_delay,
    settings$feedback_revisions)
  check_peak_window(case, settings$weeks_before_peak, settings$weeks_after_peak,
    score_truth(case, settings$score_target, settings$score_filled))
  message("Methods: ", disease_name, " | ", location, " | horizon ", horizon_id)
  runs[[location]] <- run_methods(case)
}
saveRDS(runs, file.path(result_dir, "issued_forecasts.rds"))
# Scores of every prepared location, horizon and feedback design for the trade-off pages
# (section 9), in parallel; only the scores are kept. Each design uses the settings tuned
# with it; a re-run re-runs only the evaluation weeks. The runs above are reused.
all_locations <- unique(unlist(lapply(prepared$cases, function(case)
  if (case$disease == disease_name) case$location)))
tradeoff_params <- Map(tuned_for, tradeoff_feedback$delay, tradeoff_feedback$revisions)
tradeoff_cases <- expand.grid(location = all_locations, horizon = settings$horizons,
  feedback = seq_len(nrow(tradeoff_feedback)), stringsAsFactors = FALSE)
message("Trade-off scores: ", disease_name, " | ", length(all_locations), " locations x ",
  length(settings$horizons), " horizons x ", nrow(tradeoff_feedback), " feedback designs")
all_scores <- parallel::mclapply(seq_len(nrow(tradeoff_cases)), function(i) {
  location <- tradeoff_cases$location[i]
  horizon <- tradeoff_cases$horizon[i]
  k <- tradeoff_cases$feedback[i]
  delay <- tradeoff_feedback$delay[k]
  revisions <- tradeoff_feedback$revisions[k]
  if (location %in% names(runs) && horizon == horizon_id && delay == settings$feedback_delay &&
    revisions == settings$feedback_revisions) return(runs[[location]]$scores)
  case <- get_case(prepared, disease_name, location, horizon, delay, revisions)
  run_methods(case, tradeoff_params[[k]], which(case$period == "evaluation"))$scores
}, mc.cores = cores)
failed <- vapply(all_scores, inherits, logical(1), "try-error")
if (any(failed)) stop(all_scores[[which(failed)[1]]])
all_scores <- do.call(rbind, all_scores)
rownames(all_scores) <- NULL
write_result(all_scores[all_scores$location %in% settings$location_ids, ], "summary_scores.csv")
write_result(all_scores, "summary_scores_all_locations.csv")
write_result(do.call(rbind, lapply(runs, `[[`, "scores_by_level")), "summary_scores_by_level.csv")
write_result(do.call(rbind, lapply(runs, function(run) data.frame(location = run$case$location,
  scale_method = run$scaling$method, scale_used = run$scaling$used, count_scale = run$scaling$scale,
  transform, scale_as_of = run$case$scale_as_of, reference_weeks = run$scaling$reference_weeks,
  initialization_date = run$s0$date, warmup_weeks = run$s0$warmup_weeks,
  raw_error_scale = run$s0$raw_scale, error_scale = run$s0$scale,
  sigma_window = run$sigma$window, sigma_warmup_weeks = run$median_warmup_weeks,
  constant_sigma = run$constant_sigma,
  base_sigma_levels = paste(params$median_base_levels, collapse = "-"),
  median_start = params$median_start))), "scales.csv")
write_result(do.call(rbind, lapply(runs, function(run)
  data.frame(location = run$case$location, run$rate_experts$experts))), "expert_grid.csv")
write_result(do.call(rbind, lapply(runs, function(run)
  data.frame(location = rep(run$case$location, nrow(run$restart_experts$experts)),
    run$restart_experts$experts, birth_week = run$restart_experts$birth_week))),
  "restart_experts.csv")
# The context and both expert distributions are saved at issuance, alongside the
# actual mixed weights. The gate describes declining context, not a peak probability.
switch_runs <- Filter(function(run) !is.null(run$soft_switch), runs)
if (length(switch_runs)) {
  write_result(do.call(rbind, lapply(switch_runs, function(run) data.frame(
    disease = run$case$disease, location = run$case$location,
    location_name = run$case$location_name, horizon = run$case$horizon,
    target_date = run$case$target_date, run$soft_switch$gate,
    lambda_baseline = run$soft_switch$regime_lambda[, "baseline"],
    lambda_declining = run$soft_switch$regime_lambda[, "declining"],
    # A gate-weighted summary of the two learning rates, not a third update rate.
    lambda_average = as.numeric(run$lambda$diverse_soft_switch)))), "soft_switch_gate.csv")
  write_result(do.call(rbind, lapply(switch_runs, function(run) {
    regimes <- run$soft_switch$regime_weights
    weights <- list(baseline = regimes[, , "baseline"],
      declining = regimes[, , "declining"], issued = run$weights$diverse_soft_switch[, , 1])
    do.call(rbind, lapply(names(weights), function(context) {
      values <- weights[[context]]
      data.frame(disease = run$case$disease, location = run$case$location,
        location_name = run$case$location_name, horizon = run$case$horizon,
        reference_date = rep(run$case$reference_date, ncol(values)),
        target_date = rep(run$case$target_date, ncol(values)), context,
        expert = rep(colnames(values), each = nrow(values)), weight = as.vector(values))
    }))
  })), "soft_switch_weights.csv")
}

# 7. PLOT SETTINGS ---------------------------------------------------
# Labels show each method's own value of its tuned setting (tuned_setting in R/03).
residual_label <- paste("RMS of", params$median_window)
base_width <- paste0(round(100 * diff(params$median_base_levels)), "%")
value <- function(method) par_of(method, tuned_setting[[method]])
eta <- function(method) paste0("eta=", value(method))
outer <- function(method) paste0("outer m=", value(method))
H_of <- function(method) paste0("H=", value(method))
lambda_of <- function(method) paste0("lambda=", value(method))
slow <- function(method) paste0("slow m=", value(method))
labels <- c(original = "Original",
  fixed = paste0("Fixed MultiQT (eta=", value("fixed"), "*S0)"),
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
  restart_persistent = paste0("First start only (m=", params$restart_multiplier, ")"),
  restart_uniform = paste0("Weekly starts: equal (m=", params$restart_multiplier, ")"),
  restart_weighted = paste0("Weekly starts: learned (m=", params$restart_multiplier, ", ",
    lambda_of("restart_weighted"), ")"),
  restart_grid = paste0("Starts x rate grid: learned (", lambda_of("restart_grid"), ")"),
  restart_uniform_outer = paste0("Starts equal + outer (", outer("restart_uniform_outer"), ")"),
  restart_weighted_outer = paste0("Starts learned + outer (", outer("restart_weighted_outer"), ")"),
  restart_grid_outer = paste0("Starts x rate grid + outer (", outer("restart_grid_outer"), ")"),
  diverse_uniform = paste0("Diverse experts: equal (", slow("diverse_uniform"), ")"),
  diverse_shared = paste0("Diverse experts: learned + sharing (", slow("diverse_shared"), ", ",
    "H=", par_of("diverse_shared", "H"), ")"),
  diverse_shared_common = paste0("Diverse experts: learned + sharing, common weights (",
    slow("diverse_shared_common"), ")"),
  diverse_soft_switch = paste0("Diverse experts: soft switching (strength=",
    value("diverse_soft_switch"), ")"),
  studentized = paste0("Studentized median (", residual_label, ", ", params$median_start,
    " start, ", eta("studentized"), ")"),
  studentized_base = paste0("Studentized median (base ", base_width, " width, normal start, ",
    eta("studentized_base"), ")"),
  studentized_full = paste0("Studentized full base (base ", base_width, " width, ",
    eta("studentized_full"), ")"),
  studentized_levels = paste0("Studentized full base (level scales, ", eta("studentized_levels"),
    ")"),
  studentized_level_rms = paste0("Studentized full base (level RMS of ", params$median_window,
    ", ", eta("studentized_level_rms"), ")"),
  median_constant = paste0("Median, constant sigma (", residual_label, ", ", params$median_start,
    " start, ", eta("median_constant"), ")"))
colours <- c(original = "gray40", fixed = "#0072B2", absolute_slope = "#009E73",
  oracle_target = "#CC79A7", uniform = "#E69F00", weighted = "#56B4E9",
  weighted_common = "#56B4E9", shared = "#882255", shared_common = "#882255",
  shared_outer = "#D55E00", shared_common_outer = "#D55E00",
  restart_persistent = "#0072B2", restart_uniform = "#E69F00", restart_weighted = "#009E73",
  restart_grid = "#882255", restart_uniform_outer = "#56B4E9",
  restart_weighted_outer = "#D55E00", restart_grid_outer = "#CC79A7",
  diverse_uniform = "#E69F00", diverse_shared = "#009E73", diverse_shared_common = "#009E73",
  diverse_soft_switch = "#7B3294",
  studentized = "#D55E00", studentized_base = "#009E73", median_constant = "#CC79A7",
  studentized_full = "#D55E00", studentized_levels = "#009E73", studentized_level_rms = "#CC79A7")
# Common weights across quantile levels: dashed lines and triangles.
common <- setNames(grepl("_common", names(labels)), names(labels))
common["diverse_soft_switch"] <- TRUE
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
    "diverse_shared", "diverse_shared_common", "diverse_soft_switch"),
  "Studentized median" = c("original", "fixed", "studentized", "studentized_base",
    "median_constant"),
  "Studentized full base" = c("original", "fixed", "studentized_full", "studentized_levels",
    "studentized_level_rms"))
score_units <- count_scale_label(settings$scale_method)
loss_axis_label <- paste0("Quantile loss (", score_units, ")")
mean_loss_axis_label <- paste("Mean", tolower(loss_axis_label))
evaluation_label <- prepared$settings$evaluation_label
window_titles <- c(overall = paste0(toupper(substr(evaluation_label, 1, 1)),
  substring(evaluation_label, 2)),
  before = paste0("Up to peak (weeks -", settings$weeks_before_peak - 1, " to 0)"),
  after = paste0("After peak (weeks 1 to ", settings$weeks_after_peak, ")"))
plot_title <- function(run = NULL) paste(disease_name, "| horizon", horizon_id, "|",
  if (is.null(run)) "all locations" else
    paste(run$case$location_name, paste0("(", run$case$location, ")")), "| evaluation:",
  evaluation_label, "|",
  if (transform == "log1p") "log1p(counts) updates" else "raw-count updates",
  "| feedback", feedback_label(settings$feedback_delay, settings$feedback_revisions))
# Plots leave out a method whose identical twin is drawn in the same plot.
shown_methods <- function(run, methods) {
  twins <- run$duplicate_of[names(run$duplicate_of) %in% methods]
  setdiff(methods, names(twins)[twins %in% methods])
}
truth_label <- score_truth_label(settings$score_target, settings$feedback_delay)
truth_legend <- paste0(toupper(substr(truth_label, 1, 1)), substring(truth_label, 2))
date_ticks <- function(dates) pretty(as.Date(range(dates), origin = "1970-01-01"), n = 6)
# Background of the time plots: the tuning period shaded, dotted lines at the peaks.
mark_time <- function(case, col = "gray60") {
  tuning <- case$target_date[case$period == "tuning"]
  usr <- par("usr")
  if (length(tuning)) rect(min(tuning) - 3.5, usr[3], max(tuning) + 3.5, usr[4], col = "gray93",
    border = NA)
  abline(v = as.numeric(case$peaks$peak_date), col = col, lty = 3)
}

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
  mtext(paste("Issued alpha-quantile forecasts and their weekly pinball loss in", score_units,
    "against the", truth_label, "| Dotted lines: peaks | Gray: tuning period (not scored)"),
    side = 3, outer = TRUE, line = .5)
}
dev.off()

# 9. CALIBRATION ERROR AND QUANTILE LOSS ON MATCHED WINDOWS ----------
# One page per location, after a first page with all prepared locations. Columns:
# horizons; rows: the feedback designs of tradeoff_feedback, each with its three windows.
# Every method in each panel. The loss is scored on the tracker's scale: log(counts + 1)
# for the log1p tracker, counts per 100k residents for raw counts. Each feedback design
# uses the settings tuned (R/05_tuning.R) with it.
# A method is coloured by the first group that lists it: one colour family per group,
# one shade per method; common-weight and outer variants share their method's shade.
# CE = mean over levels of |coverage - alpha| over the window's weeks. Every panel has
# its own axis ranges (panel_limits()); up to peak includes the peak week.
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
tradeoff_scores <- all_scores[all_scores$season == "all" & all_scores$method %in% tradeoff_methods, ]
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
# Short names for the panel headings: "delay d" or "re-run".
feedback_names <- with(tradeoff_feedback, ifelse(revisions,
  ifelse(delay > 0, paste("re-run, delay", delay), "re-run"), paste("delay", delay)))
panel_height <- 2.6                         # Inches per row of panels.
legend_height <- 2.6
tradeoff_title <- function(location = NULL) paste(disease_name, "|",
  if (is.null(location)) paste("all", length(all_locations), "locations") else
    paste0(tradeoff_scores$location_name[match(location, tradeoff_scores$location)],
      " (", location, ")"), "| evaluation:", evaluation_label, "|",
  if (transform == "log1p") "log1p(counts) updates, loss of log(counts + 1)" else
    paste("raw-count updates, loss in", score_units))
# One page; location NULL: all prepared locations. Every panel has its own axis ranges;
# a loss beyond the axis is drawn at its edge as ">".
draw_tradeoff_page <- function(location = NULL) {
  all_page <- is.null(location)
  page <- tradeoff_scores[if (all_page) !tradeoff_scores$method %in% hidden_everywhere else
    tradeoff_scores$location == location & !tradeoff_scores$hidden, ]
  page <- page[order(page$method == "original"), ]   # Original last, on top of the others.
  n_panels <- nrow(panel_rows) * length(settings$horizons)
  layout(rbind(matrix(seq_len(n_panels), nrow(panel_rows), byrow = TRUE), n_panels + 1L),
    heights = c(rep(panel_height, nrow(panel_rows)), legend_height))
  par(oma = c(9, 3.4, 7.8, 0), mgp = c(2.2, .6, 0))
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
  legend("center", ifelse(legend_method, labels[legend_entries], legend_entries),
    col = ifelse(legend_method, tradeoff_colours[legend_entries], NA),
    pch = ifelse(legend_method, tradeoff_shapes[legend_entries], NA),
    pt.cex = ifelse(legend_method, tradeoff_sizes[legend_entries], 1),
    pt.lwd = ifelse(legend_entries == "original", 2, 1),
    text.font = ifelse(legend_header, 2, 1), ncol = length(legend_columns), bty = "n", cex = 1.1,
    text.width = NA)                         # Column-wise widths, so every group column fits.
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
  mtext(paste("Rows: feedback design (delay d: learn from the report d weeks after the",
    "first one; re-run: re-run every week on the reports known then) and window; columns:",
    "horizon."), side = 3, outer = TRUE, line = 3.9, cex = 1.05)
  if (all_page) mtext(paste("Small: single locations; large: means over locations;",
    "axes: 1st to 99th percentile of both."), side = 3, outer = TRUE, line = 2.5, cex = 1.05)
  notes <- c("Calibration error (CE): mean absolute coverage residual after averaging weeks within each quantile level.",
    paste("Open circle: original forecasts. Filled circle: weights per quantile level; triangle: one",
      "weight vector shared by all levels; square / diamond: the same + outer MultiQT."),
    if (length(hidden_pairs)) paste("Not drawn (identical forecasts):",
      paste(hidden_pairs, collapse = "; ")),
    paste0("Each feedback design uses the settings tuned with it; the values in the legend are ",
      "those of feedback ", feedback_label(settings$feedback_delay, settings$feedback_revisions),
      " (all: output/results/tuning)."),
    # Studentized median: what sigma_t is, and the first weeks without a residual sigma_t,
    # which are the original base (e.g. most of the flu up-to-peak window).
    paste0("Studentized median methods: sigma_t = RMS of the last ", params$median_window,
      " released median residuals, or the base's ", base_width, " interval width / ",
      round(diff(qnorm(params$median_base_levels)), 3), " (a normal sd)."),
    paste("With residuals, the original base is issued until sigma_t exists; these weeks are",
      "scored as the base."),
    paste0("Studentized full base: offsets start at the base; level scales = normal sd of the ",
      "base's central interval at each level, or the RMS of each level's last ",
      params$median_window, " released residuals (weighted isotonic projection)."))
  mtext(notes, side = 1, outer = TRUE, line = .6 + 1.3 * (seq_along(notes) - 1), cex = .95)
}
cairo_pdf(file.path(figure_dir, "multiqt_peak_reasoning_tradeoff.pdf"), width = 19,
  height = nrow(panel_rows) * panel_height + legend_height + 2.5, family = "DejaVu Sans")
draw_tradeoff_page()
for (location in settings$location_ids) draw_tradeoff_page(location)
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
plot_weight_distribution <- function(run, weights, fills, title, row_label, first_column, last_row,
                                     dates = run$case$target_date) {
  times <- as.numeric(dates)
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
  abline(v = as.numeric(run$case$peaks$peak_date), col = "gray30", lty = 3)
  axis(2, at = c(0, .5, 1), las = 1, cex.axis = .75)
  if (last_row) axis.Date(1, at = date_ticks(dates), format = "%b %Y", cex.axis = .8)
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
  mtext("Target date | Weight at issuance for each target | Dotted lines: peaks",
    side = 1, outer = TRUE, line = .4, cex = .85)
  mtext("Rows alpha = 0.1-0.9: one weight vector per quantile level. Last row: one weight vector for all levels, learned from the mean pinball loss.",
    side = 1, outer = TRUE, line = 1.6, cex = .8)
  if (any(run$case$filled)) mtext(paste("Lighter bars: rounds without a CMU forecast, filled with",
    "the latest earlier CMU forecast of the same target week."), side = 1, outer = TRUE,
    line = 2.8, cex = .8)
}
dev.off()

# 11. SOFT SWITCHING: CONTEXT AND COMMON WEIGHTS ---------------------
# Issue dates align the signal with the information used for that forecast. Each
# regime learns from delayed own-horizon losses using the context saved at issuance.
if (length(switch_runs)) {
  cairo_pdf(file.path(figure_dir, "multiqt_soft_switch_timelines.pdf"), width = 15,
    height = 11, family = "DejaVu Sans")
  for (run in switch_runs) {
    signal <- run$soft_switch$gate
    regimes <- run$soft_switch$regime_weights
    dates <- run$case$reference_date
    experts <- dimnames(regimes)[[2]]
    fills <- c("gray45", hcl.colors(length(experts) - 1L, "Dark 3"))
    layout(matrix(1:5, ncol = 1), heights = c(1.1, 1, 1, 1, .45))
    par(oma = c(4.4, .5, 4.1, .5), mgp = c(2, .4, 0), cex = .9)
    par(mar = c(.4, 3.7, 2.1, .6))
    plot(dates, signal$gate, type = "n", xlim = range(dates) + c(-.5, .5) * 7,
      ylim = c(0, 1.18), xaxs = "i", xaxt = "n", yaxt = "n", xlab = "",
      ylab = "Context weight", main = "Declining context at issuance")
    axis(2, at = c(0, .5, 1), las = 1, cex.axis = .8)
    abline(h = .5, col = "gray70", lty = 3)
    abline(v = as.numeric(run$case$peaks$peak_date), col = "gray50", lty = 3)
    score_names <- c("path_score", "nowcast_score", "revision_score")
    score_colours <- c("#0072B2", "#D55E00", "#009E73")
    broken <- break_gaps(dates, as.matrix(signal[c(score_names, "gate")]))
    for (j in seq_along(score_names)) lines(broken$dates, broken[[2]][, j],
      col = score_colours[j], lwd = 1)
    lines(broken$dates, broken[[2]][, "gate"], col = colours["diverse_soft_switch"], lwd = 2)
    points(dates, signal$gate, col = colours["diverse_soft_switch"], pch = 16, cex = .3)
    legend("topleft", c("Forecast path", "Nowcast change", "Same-target revisions", "Mixed gate"),
      col = c(score_colours, colours["diverse_soft_switch"]), lty = 1, lwd = c(1, 1, 1, 2),
      ncol = 4, bty = "n", cex = .85)
    plot_weight_distribution(run, regimes[, , "baseline"], fills,
      "Baseline-context expert weights", "Weight", TRUE, FALSE, dates)
    plot_weight_distribution(run, regimes[, , "declining"], fills,
      "Declining-context expert weights", "Weight", TRUE, FALSE, dates)
    plot_weight_distribution(run, run$weights$diverse_soft_switch[, , 1], fills,
      "Issued mixture: (1 - gate) x baseline + gate x declining", "Weight", TRUE, TRUE, dates)
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", experts, fill = fills, border = NA, ncol = ceiling(length(experts) / 2),
      bty = "n", cex = .9)
    mtext(labels["diverse_soft_switch"], side = 3, outer = TRUE, line = 2.2, font = 2, cex = 1.2)
    mtext(plot_title(run), side = 3, outer = TRUE, line = .8, cex = .9)
    mtext("Issue date | Common expert weights across quantile levels | Dotted lines: hindsight peaks",
      side = 1, outer = TRUE, line = .4, cex = .85)
    mtext("The gate is a context weight, not a peak probability; missing signals contribute a neutral score of 0.5.",
      side = 1, outer = TRUE, line = 1.6, cex = .8)
    if (any(run$case$filled)) mtext("Lighter bars: filled forecast rounds.",
      side = 1, outer = TRUE, line = 2.8, cex = .8)
  }
  dev.off()
}
