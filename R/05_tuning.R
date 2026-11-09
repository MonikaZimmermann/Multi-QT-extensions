# Purpose: tune one setting per recalibration method on the tuning period of the
# design prepared by R/01_data_overview.R and save the choice for R/04_methods.R.
# Run from the project root after R/01_data_overview.R: Rscript R/05_tuning.R
# (about an hour on 18 cores). Order of the workflow: 00 -> 01 -> 05 -> 02 / 04.
# Every method runs online from the first round on, exactly as in 04; only the rows of
# the tuning period are scored (phase_weeks(..., period = "tuning")).
source("R/03_functions.R")

# 1. SETTINGS -----------------------------------------------------------------
settings <- multiqt_settings()
prepared <- readRDS("output/data/peak_aligned.rds")
design <- prepared$settings$design
diseases <- c("covid", "flu")
# Pooled over locations and these horizons. "seasons" design: horizon -1 uses the values
# tuned on 0-3 (the windowed_seasonal base has zero width there, and its 2023/24 forecasts
# used final counts); "hub_first_months": CMU's horizon -1 nowcasts are tuned as well.
tuning_horizons <- if (design == "seasons") 0:3 else settings$horizons
# Criterion per tracker, the scale it learns on (as 04 evaluates): the mean pinball loss
# of log(counts + 1) for log1p ("log", Bosse et al. 2023: relative errors) and of counts
# per 100k residents for identity ("natural"). Both are saved, with the calibration error.
tuning_score <- c(log1p = "log", identity = "natural")
criterion_of <- function(transform) ifelse(tuning_score[transform] == "log", "log_loss", "loss")
transforms <- c("log1p", "identity")       # Tuned for both trackers; 04 uses settings$transform.
# Feedback designs, each tuned separately (04 compares them): the delay in weeks after the
# first report, after which each forecast is learned once from its report current then, or
# the weekly re-run on the reports known then (feedback_protocol() in R/03_functions.R).
tuning_feedback <- data.frame(delay = c(0L, 1L, 0L), revisions = c(FALSE, FALSE, TRUE))
# One grid per tuned setting (tuned_setting in R/03_functions.R, plus the soft switch's
# threshold, width and signal weights); multiples of S0 except H, lambda, the switch
# settings and the studentized eta (units of sigma_t).
grids <- list(fixed_multiplier = c(.001, .0025, .005, .01, .02, .05, .1, .2, .4),
              slope_min_multiplier = c(.002, .004, .01, .02, .05, .1),
              median_rate = c(.003, .01, .03, .1, .3, 1),
              H = c(2L, 4L, 8L, 16L, 32L),
              restart_lambda = c(.25, .5, 1, 2.72, 7, 20),
              diverse_slow_multiplier = c(.001, .002, .004, .008, .02, .05, .1, .2),
              switch_strength = c(0, .25, .5, .75, 1),
              switch_threshold = c(.02, .05, .1, .2),
              switch_width = c(.02, .05, .1, .2),
              switch_weight_path = c(0, .5, 1),
              switch_weight_nowcast = c(0, .5, 1),
              switch_weight_revision = c(0, .5, 1),
              outer_multiplier = c(.001, .002, .004, .01, .02, .05, .1, .2))
# Stages, in order; each uses the values chosen in the earlier ones. First every
# method's own setting except the soft switch's and the outer rates; then the soft
# switch's settings one after another, on the tuned diverse common bank (strength 0 is
# that comparator); last the outer rates, on the tuned inner mixtures.
switch_settings <- c("switch_strength", "switch_threshold", "switch_width",
  "switch_weight_path", "switch_weight_nowcast", "switch_weight_revision")
stages <- c(list(setdiff(names(grids), c(switch_settings, "outer_multiplier"))),
  as.list(switch_settings), list("outer_multiplier"))
tuned_methods <- function(setting) if (setting %in% switch_settings) "diverse_soft_switch" else
  names(tuned_setting)[tuned_setting == setting]
# A value chosen on the edge of its grid is flagged in tuned_settings.csv (on_grid_edge).
cores <- max(1L, parallel::detectCores() - 2L)
result_dir <- "output/results/tuning"
figure_dir <- "output/pdf"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# 2. SCORE EVERY GRID VALUE ON THE TUNING SEASON ------------------------------
# Loss and calibration error of each method over the tuning rows of one case (NA
# without tuning rows).
score_tuning <- function(case, predictions, count_scale) {
  truth <- score_truth(case, settings$score_target, settings$score_filled)
  rows <- which(is.finite(truth) & phase_weeks(case, "overall", period = "tuning"))
  do.call(rbind, lapply(names(predictions), function(method) {
    q <- predictions[[method]] * count_scale
    on_log <- summarize_scores(score_quantiles(q, truth, case$alpha, score_scale = "log"),
      case$alpha, rows)
    natural <- summarize_scores(score_quantiles(q, truth, case$alpha, count_scale),
      case$alpha, rows)
    data.frame(method, n_weeks = length(rows), log_loss = on_log$quantile_loss,
      loss = natural$quantile_loss, coverage_error = natural$coverage_error)
  }))
}
# Every grid value of a stage's settings for every case of a disease, in parallel over
# cases, with one feedback design (a row of tuning_feedback; a re-run re-runs only the
# tuning weeks). by_method: the values chosen in earlier stages. Each case keeps one cache,
# so its scales, expert banks and mixtures are computed once for all grid values.
run_stage <- function(disease, transform, feedback, stage, by_method = list()) {
  run_settings <- settings
  run_settings$transform <- transform
  keys <- Filter(function(key) prepared$cases[[key]]$disease == disease &&
    prepared$cases[[key]]$horizon %in% tuning_horizons, names(prepared$cases))
  rows <- parallel::mclapply(keys, function(key) {
    prepared_case <- prepared$cases[[key]]
    case <- get_case(prepared, disease, prepared_case$location, prepared_case$horizon,
      feedback$delay, feedback$revisions)
    cache <- new.env()
    do.call(rbind, lapply(stage, function(setting) do.call(rbind, lapply(grids[[setting]],
      function(value) {
        methods <- tuned_methods(setting)
        params <- method_parameters()
        params$by_method <- by_method
        for (method in methods) params$by_method[[method]][[setting]] <- value
        run <- apply_methods(case, run_settings, params, methods, cache,
          weeks = which(case$period == "tuning"))
        data.frame(disease, location = case$location, horizon = case$horizon, transform,
          feedback_delay = feedback$delay, feedback_revisions = feedback$revisions,
          setting, value, score_tuning(case, run$predictions[c("original", methods)],
            run$scaling$scale))
      }))))
  }, mc.cores = cores)
  failed <- vapply(rows, inherits, logical(1), "try-error")
  if (any(failed)) stop(rows[[which(failed)[1]]])
  do.call(rbind, rows)
}
# The tuned value per method and setting: the smallest mean criterion over cases (ties:
# the smaller value, i.e. the slower rate, shorter lookback or weaker context dependence).
groups <- c("disease", "transform", "feedback_delay", "feedback_revisions")
choose_values <- function(scores) {
  scores <- scores[scores$method != "original", ]
  means <- aggregate(cbind(log_loss, loss, coverage_error) ~ disease + transform +
    feedback_delay + feedback_revisions + method + setting + value, scores, mean)
  criterion <- ifelse(criterion_of(means$transform) == "log_loss", means$log_loss, means$loss)
  means <- means[do.call(order, unname(c(means[c(groups, "method", "setting")],
    list(criterion, means$value)))), ]
  means[!duplicated(means[c(groups, "method", "setting")]), ]
}

all_scores <- list()
for (disease in diseases) for (transform in transforms) for (f in seq_len(nrow(tuning_feedback))) {
  feedback <- tuning_feedback[f, ]
  by_method <- list()
  for (stage in stages) {
    message("Tuning: ", disease, " | ", transform, " | feedback ",
      feedback_label(feedback$delay, feedback$revisions), " | ", paste(stage, collapse = ", "))
    scores <- run_stage(disease, transform, feedback, stage, by_method)
    all_scores[[length(all_scores) + 1L]] <- scores
    scores <- scores[is.finite(scores$log_loss), ]
    if (!nrow(scores)) {
      message("No tuning rows for ", disease, " (", transform, "): no further stages.")
      break
    }
    chosen <- choose_values(scores)
    for (i in seq_len(nrow(chosen)))
      by_method[[chosen$method[i]]][[chosen$setting[i]]] <- chosen$value[i]
  }
}
all_scores <- do.call(rbind, all_scores)
rownames(all_scores) <- NULL
# The feedback protocol (R/03_functions.R): tuned_parameters() uses no settings tuned under
# another one.
all_scores$feedback_protocol <- rep(feedback_protocol()$id, nrow(all_scores))
write.csv(all_scores, file.path(result_dir, paste0("tuning_scores_", design, ".csv")),
  row.names = FALSE)

# 3. CHOOSE AND SAVE ----------------------------------------------------------
# Cases without tuning rows (NA) drop out of the means; a disease without any keeps
# the defaults in 04.
scored <- all_scores[is.finite(all_scores$log_loss), ]
tuned <- choose_values(scored)
base <- aggregate(cbind(log_loss, loss, coverage_error) ~ disease + transform +
  feedback_delay + feedback_revisions, unique(scored[scored$method == "original",
    c(groups, "location", "horizon", "log_loss", "loss", "coverage_error")]), mean)
names(base)[5:7] <- paste0("original_", names(base)[5:7])
tuned <- merge(tuned, base, by = groups)
tuned$criterion <- unname(tuning_score[tuned$transform])
tuned$log_loss_change_percent <- 100 * (tuned$log_loss / tuned$original_log_loss - 1)
tuned$loss_change_percent <- 100 * (tuned$loss / tuned$original_loss - 1)
tuned$on_grid_edge <- mapply(function(setting, value) value %in% range(grids[[setting]]),
  tuned$setting, tuned$value)
tuned <- tuned[order(tuned$disease, tuned$transform, tuned$feedback_delay,
  tuned$feedback_revisions, match(tuned$method, method_names)), ]
rownames(tuned) <- NULL
tuned$feedback_protocol <- rep(feedback_protocol()$id, nrow(tuned))
write.csv(tuned, file.path(result_dir, paste0("tuned_settings_", design, ".csv")),
  row.names = FALSE)
print(tuned[c(groups, "method", "setting", "value", "on_grid_edge",
  "loss_change_percent", "log_loss_change_percent", "coverage_error",
  "original_coverage_error")], digits = 3,
  row.names = FALSE)

# 4. PLOT: CRITERION OVER EACH GRID -------------------------------------------
# One page per disease and feedback design; one panel per tuned setting; one line per
# method. Solid: log1p tracker (loss of log(counts + 1)), dashed: identity (loss per
# 100k); each tracker on its own criterion; the chosen value is circled.
means <- aggregate(cbind(log_loss, loss, coverage_error) ~ disease + transform +
  feedback_delay + feedback_revisions + method + setting + value, scored, mean)
cairo_pdf(file.path(figure_dir, "multiqt_tuning.pdf"), width = 15,
  height = 3.6 * ceiling(length(grids) / 4L) + 1.5, family = "DejaVu Sans")
for (disease in intersect(diseases, means$disease)) for (f in seq_len(nrow(tuning_feedback))) {
  feedback <- tuning_feedback[f, ]
  in_design <- function(x) x$disease == disease & x$feedback_delay == feedback$delay &
    x$feedback_revisions == feedback$revisions
  # Eight setting panels plus a separate legend strip.
  panel_rows <- ceiling(length(grids) / 4L)
  panel_ids <- matrix(seq_len(4L * panel_rows), panel_rows, 4, byrow = TRUE)
  panel_ids[panel_ids > length(grids)] <- 0L
  layout(rbind(panel_ids, rep(length(grids) + 1L, 4L)), heights = c(rep(1, panel_rows), .18))
  par(oma = c(2.5, 0, 3, 0), mar = c(4, 4.2, 3, 1), mgp = c(2.3, .6, 0))
  for (setting in names(grids)) {
    rows <- means[in_design(means) & means$setting == setting & means$method != "original", ]
    original <- base[in_design(base), ]
    original <- original[match(rows$transform, original$transform), ]
    rows$change <- 100 * ifelse(criterion_of(rows$transform) == "log_loss",
      rows$log_loss / original$original_log_loss, rows$loss / original$original_loss) - 100
    methods <- unique(rows$method)
    colours <- setNames(hcl.colors(max(3L, length(methods)), "Dark 3")[seq_along(methods)], methods)
    plot(NA, xlim = range(grids[[setting]]), ylim = range(0, rows$change, na.rm = TRUE),
      log = if (min(grids[[setting]]) > 0) "x" else "", xlab = setting,
      ylab = "Loss vs original (%)", main = setting)
    abline(h = 0, col = "gray60")
    for (method in methods) for (transform in transforms) {
      line <- rows[rows$method == method & rows$transform == transform, ]
      line <- line[order(line$value), ]
      lines(line$value, line$change, col = colours[method], lty = if (transform == "log1p") 1 else 2,
        type = "b", pch = 16, cex = .6)
      chosen <- tuned[in_design(tuned) & tuned$transform == transform &
        tuned$method == method, "value"]
      points(chosen, line$change[match(chosen, line$value)], col = colours[method], cex = 1.6)
    }
    legend("topleft", methods, col = colours, lty = 1, bty = "n", cex = .7)
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", c("log1p tracker (loss of log(counts + 1))",
    paste0("identity tracker (loss in ", count_scale_label(settings$scale_method), ")"),
    "chosen value"), lty = c(1, 2, NA),
    pch = c(NA, NA, 1), bty = "n", ncol = 3)
  mtext(paste(disease, "| design", design, "| feedback", feedback_label(feedback$delay,
    feedback$revisions), paste0("(", feedback_protocol()$id, ")"), "| tuning period",
    prepared$settings$tuning_label, "| horizons",
    paste(range(tuning_horizons), collapse = " to "), "| all locations"),
    side = 3, outer = TRUE, line = 1.2, font = 2)
  mtext(paste("Mean pinball loss relative to the original forecasts (mean over locations",
    "and horizons), on each tracker's scale"), side = 1, outer = TRUE, line = 1)
}
dev.off()
