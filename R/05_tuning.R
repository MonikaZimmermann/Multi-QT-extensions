# Purpose: tune the tracker rates on the tuning period of the data prepared by
# R/01_data_overview.R and save the choice for R/04_methods.R, each by the elbow of its loss
# and calibration-error curve (see choose_values()). Everything else is fixed
# (method_parameters()), also the look-back H of the learned mixtures.
# Run from the project root after R/01_data_overview.R: Rscript R/05_tuning.R
# (about 15 minutes on 18 cores). Order of the workflow: 00 -> 01 -> 05 -> 02 / 04 -> 06.
# Every finished stage is saved (cached() in R/03_functions.R): a restarted run, after a
# crash or a change, recomputes only the stages whose settings, data or code changed.
# Every method runs online from the first round on, exactly as in 04; only the rows of
# the tuning period are scored (phase_weeks(..., period = "tuning")).
source("R/03_functions.R")

# 1. SETTINGS -----------------------------------------------------------------
settings <- multiqt_settings()
prepared <- read_prepared(settings)       # The base model's prepared data and output folder.
base_model <- prepared$settings$base_model
diseases <- c("covid", "flu")
# Pooled over the locations except the US (settings$summary_exclude, as in the summaries of
# 02, 04 and 06; the US uses the values tuned on the others) and these horizons of each
# disease, set by its base: CMU's horizon -1 nowcasts are tuned with the others; the
# windowed_seasonal forecasts have zero width at horizon -1 once the week before is
# reported (the latest report is a predictor), so their horizon -1 uses the values tuned
# on 0-3.
tuning_horizons <- function(disease) if (prepared$settings$base[[disease]] == "windowed_seasonal")
  setdiff(settings$horizons, -1L) else settings$horizons
# Criterion per tracker, the scale it learns on (as 04 evaluates): the mean pinball loss
# of log(counts + 1) for log1p ("log", Bosse et al. 2023: relative errors) and of counts
# per 100k residents for identity ("natural"). Both are saved, with the calibration error.
tuning_score <- c(log1p = "log", identity = "natural")
criterion_of <- function(transform) ifelse(tuning_score[transform] == "log", "log_loss", "loss")
# Trackers tuned: identity for every tuned method (04 uses settings$transform), log1p only
# for fixed and leaky, whose values 04 uses for its methods on the log1p tracker
# (log_tracker_methods in R/03).
transforms <- c("log1p", "identity")
# Feedback designs, each tuned separately (04 compares them): the delay in weeks after the
# first report, after which each forecast is learned once from its report current then, or
# the weekly re-run on the reports known then (feedback_protocol() in R/03_functions.R).
tuning_feedback <- data.frame(delay = c(0L, 1L, 0L), revisions = c(FALSE, FALSE, TRUE))
# One grid per tuned setting (tuned_setting in R/03_functions.R): the tracker rates
# (rate_settings; multiples of S0, the studentized eta in units of sigma_t). Not tuned:
# H of the learned mixtures (rate grid, weekly starts, score experts; method_parameters()),
# lambda (adaptive from lambda_start = 2.72), sharing 1/(2H), the outer rate, the number of
# recent weekly starts. The median-start
# methods use the rate of their method (parameter_source in R/03).
grids <- list(fixed_multiplier = c(.001, .0025, .005, .01, .02, .05, .1, .2, .4),
              median_rate = c(.003, .01, .03, .1, .3, 1),
              leaky_multiplier = c(.001, .002, .004, .008, .02, .05, .1, .2, .3, .5, 1))
# One stage, the tracker rates. H is fixed for now: its tuning curves showed no structure.
# A later stage (e.g. "H", with a grid above and its methods in tuned_setting) would be
# tuned with the rates chosen before it (by_method).
stages <- list(rate_settings)
tuned_methods <- function(setting) names(tuned_setting)[tuned_setting == setting]
# Elbow: with the loss (criterion) and the calibration error both as % change
# from the original forecasts, the grid value with the smallest sum of the two changes, i.e.
# the largest calibration gain (%) minus loss increase (%): beyond it, 1 % more calibration
# costs more than 1 % loss. Ties: the smaller value.
# A value chosen on the edge of its grid is flagged in tuned_settings.csv (on_grid_edge).
cores <- max(1L, parallel::detectCores() - 2L)
result_dir <- file.path(settings$output_dir, "results", "tuning")
figure_dir <- file.path(settings$output_dir, "pdf")
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
  # The log1p tracker serves only standard MultiQT (fixed and leaky, for fixed_log and
  # leaky_log: log_tracker_methods in R/03); the studentized methods track counts.
  methods_of <- function(setting) if (transform == "log1p")
    intersect(tuned_methods(setting), log_tracker_methods) else tuned_methods(setting)
  stage <- Filter(function(setting) length(methods_of(setting)) > 0L, stage)
  keys <- Filter(function(key) prepared$cases[[key]]$disease == disease &&
    prepared$cases[[key]]$horizon %in% tuning_horizons(disease) &&
    !prepared$cases[[key]]$location %in% settings$summary_exclude, names(prepared$cases))
  rows <- parallel::mclapply(keys, function(key) {
    prepared_case <- prepared$cases[[key]]
    case <- get_case(prepared, disease, prepared_case$location, prepared_case$horizon,
      feedback$delay, feedback$revisions)
    cache <- new.env()
    do.call(rbind, lapply(stage, function(setting) do.call(rbind, lapply(grids[[setting]],
      function(value) {
        methods <- methods_of(setting)
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
# Means over cases of each method, setting and value, with the criterion of the tracker, that
# of the original forecasts (original_criterion), both changes in % from the original
# (loss_change, ce_change) and their sum (elbow: the smallest is chosen).
groups <- c("disease", "transform", "feedback_delay", "feedback_revisions")
original_means <- function(scores) {
  base <- aggregate(cbind(log_loss, loss, coverage_error) ~ disease + transform +
    feedback_delay + feedback_revisions, unique(scores[scores$method == "original",
      c(groups, "location", "horizon", "log_loss", "loss", "coverage_error")]), mean)
  names(base)[5:7] <- paste0("original_", names(base)[5:7])
  base
}
curve_means <- function(scores) {
  means <- aggregate(cbind(log_loss, loss, coverage_error) ~ disease + transform +
    feedback_delay + feedback_revisions + method + setting + value,
    scores[scores$method != "original", ], mean)
  means <- merge(means, original_means(scores), by = groups)
  log_scale <- criterion_of(means$transform) == "log_loss"
  means$criterion <- ifelse(log_scale, means$log_loss, means$loss)
  means$original_criterion <- ifelse(log_scale, means$original_log_loss, means$original_loss)
  means$loss_change <- 100 * (means$criterion / means$original_criterion - 1)
  means$ce_change <- 100 * (means$coverage_error / means$original_coverage_error - 1)
  means$elbow <- means$loss_change + means$ce_change
  means
}
# The tuned value per method and setting: the smallest elbow sum (ties: the smaller value).
choose_values <- function(scores) {
  means <- curve_means(scores)
  means <- means[do.call(order, unname(c(means[c(groups, "method", "setting")],
    list(means$elbow, means$value)))), ]
  means[!duplicated(means[c(groups, "method", "setting")]), ]
}

all_scores <- list()
for (disease in diseases) for (transform in transforms) for (f in seq_len(nrow(tuning_feedback))) {
  feedback <- tuning_feedback[f, ]
  by_method <- list()
  for (stage in stages) {
    message("Tuning: ", disease, " | ", transform, " | feedback ",
      feedback_label(feedback$delay, feedback$revisions), " | ", paste(stage, collapse = ", "))
    scores <- cached("05_tuning", run_stage, disease, transform, feedback, stage, by_method)
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
write.csv(all_scores, file.path(result_dir, paste0("tuning_scores_", base_model, ".csv")),
  row.names = FALSE)

# 3. CHOOSE AND SAVE ----------------------------------------------------------
# Cases without tuning rows (NA) drop out of the means; a disease without any keeps
# the defaults in 04.
scored <- all_scores[is.finite(all_scores$log_loss), ]
tuned <- choose_values(scored)
tuned <- tuned[setdiff(names(tuned), c("criterion", "original_criterion"))]
tuned$criterion <- unname(tuning_score[tuned$transform])
tuned$rule <- "elbow: smallest loss change + calibration error change (% of the original)"
tuned$log_loss_change_percent <- 100 * (tuned$log_loss / tuned$original_log_loss - 1)
tuned$loss_change_percent <- 100 * (tuned$loss / tuned$original_loss - 1)
tuned$on_grid_edge <- mapply(function(setting, value) value %in% range(grids[[setting]]),
  tuned$setting, tuned$value)
tuned <- tuned[order(tuned$disease, tuned$transform, tuned$feedback_delay,
  tuned$feedback_revisions, match(tuned$method, method_names)), ]
rownames(tuned) <- NULL
tuned$feedback_protocol <- rep(feedback_protocol()$id, nrow(tuned))
write.csv(tuned, file.path(result_dir, paste0("tuned_settings_", base_model, ".csv")),
  row.names = FALSE)
print(tuned[c(groups, "method", "setting", "value", "on_grid_edge",
  "loss_change_percent", "log_loss_change_percent", "coverage_error",
  "original_coverage_error")], digits = 3,
  row.names = FALSE)

# 4. ELBOW PLOTS: CALIBRATION ERROR AGAINST LOSS ALONG EACH GRID --------------
# One page per disease and feedback design, one panel per tuned setting: for each method and
# tracker the path of its grid values in the plane of loss and calibration error, both as %
# change from the base forecasts (the base: open square at 0, 0), the symbols growing
# with the value, as in 02's trade-off plots. Solid: log1p tracker (loss of log(counts + 1)),
# dashed: identity (loss per 100k); each tracker on its own criterion. The elbow rule
# (choose_values()) picks the value with the smallest sum of both changes: circled, with the
# dotted line of equal sum through it (slope -1), below which no grid value lies. Lower left
# is better. Only the pool the values are chosen on (the criterion); no pages per location.
panel_columns <- 3L
panel_rows <- ceiling(length(grids) / panel_columns)
value_size <- function(value, grid) .5 + 1.4 * (log(value) - log(min(grid))) /
  diff(log(range(grid)))
draw_tuning_pages <- function(scored) {
  means <- curve_means(scored)
  for (disease in intersect(diseases, means$disease)) for (f in seq_len(nrow(tuning_feedback))) {
    feedback <- tuning_feedback[f, ]
    in_design <- function(x) x$disease == disease & x$feedback_delay == feedback$delay &
      x$feedback_revisions == feedback$revisions
    panel_ids <- matrix(seq_len(panel_columns * panel_rows), panel_rows, panel_columns,
      byrow = TRUE)
    panel_ids[panel_ids > length(grids)] <- 0L
    layout(rbind(panel_ids, rep(length(grids) + 1L, panel_columns)),
      heights = c(rep(1, panel_rows), .2))
    par(oma = c(2.2, 0, 3, 0), mar = c(4, 4.2, 3.4, 1), mgp = c(2.3, .6, 0))
    for (setting in names(grids)) {
      rows <- means[in_design(means) & means$setting == setting, ]
      methods <- unique(rows$method)
      colours <- setNames(hcl.colors(max(3L, length(methods)), "Dark 3")[seq_along(methods)], methods)
      plot(NA, xlim = range(0, rows$loss_change, na.rm = TRUE),
        ylim = range(0, rows$ce_change, na.rm = TRUE), xlab = "Loss vs base (%)",
        ylab = "Calibration error vs base (%)", main = setting)
      abline(h = 0, v = 0, col = "gray80")
      points(0, 0, pch = 0, cex = 1.6, lwd = 1.5, col = "gray40")
      for (method in methods) for (transform in transforms) {
        line <- rows[rows$method == method & rows$transform == transform, ]
        if (!nrow(line)) next
        line <- line[order(line$value), ]
        lines(line$loss_change, line$ce_change, col = colours[method],
          lty = if (transform == "log1p") 1 else 2)
        points(line$loss_change, line$ce_change, col = colours[method], pch = 16,
          cex = value_size(line$value, grids[[setting]]))
        chosen <- match(tuned[in_design(tuned) & tuned$transform == transform &
          tuned$method == method & tuned$setting == setting, "value"], line$value)
        if (is.na(chosen)) next
        points(line$loss_change[chosen], line$ce_change[chosen], col = colours[method],
          cex = 1.2 + value_size(line$value[chosen], grids[[setting]]), lwd = 1.5)
        abline(a = line$elbow[chosen], b = -1, col = colours[method], lty = 3)
      }
      legend("topright", methods, col = colours, lty = 1, bty = "n", cex = .7)
      mtext(paste("symbols grow with the value:", paste(format(range(grids[[setting]])),
        collapse = " to ")), side = 3, line = .2, cex = .65)
    }
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", c("log1p tracker (loss of log(counts + 1))",
      paste0("identity tracker (loss in ", count_scale_label(settings$scale_method), ")"),
      "chosen value (circled)", "equal sum of both changes through it",
      "base forecasts (0, 0)", "symbol size: grid value (larger = larger)"),
      lty = c(1, 2, NA, 3, NA, NA), pch = c(NA, NA, 1, NA, 0, 16),
      pt.cex = c(1, 1, 1.8, 1, 1.6, 1.4), col = c(1, 1, 1, 1, "gray40", "gray60"),
      bty = "n", ncol = 3, cex = .9)
    mtext(paste(disease, "| base", base_model, "| feedback", feedback_label(feedback$delay,
      feedback$revisions), "| tuning period",
      prepared$settings$tuning_label, "| horizons",
      paste(range(tuning_horizons(disease)), collapse = " to "), "|",
      pooled_label(length(unique(scored$location[scored$disease == disease])), settings)),
      side = 3, outer = TRUE, line = 1.2, font = 2)
    mtext(paste("Calibration error against mean quantile loss, both relative to the base",
      "forecasts (means over locations and horizons) on the scale of each tracker;",
      "lower left is better"), side = 1, outer = TRUE, line = .8)
  }
}
cairo_pdf(file.path(figure_dir, "multiqt_tuning.pdf"), width = 15,
  height = 4.6 * panel_rows + 1.6, family = "DejaVu Sans")
draw_tuning_pages(scored)
dev.off()
