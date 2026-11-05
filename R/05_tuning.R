# Purpose: tune one setting per recalibration method on the tuning season (2023/24)
# and save the choice for R/04_methods.R.
# Run from the project root after R/01_data_overview.R: Rscript R/05_tuning.R
# (about 45 minutes on 18 cores). Order of the workflow: 00 -> 01 -> 05 -> 02 / 04.
# Every method runs online from the first round on, exactly as in 04; only the rows of
# the tuning season are scored (2023/24: phase_weeks(..., period = "tuning")).
source("R/03_functions.R")

# 1. SETTINGS -----------------------------------------------------------------
settings <- multiqt_settings()
diseases <- c("covid", "flu")
# Pooled over locations and these horizons. Horizon -1 uses the values tuned here: the
# windowed_seasonal base has zero width there, and its 2023/24 forecasts used final counts.
tuning_horizons <- 0:3
# Criterion: the mean pinball loss of log(counts + 1) ("log", Bosse et al. 2023:
# relative errors, comparable across locations and seasons) or of counts per 100k
# residents ("natural"). Both are saved, with the calibration error.
tuning_score <- "log"
transforms <- c("log1p", "identity")       # Tuned for both trackers; 04 uses settings$transform.
# One grid per tuned setting (tuned_setting in R/03_functions.R); multiples of S0 except
# H, lambda and the studentized eta (units of sigma_t). The outer rates come last, on
# the tuned inner mixtures.
grids <- list(fixed_multiplier = c(.001, .0025, .005, .01, .02, .05, .1, .2, .4),
              slope_min_multiplier = c(.002, .004, .01, .02, .05, .1),
              median_rate = c(.003, .01, .03, .1, .3, 1),
              H = c(2L, 4L, 8L, 16L, 32L),
              restart_lambda = c(.25, .5, 1, 2.72, 7, 20),
              diverse_slow_multiplier = c(.001, .002, .004, .008, .02, .05, .1, .2),
              outer_multiplier = c(.001, .002, .004, .01, .02, .05, .1, .2))
# A value chosen on the edge of its grid is flagged in tuned_settings.csv (on_grid_edge).
cores <- max(1L, parallel::detectCores() - 2L)
prepared <- readRDS("output/data/peak_aligned.rds")
result_dir <- "output/results/tuning"
figure_dir <- "output/pdf"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# 2. SCORE EVERY GRID VALUE ON THE TUNING SEASON ------------------------------
# Loss and calibration error of each method over the tuning rows of one case (NA
# without tuning rows, e.g. covid with the CMU hub base).
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
# One grid value of one setting, for every case of a disease, in parallel. by_method:
# the settings already tuned (the inner mixtures of the outer stage).
run_grid_value <- function(disease, transform, setting, value, by_method = list()) {
  run_settings <- settings
  run_settings$transform <- transform
  methods <- names(tuned_setting)[tuned_setting == setting]
  params <- method_parameters()
  params$by_method <- by_method
  for (method in methods) params$by_method[[method]][[setting]] <- value
  keys <- Filter(function(key) prepared$cases[[key]]$disease == disease &&
    prepared$cases[[key]]$horizon %in% tuning_horizons, names(prepared$cases))
  rows <- parallel::mclapply(keys, function(key) {
    prepared_case <- prepared$cases[[key]]
    case <- get_case(prepared, disease, prepared_case$location, prepared_case$horizon,
      run_settings$feedback_delay)
    run <- apply_methods(case, run_settings, params, methods)
    data.frame(disease, location = case$location, horizon = case$horizon, transform,
      setting, value, score_tuning(case, run$predictions[c("original", methods)],
        run$scaling$scale))
  }, mc.cores = cores)
  failed <- vapply(rows, inherits, logical(1), "try-error")
  if (any(failed)) stop(rows[[which(failed)[1]]])
  do.call(rbind, rows)
}
# The tuned value per method: the smallest mean criterion over cases (ties: the
# smaller value, i.e. the slower rate or shorter lookback).
choose_values <- function(scores) {
  scores <- scores[scores$method != "original", ]
  criterion <- if (tuning_score == "log") "log_loss" else "loss"
  means <- aggregate(cbind(log_loss, loss, coverage_error) ~ disease + transform + method +
    setting + value, scores, mean)
  means <- means[order(means$disease, means$transform, means$method, means[[criterion]],
    means$value), ]
  means[!duplicated(means[c("disease", "transform", "method")]), ]
}

all_scores <- list()
for (disease in diseases) for (transform in transforms) {
  for (setting in setdiff(names(grids), "outer_multiplier")) for (value in grids[[setting]]) {
    message("Tuning: ", disease, " | ", transform, " | ", setting, " = ", value)
    all_scores[[length(all_scores) + 1L]] <- run_grid_value(disease, transform, setting, value)
  }
  inner <- do.call(rbind, all_scores)
  inner <- inner[inner$disease == disease & inner$transform == transform & is.finite(inner$log_loss), ]
  if (!nrow(inner)) {
    message("No tuning rows for ", disease, " (", transform, "): no outer stage.")
    next
  }
  tuned_inner <- choose_values(inner)
  by_method <- list()
  for (i in seq_len(nrow(tuned_inner)))
    by_method[[tuned_inner$method[i]]][[tuned_inner$setting[i]]] <- tuned_inner$value[i]
  for (value in grids$outer_multiplier) {
    message("Tuning: ", disease, " | ", transform, " | outer_multiplier = ", value)
    all_scores[[length(all_scores) + 1L]] <- run_grid_value(disease, transform,
      "outer_multiplier", value, by_method)
  }
}
all_scores <- do.call(rbind, all_scores)
rownames(all_scores) <- NULL
write.csv(all_scores, file.path(result_dir, "tuning_scores.csv"), row.names = FALSE)

# 3. CHOOSE AND SAVE ----------------------------------------------------------
# Cases without tuning rows (NA) drop out of the means; a disease without any keeps
# the defaults in 04.
scored <- all_scores[is.finite(all_scores$log_loss), ]
tuned <- choose_values(scored)
base <- aggregate(cbind(log_loss, loss, coverage_error) ~ disease + transform,
  unique(scored[scored$method == "original", c("disease", "location", "horizon", "transform",
    "log_loss", "loss", "coverage_error")]), mean)
names(base)[3:5] <- paste0("original_", names(base)[3:5])
tuned <- merge(tuned, base, by = c("disease", "transform"))
tuned$criterion <- tuning_score
tuned$log_loss_change_percent <- 100 * (tuned$log_loss / tuned$original_log_loss - 1)
tuned$on_grid_edge <- mapply(function(setting, value) value %in% range(grids[[setting]]),
  tuned$setting, tuned$value)
tuned <- tuned[order(tuned$disease, tuned$transform, match(tuned$method, method_names)), ]
rownames(tuned) <- NULL
write.csv(tuned, file.path(result_dir, "tuned_settings.csv"), row.names = FALSE)
print(tuned[c("disease", "transform", "method", "setting", "value", "on_grid_edge",
  "log_loss_change_percent", "coverage_error", "original_coverage_error")], digits = 3,
  row.names = FALSE)

# 4. PLOT: CRITERION OVER EACH GRID -------------------------------------------
# One page per disease; one panel per tuned setting; one line per method. Solid: log1p
# tracker, dashed: identity; the chosen value is circled.
means <- aggregate(cbind(log_loss, loss, coverage_error) ~ disease + transform + method +
  setting + value, scored, mean)
criterion <- if (tuning_score == "log") "log_loss" else "loss"
criterion_label <- if (tuning_score == "log") "log(counts + 1)" else
  count_scale_label(settings$scale_method)
cairo_pdf(file.path(figure_dir, "multiqt_tuning.pdf"), width = 15, height = 9.5,
  family = "DejaVu Sans")
for (disease in intersect(diseases, means$disease)) {
  layout(rbind(matrix(1:8, 2, 4, byrow = TRUE)))
  par(oma = c(2.5, 0, 3, 0), mar = c(4, 4.2, 3, 1), mgp = c(2.3, .6, 0))
  for (setting in names(grids)) {
    rows <- means[means$disease == disease & means$setting == setting &
      means$method != "original", ]
    rows$change <- 100 * (rows[[criterion]] / base[[paste0("original_", criterion)]][
      match(paste(rows$disease, rows$transform), paste(base$disease, base$transform))] - 1)
    methods <- unique(rows$method)
    colours <- setNames(hcl.colors(max(3L, length(methods)), "Dark 3")[seq_along(methods)], methods)
    plot(NA, xlim = range(grids[[setting]]), ylim = range(0, rows$change, na.rm = TRUE),
      log = if (min(grids[[setting]]) > 0) "x" else "", xlab = setting,
      ylab = paste0("Loss vs original (%), ", tuning_score), main = setting)
    abline(h = 0, col = "gray60")
    for (method in methods) for (transform in transforms) {
      line <- rows[rows$method == method & rows$transform == transform, ]
      line <- line[order(line$value), ]
      lines(line$value, line$change, col = colours[method], lty = if (transform == "log1p") 1 else 2,
        type = "b", pch = 16, cex = .6)
      chosen <- tuned[tuned$disease == disease & tuned$transform == transform &
        tuned$method == method, "value"]
      points(chosen, line$change[match(chosen, line$value)], col = colours[method], cex = 1.6)
    }
    legend("topleft", methods, col = colours, lty = 1, bty = "n", cex = .7)
  }
  plot.new()
  legend("center", c("log1p tracker", "identity tracker", "chosen value"), lty = c(1, 2, NA),
    pch = c(NA, NA, 1), bty = "n")
  mtext(paste(disease, "| tuning season", paste(prepared$settings$tuning_seasons, collapse = ", "),
    "| horizons", paste(range(tuning_horizons), collapse = " to "), "| all locations"),
    side = 3, outer = TRUE, line = 1.2, font = 2)
  mtext(paste("Mean pinball loss of", criterion_label, "relative to the original forecasts",
    "(mean over locations and horizons)"), side = 1, outer = TRUE, line = 1)
}
dev.off()
