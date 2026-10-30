# Source R/02_multiqt_ties.R first, then this script in the same R session.
# Run both from the project root; shared settings come from script 02.
# Offsets and ensemble weights are fitted independently for each location.
# Tracker coordinates follow use_log_scale from script 02; losses/mixtures use scaled counts.

# 1. SETTINGS AND PREPARED DATA -------------------------------------
# scale_method controls loss normalization and expert weights, not tracker coordinates.
# Inherit the switch instead of maintaining a second setting that can disagree with 02.
tracker_scale <- if (use_log_scale) "log1p_counts" else "counts"
to_tracker <- if (use_log_scale) log1p else identity
from_tracker <- if (use_log_scale) expm1 else identity
tracker_method <- if (use_log_scale) "log-MultiQT" else "count-MultiQT"
prediction_scale <- "counts / count_scale"
combined_weeks <- seq.int(-weeks_before_peak, weeks_after_peak)
# Baseline eta = reference_multiplier * S0, using the same frozen scale as the bank.
reference_multiplier <- .05
# Other rates and caps below are in tracker units (log1p counts or raw counts).
outer_rate <- .005
slope_rate_min <- .005
slope_rate_max <- .1                       # Tracker-unit cap; independent of the multiplier grid.
# Shared bank settings come from script 02; edit the grid there for both scripts.
# Defaults also support sessions that supply the selectors and use_log_scale directly.
H <- 8L                                  # Fixed lookback for expert-weight updates.
plot_levels <- (1:9) / 10
lambda_start <- 2.72
# Weekly restart comparison: each copy votes from birth, initially forecasting the base.
restart_multiplier <- .05
restart_multipliers <- c(.01, .05, .2)       # Smaller grid at each new start.
restart_weight_rate <- 2.72                 # Fixed meta-rate on scaled-count pinball loss.
restart_share <- 1 / (2 * H)                # Recovery mass within eligible historical issuers.
prepared <- readRDS("output/data/peak_aligned.rds")
selected_cases <- Filter(function(case) case$disease == disease_name &&
  case$horizon == horizon_id && if (combined) case$location != "US" else
    case$location == location_id, prepared$cases)
if (!length(selected_cases)) stop("No prepared cases match disease, location and horizon.")

# A common set of peak-relative weeks gives every retained location the same weight.
# Individual mode retains the full history, even for locations missing a peak-window week.
location_support <- do.call(rbind, lapply(names(selected_cases), function(case_id) {
  case <- selected_cases[[case_id]]
  indices <- match(combined_weeks, case$relative_week)
  observed <- !is.na(indices)
  observed[observed] <- is.finite(case$truth_latest[indices[observed]]) &
    rowSums(is.finite(case$base[indices[observed], , drop = FALSE])) == length(case$alpha) &
    case$within_winter[indices[observed]]
  data.frame(case_id, location = case$location, location_name = case$location_name,
    peak_date = case$peak_date, boundary_peak = case$peak_boundary,
    expected_weeks = length(combined_weeks),
    observed_weeks = sum(observed), complete = all(observed),
    missing_relative_weeks = paste(combined_weeks[!observed], collapse = ","))
}))
if (combined) {
  # A maximum on the winter edge is likely an earlier wave's tail, not a winter peak.
  selected_cases <- selected_cases[location_support$case_id[location_support$complete &
    !location_support$boundary_peak]]
  if (!length(selected_cases))
    stop("No complete combined windows. Reduce weeks_before_peak/weeks_after_peak or inspect an individual location.")
}
n_locations <- length(selected_cases)
output_dir <- if (combined) file.path("output/results/combined",
  paste0(disease_name, "_h", horizon_id), "03") else "output/results/peak_reasoning"
figure_dir <- if (combined) file.path("output/pdf/combined",
  paste0(disease_name, "_h", horizon_id), "03") else "output/pdf"
figure_file <- file.path(figure_dir, "multiqt_peak_reasoning.pdf")
overall_figure_file <- file.path(figure_dir, "multiqt_peak_reasoning_overall.pdf")
tradeoff_figure_file <- file.path(figure_dir, "multiqt_peak_reasoning_tradeoff.pdf")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write_result <- function(rows, filename)
  write.csv(rows, file.path(output_dir, filename), row.names = FALSE)
write_result(location_support, "location_support.csv")

# Restart copies share the frozen error scale and vote immediately at birth.
# There is no separate base expert: each untrained copy forecasts the base.
# Births are actual issue dates; forecasts before birth stay NA.
run_restart_experts <- function(base_counts, feedback_counts, feedback_by_week,
                                reference_date, alpha, multipliers, rate_scale,
                                warmup_weeks, count_scale = 1) {
  n <- nrow(base_counts); a <- length(alpha)
  stopifnot(n >= 1L, ncol(base_counts) == a, length(feedback_counts) == n,
    length(feedback_by_week) == n, length(reference_date) == n,
    all(is.finite(base_counts)), all(base_counts >= 0),
    is.finite(count_scale), count_scale > 0,
    length(multipliers) > 0, all(is.finite(multipliers)), all(multipliers > 0),
    !anyDuplicated(multipliers), warmup_weeks >= 0, warmup_weeks <= n)
  released <- unlist(feedback_by_week, use.names = FALSE)
  stopifnot(!anyDuplicated(released), all(is.finite(feedback_counts[released])),
    all(feedback_counts[released] >= 0))
  for (week in seq_len(n)) stopifnot(all(feedback_by_week[[week]] < week),
    all(feedback_by_week[[week]] >= 1))
  births <- if (warmup_weeks < n) seq.int(warmup_weeks + 1L, n) else integer()
  if (length(births)) stopifnot(is.finite(rate_scale), rate_scale > 0)
  birth_week <- rep(births, each = length(multipliers))
  multiplier <- rep(multipliers, times = length(births))
  type <- rep("restart", length(birth_week))
  expert_names <- if (length(birth_week))
    paste0("start=", reference_date[birth_week], ";m=", multiplier) else character()
  k <- length(birth_week)
  forecasts <- array(NA_real_, c(n, k, a),
    dimnames = list(as.character(reference_date), expert_names, as.character(alpha)))
  own_feedback <- matrix(0L, n, k, dimnames = dimnames(forecasts)[1:2])
  admission_week <- birth_week
  base_tracker <- to_tracker(base_counts); y_tracker <- to_tracker(feedback_counts)
  for (expert in seq_len(k)) {
    birth <- birth_week[expert]; eta <- multiplier[expert] * rate_scale
    offsets <- numeric(a); q_tracker <- matrix(NA_real_, n, a); seen <- 0L
    for (week in seq.int(birth, n)) {
      # Ignore delayed reports for rows this expert never forecast.
      past <- feedback_by_week[[week]]
      past <- past[past >= birth]
      for (row in past) offsets <- offsets + eta * (alpha - (y_tracker[row] <= q_tracker[row, ]))
      seen <- seen + length(past)
      own_feedback[week, expert] <- seen
      q_tracker[week, ] <- isoreg(base_tracker[week, ] + offsets)$yf
    }
    issued <- seq.int(birth, n)
    q <- from_tracker(q_tracker[issued, , drop = FALSE])
    unchanged <- q_tracker[issued, , drop = FALSE] == base_tracker[issued, , drop = FALSE]
    q[unchanged] <- base_counts[issued, , drop = FALSE][unchanged]
    stopifnot(all(is.finite(q)), !use_log_scale || all(q > -1))
    forecasts[issued, expert, ] <- q / count_scale
  }
  experts <- data.frame(expert = expert_names, type, birth_week,
    birth_date = reference_date[birth_week], multiplier, eta = multiplier * rate_scale,
    admission_week, admission_date = reference_date[admission_week])
  list(forecasts = forecasts, experts = experts, feedback_counts = own_feedback)
}

# A practical delayed specialist-style mixture. This is a separate weighting
# experiment, not a claim of the published sleeping-expert regret guarantee.
# Losses are count-scale pinball losses, with a fixed meta-rate for clarity.
run_restart_weights <- function(bank, base, feedback_observed, feedback_by_week,
                                alpha, learn = TRUE, weight_rate = 2.72, share = .125) {
  n <- nrow(base); a <- length(alpha); k <- nrow(bank$experts)
  stopifnot(identical(dim(bank$forecasts), c(n, k, a)),
    identical(dim(bank$feedback_counts), c(n, k)),
    all(bank$experts$admission_week == bank$experts$birth_week),
    is.finite(weight_rate), weight_rate >= 0, is.finite(share), share > 0, share <= 1)
  p <- matrix(0, k, a)
  weights <- array(0, c(n, k, a), dimnames = dimnames(bank$forecasts))
  q <- base
  diagnostics <- data.frame(week = seq_len(n), n_started = 0L, n_eligible = 0L,
    n_waiting = 0L, n_admitted = 0L, n_weight_updates = 0L,
    mean_age_weeks = NA_real_, effective_experts = NA_real_, largest_weight = NA_real_,
    n_untrained = 0L, untrained_weight = 0)
  for (week in seq_len(n)) {
    incumbents <- which(bank$experts$birth_week < week)
    if (learn) for (past in feedback_by_week[[week]]) {
      # Every historical issuer is eligible, including on its first released loss.
      eligible <- incumbents[bank$experts$birth_week[incumbents] <= past]
      if (length(eligible) < 2L) next
      past_q <- matrix(bank$forecasts[past, eligible, ], nrow = length(eligible), ncol = a)
      error <- feedback_observed[past] - past_q
      loss <- sweep(pmax(error, 0), 2, alpha, `*`) +
        sweep(pmax(-error, 0), 2, 1 - alpha, `*`)
      stopifnot(all(is.finite(loss)))
      # Normalize within historical issuers, preserving all other experts' mass.
      mass <- colSums(p[eligible, , drop = FALSE])
      conditional <- sweep(p[eligible, , drop = FALSE], 2, mass, `/`)
      loss <- sweep(loss, 2, apply(loss, 2, min), `-`)
      log_p <- log(conditional) - weight_rate * loss
      log_p <- sweep(log_p, 2, apply(log_p, 2, max), `-`)
      conditional <- exp(log_p)
      conditional <- sweep(conditional, 2, colSums(conditional), `/`)
      conditional <- (1 - share) * conditional + share / length(eligible)
      p[eligible, ] <- sweep(conditional, 2, mass, `*`)
      diagnostics$n_weight_updates[week] <- diagnostics$n_weight_updates[week] + 1L
    }
    newcomers <- which(bank$experts$birth_week == week)
    active <- c(incumbents, newcomers)
    if (learn && length(newcomers)) {
      # Each newborn gets 1/K among K active copies, before issuing its first forecast.
      # Preserve the incumbents' relative weights after processing released losses.
      new_mass <- length(newcomers) / length(active)
      p <- p * (1 - new_mass)
      p[newcomers, ] <- new_mass / length(newcomers)
    }
    # Before the first birth (or without initialization), retain the exact base.
    if (!learn && length(active)) p[active, ] <- 1 / length(active)
    if (length(active)) {
      stopifnot(all(is.finite(p)), max(abs(colSums(p) - 1)) < 1e-10)
      current <- matrix(bank$forecasts[week, active, ], nrow = length(active), ncol = a)
      # Weight first, then repair crossings. Preserve an all-base mixture exactly.
      q[week, ] <- if (all(current == rep(base[week, ], each = length(active)))) base[week, ] else
        isoreg(colSums(current * p[active, , drop = FALSE]))$yf
      issue_date <- as.Date(dimnames(bank$forecasts)[[1]][week])
      age <- as.numeric(issue_date - bank$experts$birth_date[active]) / 7
      diagnostics$mean_age_weeks[week] <- sum(age * rowMeans(p[active, , drop = FALSE]))
      diagnostics$effective_experts[week] <- mean(1 / colSums(p^2))
      diagnostics$largest_weight[week] <- mean(apply(p, 2, max))
      untrained <- active[bank$feedback_counts[week, active] == 0L]
      diagnostics$n_untrained[week] <- length(untrained)
      diagnostics$untrained_weight[week] <- mean(colSums(p[untrained, , drop = FALSE]))
    }
    weights[week, , ] <- p
    diagnostics$n_started[week] <- sum(bank$experts$birth_week <= week)
    diagnostics$n_eligible[week] <- length(active)
    diagnostics$n_waiting[week] <- diagnostics$n_started[week] - length(active)
    diagnostics$n_admitted[week] <- length(newcomers)
  }
  list(q = q, weights = weights, diagnostics = diagnostics)
}

# 2. RUN EACH LOCATION INDEPENDENTLY --------------------------------
# Peaks are the largest latest-reported count in the prepared winter window.
# Only the target oracle uses that hindsight date to choose forecasts.
# Alignment changes evaluation coordinates, never when feedback becomes available.
# Every location starts fresh: offsets, expert weights, lambda and loss history.
extension_runs <- list()
error_rows <- list()
scale_rows <- list()
expert_grid_rows <- list()
restart_expert_rows <- restart_diagnostic_rows <- list()
for (case_id in names(selected_cases)) {
  case <- selected_cases[[case_id]]
  alpha <- case$alpha
  scale_history <- case$scale_history[is.finite(case$scale_history)]
  count_scale <- 1
  scale_used <- scale_method
  if (scale_method != "none" && length(scale_history)) {
    count_scale <- if (scale_method == "median") median(scale_history) else max(scale_history)
    if (count_scale <= 0) {
      count_scale <- max(scale_history)
      scale_used <- "max"
    }
  }
  if (!is.finite(count_scale) || count_scale <= 0 || !length(scale_history)) {
    count_scale <- 1
    scale_used <- "none"
  }
  if (any(!is.finite(case$base)) || any(case$base < 0))
    stop("Base counts must be finite and nonnegative for case ", case_id, ".")
  available_feedback <- !is.na(case$feedback_observed)
  if (any(!is.finite(case$feedback_observed[available_feedback])) ||
      any(case$feedback_observed[available_feedback] < 0))
    stop("Available feedback counts must be finite and nonnegative for case ", case_id, ".")
  feedback_indices <- unlist(case$feedback_by_week, use.names = FALSE)
  if (any(!is.finite(case$feedback_observed[feedback_indices])))
    stop("Scheduled feedback must have finite observations for case ", case_id, ".")
  base_tracker <- to_tracker(case$base)
  feedback_tracker <- to_tracker(case$feedback_observed)
  # Keep stored predictions and the scoring/weighting inputs in counts / count_scale.
  base <- case$base / count_scale
  feedback_observed <- case$feedback_observed / count_scale
  truth_latest <- case$truth_latest / count_scale  # Scoring target; updates use feedback only.
  feedback_by_week <- case$feedback_by_week
  dates <- case$target_date
  peak_date <- case$peak_date
  stopifnot(!is.na(peak_date), all(is.finite(base)))
  n_weeks <- nrow(base)
  n_levels <- length(alpha)
  observed_weeks <- which(is.finite(truth_latest))
  # E01 scales rates by recent base errors. Freeze that scale once in tracker
  # coordinates for both the fixed baseline and the constant-rate expert bank.
  # https://delphi-forecasting-reports.netlify.app/calibration_experiments/e01_eta_settings
  # Only released rows enter the first 20 rounds; pool their errors over all levels.
  initialization_week <- match(TRUE, cumsum(lengths(feedback_by_week)) >= expert_warmup_rounds)
  initialized <- !is.na(initialization_week)
  expert_warmup_weeks <- if (initialized) initialization_week else n_weeks
  released <- Reduce(c, feedback_by_week[seq_len(expert_warmup_weeks)], init = integer())
  warmup_rows <- head(released, expert_warmup_rounds)
  expert_scale <- raw_expert_scale <- NA_real_
  if (initialized) {
    residuals <- sweep(base_tracker[warmup_rows, , drop = FALSE], 1,
      feedback_tracker[warmup_rows], "-")
    raw_expert_scale <- as.numeric(quantile(abs(residuals), .9, names = FALSE))
    expert_scale <- max(expert_scale_floor, raw_expert_scale)
  }
  # Consume the whole freeze-date batch as warmup, even if it contains >20 rows.
  # With fewer than 20 releases, the fixed baseline and all experts retain the base.
  reference_rate <- reference_multiplier * expert_scale
  expert_rates <- expert_multipliers * expert_scale
  expert_rates[1] <- 0
  expert_names <- paste0("m=", expert_multipliers)
  expert_initialization <- list(initialized = initialized,
    week = if (initialized) expert_warmup_weeks else NA_integer_,
    date = if (initialized) case$reference_date[expert_warmup_weeks] else as.Date(NA),
    rows = warmup_rows, consumed_rows = released,
    raw_scale = raw_expert_scale, scale = expert_scale,
    required_rounds = expert_warmup_rounds, scale_floor = expert_scale_floor)

  # 3. ONE TRACKER: FIXED RATE OR A RATE BASED ON ABSOLUTE SLOPE --------
  run_multiqt <- function(rates, base_forecasts_tracker = base_tracker, base_forecasts = base,
                         warmup_weeks = 0L) {
    if (length(rates) == 1) rates <- rep(rates, n_weeks)
    stopifnot(length(rates) == n_weeks, all(is.finite(rates)), all(rates >= 0),
              all(is.finite(base_forecasts_tracker)))
    offsets <- numeric(n_levels)
    q_tracker <- base_forecasts_tracker
    for (week in seq_len(n_weeks)) {
      # Preserve the exact base during initialization; do not replay these reports later.
      if (week <= warmup_weeks) next
      for (past_week in feedback_by_week[[week]]) {
        covered <- feedback_tracker[past_week] <= q_tracker[past_week, ]
        offsets <- offsets + rates[week] * (alpha - covered)
      }
      q_tracker[week, ] <- isoreg(base_forecasts_tracker[week, ] + offsets)$yf
    }
    # Do not clip negative predictions: that would change the issued coverage indicators.
    q <- from_tracker(q_tracker) / count_scale
    # Keep unchanged bases exact, including eta=0, so transformation roundoff cannot flip ties.
    unchanged <- q_tracker == base_forecasts_tracker
    q[unchanged] <- base_forecasts[unchanged]
    q
  }
  fixed <- if (initialized) run_multiqt(reference_rate, warmup_weeks = expert_warmup_weeks) else base
  # Use only the current base forecast's weekly slope in tracker units, never future counts.
  slope <- c(0, diff(base_tracker[, match(0.5, alpha)]) / (as.numeric(diff(case$reference_date)) / 7))
  slope_rates <- pmin(slope_rate_max, slope_rate_min + abs(slope))
  adaptive <- run_multiqt(slope_rates)
  # Larger rates can unwind offsets faster, but can also amplify stale feedback.

  # 4. WHAT IF WE KNEW THE TRUE PEAK? ----------------------------------
  predictions <- list(original = base, fixed = fixed, absolute_slope = adaptive)
  # Hindsight benchmark: bypass correction for postpeak targets, even if issued prepeak.
  # Setting eta=0 instead would only freeze the accumulated offset.
  predictions$oracle_target <- fixed
  predictions$oracle_target[dates > peak_date, ] <- base[dates > peak_date, ]

  # 5. SEVERAL TRACKERS: AVERAGE THEIR FORECASTS ------------------------
  expert_grid_rows[[case_id]] <- data.frame(case_id, disease = disease_name,
    location = case$location, horizon = horizon_id, tracker_scale, expert = expert_names,
    multiplier = expert_multipliers, eta = expert_rates,
    initialized, initialization_date = expert_initialization$date,
    initialization_rounds = length(warmup_rows), consumed_rounds = length(released),
    raw_error_scale = raw_expert_scale, error_scale = expert_scale,
    scale_floor = expert_scale_floor)
  # Each expert learns from its own issued forecasts and keeps its own offsets.
  experts <- setNames(lapply(expert_rates, function(eta)
    if (!initialized || eta == 0) base else
      run_multiqt(eta, warmup_weeks = expert_warmup_weeks)), expert_names)
  n_experts <- length(experts)
  # Both weighting methods use the same forecasts and count-scale pinball losses.
  # Precomputation does not release losses early: only feedback_by_week does that.
  expert_forecasts <- lapply(seq_len(n_weeks), function(week)
    matrix(vapply(experts, function(expert) expert[week, ], numeric(n_levels)), nrow = n_levels))
  expert_losses <- lapply(seq_len(n_weeks), function(week) {
    error <- feedback_observed[week] - expert_forecasts[[week]]
    t(pmax(error, 0) * alpha + pmax(-error, 0) * (1 - alpha))
  })
  # Each quantile learns its own weights over this location's frozen rate grid.
  run_weights <- function(sharing = TRUE) {
    p <- matrix(1 / n_experts, n_experts, n_levels)
    log_weights <- log(p)
    lambda <- rep(lambda_start, n_levels)
    share <- if (sharing) 1 / (2 * H) else 0
    moments <- matrix(0, n_weeks, n_levels)
    n_feedback <- 0L
    issued_weights <- array(NA_real_, c(n_weeks, n_experts, n_levels),
      dimnames = list(as.character(case$reference_date), expert_names, as.character(alpha)))
    q <- base
    for (week in seq_len(n_weeks)) {
      if (week <= expert_warmup_weeks) {
        issued_weights[week, , ] <- p
        next
      }
      for (past_week in feedback_by_week[[week]]) {
        # Use the PREVIOUS H loss moments, excluding the feedback being processed now.
        if (n_feedback >= H) {
          recent <- seq.int(n_feedback - H + 1L, n_feedback)
          denominator <- colSums(moments[recent, , drop = FALSE])
          positive <- denominator > 0
          lambda[positive] <- sqrt((log(2 * n_experts * H) + 1) / denominator[positive])
          # Tiny positive moments can overflow the ratio before its square root.
          overflow <- positive & !is.finite(lambda)
          lambda[overflow] <- sqrt(log(2 * n_experts * H) + 1) / sqrt(denominator[overflow])
          # If recent losses are all zero, keep the previous finite learning rate.
        }
        loss <- expert_losses[[past_week]]
        # Paper's second moment: average squared expert losses using pre-update weights.
        n_feedback <- n_feedback + 1L
        moments[n_feedback, ] <- colSums(p * loss^2)
        # Exponential weights, stabilized and normalized separately for each quantile.
        log_weights <- log_weights - sweep(loss, 2, lambda, `*`)
        log_weights <- sweep(log_weights, 2, apply(log_weights, 2, max), `-`)
        p <- exp(log_weights)
        p <- sweep(p, 2, colSums(p), `/`)
        if (sharing) {
          p <- (1 - share) * p + share / n_experts
          log_weights <- log(p)
        }
        # Without sharing, retain log_weights even if exp() underflows.
      }
      issued_weights[week, , ] <- p
      # Arithmetic mixture in scaled counts, followed by repair of quantile crossings.
      q[week, ] <- rowSums(expert_forecasts[[week]] * t(p))
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
  weighted <- run_weights(sharing = FALSE)
  # Each outer tracker has its own offsets and learns from its OWN final forecasts.
  # Its state stays continuous when new restart copies join the underlying mixture.
  run_outer <- function(mixture) {
    counts <- mixture * count_scale
    if (any(!is.finite(counts)) || (use_log_scale && any(counts <= -1)))
      stop("Mixture is outside the ", tracker_scale, " domain for case ", case_id, ".")
    run_multiqt(outer_rate, base_forecasts_tracker = to_tracker(counts),
      base_forecasts = mixture, warmup_weeks = expert_warmup_weeks)
  }
  predictions$uniform <- Reduce(`+`, experts) / n_experts
  predictions$uniform[seq_len(expert_warmup_weeks), ] <- base[seq_len(expert_warmup_weeks), ]
  predictions$weighted <- weighted$q
  predictions$shared <- shared$q
  predictions$shared_outer <- run_outer(shared$q)

  # Weekly starts are a separate experiment from the persistent rate-grid bank.
  restart_bank <- run_restart_experts(case$base, case$feedback_observed,
    feedback_by_week, case$reference_date, alpha, restart_multipliers, expert_scale,
    expert_warmup_weeks, count_scale)
  keep <- which(restart_bank$experts$multiplier == restart_multiplier)
  fixed_restart_bank <- list(forecasts = restart_bank$forecasts[, keep, , drop = FALSE],
    experts = restart_bank$experts[keep, , drop = FALSE],
    feedback_counts = restart_bank$feedback_counts[, keep, drop = FALSE])
  restart_equal <- run_restart_weights(fixed_restart_bank, base, feedback_observed,
    feedback_by_week, alpha, learn = FALSE)
  restart_learned <- run_restart_weights(fixed_restart_bank, base, feedback_observed,
    feedback_by_week, alpha, weight_rate = restart_weight_rate, share = restart_share)
  restart_joint <- run_restart_weights(restart_bank, base, feedback_observed,
    feedback_by_week, alpha, weight_rate = restart_weight_rate, share = restart_share)
  # The first-start control has exactly the same birth rule and own-feedback mask.
  # Both this control and all ensemble copies issue forecasts immediately at birth.
  predictions$restart_persistent <- base
  copies <- which(fixed_restart_bank$experts$type == "restart")
  if (length(copies)) {
    first_copy <- copies[1]
    first <- fixed_restart_bank$experts$birth_week[first_copy]
    issued <- seq.int(first, n_weeks)
    predictions$restart_persistent[issued, ] <- fixed_restart_bank$forecasts[issued, first_copy, ]
  }
  predictions$restart_uniform <- restart_equal$q
  predictions$restart_weighted <- restart_learned$q
  predictions$restart_grid <- restart_joint$q
  for (method in c("restart_uniform", "restart_weighted", "restart_grid"))
    predictions[[paste0(method, "_outer")]] <- run_outer(predictions[[method]])
  restart_expert_rows[[case_id]] <- data.frame(case_id = rep(case_id, nrow(restart_bank$experts)),
    location = rep(case$location, nrow(restart_bank$experts)), restart_bank$experts)
  restart_diagnostics <- list(restart_uniform = restart_equal$diagnostics,
    restart_weighted = restart_learned$diagnostics, restart_grid = restart_joint$diagnostics)
  restart_diagnostic_rows[[case_id]] <- do.call(rbind, lapply(names(restart_diagnostics), function(method)
    data.frame(case_id, location = case$location, method, reference_date = case$reference_date,
      restart_diagnostics[[method]])))

  # 6. KEEP THE ISSUED FORECASTS AND SIGNED WEEKLY DIAGNOSTICS --------
  # Score against the latest vintage; the online updates above saw only first-use feedback.
  # Counts are already scaled: do not divide quantile loss a second time.
  extension_runs[[case_id]] <- list(location = case$location, base = base, alpha = alpha,
    base_tracker = base_tracker, feedback_tracker = feedback_tracker,
    feedback_observed = feedback_observed, truth_latest = truth_latest,
    feedback_by_week = feedback_by_week,
    reference_date = case$reference_date, target_date = dates, peak_date = peak_date,
    relative_week = case$relative_week, predictions = predictions,
    weights = list(weighted = weighted$weights, shared = shared$weights,
      restart_uniform = restart_equal$weights, restart_weighted = restart_learned$weights,
      restart_grid = restart_joint$weights),
    restart = list(includes_base = FALSE, admission = "at_birth",
      fixed_multiplier = restart_multiplier, multipliers = restart_multipliers,
      weight_rate = restart_weight_rate, share = restart_share,
      experts = restart_bank$experts, forecasts = restart_bank$forecasts,
      feedback_counts = restart_bank$feedback_counts, diagnostics = restart_diagnostics),
    experts = experts, expert_multipliers = expert_multipliers, expert_rates = expert_rates,
    expert_initialization = expert_initialization,
    reference_multiplier = reference_multiplier, reference_rate = reference_rate,
    slope_rates = slope_rates, count_scale = count_scale, scale_used = scale_used,
    tracker_scale = tracker_scale, prediction_scale = prediction_scale)
  scale_rows[[case_id]] <- data.frame(location = case$location, scale_method, scale_used,
    count_scale, tracker_scale, prediction_scale,
    scale_as_of = case$scale_as_of, reference_weeks = length(scale_history))
  for (j in seq_along(alpha)) {
    for (method in names(predictions)) {
      q <- predictions[[method]][observed_weeks, j]
      covered <- truth_latest[observed_weeks] <= q
      error <- truth_latest[observed_weeks] - q
      error_rows[[paste(case_id, j, method)]] <- data.frame(location = case$location,
        target_date = dates[observed_weeks], relative_week = case$relative_week[observed_weeks],
        within_winter = case$within_winter[observed_weeks],
        alpha = alpha[j], method, calibration_residual = covered - alpha[j],
        quantile_error = pmax(error, 0) * alpha[j] + pmax(-error, 0) * (1 - alpha[j]))
    }
  }
}
location_errors <- do.call(rbind, error_rows)
rownames(location_errors) <- NULL
write_result(location_errors, "errors_by_location.csv")
scale_table <- do.call(rbind, scale_rows)
write_result(scale_table, "scales.csv")
write_result(do.call(rbind, expert_grid_rows), "expert_grid.csv")
write_result(do.call(rbind, restart_expert_rows), "restart_experts.csv")
write_result(do.call(rbind, restart_diagnostic_rows), "restart_diagnostics.csv")
write_result(data.frame(tracker_scale, fixed_multiplier = restart_multiplier,
  grid_multipliers = paste(restart_multipliers, collapse = ","),
  weight_rate = restart_weight_rate, share = restart_share,
  outer_rate, base_expert = "No permanent base; newborns forecast the base; exact base fallback before births",
  birth = "Restart copies at every issue after frozen-scale warmup",
  admission = "Immediate at birth, before the first forecast",
  loss_eligibility = "All historical issuers, including each expert's first released loss",
  new_mass = "Each newborn receives 1/K among all K active copies; incumbents retain relative weights"),
  "restart_settings.csv")
saveRDS(extension_runs, file.path(output_dir, "issued_forecasts.rds"))

# 7. AVERAGE LOCATIONS ONLY AFTER CALCULATING THEIR OWN ERRORS -------
# Per-quantile weekly curves and their overall average both remain signed.
# These weekly diagnostics are not time-window calibration-error scores.
# CE = mean over locations and quantiles of |mean residual over that location's weeks|.
# This prevents overcoverage in one location cancelling undercoverage elsewhere.
summarize_errors <- function(rows) {
  local <- aggregate(cbind(calibration_residual, quantile_error) ~ location + alpha + method, rows, mean)
  local$calibration_error <- abs(local$calibration_residual)
  levels <- aggregate(cbind(calibration_residual, calibration_error, quantile_error) ~
    alpha + method, local, mean)
  # Pooled CE takes the absolute value AFTER averaging locations, allowing cancellation.
  levels$pooled_calibration_error <- abs(levels$calibration_residual)
  scores <- aggregate(cbind(calibration_residual, calibration_error, quantile_error) ~ method, local, mean)
  pooled <- aggregate(pooled_calibration_error ~ method, levels, mean)
  list(by_location = local, by_level = levels, by_method = merge(scores, pooled, by = "method"))
}
evaluation_errors <- if (combined)
  location_errors[location_errors$relative_week %in% combined_weeks, ] else location_errors
whole_scores <- summarize_errors(evaluation_errors)
local_time_averages <- whole_scores$by_location
time_averages <- whole_scores$by_level
overall_time_average <- aggregate(cbind(calibration_residual, calibration_error,
  pooled_calibration_error, quantile_error) ~ method, time_averages, mean)
write_result(local_time_averages, "time_averages_by_location.csv")
write_result(rbind(time_averages, data.frame(alpha = NA_real_, overall_time_average)), "time_averages.csv")

# The complete combined cohort contributes once at every relative week.
# In individual mode retain calendar dates and the original full-history view.
if (combined) {
  errors <- aggregate(cbind(calibration_residual, quantile_error) ~
    relative_week + alpha + method, evaluation_errors, mean)
  errors$plot_time <- errors$relative_week
} else {
  errors <- evaluation_errors
  errors$plot_time <- as.numeric(errors$target_date)
}
# Signed mean over locations and quantiles: the overall direction of miscoverage.
# Opposite gaps (e.g. too-narrow intervals) can cancel; see the per-quantile panels.
overall <- aggregate(cbind(mean_calibration_residual = calibration_residual,
  quantile_error) ~ plot_time + method, errors, mean)
write_result(errors, "errors_over_time.csv")
write_result(overall, "errors_average_quantiles.csv")

# Pre ends with the peak week (0); post starts one week later. Whole means full
# history for one location, or the common aligned window for the combined cohort.
summary_rows <- list()
summary_by_location <- list()
for (window in c("whole", "pre", "post")) {
  rows <- if (window == "whole") evaluation_errors else location_errors[
    location_errors$within_winter & location_errors$relative_week %in%
      if (window == "pre") -weeks_before_peak:0 else seq_len(weeks_after_peak), ]
  if (!nrow(rows)) {
    summary_rows[[window]] <- data.frame(method = names(predictions),
      calibration_residual = NA_real_, calibration_error = NA_real_, quantile_error = NA_real_,
      pooled_calibration_error = NA_real_, window, n_locations = 0L, min_weeks = 0L, max_weeks = 0L)
    next
  }
  window_scores <- if (window == "whole") whole_scores else summarize_errors(rows)
  local <- window_scores$by_location
  local$window <- window
  summary_by_location[[window]] <- local
  score <- window_scores$by_method
  support <- aggregate(target_date ~ location, unique(rows[c("location", "target_date")]), length)
  score$window <- window
  score$n_locations <- nrow(support)
  score$min_weeks <- min(support$target_date)
  score$max_weeks <- max(support$target_date)
  summary_rows[[window]] <- score
}
summary_scores <- do.call(rbind, summary_rows)
write_result(summary_scores, "summary_scores.csv")
write_result(do.call(rbind, summary_by_location), "summary_scores_by_location_quantile.csv")

# 8. WEEKLY COVERAGE RESIDUALS AND QUANTILE LOSS ---------------------
labels <- c(original = "Original",
  fixed = paste0("Fixed ", tracker_method, " (eta=", reference_multiplier, "*S0)"),
  absolute_slope = if (use_log_scale) "Log-slope rate" else "Count-slope rate",
  oracle_target = "Target oracle (hindsight)",
  uniform = "Rate grid: equal", weighted = "Rate grid: learned, no sharing",
  shared = paste0("Rate grid: learned + sharing (H=", H, ")"),
  shared_outer = paste0("Rate grid: sharing + outer (eta=", outer_rate, ")"),
  restart_persistent = paste0("First start only (m=", restart_multiplier, ")"),
  restart_uniform = paste0("Weekly starts: equal (m=", restart_multiplier, ")"),
  restart_weighted = paste0("Weekly starts: learned (m=", restart_multiplier, ")"),
  restart_grid = "Starts x rate grid: learned",
  restart_uniform_outer = "Starts equal + outer",
  restart_weighted_outer = "Starts learned + outer",
  restart_grid_outer = "Starts x rate grid + outer")
colours <- c(original = "gray40", fixed = "#0072B2", absolute_slope = "#009E73",
  oracle_target = "#CC79A7", uniform = "#E69F00", weighted = "#56B4E9", shared = "#882255",
  shared_outer = "#D55E00",
  restart_persistent = "#0072B2", restart_uniform = "#E69F00", restart_weighted = "#009E73",
  restart_grid = "#882255", restart_uniform_outer = "#56B4E9",
  restart_weighted_outer = "#D55E00", restart_grid_outer = "#CC79A7")
method_names <- names(labels)
method_groups <- list(
  "Frozen rate-grid comparison" = c("original", "fixed", "absolute_slope", "oracle_target",
    "uniform", "weighted", "shared", "shared_outer"),
  "Weekly starts: immediate admission" = c("original", "restart_persistent", "restart_uniform",
    "restart_weighted", "restart_grid", "restart_uniform_outer", "restart_weighted_outer",
    "restart_grid_outer"))
loss_units <- if (all(scale_table$scale_used == "none")) "counts" else
  if (all(scale_table$scale_used == "max")) "counts / max" else
    if (all(scale_table$scale_used == "median")) "counts / median" else "counts / local scale"
loss_axis_label <- paste0("Quantile loss (", loss_units, ")")
mean_loss_axis_label <- paste("Mean", tolower(loss_axis_label))
time_axis_label <- if (combined) "Weeks from winter peak (peak = 0)" else "Target date"
plot_levels <- plot_levels[plot_levels %in% alpha]
peak_marker <- if (combined) 0 else as.numeric(peak_date)
ticks <- if (combined) combined_weeks else
  pretty(as.Date(range(errors$plot_time), origin = "1970-01-01"), n = 3)
mean_limit <- 1.15 * max(.01, abs(time_averages$calibration_residual))
mean_limits <- list(calibration_residual = c(-mean_limit, mean_limit),
                    quantile_error = c(0, 1.15 * max(1e-8, time_averages$quantile_error)))
cohort_label <- if (combined) paste(n_locations, "locations with all", length(combined_weeks),
  "aligned weeks; US excluded") else paste(case$location_name, "| full history")
plot_title <- paste(disease_name, "| horizon", horizon_id, "|", cohort_label,
  if (use_log_scale) "| log1p(counts) updates" else "| raw-count updates")
plot_errors <- function(rows, title, calibration_metric = "calibration_residual") {
  rows <- rows[rows$method %in% method_names, ]
  for (metric in c(calibration_metric, "quantile_error")) {
    residual <- metric != "quantile_error"
    # Axes follow the plotted values, like the loss panels; residual axes keep 0 visible.
    limits <- if (residual) 1.05 * range(-.01, .01, rows[[metric]]) else
      c(0, 1.05 * max(1e-8, rows[[metric]]))
    par(mar = if (residual) c(.5, 5.2, 2, 9) else c(3, 5.2, .5, 9))
    axis_label <- if (residual) "Coverage - alpha" else loss_axis_label
    plot(rows$plot_time, rows[[metric]], type = "n", xaxt = "n", ylim = limits,
      xlab = "", ylab = paste(strwrap(axis_label, width = if (residual) 14 else 20), collapse = "\n"),
      main = if (residual) title else "")
    if (residual) abline(h = 0, col = "gray85") else {
      if (combined) axis(1, at = ticks) else axis.Date(1, at = ticks, format = "%b %Y")
    }
    abline(v = peak_marker, col = "gray60", lty = 3)
    for (method in method_names) {
      line <- rows[rows$method == method, ]
      line <- line[order(line$plot_time), ]
      lines(line$plot_time, line[[metric]], col = colours[method], lwd = 1.5)
    }
    means <- vapply(method_names, function(method) mean(rows[[metric]][rows$method == method]), numeric(1))
    # Right-hand scale: shared across quantile panels; the overall panel uses its own means.
    mean_range <- if (metric == "mean_calibration_residual")
      c(-1, 1) * 1.15 * max(.01, abs(means)) else mean_limits[[metric]]
    mean_ticks <- pretty(mean_range, n = 3)
    mean_ticks <- mean_ticks[mean_ticks >= mean_range[1] & mean_ticks <= mean_range[2]]
    scale <- diff(limits) / diff(mean_range)
    offset <- limits[1] - scale * mean_range[1]
    edge <- par("usr")[2]
    inch <- diff(par("usr")[1:2]) / par("pin")[1]
    rect(edge + .04 * inch, limits[1], edge + .68 * inch, limits[2], col = "gray97", border = NA, xpd = NA)
    segments(edge + .04 * inch, offset, edge + .68 * inch, offset, col = "gray80", xpd = NA)
    points(edge + seq(.12, .62, length.out = length(means)) * inch,
      offset + scale * means, col = colours[names(means)], pch = 16, cex = .85, xpd = NA)
    axis(4, at = offset + scale * mean_ticks, labels = mean_ticks,
      pos = edge + .75 * inch, las = 1, cex.axis = .8, xpd = NA)
    # The right-hand means use the same units as the left axis.
    mean_axis_label <- if (residual) "Mean coverage - alpha" else "Mean quantile loss"
    mtext(paste(strwrap(mean_axis_label, width = 16), collapse = "\n"),
      side = 4, line = 8, cex = .8)
  }
}
panels <- matrix(0L, nrow = 2 * ceiling(length(plot_levels) / 3), ncol = 3)
for (j in seq_along(plot_levels)) {
  row <- 2 * ((j - 1L) %/% 3) + 1L
  column <- (j - 1L) %% 3 + 1L
  panels[row:(row + 1L), column] <- c(2 * j - 1L, 2 * j)
}
last_panel <- 2 * length(plot_levels)
panels <- rbind(panels, rep(last_panel + 1L, 3), rep(last_panel + 2L, 3), rep(last_panel + 3L, 3))
cairo_pdf(figure_file, width = 16, height = 16, family = "DejaVu Sans")
for (group_name in names(method_groups)) {
  method_names <- method_groups[[group_name]]
  layout(panels, heights = c(rep(1, nrow(panels) - 1L), .7))
  par(oma = c(2, 0, 4, 0), mgp = c(2.3, .6, 0), cex = 1)
  for (level in plot_levels) plot_errors(errors[errors$alpha == level, ], paste("alpha =", level))
  plot_errors(overall, "Mean over quantile levels",
    calibration_metric = "mean_calibration_residual")
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", labels[method_names], col = colours[method_names], lty = 1, lwd = 2, ncol = 3, bty = "n")
  mtext(paste(plot_title, "|", group_name), side = 3, outer = TRUE, line = 2, font = 2, cex = 1.15)
  mtext("Coverage residual = coverage - alpha | Residuals stay signed; opposite residuals can cancel in the mean over quantile levels",
    side = 3, outer = TRUE, line = .5)
  mtext(paste(time_axis_label, "| Dotted line: winter peak"),
    side = 1, outer = TRUE, line = .5)
}
dev.off()

cairo_pdf(overall_figure_file, width = 12, height = 8, family = "DejaVu Sans")
for (group_name in names(method_groups)) {
  method_names <- method_groups[[group_name]]
  layout(matrix(1:3, ncol = 1), heights = c(1, 1, .25))
  par(oma = c(2, 0, 3, 0), mgp = c(2.5, .6, 0), cex = .9)
  plot_errors(overall, "Mean over quantile levels",
    calibration_metric = "mean_calibration_residual")
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", labels[method_names], col = colours[method_names], lty = 1, lwd = 2, ncol = 2, bty = "n", cex = .85)
  mtext(paste(plot_title, "|", group_name), side = 3, outer = TRUE, line = 1, font = 2, cex = .85)
  mtext(paste(time_axis_label, "| Coverage residuals stay signed; opposite residuals can cancel"),
    side = 1, outer = TRUE, line = .5, cex = .85)
}
dev.off()

# 9. CALIBRATION ERROR AND QUANTILE LOSS ON MATCHED WINDOWS ----------
# Do not average data or signed gaps before the absolute value. Every location
# and quantile contributes one calibration gap; scale is already applied above.
cairo_pdf(tradeoff_figure_file, width = 14, height = 5.5, family = "DejaVu Sans")
for (group_name in names(method_groups)) {
  method_names <- method_groups[[group_name]]
  layout(rbind(1:3, c(4, 4, 4)), heights = c(1, .22))
  par(oma = c(3, 0, 3, 0), mgp = c(2.6, .7, 0))
  window_titles <- c(whole = if (combined) "Whole aligned window" else "Full observed history",
                     pre = "Up to winter peak", post = "After winter peak")
  xlim <- c(0, 1.1 * max(1e-8, summary_scores$quantile_error, na.rm = TRUE))
  ylim <- c(0, 1.1 * max(.01, summary_scores$calibration_error, na.rm = TRUE))
  for (window in c("whole", "pre", "post")) {
    rows <- summary_scores[summary_scores$window == window & summary_scores$method %in% method_names, ]
    par(mar = c(4, 4, 4, 1))
    plot(NA, xlim = xlim, ylim = ylim, xlab = mean_loss_axis_label, ylab = "Calibration error",
      main = "")
    mtext(window_titles[window], side = 3, line = 1.8, font = 2, cex = .95)
    grid(col = "gray90")
    if (nrow(rows)) {
      support <- if (rows$min_weeks[1] == rows$max_weeks[1]) rows$min_weeks[1] else
        paste0(rows$min_weeks[1], "-", rows$max_weeks[1])
      mtext(paste(rows$n_locations[1], "locations;", support, "weeks each"), side = 3, line = .25, cex = .8)
      points(rows$quantile_error, rows$calibration_error, col = colours[rows$method], pch = 16, cex = 1.2)
      if (combined) {
        # Open circles: pooled CE of the same method, directly below its local CE.
        segments(rows$quantile_error, rows$pooled_calibration_error, rows$quantile_error,
          rows$calibration_error, col = adjustcolor(colours[rows$method], .4))
        points(rows$quantile_error, rows$pooled_calibration_error, col = colours[rows$method], pch = 1, cex = 1.2)
      }
      if (!any(is.finite(rows$calibration_error)))
        text(mean(xlim), mean(ylim), "No observations inside winter", cex = .8)
    }
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", labels[method_names], col = colours[method_names], pch = 16, ncol = 3, bty = "n", cex = .8)
  mtext(paste(plot_title, "|", group_name), side = 3, outer = TRUE, line = 1, font = 2, cex = .9)
  mtext("Calibration error (CE): mean absolute coverage residual after averaging weeks within each location and quantile level.",
    side = 1, outer = TRUE, line = .5, cex = .8)
  mtext(paste(if (combined) "Filled: per-location CE; open: pooled CE (average locations before taking absolute values).",
    "Up to winter peak includes the peak week."), side = 1, outer = TRUE, line = 1.8, cex = .8)
}
dev.off()

# 10. ISSUED ENSEMBLE WEIGHTS OVER TIME -----------------------------
# These are the weights used to mix forecasts before crossing repair or an outer
# correction. The outer methods keep their underlying ensemble's weights.
weight_figure_file <- file.path(figure_dir, "multiqt_weight_timelines.pdf")
weight_times <- if (combined) combined_weeks else as.numeric(extension_runs[[1]]$target_date)
weight_step <- if (combined) 1 else 7
age_breaks <- c(0, 1, 4, 8, 16, Inf)
age_labels <- c("<1 week", "1-3 weeks", "4-7 weeks", "8-15 weeks", "16+ weeks")
age_colours <- c("#F1C453", "#E79B4F", "#C75643", "#864C7D", "#3D4E78")
weight_multipliers <- sort(unique(c(expert_multipliers, restart_multipliers)))
weight_colours <- setNames(c("gray45", hcl.colors(length(weight_multipliers) - 1L, "Viridis")),
  as.character(weight_multipliers))
rate_colours <- weight_colours[as.character(expert_multipliers)]
restart_rate_colours <- weight_colours[as.character(restart_multipliers)]
fallback_colour <- "gray90"
weight_pages <- list(
  "Frozen rate-grid weights" = data.frame(
    method = c("uniform", "weighted", "shared"), distribution = "multiplier",
    title = c("Equal weights", "Learned, no sharing", "Learned + sharing")))

# Bin within each location before averaging: restart identities and start dates
# differ between locations, whereas rate multipliers and calendar ages are comparable.
weight_distribution <- function(run, method, distribution, level) {
  weeks <- if (combined) match(combined_weeks, run$relative_week) else seq_along(run$target_date)
  level_index <- match(level, run$alpha)
  restart <- startsWith(method, "restart_")
  if (method == "uniform") {
    p <- matrix(1 / length(run$expert_multipliers), length(weeks), length(run$expert_multipliers))
    multipliers <- run$expert_multipliers
  } else {
    weights <- run$weights[[method]]
    p <- matrix(weights[weeks, , level_index, drop = FALSE], nrow = length(weeks), ncol = dim(weights)[2])
    if (restart) {
      experts <- run$restart$experts[match(dimnames(weights)[[2]], run$restart$experts$expert), , drop = FALSE]
      stopifnot(!anyNA(experts$expert))
      multipliers <- experts$multiplier
    } else multipliers <- run$expert_multipliers
  }
  components <- if (distribution == "age") age_labels else
    as.character(if (restart) restart_multipliers else expert_multipliers)
  mass <- matrix(0, length(weeks), length(components), dimnames = list(NULL, components))
  if (distribution == "age") {
    ages <- outer(as.numeric(run$reference_date[weeks]), as.numeric(experts$birth_date), `-`) / 7
    for (j in seq_along(components))
      mass[, j] <- rowSums(p * (ages >= age_breaks[j] & ages < age_breaks[j + 1L]))
  } else {
    for (j in seq_along(components))
      mass[, j] <- rowSums(p[, multipliers == as.numeric(components[j]), drop = FALSE])
  }
  # Before the first restart, no expert has weight and the original forecast is
  # issued instead. Keep that fallback visible, including partial combined cohorts.
  total <- rowSums(p)
  mass <- cbind(mass, "Base fallback" = pmax(0, 1 - total))
  stopifnot(all(is.finite(mass)), all(mass >= 0), max(abs(rowSums(mass) - 1)) < 1e-8)
  list(mass = mass, active = total > 0)
}

plot_weight_distribution <- function(mass, fills, title, level, first_column, last_row) {
  par(mar = c(if (last_row) 2.4 else .3, if (first_column) 3.7 else 2.2, 1.5, .6))
  plot(NA, xlim = range(weight_times) + c(-.5, .5) * weight_step, ylim = c(0, 1),
    xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "",
    ylab = if (first_column) paste("alpha =", level) else "", main = title)
  bottom <- numeric(nrow(mass))
  for (j in seq_len(ncol(mass))) {
    top <- bottom + mass[, j]
    rect(weight_times - .45 * weight_step, bottom, weight_times + .45 * weight_step,
      top, col = fills[j], border = NA)
    bottom <- top
  }
  abline(h = .5, col = adjustcolor("white", .6), lty = 3)
  abline(v = peak_marker, col = "gray30", lty = 3)
  axis(2, at = c(0, .5, 1), las = 1, cex.axis = .75)
  if (last_row) {
    if (combined) axis(1, at = ticks, cex.axis = .8) else
      axis.Date(1, at = ticks, format = "%b %Y", cex.axis = .8)
  }
  box(col = "gray70")
}

weight_rows <- list()
cairo_pdf(weight_figure_file, width = 18, height = max(8, 1.45 * length(plot_levels) + 3.7), family = "DejaVu Sans")
for (page in names(weight_pages)) {
  specs <- weight_pages[[page]]
  panel_ids <- matrix(seq_len(length(plot_levels) * nrow(specs)), ncol = nrow(specs), byrow = TRUE)
  layout(rbind(panel_ids, rep(max(panel_ids) + 1L, nrow(specs))),
    heights = c(rep(1, length(plot_levels)), .8))
  par(oma = c(4.8, .8, 5.2, .8), mgp = c(2, .4, 0), cex = .9)
  for (level in plot_levels) for (j in seq_len(nrow(specs))) {
    method <- specs$method[j]; distribution <- specs$distribution[j]
    local <- lapply(extension_runs, weight_distribution, method = method,
      distribution = distribution, level = level)
    mass <- Reduce(`+`, lapply(local, `[[`, "mass")) / n_locations
    active_locations <- Reduce(`+`, lapply(local, function(x) as.integer(x$active)))
    fills <- if (distribution == "age") age_colours else if (startsWith(method, "restart_"))
      restart_rate_colours else rate_colours
    plot_weight_distribution(mass, c(fills, fallback_colour),
      if (level == plot_levels[1]) specs$title[j] else "", level,
      first_column = j == 1L, last_row = level == tail(plot_levels, 1))
    weight_rows[[paste(method, distribution, level)]] <- data.frame(
      plot_time = rep(weight_times, ncol(mass)),
      target_date = if (combined) as.Date(NA) else rep(extension_runs[[1]]$target_date, ncol(mass)),
      relative_week = if (combined) rep(combined_weeks, ncol(mass)) else
        rep(extension_runs[[1]]$relative_week, ncol(mass)),
      alpha = level, method, distribution, component = rep(colnames(mass), each = nrow(mass)),
      weight = as.vector(mass), n_locations, active_locations = rep(active_locations, ncol(mass)))
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  if (page == "Frozen rate-grid weights") {
    legend("center", paste("m =", names(rate_colours)), fill = rate_colours,
      border = NA, ncol = 6, bty = "n", cex = .95)
    note <- paste0("")
  } else {
    legend("top", c(paste("Age:", age_labels), "Base fallback"),
      fill = c(age_colours, fallback_colour), border = NA, ncol = 6, bty = "n", cex = .9)
    legend("bottom", paste("m =", names(restart_rate_colours)), fill = restart_rate_colours,
      border = NA, ncol = length(restart_rate_colours), bty = "n", cex = .9)
    note <- "Age is calendar weeks since start. Age bands sum weights over copies; grid panels show separate age/rate marginals. Outer corrections retain these weights."
  }
  mtext(page, side = 3, outer = TRUE, line = 2.4, font = 2, cex = 1.3)
  mtext(plot_title, side = 3, outer = TRUE, line = .9, cex = .9)
  mtext(paste(time_axis_label, "| Weight at issuance for each target | Dotted line: winter peak"),
    side = 1, outer = TRUE, line = .4, cex = .85)
  mtext(note, side = 1, outer = TRUE, line = 1.6, cex = .8)
  # mtext(if (page == "Frozen rate-grid weights") {
  #   if (combined) "Each bar is the equal-location mean of issued weight distributions; total weight = 1." else
  #     "Each bar shows issued expert weights; total weight = 1. White gaps have no issued forecasts."
  # } else if (combined)
  #   "Equal-location average of weight distributions. Gray mass is the fraction of locations still using the base fallback." else
  #   "Gray bars indicate base fallback before any copies start. White gaps have no issued forecasts.",
  #   side = 1, outer = TRUE, line = 2.7, cex = .8)
}
dev.off()
write_result(do.call(rbind, weight_rows), "weight_distributions.csv")
