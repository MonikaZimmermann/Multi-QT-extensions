# Functions for scripts 02 and 04. Sourcing this file only defines functions:
#   source("R/03_functions.R")
# Run R/01_data_overview.R once to prepare <output_dir>/data/peak_aligned.rds.
# Forecasts are tracked in tracker units (log1p(counts) or counts; Delphi's production
# trackers always on the sqrt scale); predictions
# and scores use counts / count_scale (per 100k residents by default); the expert
# weights learn from losses relative to the base's (weighting_loss_scale()).

# 1. SHARED SETTINGS -----------------------------------------------------
# Scripts 01, 02, 04, 05 and 06 read these settings, so a change here applies to all.
# base_model: the base forecasts of both diseases, "cmu_submission" (the CMU-TimeSeries hub
# submissions) or "windowed_seasonal" (our forecasts of R/00_base_forecasts.R, the longer
# history); R/01_data_overview.R sets the dates from it. Each base has its own output folder
# (output_dir: data, results, pdf and cache): output/ and output/windowed_seasonal/. The
# default is cmu_submission; the environment variable MULTIQT_BASE_MODEL chooses the base
# of one run without changing the default, e.g.
#   MULTIQT_BASE_MODEL=windowed_seasonal Rscript R/01_data_overview.R
# horizons: the horizons run, tuned and evaluated. windowed_seasonal leaves out horizon -1:
# its forecasts there have zero width once the week before is reported (the latest report
# is a predictor). R/01 prepares horizon -1 for both bases: the methods read the horizon -1
# forecasts as the start of each forecast path.
multiqt_settings <- function() {
  base_model <- Sys.getenv("MULTIQT_BASE_MODEL", "cmu_submission")
  stopifnot(base_model %in% c("cmu_submission", "windowed_seasonal"))
  list(
  base_model = base_model,
  output_dir = if (base_model == "cmu_submission") "output" else file.path("output", base_model),
  location_ids = c("US", "06", "48", "36", "42"),  # US, California, Texas, New York, Pennsylvania.
  # The pooled summaries (02, 04, 06) and the tuning (05): the 10 most populous states
  # without Florida (California, Texas, New York, Pennsylvania, Illinois, Ohio, Georgia,
  # North Carolina, Michigan, New Jersey). The US, a sum of the states, is never pooled.
  pooled_locations = c("06", "48", "36", "42", "17", "39", "13", "37", "26", "34"),
  horizons = if (base_model == "windowed_seasonal") 0:3 else -1:3,
  transform = "identity",           # MultiQT tracks raw counts (log1p: fixed_log, leaky_log).
  scale_method = "population",      # Scores/plots: per 100k residents ("population"), or "max", "median" or "none".
  expert_multipliers = c(0, .01, .015, .02, .03, .05, .075, .1, .15, .2, .3, .45),
  warmup_rounds = 20L,              # First released rows that freeze the error scale S0.
  replay_warmup = TRUE,             # TRUE: learn from the first week on, S0 from those rows in hindsight.
  feedback_delay = 0L,              # Wait this many weeks after the first report, then learn once from the current one.
  feedback_revisions = FALSE,       # TRUE: re-run each week on the reports known then (Delphi).
  feedback_corrected = FALSE,       # TRUE: young reports as corrected labels (label_correction).
  # The feedback designs of 02, 04 and 06: the three designs, then their corrected versions
  # (notes, "Revised reports as labels: corrected labels").
  feedback_designs = data.frame(delay = c(0L, 1L, 0L, 0L, 1L, 0L),
    revisions = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE),
    corrected = rep(c(FALSE, TRUE), each = 3L)),
  # Corrected labels (correct_labels()): a report younger than `maturity` weeks is replaced by
  # plausible counts from the revisions of the last `history` rows at least that old; the
  # location's own rows once it has `own_rows` of them, else the pooled locations'; fewer
  # than `min_rows`: the report itself.
  label_correction = list(maturity = 2L, history = 52L, own_rows = 13L, min_rows = 8L),
  score_target = "latest",          # Scores against the "latest" report or the "feedback" learned from.
  score_filled = TRUE,              # Score filled rounds (fill_gaps_up_to in R/01_data_overview.R).
  scale_floor = 1e-6,               # Guard on S0, sigma_t and the level scales, not each rate.
  # The methods' own settings (rates, H, lambda) are in method_parameters() (section 7).
  # Window lengths in weeks. No forecasts were issued 2025-10-25 to 2025-11-15, so
  # COVID after-peak windows of September 2025 peaks miss some weeks.
  weeks_before_peak = 8L,           # before = weeks -7..0: 8 weeks, including the peak week.
  weeks_after_peak = 8L,            # after = weeks 1..8: 8 weeks.
  # Phase probabilities of phase_angle, phase_hidden and the angle controls: R/05's file of this
  # base, or another complete file of the same structure (MULTIQT_PHASE_PROBABILITIES, e.g. a
  # candidate of agent/outputs/phase_probability_regularization/); phase_inputs() records its
  # path and checksum in every case, so a changed source changes R/04's cache keys.
  phase_probabilities = Sys.getenv("MULTIQT_PHASE_PROBABILITIES",
    file.path(if (base_model == "cmu_submission") "output" else file.path("output", base_model),
      "data", "phase_probabilities.rds")))
}
# The prepared data of the base model of `settings` (R/01_data_overview.R), checked to be
# of that base.
prepared_path <- function(settings) file.path(settings$output_dir, "data", "peak_aligned.rds")
read_prepared <- function(settings) {
  path <- prepared_path(settings)
  if (!file.exists(path)) stop("No ", path, ": run Rscript R/01_data_overview.R.")
  prepared <- readRDS(path)
  if (!identical(prepared$settings$base_model, settings$base_model))
    stop(path, " is not prepared for base model ", settings$base_model,
      ". Rerun Rscript R/01_data_overview.R.")
  prepared
}

# 2. DELAYED FEEDBACK, CASES, SCALES AND WINDOWS ----------------------------
# Reviewed exclusions, not an automatic outlier rule: one row per confirmed erroneous
# report version (disease, location, target week, report date = version, reason). These
# five were identified from later corrections in nhsn_vintages_epidata.parquet, so they
# are labelled "retrospective cleaning": the error was not known when the report was
# released ("real-time cleaning" would be an error identifiable then). A large revision
# alone only flags a report for review (revision_review_flags()); it is excluded only
# once its row is added here. Keep the original value and correction evidence; never
# exclude a whole target week. Any change here changes report_filter_id(), so results
# made with another table are rejected (check_feedback_protocol()).
reviewed_report_exclusions <- function() data.frame(
  disease = c("covid", "covid", "covid", "flu", "flu"),
  location = c("42", "01", "01", "01", "01"),
  target_end_date = as.Date(c("2024-11-16", "2026-01-24", "2026-02-14", "2026-02-14",
    "2026-01-24")),
  version = as.Date(c("2024-11-20", "2026-01-28", "2026-02-18", "2026-02-18", "2026-01-28")),
  value = c(6124, 1, 4, 6, 8),
  correction_version = as.Date(c("2024-11-21", "2026-02-04", "2026-02-20", "2026-02-20",
    "2026-02-04")),
  correction_value = c(1005, 105, 130, 425, 248),
  reason = c("Erroneous high Pennsylvania COVID report",
    "Incomplete Alabama COVID report", "Incomplete Alabama COVID report",
    "Incomplete Alabama flu report", "Incomplete Alabama flu report"),
  evidence = rep("Reviewed subsequent correction in nhsn_vintages_epidata.parquet", 5),
  cleaning = "retrospective cleaning")

# The table is part of the protocol and own-base cache identity, including its evidence.
report_filter_id <- function() paste0("reviewed-", unname(tools::md5sum(
  bytes = serialize(reviewed_report_exclusions(), NULL, version = 2))))

# Mask exact reviewed versions, retaining every row and its raw count for auditing.
# Selecting the latest version from this table therefore keeps an excluded latest
# version missing, rather than falling back to an older acceptable report.
filter_reports <- function(reports, exclusions = reviewed_report_exclusions()) {
  keys <- c("disease", "location", "target_end_date", "version")
  stopifnot(all(c(keys, "value") %in% names(reports)),
    all(c(keys, "value", "reason", "cleaning") %in% names(exclusions)),
    !anyNA(exclusions[keys]), !anyDuplicated(exclusions[keys]),
    all(valid_report(exclusions$value)), all(nzchar(exclusions$reason)),
    all(exclusions$cleaning %in% c("retrospective cleaning", "real-time cleaning")))
  key <- function(x) do.call(paste, c(x[keys], sep = "|"))
  i <- match(key(reports), key(exclusions))
  if (is.null(reports$raw_value)) reports$raw_value <- reports$value
  matched <- !is.na(i)
  if (any(matched & (is.na(reports$raw_value) |
      reports$raw_value != exclusions$value[i]), na.rm = TRUE))
    stop("A reviewed report has changed in the archive; review its exclusion before use.",
      call. = FALSE)
  reports$report_excluded <- matched
  reports$value <- reports$raw_value
  reports$value[matched] <- NA_real_
  reports
}

# Large revisions, for review only: successive stored versions a -> b of the same target
# week (raw counts, both valid) with max((b + 1) / (a + 1), (a + 1) / (b + 1)) > ratio and
# |b - a| >= min_diff admissions. A flag does not say which version is wrong (the newer
# one is often the correction) and never excludes a report: only a reviewed row in
# reviewed_report_exclusions() does. status: "excluded after review" when either version
# of the pair is excluded there (with its reason and cleaning label), else "flagged for
# review". reports: filter_reports() output (raw_value, report_excluded).
revision_review_flags <- function(reports, exclusions = reviewed_report_exclusions(),
                                  ratio = 2, min_diff = 10) {
  r <- reports[order(reports$disease, reports$location, reports$target_end_date,
    reports$version), , drop = FALSE]
  n <- nrow(r)
  i <- which(c(FALSE, r$disease[-1] == r$disease[-n] & r$location[-1] == r$location[-n] &
    r$target_end_date[-1] == r$target_end_date[-n]))
  a <- r$raw_value[i - 1L]
  b <- r$raw_value[i]
  size <- pmax((b + 1) / (a + 1), (a + 1) / (b + 1))
  flagged <- valid_report(a) & valid_report(b) & size > ratio & abs(b - a) >= min_diff
  i <- i[flagged]
  flags <- data.frame(disease = r$disease[i], location = r$location[i],
    target_end_date = r$target_end_date[i], previous_version = r$version[i - 1L],
    previous_value = r$raw_value[i - 1L], version = r$version[i], value = r$raw_value[i],
    ratio = size[flagged], difference = r$raw_value[i] - r$raw_value[i - 1L])
  key <- function(x, version) paste(x$disease, x$location, x$target_end_date, version)
  e <- match(key(flags, flags$previous_version), key(exclusions, exclusions$version))
  later <- match(key(flags, flags$version), key(exclusions, exclusions$version))
  e[is.na(e)] <- later[is.na(e)]
  flags$status <- ifelse(is.na(e), "flagged for review", "excluded after review")
  flags$excluded_version <- exclusions$version[e]
  flags$reason <- exclusions$reason[e]
  flags$cleaning <- exclusions$cleaning[e]
  rownames(flags) <- NULL
  flags
}

# The delayed-feedback protocol of notes/notes.tex, Section 2.1.1 (wait for the delay, then
# learn once from the report current at the update). R/01_data_overview.R builds every
# case's feedback with feedback_schedule() and stores this protocol with the prepared data;
# R/04 writes its id into the scores. get_case() and R/06 stop on any other id, so results
# of different feedback
# definitions are never combined. Change the id when the protocol changes. The
# initialization is a separate convention (error_scale(), replay_warmup).
feedback_protocol <- function() list(
  id = paste0("latest-at-use-v3-", report_filter_id()),
  report_filter = report_filter_id(),
  exclusions = "Exact reviewed report versions masked; retrospective cleaning using later corrections.",
  cutoff = paste("Reference cutoff of a round: the base's actual generation date where known",
    "(windowed_seasonal; CMU hub rounds that Delphi's production schedule ran late),",
    "otherwise the hub's Wednesday deadline (reference date - 3, also for filled rounds).",
    "Reports have dates only: a report dated on the cutoff day is available (same-day",
    "convention)."),
  vintage = paste("Delay d: each forecast is learned once, at the first later round whose",
    "cutoff is at least its first acceptable report + 7d days and whose current report (the",
    "latest by that cutoff) is valid, from that report; later revisions are ignored. An",
    "invalid current report waits for the next round; an older report is never used",
    "instead. Re-run: the same rows in the same rounds, each with its report at the re-run",
    "round's cutoff; invalid ones are left out. Reviewed excluded versions count as",
    "missing for feedback, training, initialization, scoring and peak detection."))
# Stop unless `found` (the protocol id of prepared data, tuned settings or scores, read
# from `what`) is the current one.
check_feedback_protocol <- function(found, what, rerun) {
  found <- unique(found)
  if (!length(found)) found <- "none (made before the protocol was versioned)"
  if (!identical(found, feedback_protocol()$id)) stop(what, " uses feedback protocol ",
    paste(found, collapse = ", "), ", not ", feedback_protocol()$id, ". Rerun ", rerun, ".",
    call. = FALSE)
  invisible(TRUE)
}
# A report is valid if it is a finite, nonnegative count; an NA in a stored version is an
# explicitly missing report.
valid_report <- function(x) is.finite(x) & x >= 0

# The report of each of `rows` that is current at `cutoff` (one date): the latest version
# published by then. versions: every stored report of a case, columns row, version and
# value, sorted by row and version. Value NA: no report yet, or a missing latest report
# (never replaced by an older valid one).
report_snapshot <- function(versions, rows, cutoff) {
  stopifnot(length(cutoff) == 1L, !is.na(cutoff))
  known <- versions[versions$version <= cutoff & versions$row %in% rows, , drop = FALSE]
  known <- known[!duplicated(known$row, fromLast = TRUE), , drop = FALSE]
  i <- match(rows, known$row)
  data.frame(row = rows, version = known$version[i], value = known$value[i])
}

# The delayed feedback of one case with a delay of `delay` weeks (notes/notes.tex, Section
# 2.1.1). versions: every stored report of each row's target week (row, version, value and
# optionally source; sorted by row and version); cutoff: the issue cutoff of each row's
# round, non-decreasing. For each row s:
#   first_report  nu_s, the date of its first valid report;
#   maturity      m_s = nu_s + 7 * delay days: calendar weeks after the first report,
#                 not forecast rounds, so neither a missing round nor a revision restarts
#                 this clock;
#   consumed_row  tau(s): the first later row t > s whose cutoff is at least m_s and whose
#                 current report there (the latest version on or before the cutoff, a
#                 report of the cutoff day included) is valid (NA: none, no feedback). An
#                 invalid current report waits for the next round (invalid_cutoffs counts
#                 these rounds); an older valid report is never used instead. The strict
#                 t > s keeps a forecast's own error out of it, also at horizon -1, whose
#                 first report is often known at its own issue;
#   selected      that current report, learned once at tau(s); later revisions are ignored.
# Returns observed (the selected values of the consumed rows, NA otherwise), date (the
# reference date of each row's consumption round), by_week (the rows consumed at each row,
# before its forecast is issued) and audit.
feedback_schedule <- function(versions, cutoff, delay, reference_date = NULL) {
  n <- length(cutoff)
  stopifnot(length(delay) == 1L, is.finite(delay), delay >= 0, !anyNA(cutoff),
    !is.unsorted(cutoff), all(versions$row %in% seq_len(n)), !anyNA(versions$version),
    identical(order(versions$row, versions$version), seq_len(nrow(versions))),
    !anyDuplicated(versions[c("row", "version")]),
    is.null(reference_date) || length(reference_date) == n)
  source <- if (is.null(versions$source)) rep("archive", nrow(versions)) else versions$source
  first_report <- maturity <- consumed_cutoff <- selected_version <- rep(as.Date(NA), n)
  selected_value <- rep(NA_real_, n)
  first_report_source <- selected_source <- rep(NA_character_, n)
  consumed_row <- invalid_cutoffs <- rep(NA_integer_, n)
  by_row <- split(seq_len(nrow(versions)), factor(versions$row, seq_len(n)))
  for (s in seq_len(n)) {
    k <- by_row[[s]]
    ok <- valid_report(versions$value[k])
    first <- match(TRUE, ok)
    if (is.na(first)) next                        # Never a valid report: no feedback.
    v <- versions$version[k]
    first_report[s] <- v[first]
    first_report_source[s] <- source[k[first]]
    maturity[s] <- v[first] + 7 * delay
    # Eligible rounds: later rows whose cutoff is at least the maturity, each with its
    # current report (the latest version by its cutoff; at least the first valid one).
    eligible <- which(seq_len(n) > s & cutoff >= maturity[s])
    current <- findInterval(as.numeric(cutoff[eligible]), as.numeric(v))
    use <- match(TRUE, ok[current])
    invalid_cutoffs[s] <- if (is.na(use)) length(eligible) else use - 1L
    if (is.na(use)) next                          # No later round with a valid report.
    consumed_row[s] <- eligible[use]
    consumed_cutoff[s] <- cutoff[eligible[use]]
    selected_version[s] <- v[current[use]]
    selected_value[s] <- versions$value[k[current[use]]]
    selected_source[s] <- source[k[current[use]]]
  }
  consumed <- which(!is.na(consumed_row))
  by_week <- lapply(seq_len(n), function(t) which(consumed_row == t))
  stopifnot(all(consumed_row[consumed] > consumed),
    all(maturity[consumed] <= consumed_cutoff[consumed]),
    all(selected_version[consumed] <= consumed_cutoff[consumed]),
    all(valid_report(selected_value[consumed])), !anyDuplicated(unlist(by_week)))
  observed <- rep(NA_real_, n)
  observed[consumed] <- selected_value[consumed]
  date <- if (is.null(reference_date)) rep(as.Date(NA), n) else reference_date[consumed_row]
  list(observed = observed, date = date, by_week = by_week,
    audit = data.frame(row = seq_len(n), issue_cutoff = cutoff, delay = as.integer(delay),
      first_report, first_report_source, maturity, consumed_row, consumed_date = date,
      consumed_cutoff, invalid_cutoffs, selected_version, selected_value, selected_source))
}
# Every case's feedback schedule in one table (one row per case row and delay), with the
# reference and target dates, for audits: feedback_audit(read_prepared(multiqt_settings())).
feedback_audit <- function(prepared) do.call(rbind, lapply(prepared$cases, function(case)
  do.call(rbind, lapply(case$feedback_by_delay, function(feedback) data.frame(
    disease = case$disease, location = case$location, horizon = case$horizon,
    reference_date = case$reference_date, target_date = case$target_date,
    filled = case$filled, cutoff_source = case$cutoff_source, feedback$audit)))))

# Prepared cases are named disease_location_horizon, e.g. covid_06_0. feedback_delay
# selects the stateful feedback of feedback_schedule(): each forecast is learned once, at
# the first later round at least feedback_delay weeks after its first valid report whose
# current report is valid, from that report.
# revisions = TRUE: the weekly re-run design (Delphi's "exact" re-run). At every issue week
# each method is re-run from the first week, with the same rows in the same rounds but each
# with its report at that week's cutoff (rerun_methods(), case_as_of()).
get_case <- function(prepared, disease, location, horizon, feedback_delay = 0L,
                     revisions = FALSE, corrected = FALSE) {
  check_feedback_protocol(prepared$settings$feedback_protocol$id,
    "output/data/peak_aligned.rds", "Rscript R/01_data_overview.R")
  case <- prepared$cases[[paste(disease, location, horizon, sep = "_")]]
  if (is.null(case)) stop("No prepared case ", disease, " ", location, " h", horizon,
    ". Rerun Rscript R/01_data_overview.R.")
  feedback <- case$feedback_by_delay[[as.character(feedback_delay)]]
  if (is.null(feedback)) stop("No feedback with a ", feedback_delay, "-week delay. Add it to ",
    "feedback_delays in R/01_data_overview.R and rerun it.")
  case$feedback_protocol <- prepared$settings$feedback_protocol$id
  case$feedback_delay <- feedback_delay
  case$feedback_revised <- revisions
  case$feedback_observed <- feedback$observed
  case$feedback_date <- feedback$date
  case$feedback_by_week <- feedback$by_week
  # Cross-horizon inputs for the soft-switch method. Keep the workflow's selected
  # base and its original availability/filling flags; no outcomes enter this panel.
  case$context_forecasts <- do.call(rbind, lapply(-1:3, function(h) {
    other <- prepared$cases[[paste(disease, location, h, sep = "_")]]
    if (is.null(other) || !identical(other$base_name, case$base_name)) return(NULL)
    mid <- match(.5, other$alpha)
    if (is.na(mid)) return(NULL)
    data.frame(reference_date = other$reference_date, target_date = other$target_date,
      as_of = other$as_of, horizon = h, median = other$base[, mid], filled = other$filled)
  }))
  # corrected = TRUE: every consumed row whose report is young at its consumption cutoff
  # learns plausible counts instead (correct_labels()); the revisions come from the case
  # itself and from the pooled locations' cases of the same disease and horizon.
  case$feedback_corrected <- corrected
  if (corrected) {
    settings <- multiqt_settings()
    case$label_correction <- settings$label_correction
    case$label_reference <- list(own = report_history(case),
      pool = Filter(Negate(is.null), lapply(settings$pooled_locations, function(l) {
        other <- prepared$cases[[paste(disease, l, horizon, sep = "_")]]
        if (!is.null(other)) report_history(other)
      })))
    consumed <- unlist(feedback$by_week)
    labels <- correct_labels(case, case$feedback_observed, consumed,
      feedback$audit$consumed_cutoff[consumed])
    case$feedback_observed <- labels$observed
    case$feedback_samples <- labels$samples
  }
  case$phase <- phase_inputs(case)
  case
}
# Phase inputs of a case (phase_angle, phase_hidden; R/05_phase_probabilities.R), or NULL for
# locations without them: p_bias, per row the probability that the report falls below the base
# median (known at the row's cutoff); P, per row the hidden-phase probabilities (rise, overshoot,
# tail) of its target week at its cutoff; W[s, t, ], the probabilities of row s's target week as
# known at the cutoff of issue week t (updated once reported; the weights when row s is learned
# at week t). The file is read once per session.
# source, md5: the file read (multiqt_settings()$phase_probabilities) and its checksum.
phase_cache <- new.env()
phase_inputs <- function(case) {
  file <- multiqt_settings()$phase_probabilities
  if (!file.exists(file)) return(NULL)
  stamp <- paste(file, file.mtime(file))
  if (!identical(phase_cache$stamp, stamp)) {
    x <- readRDS(file)
    phase_cache$bias <- split(x$bias, paste(x$bias$disease, x$bias$location, x$bias$horizon))
    phase_cache$hidden <- split(x$hidden, paste(x$hidden$disease, x$hidden$location))
    phase_cache$md5 <- unname(tools::md5sum(file))
    phase_cache$stamp <- stamp
  }
  bias <- phase_cache$bias[[paste(case$disease, case$location, case$horizon)]]
  hidden <- phase_cache$hidden[[paste(case$disease, case$location)]]
  if (is.null(bias) || is.null(hidden)) return(NULL)
  n <- nrow(case$base)
  stopifnot(nrow(bias) == n)
  classes <- c("rise", "over", "tail")
  lookup <- function(cutoff, week) {
    i <- match(paste(cutoff, week), paste(hidden$cutoff, hidden$week))
    out <- as.matrix(hidden[i, classes]); out[is.na(i), ] <- 1 / 3
    pmin(pmax(out, .005), 1)
  }
  pairs <- expand.grid(s = seq_len(n), t = seq_len(n))
  list(p_bias = bias$p_bias[order(bias$row)],
    P = lookup(case$as_of, case$target_date),
    W = array(lookup(case$as_of[pairs$t], case$target_date[pairs$s]), c(n, n, 3)),
    source = file, md5 = phase_cache$md5)
}
# Label of a feedback design.
feedback_label <- function(delay, revisions = FALSE, corrected = FALSE)
  paste0(ifelse(revisions, ifelse(delay > 0, paste("re-run, delay", delay), "re-run"),
    paste("delay", delay)), ifelse(corrected, ", corrected", ""))

# Corrected labels (notes, "Revised reports as labels: corrected labels"). A report y_s of
# age a (days since the row's first valid report) younger than label_correction$maturity
# weeks is replaced by plausible counts (1 + y_s) exp(d_j) - 1, with d_j = log((1 + y_j now)
# / (1 + y_j at age a)) the revisions of the last `history` rows j whose first report is at
# least `maturity` weeks old at the cutoff: everything known at the cutoff. The location's
# own rows once it has `own_rows` of them, else those of the pooled locations (same disease
# and horizon); fewer than `min_rows`: the report itself. Trackers learn the share of
# plausible counts at or below the issued forecast (covered_share()), mixtures their mean
# quantile loss; S0, scales and rates use their median, the point nowcast.
# Every version of each row of a case and its first valid report.
report_history <- function(case) {
  v <- case$report_versions[order(case$report_versions$row, case$report_versions$version), ]
  rows <- split(v[c("version", "value")], factor(v$row, seq_len(nrow(case$base))))
  first <- vapply(rows, function(r) {
    ok <- valid_report(r$value)
    if (any(ok)) as.numeric(r$version[which(ok)[1]]) else NA_real_
  }, numeric(1))
  list(rows = rows, first = as.Date(first, origin = "1970-01-01"))
}
# The report of row j current at `date`: NA before its first version or if the current
# version is invalid (never an older one, as in feedback_schedule()).
report_at <- function(history, j, date) {
  r <- history$rows[[j]]
  k <- findInterval(as.numeric(date), as.numeric(r$version))
  if (k == 0L || !valid_report(r$value[k])) NA_real_ else r$value[k]
}
# The revisions d_j of a report of age `age` days at `cutoff` (log1p units).
revision_sample <- function(history, cutoff, age, maturity, rows) {
  j <- tail(which(!is.na(history$first) & history$first + 7 * maturity <= cutoff), rows)
  young <- vapply(j, function(i) report_at(history, i, history$first[i] + age), numeric(1))
  now <- vapply(j, function(i) report_at(history, i, cutoff), numeric(1))
  ok <- is.finite(young) & is.finite(now)
  log1p(now[ok]) - log1p(young[ok])
}
# The plausible counts of report y of age `age` at `cutoff`, or NULL (keep the report).
plausible_counts <- function(reference, y, cutoff, age, settings) {
  d <- revision_sample(reference$own, cutoff, age, settings$maturity, settings$history)
  if (length(d) < settings$own_rows) d <- unlist(lapply(reference$pool, revision_sample,
    cutoff, age, settings$maturity, settings$history))
  if (length(d) < settings$min_rows) return(NULL)
  pmax(expm1(log1p(y) + d), 0)
}
# Corrected labels of the learned `rows` (reports in `observed`, learned at `cutoffs`):
# observed with the median plausible count of each corrected row, and samples, a list over
# all rows with the plausible counts of the corrected rows (NULL for the others).
correct_labels <- function(case, observed, rows, cutoffs) {
  samples <- vector("list", length(observed))
  for (i in seq_along(rows)) {
    s <- rows[i]
    age <- as.numeric(cutoffs[i] - case$label_reference$own$first[s])
    if (!is.finite(observed[s]) || !is.finite(age) ||
      age >= 7 * case$label_correction$maturity) next
    x <- plausible_counts(case$label_reference, observed[s], cutoffs[i], age,
      case$label_correction)
    if (is.null(x)) next
    samples[s] <- list(x)
    observed[s] <- median(x)
  }
  list(observed = observed, samples = samples)
}
# The coverage a tracker learns from row s at its issued forecast q (one value per level):
# the indicator of the report, or the share of its plausible counts at or below q.
covered_share <- function(y, samples, q) if (is.null(samples)) y <= q else
  vapply(q, function(z) mean(samples <= z), numeric(1))
# Plausible counts in other units (divided by `divisor`); NULL stays NULL.
scale_samples <- function(samples, divisor) if (is.null(samples)) NULL else
  lapply(samples, function(x) if (!is.null(x)) x / divisor)

# The case as re-run at issue week `week` (re-run design): its first `week` rows. The rows
# consumed by then keep their rounds of the stateful schedule (every delay; later
# revisions never move them), but each is learned from its report current at this week's
# cutoff (case$as_of[week]); a row whose current report is invalid is left out of this
# replay (rerun_omitted), not replaced by an older or later report. The feedback of every
# delay is replayed alike. S0 and the weighting loss scale stay frozen
# (frozen_s0, frozen_loss_scale: those of the stateful design with the same delay, see
# rerun_methods()).
case_as_of <- function(case, week, frozen_s0 = NULL, frozen_loss_scale = NULL) {
  rows <- seq_len(week)
  current <- report_snapshot(case$report_versions, rows, case$as_of[week])$value
  replay <- function(feedback) {
    by_week <- feedback$by_week[rows]
    released <- unlist(by_week)
    omitted <- released[!valid_report(current[released])]
    observed <- rep(NA_real_, week)
    kept <- setdiff(released, omitted)
    observed[kept] <- current[kept]
    list(observed = observed, date = feedback$date[rows],
      by_week = lapply(by_week, setdiff, omitted), omitted = omitted)
  }
  case$feedback_by_delay <- lapply(case$feedback_by_delay, replay)
  own <- case$feedback_by_delay[[as.character(case$feedback_delay)]]
  # Corrected labels: every replayed row whose report is young at this week's cutoff.
  samples <- NULL
  if (isTRUE(case$feedback_corrected)) {
    released <- unlist(own$by_week)
    labels <- correct_labels(case, own$observed, released,
      rep(case$as_of[week], length(released)))
    own$observed <- labels$observed
    samples <- labels$samples
  }
  for (name in c("reference_date", "target_date", "truth_finalized", "filled", "filled_from",
    "truth_latest", "as_of", "cutoff_source", "period", "peak_index", "peak_date",
    "relative_week")) case[[name]] <- case[[name]][rows]
  case$base <- case$base[rows, , drop = FALSE]
  # Phase inputs: the re-run learns every row with its phase probabilities as known at this
  # week's cutoff.
  if (!is.null(case$phase)) {
    W <- case$phase$W[rows, rows, , drop = FALSE]
    for (w in rows) W[, w, ] <- case$phase$W[rows, week, ]
    case$phase <- list(p_bias = case$phase$p_bias[rows],
      P = case$phase$P[rows, , drop = FALSE], W = W, source = case$phase$source,
      md5 = case$phase$md5)
  }
  case$feedback_observed <- own$observed
  case$feedback_samples <- samples
  case$feedback_date <- own$date
  case$feedback_by_week <- own$by_week
  case$rerun_omitted <- own$omitted
  case$feedback_revised <- FALSE
  case$frozen_s0 <- frozen_s0
  case$frozen_loss_scale <- frozen_loss_scale
  case
}

# The log1p median of the base forecast for each (reference date, horizon) pair, NA unless
# its exact reference and target dates exist, it was available by `cutoff` and it is not
# filled. history: case$context_forecasts (get_case()).
forecast_lookup <- function(history) {
  stopifnot(!anyDuplicated(history[c("reference_date", "horizon")]),
    all(history$target_date == history$reference_date + 7 * history$horizon))
  keys <- paste(history$reference_date, history$horizon)
  function(reference, horizon, cutoff) {
    j <- match(paste(reference, horizon), keys)
    ok <- !is.na(j) & !is.na(history$as_of[j]) & !is.na(cutoff) &
      history$as_of[j] <= cutoff & !is.na(history$filled[j]) & !history$filled[j] &
      is.finite(history$median[j]) & history$median[j] >= 0
    value <- rep(NA_real_, length(j))
    value[ok] <- log1p(history$median[j[ok]])
    value
  }
}
# Plotting-only signals from the former transition_gate() (commit db1262d), in
# log1p(counts). Path: largest adjacent decrease in the current forecast curve;
# nowcast: max(-growth, min(previous growth, previous growth - growth)); revision:
# min(-median revisions, -median paired revision changes). Each revision compares
# forecasts for the same target. Prior revisions use the prior reference-week cutoff.
# Indicators exceed log1p(threshold); missing evidence stays NA, not FALSE.
forecast_path_indicators <- function(case, threshold = .05) {
  stopifnot(length(threshold) == 1L, is.finite(threshold), threshold >= 0)
  dates <- case$reference_date
  signals <- matrix(NA_real_, length(dates), 3L,
    dimnames = list(NULL, c("path", "nowcast", "revision")))
  history <- case$context_forecasts
  if (!is.null(history) && nrow(history)) {
    lookup <- forecast_lookup(history)
    keys <- paste(history$reference_date, history$horizon)
    for (i in seq_along(dates)) {
      t <- dates[i]; cutoff <- case$as_of[i]
      current <- lookup(t, -1:3, cutoff)
      # Horizons 0:3 are required; include the -1 to 0 edge only when available.
      if (all(is.finite(current[2:5]))) signals[i, "path"] <-
        max(-diff(current), na.rm = TRUE)
      previous_row <- match(paste(t - 7, -1), keys)
      previous_cutoff <- if (is.na(previous_row)) as.Date(NA) else history$as_of[previous_row]
      if (is.na(cutoff) || is.na(previous_cutoff) || previous_cutoff > cutoff) next
      previous <- lookup(t - 7, -1, previous_cutoff)
      before_previous <- lookup(t - 14, -1, previous_cutoff)
      growth <- current[1] - previous
      previous_growth <- previous - before_previous
      if (is.finite(growth) && is.finite(previous_growth))
        signals[i, "nowcast"] <- max(-growth, min(previous_growth, previous_growth - growth))
      revision <- current[1] - lookup(t - 7 * (1:4), 0:3, cutoff)
      previous_revision <- previous - lookup(t - 7 * (2:5), 0:3, previous_cutoff)
      common <- is.finite(revision) & is.finite(previous_revision)
      if (sum(common) >= 2L) signals[i, "revision"] <- min(-median(revision[common]),
        -median(revision[common] - previous_revision[common]))
    }
  }
  data.frame(reference_date = dates, signals,
    path_indicator = signals[, "path"] > log1p(threshold),
    nowcast_indicator = signals[, "nowcast"] > log1p(threshold),
    revision_indicator = signals[, "revision"] > log1p(threshold))
}

# Tracker units: "log1p", "identity" (counts) or "sqrt" (Delphi's REF-op; negative
# tracker values map back to 0, as in Delphi's calibration code).
to_tracker <- function(x, transform) {
  stopifnot(transform %in% c("log1p", "identity", "sqrt"))
  switch(transform, log1p = log1p(x), identity = x, sqrt = sqrt(x))
}
from_tracker <- function(x, transform) {
  stopifnot(transform %in% c("log1p", "identity", "sqrt"))
  switch(transform, log1p = expm1(x), identity = x, sqrt = pmax(x, 0)^2)
}

# Count scale for scores and plots; MultiQT never tracks it, and the expert weights
# divide it out (weighting_loss_scale()). "population": counts per 100k residents
# (2019 census, stored by R/01_data_overview.R). "max" and "median": of the reports
# published before the location's first forecast. "none": counts.
choose_count_scale <- function(case, method = "population") {
  if (method == "population") {
    if (is.null(case$population)) stop("No population in the prepared data. ",
      "Rerun Rscript R/01_data_overview.R.")
    return(list(scale = case$population / 1e5, used = method, method = method,
      reference_weeks = NA_integer_))
  }
  history <- case$scale_history[is.finite(case$scale_history)]
  scale <- 1
  used <- method
  if (method != "none" && length(history)) {
    scale <- if (method == "median") median(history) else max(history)
    if (scale <= 0) {
      scale <- max(history)
      used <- "max"
    }
  }
  if (!is.finite(scale) || scale <= 0 || !length(history)) {
    scale <- 1
    used <- "none"
  }
  list(scale = scale, used = used, method = method, reference_weeks = length(history))
}
count_scale_label <- function(method) switch(method,
  population = "counts per 100k residents",
  none = "counts",
  paste0("counts / historical ", method))
# The pooled summaries of 02, 04 and 06 and the tuning of 05 use settings$pooled_locations
# only. pooled_rows() keeps their rows of a score table; pooled_label() names the pooled
# set in titles.
pooled_rows <- function(rows, settings) rows[rows$location %in% settings$pooled_locations, ,
  drop = FALSE]
pooled_label <- function(n, settings) paste(n, "largest locations without Florida")

# Frozen error scale S0: Q90 of the absolute base errors in tracker units, pooled
# over all quantile levels of the first `rounds` released rows. eta = m * S0.
# https://delphi-forecasting-reports.netlify.app/calibration_experiments/e01_eta_settings
# The whole release batch of the freeze week is consumed as warmup (no replay), but
# only its first `rounds` rows enter S0. With fewer releases every week is warmup.
# replay = TRUE keeps S0 but has no warmup weeks: MultiQT learns from the first
# week on, so S0 uses reports released later (a retrospective initialization, declared
# separately from the feedback protocol: forecasts before the 20th release depend on it).
# The released rows and values are those of the case's stateful schedule (the reports
# learned from); the re-run design keeps this S0 of its stateful design (rerun_methods()).
error_scale <- function(case, rounds = 20L, floor = 1e-6, transform = "log1p",
                        replay = FALSE) {
  feedback_by_week <- case$feedback_by_week
  week <- match(TRUE, cumsum(lengths(feedback_by_week)) >= rounds)
  initialized <- !is.na(week)
  freeze_week <- if (initialized) week else length(feedback_by_week)
  released <- Reduce(c, feedback_by_week[seq_len(freeze_week)], init = integer())
  rows <- head(released, rounds)
  raw_scale <- scale <- NA_real_
  if (initialized) {
    raw_scale <- abs_error_q90(case, rows, transform)
    scale <- max(floor, raw_scale)
  }
  warmup_weeks <- if (replay && initialized) 0L else freeze_week
  list(initialized = initialized, warmup_weeks = warmup_weeks,
    date = if (initialized) case$reference_date[freeze_week] else as.Date(NA),
    rows = rows, consumed_rows = if (warmup_weeks > 0L) released else integer(),
    raw_scale = raw_scale, scale = scale, rounds = rounds, floor = floor)
}
# Q90 of the absolute base errors in tracker units, pooled over all quantile levels of
# the given rows, against the reports learned from: S0 above, Delphi's adaptive rate below.
abs_error_q90 <- function(case, rows, transform) {
  residuals <- sweep(to_tracker(case$base[rows, , drop = FALSE], transform), 1,
    to_tracker(case$feedback_observed[rows], transform), "-")
  as.numeric(quantile(abs(residuals), .9, names = FALSE))
}

# Delphi's adaptive rate ("adaptive+" in their calibration code, used by REF-op): at each
# issue week, mult * abs_error_q90() of the last `rounds` rows learned by then (Inf: all
# of them), at least floor. The statistic of S0, but rolling with the feedback instead of
# frozen, so it needs no warmup and no hindsight. 0 before the first learned row (no
# update then). In the re-run design it uses the reports current at the issue week.
adaptive_rate <- function(case, mult, rounds = 20L, floor = 0, transform = "sqrt") {
  released <- integer()
  eta <- numeric(nrow(case$base))
  for (week in seq_along(eta)) {
    released <- c(released, case$feedback_by_week[[week]])
    if (length(released)) eta[week] <- max(floor, mult * abs_error_q90(case,
      if (is.finite(rounds)) tail(released, rounds) else released, transform))
  }
  eta
}

# The adaptive rate of every tracker (notes, "Learning rates"): Delphi's rule for any
# score, mult times the 90% quantile of the absolute scores |S_s^alpha(y_s)| of the last
# `rounds` learned rows, pooled over the levels. score(rows) returns the scores [row x level]
# of learned rows, each with the report learned and the scale of its own cutoff; rows
# without a score (before a scale exists) are left out. replay: until `rounds` rows are
# learned, the first `rounds` learned rows in hindsight (the replay of the first 20 weeks,
# as for S0), so the rate never rests on fewer rows; FALSE: the rows learned so far (0
# before the first). In the re-run design the reports current at the issue week.
rolling_rate <- function(case, score, mult = 1, rounds = 20L, replay = TRUE, floor = 0) {
  released <- unlist(case$feedback_by_week)
  n_released <- cumsum(lengths(case$feedback_by_week))
  eta <- numeric(length(n_released))
  if (!length(released)) return(eta)
  scores <- abs(score(released))
  usable <- which(apply(is.finite(scores), 1, all))      # Positions in release order.
  for (week in seq_along(eta)) {
    k <- sum(usable <= n_released[week])
    use <- if (k >= rounds) usable[(k - rounds + 1L):k] else if (replay)
      head(usable, rounds) else usable[seq_len(k)]
    if (length(use)) eta[week] <- max(floor, mult * quantile(scores[use, ], .9,
      names = FALSE))
  }
  eta
}
# The scores of learned rows for rolling_rate(): the residual in tracker units relative to the
# base (or another centre [week x level]), divided by sigma (one scale per week or one per
# week and level) for a studentized tracker.
tracker_score <- function(case, transform, sigma = NULL, centre = case$base) function(rows) {
  r <- to_tracker(case$feedback_observed[rows], transform) -
    to_tracker(centre[rows, , drop = FALSE], transform)
  if (is.null(sigma)) r else r / (if (is.matrix(sigma)) sigma[rows, , drop = FALSE] else
    sigma[rows])
}
# The scores of the PIT tracker: PIT of the report under the base minus alpha.
pit_score <- function(case) function(rows) t(vapply(rows, function(s)
  base_pit(case$feedback_observed[s], case$base[s, ], case$alpha) - case$alpha,
  numeric(length(case$alpha))))

# Reference loss that makes the expert weights unit-free: the base's mean pinball loss
# over the S0 rows (error_scale()), against the reports learned from, in counts /
# count_scale like the weighting losses. exp_weights() divides its losses by it, so
# lambda_start and fixed lambdas act alike under every count scale and location.
# Like S0, it uses the S0 rows in hindsight when the warmup is replayed. 1 without S0
# (then every week is warmup and the weights never update).
weighting_loss_scale <- function(case, scale, count_scale = 1, floor = 1e-6) {
  if (!scale$initialized) return(1)
  rows <- scale$rows
  error <- (case$feedback_observed[rows] - case$base[rows, , drop = FALSE]) / count_scale
  loss <- pmax(error, 0) * rep(case$alpha, each = length(rows)) +
    pmax(-error, 0) * rep(1 - case$alpha, each = length(rows))
  max(floor, mean(loss))
}

# Scale of recent residuals r: their root mean square, or the median absolute residual
# divided by 0.6745 (a normal sd that one outlier cannot inflate).
residual_scale <- function(r, estimator = "rms") switch(estimator,
  rms = sqrt(mean(r^2)),
  median = median(abs(r)) / qnorm(.75),
  stop("Unknown residual scale estimator ", estimator))

# Scale sigma_t of the median residuals y - b (tracker units) for the studentized
# median MultiQT: one estimate of their standard deviation per issue week.
#   "residuals" root mean square of the last `window` median residuals released by
#               issue week t (only released reports, no hindsight). NA until min_rows
#               are released. A small window adapts quickly but is noisier. Reports
#               arrive weeks after issue, so this scale still reacts to a wave late.
#   "base"      the sd implied by the issued base: the width between its quantile `levels`
#               divided by that of a standard normal, e.g. (q.9 - q.1) / 2.563 for the 80%
#               interval. No feedback lag; a point-mass base gives the floor.
# warmup_weeks: the weeks before the first sigma_t; the method issues the original base there.
# estimator (residuals): "rms", or "median": median |residual| / 0.6745 (the normal sd), which
# one outlying residual (e.g. from a broken base forecast) cannot inflate.
median_sigma <- function(case, method = "residuals", window = 5L, min_rows = 5L,
                         floor = 1e-6, transform = "log1p", levels = c(.1, .9),
                         estimator = "rms") {
  stopifnot(window >= 1L, min_rows >= 1L)   # Otherwise a later report enters sigma_t.
  base_level <- function(alpha) to_tracker(case$base[, match(alpha, case$alpha)], transform)
  scale <- switch(method,
    residuals = {
      released <- unlist(case$feedback_by_week)              # Rows in release order.
      residual <- to_tracker(case$feedback_observed[released], transform) -
        base_level(.5)[released]
      n_released <- cumsum(lengths(case$feedback_by_week))   # By each issue week.
      vapply(n_released, function(k) if (k < min_rows) NA_real_ else
        max(floor, residual_scale(residual[max(1L, k - window + 1L):k], estimator)), numeric(1))
    },
    base = {
      stopifnot(length(levels) == 2L, levels[1] < levels[2], all(levels %in% case$alpha))
      pmax(floor, (base_level(levels[2]) - base_level(levels[1])) / diff(qnorm(levels)))
    },
    stop("Unknown median sigma ", method))
  first <- match(TRUE, is.finite(scale))
  list(scale = scale, warmup_weeks = if (is.na(first)) length(scale) else first - 1L,
    method = method, window = window, min_rows = min_rows, floor = floor)
}

# Quantile-specific scales sigma_t^alpha from recent residuals (tracker units): for every
# level, the RMS (or robust scale, estimator as in median_sigma()) of y_s - b_s^alpha over
# the last `window` rows s released by issue week t (no hindsight). Rows are NA until
# min_rows are released, like median_sigma("residuals").
# Returns a matrix [week x level]. These scales often leave the inward-flow band (notes,
# "Common and constant scales"); studentized_level_rms uses them as they are.
level_residual_sigma <- function(case, window = 5L, min_rows = 5L, floor = 1e-6,
                                 transform = "log1p", estimator = "rms") {
  stopifnot(window >= 1L, min_rows >= 1L)   # Otherwise a later report enters sigma_t.
  released <- unlist(case$feedback_by_week)                # Rows in release order.
  residual <- to_tracker(case$feedback_observed[released], transform) -
    to_tracker(case$base[released, , drop = FALSE], transform)
  n_released <- cumsum(lengths(case$feedback_by_week))     # By each issue week.
  scale <- t(vapply(n_released, function(k) if (k < min_rows) rep(NA_real_, ncol(residual)) else
    pmax(apply(residual[max(1L, k - window + 1L):k, , drop = FALSE], 2, residual_scale,
      estimator), floor),
    numeric(ncol(residual))))
  dimnames(scale) <- dimnames(case$base)
  scale
}

# Scales of the studentized methods chosen for the most stationary scores (notes, "Which
# scale of the residual is most stationary?"), each a function of the forecast and of
# reports known at the cutoff, in counts (identity tracker).
# power_sigma(): the forecast median to the power p, (1 + b_t^.5)^p, one scale per week.
power_sigma <- function(case, p = .5) (1 + case$base[, match(.5, case$alpha)])^p
# spread_sigma(): per level, the derivative scale of a symmetric Box-Cox transform centred
# at the forecast median and scaled by half the forecast's interquartile range s_t,
# sigma_t^alpha = s_t (1 + |b_t^alpha - b_t^.5| / s_t)^(1 - lambda): U-shaped around the
# median. Matrix [week x level] (multiqt() then uses the weighted projection).
spread_sigma <- function(case, lambda = .4) {
  B <- case$base
  s <- half_iqr_sigma(case)
  s * (1 + abs(B - B[, match(.5, case$alpha)]) / s)^(1 - lambda)
}
# half_iqr_sigma(): half the forecast's interquartile range, s_t = (b_t^.75 - b_t^.25) / 2,
# at least one count, one scale per week (spread_sigma() at lambda = 1). Since 20 November
# 2026 the scale of width_flip's tracker (the base width before).
half_iqr_sigma <- function(case)
  pmax(1, (case$base[, match(.75, case$alpha)] - case$base[, match(.25, case$alpha)]) / 2)
# A scale in counts comparable with S0: divided by its median over the S0 rows (reference)
# and multiplied by S0, so a rate m in units of the scale acts like m * S0 in a typical week
# (the rate grid and the restart rate carry over). Each week first stays within a factor
# `guard` of the median of the scale's values so far. reference: given in the re-run design,
# frozen from the stateful design like S0 (rerun_methods()). Returns the scale (vector or
# matrix), its warmup (the weeks before the first scale, all weeks without S0) and reference.
normalized_scale <- function(sigma, scale, guard = 20, reference = NULL) {
  S <- as.matrix(sigma)
  S[!is.finite(S) | S <= 0] <- NA
  centre <- vapply(seq_len(nrow(S)), function(t)
    suppressWarnings(median(S[seq_len(t), ], na.rm = TRUE)), numeric(1))
  S <- pmin(pmax(S, centre / guard), centre * guard)
  if (is.null(reference)) reference <- median(S[intersect(scale$rows, seq_len(nrow(S))), ],
    na.rm = TRUE)
  if (!isTRUE(scale$initialized) || !is.finite(reference))
    return(list(scale = sigma, warmup_weeks = nrow(S), reference = NA_real_))
  S <- S / reference * scale$scale
  first <- match(TRUE, apply(is.finite(S), 1, all))
  list(scale = if (is.matrix(sigma)) S else S[, 1],
    warmup_weeks = if (is.na(first)) nrow(S) else first - 1L, reference = reference)
}
# The raw scale of one of the scaled methods ("power", "spread"; params:
# method_parameters()).
studentized_scale <- function(case, name, params) switch(name,
  power = power_sigma(case, params$power_p),
  spread = spread_sigma(case, params$spread_lambda))

# Peak-relative windows, as a logical vector over the case's weeks, within the given
# period(s) of the target week: "evaluation", "tuning" or "other" (the dates in
# R/01_data_overview.R). A peak is a local maximum of the latest-reported counts
# (find_peaks() in R/01_data_overview.R; hindsight label only); each week counts from
# its nearest peak. before and after are window lengths:
# before = weeks -(before - 1)..0, including the peak week, and after = weeks
# 1..after; overall = all weeks of the period(s).
phase_weeks <- function(case, phase, before = 8L, after = 8L, period = "evaluation") {
  week <- case$relative_week
  in_period <- case$period %in% period
  switch(phase,
    before = in_period & week > -before & week <= 0,
    after = in_period & week >= 1 & week <= after,
    overall = in_period,
    stop("Unknown phase ", phase))
}

# Missing weeks shorten the before/after windows: report them, but keep the peak.
check_peak_window <- function(case, before = 8L, after = 8L, truth = case$truth_latest,
                              period = "evaluation") {
  label <- paste0(case$disease, " | ", case$location_name, " (", case$location,
    ") | horizon ", case$horizon)
  weeks <- seq.int(1L - before, after)
  ok_all <- TRUE
  for (k in which(case$peaks$period %in% period)) {
    peak <- case$peaks[k, ]
    peak_label <- paste0(label, " | peak ", peak$peak_date)
    rows <- which(case$peak_index == k)
    i <- rows[match(weeks, case$relative_week[rows])]
    ok <- !is.na(i)
    ok[ok] <- case$period[i[ok]] %in% period & is.finite(truth[i[ok]])
    if (!all(ok))
      warning(peak_label, ": no scored target in peak-relative week(s) ",
        paste(weeks[!ok], collapse = ","), "; peak kept.", call. = FALSE, immediate. = TRUE)
    ok_all <- ok_all && all(ok)
  }
  invisible(ok_all)
}

# 3. SCORING --------------------------------------------------------------
# Truth for coverage, coverage residuals and quantile losses. "latest": the latest
# stored report of each target week (the archive's last version, prepared$settings$
# truth_version), the same truth for every feedback design. "feedback": the report the
# stateful updates learned from at the case's feedback delay, current at the update (NA
# for rows never consumed; in the re-run design, too, the stateful design's report).
# Calibration against these learning vintages is not calibration against the latest counts.
# filled = FALSE leaves filled rounds unscored; MultiQT still issues and learns there.
score_truth <- function(case, target = "latest", filled = TRUE) {
  truth <- switch(target,
    latest = case$truth_latest,
    feedback = case$feedback_observed,
    stop("Unknown score target ", target))
  if (!filled) truth[case$filled] <- NA
  truth
}
score_truth_label <- function(target, delay = 0L) switch(target,
  latest = "latest reported counts",
  feedback = paste0("the reports learned from (delay ", delay, ": current at the update)"),
  stop("Unknown score target ", target))

# Coverage and pinball loss of quantile forecasts q [week x level] against one
# truth per week. Pass counts with their count_scale to compare coverage on the
# count scale (exact integer ties stay ties); the loss is in counts / count_scale.
# The pinball loss penalizes underprediction by alpha, overprediction by 1 - alpha.
# score_scale "log": the loss of log(counts + 1) instead (q and truth in counts;
# count_scale unused), which scores relative errors and growth rates (Bosse et al.
# 2023, PLOS Comput Biol 19(8): e1011393). Negative quantiles count as 0 there.
# Quantiles commute with the monotone log, so the score stays proper; coverage is
# the same on both scales.
score_quantiles <- function(q, truth, alpha, count_scale = 1, score_scale = "natural") {
  stopifnot(ncol(q) == length(alpha), length(truth) == nrow(q))
  error <- switch(score_scale,
    natural = truth / count_scale - q / count_scale,
    log = log1p(truth) - log1p(pmax(q, 0)),
    stop("Unknown score scale ", score_scale))
  list(covered = truth <= q,
    loss = pmax(error, 0) * rep(alpha, each = nrow(q)) +
      pmax(-error, 0) * rep(1 - alpha, each = nrow(q)))
}

# Average scores over the given weeks. Calibration error (CE) is checked per level,
# |coverage - alpha|, and then averaged over levels; every level has equal weight. Its two
# signed parts: bias, the mean over levels of coverage - alpha (positive: the forecasts lie
# too high, the outcome falls below them too often), and dispersion, the mean of
# coverage - alpha over the levels above the median minus that over the levels below it
# (positive: the intervals are too wide, negative: too narrow).
summarize_scores <- function(scores, alpha, weeks) {
  if (!length(weeks)) {
    na <- rep(NA_real_, length(alpha))
    return(list(by_level = data.frame(alpha, coverage = na, coverage_error = na,
      quantile_loss = na), coverage_error = NA_real_, bias = NA_real_,
      dispersion = NA_real_, quantile_loss = NA_real_, n_weeks = 0L))
  }
  coverage <- colMeans(scores$covered[weeks, , drop = FALSE])
  loss <- colMeans(scores$loss[weeks, , drop = FALSE])
  by_level <- data.frame(alpha, coverage, coverage_error = abs(coverage - alpha),
    quantile_loss = loss)
  gap <- coverage - alpha
  list(by_level = by_level, coverage_error = mean(by_level$coverage_error),
    bias = mean(gap), dispersion = mean(gap[alpha > .5]) - mean(gap[alpha < .5]),
    quantile_loss = mean(loss), n_weeks = length(weeks))
}

# Plots: an NA row after every gap of more than max_days between consecutive dates, so
# lines and polygons break there (e.g. May to November 2024, between the seasons).
# Returns the dates and each further vector or matrix with the same rows.
break_gaps <- function(dates, ..., max_days = 28) {
  gap <- which(diff(as.numeric(dates)) > max_days)
  index <- sort(c(seq_along(dates), gap + .5))
  rows <- ifelse(index %% 1 == 0, index, NA)
  c(list(dates = dates[rows]), lapply(list(...), function(x)
    if (is.matrix(x)) x[rows, , drop = FALSE] else x[rows]))
}

# 4. MULTIQT: STANDARD AND STUDENTIZED -----------------------------------
# Weighted isotonic regression (pool adjacent violators): the increasing q that
# minimizes sum(w * (q - v)^2). Equal weights give isoreg().
weighted_isoreg <- function(v, w) {
  value <- weight <- numeric(length(v))
  size <- integer(length(v))
  k <- 0L
  for (i in seq_along(v)) {
    k <- k + 1L
    value[k] <- v[i]; weight[k] <- w[i]; size[k] <- 1L
    while (k > 1L && value[k - 1L] > value[k]) {         # Pool the violating neighbours.
      total <- weight[k - 1L] + weight[k]
      value[k - 1L] <- (weight[k - 1L] * value[k - 1L] + weight[k] * value[k]) / total
      weight[k - 1L] <- total; size[k - 1L] <- size[k - 1L] + size[k]
      k <- k - 1L
    }
  }
  rep(value[seq_len(k)], size[seq_len(k)])
}

# base: issued quantiles [week x level] in counts; feedback: the reported counts learned
# from (the case's feedback schedule); feedback_by_week[[w]]: the earlier rows consumed at
# issue week w, before its forecast is issued. Each consumed row is compared with ITS
# historical issued forecast, and the steps of a batch are summed (not averaged), all with
# eta[w], the rate of the update week. eta: one rate, or one per week, in tracker units.
# Only the issued forecasts are projected; the offsets stay unprojected. The first
# warmup_weeks forecasts keep the base and their reports are never replayed.
# Optional studentized offsets (sigma = 1, start = 0: standard MultiQT): week w issues
# base + sigma[w] * offsets. sigma: one value, or one per week, in tracker units; then
# offsets and eta are in units of sigma, and a change of sigma rescales the whole
# accumulated offset. A matrix sigma [week x level] gives each level its own scale; the
# projection is then weighted
# isotonic regression with weights 1 / sigma^2 (Euclidean in units of sigma).
# start: the offsets before the first update, one value or one per level.
# Optional leak and slow term (Delphi's qt_track(): fast_decay and lr_slow; fast_decay = 0,
# eta_slow = 0: standard MultiQT). The offsets are multiplied by 1 - fast_decay at every
# issue week after the first, before its update, so a correction that new reports do not
# refresh fades back to the base (Leaky MultiQT, REF-op's fast term). Week w issues
# base + sigma[w] * (offsets + slow); both terms take the same steps, the slow one with
# eta_slow (one rate, or one per week), and the slow term never decays.
# Returns the adjusted forecasts in counts.
multiqt <- function(base, feedback, feedback_by_week, alpha, eta,
                    warmup_weeks = 0L, transform = "log1p", sigma = 1, start = 0,
                    eta_slow = 0, fast_decay = 0, feedback_samples = NULL) {
  n <- nrow(base)
  # Only warmup weeks (e.g. S0 never initialized), or nothing to learn and start = 0:
  # keep the exact base.
  if (warmup_weeks >= n || (isTRUE(all(eta == 0)) && isTRUE(all(eta_slow == 0)) &&
    isTRUE(all(start == 0))))
    return(base)
  if (length(eta) == 1L) eta <- rep(eta, n)
  if (length(eta_slow) == 1L) eta_slow <- rep(eta_slow, n)
  if (length(sigma) == 1L) sigma <- rep(sigma, n)
  by_level <- is.matrix(sigma)              # One scale per level: weighted projection.
  base_tracker <- to_tracker(base, transform)
  feedback_tracker <- to_tracker(feedback, transform)
  warmup <- seq_len(n) <= warmup_weeks      # sigma may be missing in warmup weeks.
  issued_sigma <- if (by_level) sigma[!warmup, , drop = FALSE] else sigma[!warmup]
  stopifnot(ncol(base) == length(alpha), length(feedback) == n,
    length(feedback_by_week) == n, length(eta) == n, all(is.finite(eta)), all(eta >= 0),
    length(eta_slow) == n, all(is.finite(eta_slow)), all(eta_slow >= 0),
    length(fast_decay) == 1L, fast_decay >= 0, fast_decay < 1,
    NROW(sigma) == n, !by_level || ncol(sigma) == length(alpha),
    all(is.finite(issued_sigma) & issued_sigma > 0),
    length(start) %in% c(1L, length(alpha)), all(is.finite(start)),
    all(is.finite(base_tracker)), all(is.finite(feedback_tracker[unlist(feedback_by_week)])))
  offsets <- rep_len(start, length(alpha))
  slow <- rep(0, length(alpha))
  q_tracker <- base_tracker
  for (week in seq_len(n)) {
    if (week <= warmup_weeks) next
    if (fast_decay > 0 && week > warmup_weeks + 1L) offsets <- offsets * (1 - fast_decay)
    for (past_week in feedback_by_week[[week]]) {
      # Compare with the HISTORICAL ISSUED forecast in tracker units.
      # At alpha=.9, eta=.01, a miss adds .009; coverage subtracts .001.
      # Corrected labels: the share of the plausible counts (covered_share()).
      plausible <- feedback_samples[[past_week]]
      covered <- if (is.null(plausible)) feedback_tracker[past_week] <= q_tracker[past_week, ] else
        covered_share(NULL, to_tracker(plausible, transform), q_tracker[past_week, ])
      offsets <- offsets + eta[week] * (alpha - covered)
      slow <- slow + eta_slow[week] * (alpha - covered)
    }
    # Project only the issued forecasts, NOT the offsets ("lazy" MultiQT).
    # isoreg restores quantile order: e.g., (12,10) becomes (11,11).
    # Weights (min sigma / sigma)^2: the same projection as 1 / sigma^2, without overflow.
    q_tracker[week, ] <- if (by_level) weighted_isoreg(base_tracker[week, ] +
      sigma[week, ] * (offsets + slow), (min(sigma[week, ]) / sigma[week, ])^2) else
      isoreg(base_tracker[week, ] + sigma[week] * (offsets + slow))$yf
  }
  # Do not clip negative predictions: that would change the issued coverage indicators.
  # (The sqrt tracker maps negative values back to 0 counts, as Delphi's code does; its
  # coverage indicators above stay those of the unclipped tracker values.)
  q <- from_tracker(q_tracker, transform)
  # Keep unchanged bases exact, so transformation roundoff cannot flip ties.
  unchanged <- q_tracker == base_tracker
  q[unchanged] <- base[unchanged]
  q
}

# Studentized median MultiQT: the base is reduced to its median b_t, and week t issues
# b_t + sigma_t * offsets at every level (tracker units; isoreg restores the order).
# Coverage and update are those of multiqt(): 1{y_s <= q_s} = 1{(y_s - b_s) / sigma_s <=
# issued offset_s}, so offsets and eta are in units of sigma. sigma: one value, or
# sigma_t per week (median_sigma()). start: the offsets before the first update, e.g.
# qnorm(alpha), a normal forecast with sd sigma_t, or qt(alpha, window) for sigma_t from
# the RMS of `window` residuals. The first warmup_weeks forecasts keep the original base,
# like every method; reports on them released later are compared with it.
# Returns the adjusted forecasts in counts.
median_multiqt <- function(case, eta, sigma, start, warmup_weeks = 0L, transform = "log1p") {
  median_base <- case$base
  median_base[] <- case$base[, match(.5, case$alpha)]   # Every level: the base median.
  warmup <- seq_len(nrow(median_base)) <= warmup_weeks
  median_base[warmup, ] <- case$base[warmup, ]
  multiqt(median_base, case$feedback_observed, case$feedback_by_week, case$alpha, eta,
    warmup_weeks, transform, sigma = sigma, start = start,
    feedback_samples = case$feedback_samples)
}

# The base's quantile function at probability levels p (one row b of base quantiles at the
# levels alpha): linear in qnorm(level) between the base quantiles and, beyond the outer
# levels, along the line through the two outer quantiles; at least 0. At p = alpha it
# returns the (ordered) base quantiles.
base_quantile_function <- function(b, alpha, p) {
  z <- qnorm(alpha); m <- length(z); b <- cummax(b)
  zp <- qnorm(p)
  q <- approx(z, b, xout = zp, rule = 2)$y
  lower <- zp < z[1]; upper <- zp > z[m]
  q[lower] <- b[1] + (zp[lower] - z[1]) * (b[2] - b[1]) / (z[2] - z[1])
  q[upper] <- b[m] + (zp[upper] - z[m]) * (b[m] - b[m - 1]) / (z[m] - z[m - 1])
  pmax(q, 0)
}
# The PIT of y under the base (its inverse): the level of base_quantile_function() at
# which the base reaches y, within [floor, 1 - floor].
base_pit <- function(y, b, alpha, floor = 1e-4) {
  z <- qnorm(alpha); m <- length(z); b <- cummax(b)
  if (!is.finite(y) || b[m] <= b[1]) return(NA_real_)
  slope <- function(i, j) (b[j] - b[i]) / (z[j] - z[i])
  zy <- if (y < b[1]) {
    if (slope(1, 2) > 0) z[1] + (y - b[1]) / slope(1, 2) else -Inf
  } else if (y > b[m]) {
    if (slope(m - 1, m) > 0) z[m] + (y - b[m]) / slope(m - 1, m) else Inf
  } else approx(b, z, xout = y, ties = mean)$y
  min(max(pnorm(zy), floor), 1 - floor)
}
# Studentized MultiQT on the PIT of the base (quantile-level tracking, as ACI does): one
# offset per level on the probability scale; week t issues the base quantile function at
# the levels alpha + offsets, ordered (isoreg) and kept within [floor, 1 - floor], and
# each released row moves the offsets by eta * (alpha - covered), with covered = 1{y_s <=
# issued q_s}. The PIT is the same for any increasing transform of the counts, so there
# is no tracker transform and no S0: eta is in probability units. The first warmup_weeks
# forecasts keep the base, as in multiqt(). Returns the forecasts in counts.
pit_multiqt <- function(base, feedback, feedback_by_week, alpha, eta, warmup_weeks = 0L,
                        floor = 1e-4, feedback_samples = NULL) {
  n <- nrow(base)
  if (length(eta) == 1L) eta <- rep(eta, n)
  if (warmup_weeks >= n || all(eta == 0)) return(base)
  stopifnot(ncol(base) == length(alpha), length(feedback) == n,
    length(feedback_by_week) == n, length(eta) == n, all(is.finite(eta)), all(eta >= 0),
    all(is.finite(base)), all(is.finite(feedback[unlist(feedback_by_week)])))
  offsets <- rep(0, length(alpha))
  q <- base
  for (week in seq_len(n)) {
    if (week <= warmup_weeks) next
    for (past_week in feedback_by_week[[week]])
      offsets <- offsets + eta[week] * (alpha - covered_share(feedback[past_week],
        feedback_samples[[past_week]], q[past_week, ]))
    levels <- pmin(pmax(isoreg(alpha + offsets)$yf, floor), 1 - floor)
    # Unchanged levels keep the exact base, so roundoff cannot flip ties.
    if (any(levels != alpha)) q[week, ] <- base_quantile_function(base[week, ], alpha, levels)
  }
  q
}

# Two-phase MultiQT (phase_power, phase_oracle; notes, "A phase probability in the
# score"): two hidden states, one for the rise and one for the decline of a wave, issued as
# b_t + sigma_t ((1 - p_t) offsets_rise + p_t offsets_decline) with p_t the probability that
# week t lies in a decline (known at its cutoff, or the hindsight phase for the oracle). A
# report on row s gives the subgradient of the studentized quantile loss of the issued q_s
# with respect to each state, (1 - p_s) (alpha - covered) and p_s (alpha - covered): each
# state learns from the weeks of its phase. eta per week, sigma per week or per week and
# level (weighted projection); counts (identity tracker). The first warmup_weeks keep the
# base.
phase_multiqt <- function(base, feedback, feedback_by_week, alpha, eta, p, sigma = 1,
                          warmup_weeks = 0L, feedback_samples = NULL) {
  n <- nrow(base)
  if (length(eta) == 1L) eta <- rep(eta, n)
  if (length(sigma) == 1L) sigma <- rep(sigma, n)
  by_level <- is.matrix(sigma)
  stopifnot(length(eta) == n, all(is.finite(eta)), length(p) == n,
    all(is.finite(p) & p >= 0 & p <= 1), NROW(sigma) == n)
  rise <- decline <- rep(0, length(alpha))
  q <- base
  for (week in seq_len(n)) {
    if (week <= warmup_weeks) next
    for (past_week in feedback_by_week[[week]]) {
      step <- eta[week] * (alpha - covered_share(feedback[past_week],
        feedback_samples[[past_week]], q[past_week, ]))
      rise <- rise + (1 - p[past_week]) * step
      decline <- decline + p[past_week] * step
    }
    offsets <- (1 - p[week]) * rise + p[week] * decline
    q[week, ] <- if (by_level) weighted_isoreg(base[week, ] + sigma[week, ] * offsets,
      (min(sigma[week, ]) / sigma[week, ])^2) else isoreg(base[week, ] + sigma[week] * offsets)$yf
  }
  q
}
# The probability that week t lies in the decline of a wave, learned from the recent errors
# of the base: the mean of 1 - PIT (base_pit()) of the reports of the last `rows` rows
# learned by the cutoff, near 1 when the reports fall below the base (it forecasts too high,
# as after a peak), near 0 when they rise above it; 1/2 before the first report. A PIT the
# base cannot give (e.g. a flat nowcast) counts as 1/2, no information.
decline_probability <- function(case, rows = 4L) {
  released <- unlist(case$feedback_by_week)
  n_released <- cumsum(lengths(case$feedback_by_week))
  pit <- vapply(released, function(s) base_pit(case$feedback_observed[s], case$base[s, ],
    case$alpha), numeric(1))
  pit[!is.finite(pit)] <- .5
  vapply(n_released, function(k) if (k == 0) .5 else mean(1 - pit[max(1L, k - rows + 1L):k]),
    numeric(1))
}

# Error scale from the base's own revisions, known at the cutoff and without feedback lag:
# the forecasts of one target week at successive cutoffs disagree roughly as much as the
# base errs. For row t: the median absolute difference, in tracker units, between a base
# median of this cutoff and the medians the base issued for the same target week at earlier
# cutoffs (horizons -1 to 3, so up to 4 weeks earlier). type "cutoff": every target of this
# cutoff (horizons -1 to 2 have earlier forecasts), one scale for all horizons; "target":
# only row t's own target (none at the longest horizon). Filled forecasts and forecasts not
# available by the cutoff are left out; a row without pairs keeps the last scale.
# history: case$context_forecasts (get_case()). A candidate scale of R/07, not a method.
revision_sigma <- function(case, type = "cutoff", floor = 1e-6, transform = "identity") {
  stopifnot(type %in% c("cutoff", "target"))
  history <- case$context_forecasts
  history <- history[!is.na(history$filled) & !history$filled & !is.na(history$as_of) &
    is.finite(history$median) & history$median >= 0, ]
  value <- to_tracker(history$median, transform)
  scale <- rep(NA_real_, length(case$reference_date))
  for (t in seq_along(scale)) {
    known <- history$as_of <= case$as_of[t]
    now <- which(known & history$reference_date == case$reference_date[t])
    if (type == "target") now <- now[history$target_date[now] == case$target_date[t]]
    earlier <- which(known & history$reference_date < case$reference_date[t] &
      history$target_date %in% history$target_date[now])
    if (!length(earlier)) next
    current <- value[now][match(history$target_date[earlier], history$target_date[now])]
    scale[t] <- max(floor, median(abs(value[earlier] - current)))
  }
  for (t in seq_along(scale)[-1]) if (!is.finite(scale[t])) scale[t] <- scale[t - 1L]
  first <- match(TRUE, is.finite(scale))
  list(scale = scale, warmup_weeks = if (is.na(first)) length(scale) else first - 1L,
    type = type)
}

# Scores tested for a stable alpha-quantile (R/07_score_selection.R; notes, "Scores and
# scales tested"). For one case, each candidate turns an outcome y (counts) into a score
# S_t^alpha(y) relative to the base quantile; the alpha-quantile of the realized scores is
# the threshold a tracker on that score needs, constant over time if the score is
# quantile-stable. Three kinds:
#   scale      (y - centre_t^alpha) / sigma_t, sigma one scale per week (vector) or one per
#              week and level (matrix); the centre is the base, or for the median rebuild
#              b_t^.5 + sigma_t qnorm(alpha);
#   transform  psi_t(y) - psi_t(b_t^alpha), psi_t increasing and fixed at the cutoff:
#              f(x, t) and its inverse inv(v, t);
#   pit        F_t(y) - alpha, the PIT of the base (base_pit()).
# Scales from reports use those known at each cutoff (the case's feedback schedule); a scale
# stays within a factor `guard` of the median of its values so far (it binds only where a
# scale explodes). The scales of the base width, the recent errors and
# the revisions keep the bounds sigma_bounds * S0 of their methods. origin: where the idea
# came from; kind: transform, common scale (one per week) or per-level scale; methods: the
# methods that use it (some paused). The spread scale runs over five lambdas; lambda = 1 is
# half the base IQR (since 20 November 2026).
candidate_row <- function(id, label, origin, kind, methods = "")
  data.frame(id = id, label = label, origin = origin, kind = kind, methods = methods)
score_candidates <- rbind(
  candidate_row("identity", "Identity", "Trackers in use", "transform", "fixed, leaky"),
  candidate_row("sqrt", "Square root", "Trackers in use", "transform", "ref_op, sqrt_constant"),
  candidate_row("log1p", "log1p", "Trackers in use", "transform", "fixed_log, leaky_log"),
  candidate_row("boxcox_25", "Box-Cox 0.25 on counts", "Level of the forecast", "transform"),
  candidate_row("boxcox_75", "Box-Cox 0.75 on counts", "Level of the forecast", "transform"),
  candidate_row("power_m100", "(1 + median)^-1", "Level of the forecast", "common scale"),
  candidate_row("power_25", "(1 + median)^0.25", "Level of the forecast", "common scale"),
  candidate_row("power_50", "(1 + median)^0.5", "Level of the forecast", "common scale",
    "studentized_power"),
  candidate_row("power_75", "(1 + median)^0.75", "Level of the forecast", "common scale"),
  candidate_row("power_100", "1 + median", "Level of the forecast", "common scale"),
  candidate_row("power_150", "(1 + median)^1.5", "Level of the forecast", "common scale"),
  candidate_row("power_200", "(1 + median)^2", "Level of the forecast", "common scale"),
  candidate_row("power_level_75", "(1 + b)^0.75 per level", "Level of the forecast",
    "per-level scale"),
  candidate_row("level_rms", "Level RMS of recent errors", "Recent errors", "per-level scale",
    "studentized_level_rms, level_flip"),
  candidate_row("median_rms", "Median rebuild, recent errors", "Recent errors", "common scale",
    "studentized"),
  candidate_row("base_width", "Base width", "Spread of the forecast", "common scale",
    "studentized_full, studentized_base"),
  candidate_row("spread_0", "Spread per level, \u03bb = 0", "Spread of the forecast",
    "per-level scale"),
  candidate_row("spread_20", "Spread per level, \u03bb = 0.2", "Spread of the forecast",
    "per-level scale"),
  candidate_row("spread", "Spread per level, \u03bb = 0.4", "Spread of the forecast",
    "per-level scale", "studentized_spread"),
  candidate_row("spread_60", "Spread per level, \u03bb = 0.6", "Spread of the forecast",
    "per-level scale"),
  candidate_row("spread_80", "Spread per level, \u03bb = 0.8", "Spread of the forecast",
    "per-level scale"),
  candidate_row("half_iqr", "Half the base IQR (\u03bb = 1)", "Spread of the forecast",
    "common scale", "width_flip"),
  candidate_row("pit", "PIT of the base", "Spread of the forecast", "transform",
    "studentized_pit"),
  candidate_row("revisions_cutoff", "Revisions, every target", "Revisions of the forecast",
    "common scale"),
  candidate_row("revisions_target", "Revisions, own target", "Revisions of the forecast",
    "common scale"),
  candidate_row("window_asinh", "asinh, window, exact", "Uniejewski (2025)", "transform"),
  candidate_row("window_asinh_median", "asinh, window, at the median", "Uniejewski (2025)",
    "common scale"),
  candidate_row("window_boxcox_median", "Box-Cox 0, window, at the median",
    "Uniejewski (2025)", "common scale"),
  candidate_row("window_mlog_median", "Mirror log 0.3, window, at the median",
    "Uniejewski (2025)", "common scale"),
  candidate_row("centred_asinh", "asinh, centred, exact", "Uniejewski (2025)", "transform"),
  candidate_row("centred_mlog", "Mirror log 1/3, centred, exact", "Uniejewski (2025)",
    "transform"))

# Transforms T(z) of Uniejewski (2025), their inverse and derivative: the parametrized asinh
# (slope c at the origin), the symmetric Box-Cox and the mirror log.
asinh_transform <- function(c) {
  a <- sqrt(1 / c^2 - 1)
  list(f = function(z) sign(z) * (asinh(abs(z) + a) - asinh(a)),
    inv = function(v) sign(v) * (sinh(abs(v) + asinh(a)) - a),
    d = function(z) 1 / sqrt(1 + (abs(z) + a)^2))
}
boxcox_transform <- function(lambda) list(
  f = function(z) if (lambda == 0) sign(z) * log1p(abs(z)) else
    sign(z) * ((abs(z) + 1)^lambda - 1) / lambda,
  d = function(z) (abs(z) + 1)^(lambda - 1))
mlog_transform <- function(c) list(f = function(z) sign(z) * log1p(c * abs(z)),
  inv = function(v) sign(v) * expm1(abs(v)) / c, d = function(z) c / (1 + c * abs(z)))
# Box-Cox on counts, psi(x) = ((1 + x)^lambda - 1) / lambda (lambda = 0: log1p), and inverse.
boxcox_counts <- function(lambda) list(
  f = function(x, t) if (lambda == 0) log1p(x) else ((1 + x)^lambda - 1) / lambda,
  inv = function(v, t) if (lambda == 0) expm1(v) else pmax(1 + lambda * v, 0)^(1 / lambda) - 1)

# The candidates of one case (score_candidates$id): a list of score objects. s0: the case's
# error_scale() (for the bounds); params: method_parameters().
candidate_scores <- function(case, s0, params = method_parameters(), guard = 20,
                             window = 26L, min_reports = 8L) {
  B <- case$base; a <- case$alpha
  m <- B[, match(.5, a)]
  n <- nrow(B)
  bound <- function(scale) if (!isTRUE(s0$initialized)) scale else
    pmin(pmax(scale, params$sigma_bounds[1] * s0$scale), params$sigma_bounds[2] * s0$scale)
  guarded <- function(sigma) {
    S <- as.matrix(sigma)
    S[!is.finite(S) | S <= 0] <- NA
    centre <- vapply(seq_len(nrow(S)), function(t)
      suppressWarnings(median(S[seq_len(t), ], na.rm = TRUE)), numeric(1))
    S <- pmin(pmax(S, centre / guard), centre * guard)
    if (is.matrix(sigma)) S else S[, 1]
  }
  scaled <- function(sigma, centre = B) list(kind = "scale", sigma = guarded(sigma),
    centre = centre)
  transform <- function(f, inv) list(kind = "transform", f = f, inv = inv)
  # Median and median absolute deviation of the last `window` reports known at each cutoff.
  released <- unlist(case$feedback_by_week)
  n_released <- cumsum(lengths(case$feedback_by_week))
  reports <- lapply(n_released, function(k) if (k < min_reports) NULL else {
    x <- case$feedback_observed[released[max(1L, k - window + 1L):k]]
    x[is.finite(x)]
  })
  m_w <- vapply(reports, function(x) if (length(x)) median(x) else NA_real_, numeric(1))
  s_w <- vapply(reports, function(x) if (length(x)) max(1, median(abs(x - median(x)))) else
    NA_real_, numeric(1))
  s_c <- half_iqr_sigma(case)
  # Standardized transform T((x - centre_t) / spread_t) as a score.
  standardized <- function(tr, centre, spread) transform(
    function(x, t) tr$f((x - centre[t]) / spread[t]),
    function(v, t) centre[t] + spread[t] * tr$inv(v))
  at_median <- function(tr) scaled(s_w / tr$d((m - m_w) / s_w))
  median_scale <- guarded(bound(median_sigma(case, "residuals", params$median_window,
    params$median_min_rows, 1e-6, "identity", estimator = params$residual_estimator)$scale))
  list(
    identity = scaled(rep(1, n)),
    sqrt = transform(function(x, t) sqrt(pmax(x, 0)), function(v, t) pmax(v, 0)^2),
    log1p = transform(boxcox_counts(0)$f, boxcox_counts(0)$inv),
    boxcox_25 = transform(boxcox_counts(.25)$f, boxcox_counts(.25)$inv),
    boxcox_75 = transform(boxcox_counts(.75)$f, boxcox_counts(.75)$inv),
    power_m100 = scaled((1 + m)^-1), power_25 = scaled((1 + m)^.25),
    power_50 = scaled(power_sigma(case, .5)), power_75 = scaled((1 + m)^.75),
    power_100 = scaled(1 + m), power_150 = scaled((1 + m)^1.5), power_200 = scaled((1 + m)^2),
    power_level_75 = scaled((1 + B)^.75),
    level_rms = scaled(bound(level_residual_sigma(case, params$median_window,
      params$median_min_rows, 1e-6, "identity", params$residual_estimator))),
    median_rms = list(kind = "scale", sigma = median_scale,
      centre = m + outer(median_scale, qnorm(a))),
    base_width = scaled(bound(median_sigma(case, "base", floor = 1e-6, transform = "identity",
      levels = params$median_base_levels)$scale)),
    half_iqr = scaled(s_c),
    spread_0 = scaled(spread_sigma(case, 0)), spread_20 = scaled(spread_sigma(case, .2)),
    spread = scaled(spread_sigma(case, params$spread_lambda)),
    spread_60 = scaled(spread_sigma(case, .6)), spread_80 = scaled(spread_sigma(case, .8)),
    pit = list(kind = "pit"),
    revisions_cutoff = scaled(bound(revision_sigma(case, "cutoff")$scale)),
    revisions_target = scaled(bound(revision_sigma(case, "target")$scale)),
    window_asinh = standardized(asinh_transform(1), m_w, s_w),
    window_asinh_median = at_median(asinh_transform(1)),
    window_boxcox_median = at_median(boxcox_transform(0)),
    window_mlog_median = at_median(mlog_transform(.3)),
    centred_asinh = standardized(asinh_transform(1), m, s_c),
    centred_mlog = standardized(mlog_transform(1 / 3), m, s_c))
}
# Scores [rows x level] of outcomes y (one per row of rows) under one candidate, each row
# with the scale or transform of its own cutoff.
candidate_score <- function(candidate, case, y, rows) {
  B <- case$base[rows, , drop = FALSE]
  switch(candidate$kind,
    scale = (y - candidate$centre[rows, , drop = FALSE]) / (if (is.matrix(candidate$sigma))
      candidate$sigma[rows, , drop = FALSE] else candidate$sigma[rows]),
    transform = t(vapply(seq_along(rows), function(i) candidate$f(y[i], rows[i]) -
      candidate$f(B[i, ], rows[i]), numeric(ncol(B)))),
    pit = t(vapply(seq_along(rows), function(i) base_pit(y[i], B[i, ], case$alpha) -
      case$alpha, numeric(ncol(B)))))
}
# The forecast of row t with thresholds u (one per level) on a candidate's score: the
# inversion of S_t^alpha(y) <= u^alpha, levels sorted.
candidate_forecast <- function(candidate, case, u, t) {
  b <- case$base[t, ]
  sort(na.last = TRUE, switch(candidate$kind,
    scale = candidate$centre[t, ] + (if (is.matrix(candidate$sigma)) candidate$sigma[t, ] else
      candidate$sigma[t]) * u,
    transform = candidate$inv(candidate$f(b, t) + u, t),
    pit = base_quantile_function(b, case$alpha, pmin(pmax(case$alpha + u, 1e-4), 1 - 1e-4))))
}

# 5. EXPERT CREATORS ------------------------------------------------------
# All creators return the same expert object:
#   forecasts  [week x expert x level] in counts / count_scale, NA before birth;
#   birth_week the first issue week of each expert;
#   experts    metadata: expert name, multiplier or leak (and start date).
# Rates are multiples of a rolling scale per week (rolling_rate()).

# One MultiQT expert per rate multiplier. All experts issue from the first week, keep the
# base during the warmup and learn from their own issued forecasts. rate_unit: the rolling
# scale of the tracker's scores per week (rolling_rate() with mult 1), so expert j learns at
# multipliers[j] * rate_unit (Delphi's adaptive rate: multiplier .03). sigma (optional): a
# studentized scale in counts (normalized_scale()); rate_unit is then in units of sigma.
make_rate_experts <- function(case, multipliers, rate_unit, count_scale = 1,
                              transform = "log1p", sigma = NULL, warmup_weeks = 0L) {
  n <- nrow(case$base)
  names <- paste0("m=", multipliers)
  forecasts <- array(NA_real_, c(n, length(multipliers), length(case$alpha)),
    dimnames = list(as.character(case$reference_date), names, as.character(case$alpha)))
  for (j in seq_along(multipliers))
    forecasts[, j, ] <- multiqt(case$base, case$feedback_observed, case$feedback_by_week,
      case$alpha, multipliers[j] * rate_unit, warmup_weeks, transform,
      sigma = if (is.null(sigma)) 1 else sigma,
      feedback_samples = case$feedback_samples) / count_scale
  list(forecasts = forecasts, birth_week = rep(1L, length(multipliers)),
    experts = data.frame(expert = names, multiplier = multipliers))
}

# One fresh MultiQT copy per multiplier at every issue week after the warmup. A copy
# starts at the base with zero offsets and learns only from forecasts issued since its
# start, at multiplier * rate_unit of those weeks (make_rate_experts()).
make_restart_experts <- function(case, multipliers, rate_unit, count_scale = 1,
                                 transform = "log1p", sigma = NULL, warmup_weeks = 0L) {
  n <- nrow(case$base)
  births <- if (warmup_weeks < n) seq.int(warmup_weeks + 1L, n) else integer()
  birth_week <- rep(births, each = length(multipliers))
  multiplier <- rep(multipliers, times = length(births))
  names <- if (length(birth_week))
    paste0("start=", case$reference_date[birth_week], ";m=", multiplier) else character()
  forecasts <- array(NA_real_, c(n, length(birth_week), length(case$alpha)),
    dimnames = list(as.character(case$reference_date), names, as.character(case$alpha)))
  for (j in seq_along(birth_week)) {
    issued <- seq.int(birth_week[j], n)
    # Renumber the release schedule from the copy's birth; drop rows before it.
    schedule <- lapply(case$feedback_by_week[issued], function(rows)
      rows[rows >= issued[1]] - issued[1] + 1L)
    copy_sigma <- if (is.null(sigma)) 1 else if (is.matrix(sigma))
      sigma[issued, , drop = FALSE] else sigma[issued]
    forecasts[issued, j, ] <- multiqt(case$base[issued, , drop = FALSE],
      case$feedback_observed[issued], schedule, case$alpha, multiplier[j] * rate_unit[issued],
      transform = transform, sigma = copy_sigma,
      feedback_samples = case$feedback_samples[issued]) / count_scale
  }
  list(forecasts = forecasts, birth_week = birth_week,
    experts = data.frame(expert = names, multiplier,
      start_date = case$reference_date[birth_week]))
}

# Leaky MultiQT at the rate eta (per week) for every leak of `leaks`: the leak bank of
# leaky_select and leaky_weighted, the same expert object as make_rate_experts().
make_leak_experts <- function(case, leaks, eta, count_scale = 1, transform = "log1p") {
  n <- nrow(case$base)
  names <- paste0("leak=", leaks)
  forecasts <- array(NA_real_, c(n, length(leaks), length(case$alpha)),
    dimnames = list(as.character(case$reference_date), names, as.character(case$alpha)))
  for (j in seq_along(leaks))
    forecasts[, j, ] <- multiqt(case$base, case$feedback_observed, case$feedback_by_week,
      case$alpha, eta, 0L, transform, fast_decay = leaks[j],
      feedback_samples = case$feedback_samples) / count_scale
  list(forecasts = forecasts, birth_week = rep(1L, length(leaks)),
    experts = data.frame(expert = names, leak = leaks))
}

# Prequential choice among experts that issue every week (leaky_select): at each issue
# week the expert whose issued forecasts had the smallest quantile loss, summed over the
# levels, over the last `window` rows learned by then (each row's report as learned);
# expert `first` before any row is learned; ties: the first. Returns the issued forecasts
# [week x level] (counts / count_scale) and the chosen expert per week.
select_expert <- function(experts, case, window, count_scale = 1, first = 1L) {
  f <- experts$forecasts
  released <- unlist(case$feedback_by_week)
  n_released <- cumsum(lengths(case$feedback_by_week))
  y <- case$feedback_observed / count_scale
  plausible <- scale_samples(case$feedback_samples, count_scale)   # Corrected labels.
  loss <- matrix(0, length(released), dim(f)[2])
  for (i in seq_along(released)) for (j in seq_len(dim(f)[2])) {
    pinball <- function(y) {
      e <- y - f[released[i], j, ]
      sum(pmax(e, 0) * case$alpha + pmax(-e, 0) * (1 - case$alpha))
    }
    x <- plausible[[released[i]]]
    loss[i, j] <- if (is.null(x)) pinball(y[released[i]]) else mean(vapply(x, pinball, 0))
  }
  chosen <- vapply(n_released, function(k) if (k == 0) as.integer(first) else
    which.min(colSums(loss[max(1L, k - window + 1L):k, , drop = FALSE])), integer(1))
  q <- t(vapply(seq_along(chosen), function(t) f[t, chosen[t], ], numeric(dim(f)[3])))
  dimnames(q) <- dimnames(f)[c(1, 3)]
  list(q = q, chosen = chosen)
}

# Polar MultiQT (polar; notes, "A shared angle for the sign"): the offset of level alpha is
# theta^alpha = m^alpha * cos(phi), a signed length per level and one angle phi in
# [eps, pi - eps] shared by all levels, so that the sign of the whole correction can turn
# at once. Both learn by subgradient steps of the quantile loss through the coverage
# error g = covered - alpha: m^alpha at the tracker's rate eta_t (counts), phi at
# angle_rate (radians per unit of m / unit_t, unit_t the rolling scale of rolling_rate()):
#   m^alpha <- m^alpha - eta_t g^alpha cos(phi),
#   phi <- phi + angle_rate * sum_alpha g^alpha m^alpha / unit_t * sin(phi).
# The issued forecasts are projected (isoreg), the state is not, as in multiqt(); with phi
# fixed at 0 (eps = 0, angle_rate = 0) it is Fixed MultiQT. Returns the forecasts in counts.
polar_multiqt <- function(base, feedback, feedback_by_week, alpha, eta, unit, angle_rate,
                          eps = .1, feedback_samples = NULL) {
  n <- nrow(base); a <- length(alpha)
  stopifnot(length(eta) == n, length(unit) == n, angle_rate >= 0, eps >= 0, eps < pi / 2)
  m <- rep(0, a); phi <- eps; q <- base
  for (week in seq_len(n)) {
    for (past_week in feedback_by_week[[week]]) {
      g <- covered_share(feedback[past_week], feedback_samples[[past_week]],
        q[past_week, ]) - alpha
      push <- sum(g * m) / max(unit[week], 1e-9) * sin(phi)
      m <- m - eta[week] * g * cos(phi)
      phi <- min(max(phi + angle_rate * push, eps), pi - eps)
    }
    if (any(m != 0)) q[week, ] <- isoreg(base[week, ] + m * cos(phi))$yf
  }
  q
}

# Phase-weighted MultiQT (phase_angle, phase_hidden; notes, "One family" and Proposition
# "Phase-weighted coverage"): per level the correction sum_j (P_t V)_j z_j of states z_j,
# issued as b_t + sigma_t * correction (isoreg; the states are not projected); P (n x K): the
# phase probabilities of each row known at its cutoff; V (K x J) maps the states to the
# phases; rate (J): their rate factors. A report on row s learned at issue week t moves z_j by
# eta_t rate_j (W_ts V)_j (alpha - covered), W_ts = weights(t, s) the phase probabilities of
# row s as known then. phase_angle: P = (1 - p, p), V = (shared, flip) = (1 1; 1 -1): a slow
# per-level offset and one correction whose sign follows 1 - 2p. phase_hidden: three free
# phases, V = I.
# diagnostics = TRUE: also the issued components (angle-controls diagnostics, 20 Nov 2026),
# in the units of base: states [week x state x level], the states Z as issued at each week
# (after that week's updates; 0 in the warmup); multiplier [week x state], (P_t V)_j; sigma,
# eta; components [week x state x level] = sigma_t (P_t V)_j Z_j; raw = base + their sum (the
# forecast before the ordering projection); issued = the forecast; adjustment = issued - raw,
# the projection's change, attributed to no state. The forecast is computed exactly as without
# diagnostics. Also (hidden-phase diagnostics, 22 Nov 2026): P, the issuance probabilities
# given; and per learning week t, over the rows s learned that week (F_t): learn_rows = |F_t|,
# learn_issued[t, k] = sum_s P[s, k] (the probabilities the rows were issued with) and
# learn_weights[t, k] = sum_s weights(t, s)[k] (the weights they are learned with).
phase_weighted_multiqt <- function(base, feedback, feedback_by_week, alpha, eta, sigma, P,
                                   weights, V = diag(ncol(P)), rate = rep(1, ncol(V)),
                                   warmup_weeks = 0L, feedback_samples = NULL,
                                   diagnostics = FALSE) {
  n <- nrow(base)
  stopifnot(nrow(P) == n, nrow(V) == ncol(P), length(rate) == ncol(V), length(eta) == n)
  Z <- matrix(0, ncol(V), length(alpha))
  q <- base
  if (diagnostics) {
    states <- array(0, c(n, ncol(V), length(alpha)), list(NULL, colnames(V), alpha))
    multiplier <- P %*% V
    learn_rows <- integer(n)
    learn_issued <- learn_weights <- matrix(0, n, ncol(P), dimnames = list(NULL, colnames(P)))
  }
  for (week in seq_len(n)) {
    if (week <= warmup_weeks) next
    for (s in feedback_by_week[[week]]) {
      w <- weights(week, s)
      step <- eta[week] * (alpha - covered_share(feedback[s], feedback_samples[[s]], q[s, ]))
      Z <- Z + outer(rate * drop(w %*% V), step)
      if (diagnostics) {
        learn_rows[week] <- learn_rows[week] + 1L
        learn_issued[week, ] <- learn_issued[week, ] + P[s, ]
        learn_weights[week, ] <- learn_weights[week, ] + w
      }
    }
    q[week, ] <- isoreg(base[week, ] + sigma[week] * drop((P[week, ] %*% V) %*% Z))$yf
    if (diagnostics) states[week, , ] <- Z
  }
  if (!diagnostics) return(q)
  components <- states * array(multiplier * sigma, dim(states))
  raw <- base + apply(components, c(1, 3), sum)
  list(q = q, diagnostics = list(states = states, multiplier = multiplier, sigma = sigma,
    eta = eta, components = components, raw = raw, issued = q, adjustment = q - raw,
    warmup = seq_len(n) <= warmup_weeks, P = P, learn_rows = learn_rows,
    learn_issued = learn_issued, learn_weights = learn_weights))
}

# Sign-flip experts (power_flip, level_flip, width_flip): the base, a tracker's forecasts
# and its mirror, the base minus the tracker's correction, 2 b_t - q_t, levels sorted. All
# [week x level] in counts / count_scale; experts named base, tracker and mirror.
flip_bank <- function(base, tracker, case) {
  mirror <- t(apply(2 * base - tracker, 1, sort))
  forecasts <- aperm(array(c(base, tracker, mirror), c(nrow(base), ncol(base), 3L)),
    c(1, 3, 2))
  dimnames(forecasts) <- list(as.character(case$reference_date), c("base", "tracker",
    "mirror"), as.character(case$alpha))
  list(forecasts = forecasts, birth_week = rep(1L, 3L),
    experts = data.frame(expert = c("base", "tracker", "mirror")))
}

# 6. EXPONENTIAL WEIGHTING OF EXPERTS --------------------------------------
# Gibbs and Candes (2024), Algorithms 1-2, motivate exponential weighting and sharing.
# https://jmlr.org/papers/v25/22-1218.html
# https://github.com/isgibbs/DtACI/blob/main/DtACI.R
# experts: object from make_rate_experts() or make_restart_experts().
# base, feedback: counts / count_scale; the weights learn from pinball losses in
#   these units divided by loss_scale (weighting_loss_scale(): the base's mean loss
#   over the S0 rows), released by feedback_by_week like the MultiQT updates. The
#   division makes lambda unit-free: the count scale does not change the weights.
# lambda: "learned" starts at lambda_start and, after H updates, uses
#   sqrt((log(2 K H) + 1) / sum of the previous H moments); a number fixes lambda; 0 gives
#   the uniform average of the active experts. The moment is the paper's second moment of
#   the experts' losses under the pre-update weights; centered = TRUE: their weighted
#   variance instead (the score experts: the weights ignore a loss shared by all experts,
#   e.g. the base error that recalibrations of one base share, so the variance sets lambda
#   for the differences).
# sharing: after each update, share (default 1/(2H)) of the mass is spread uniformly.
# per_quantile: one weight vector per level. FALSE: one weight vector for all levels,
#   learned from the mean pinball loss over levels, with scalar lambda and moment.
#   A common-weight average of ordered quantiles stays ordered, so it is not repaired.
# Later starts: an expert joins at its birth week with 1/K of the weight among the K
#   active experts; the others keep their relative weights. A released loss of row t
#   updates only experts that issued row t; the others keep their mass.
# keep: only the keep youngest experts stay active (Inf: all). When a younger one joins, the
#   oldest retires for good; its mass goes to the others in proportion to their weights.
# The first warmup_weeks forecasts keep the base and their reports are not replayed.
# Returns the mixture q [week x level] in counts / count_scale, the issued weights
# [week x expert x level, or x 1 if common], the lambda used at each issue week and the
# experts' forecasts it weighs (experts$forecasts).
exp_weights <- function(experts, base, feedback, feedback_by_week, alpha,
                        lambda = "learned", centered = FALSE, sharing = TRUE, H = 8L,
                        share = 1 / (2 * H), per_quantile = TRUE, lambda_start = 2.72,
                        warmup_weeks = 0L, loss_scale = 1, keep = Inf,
                        feedback_samples = NULL) {
  forecasts <- experts$forecasts
  birth <- experts$birth_week
  n <- nrow(base); k <- length(birth); a <- length(alpha)
  learned <- identical(lambda, "learned")
  stopifnot(identical(dim(forecasts), c(n, k, a)), length(feedback) == n,
    length(feedback_by_week) == n, learned || (is.numeric(lambda) &&
      length(lambda) == 1L && is.finite(lambda) && lambda >= 0),
    !sharing || (is.finite(share) && share > 0 && share <= 1),
    length(loss_scale) == 1L, is.finite(loss_scale), loss_scale > 0,
    length(keep) == 1L, keep >= 1)
  learn <- learned || lambda > 0
  columns <- if (per_quantile) a else 1L     # One weight vector per level, or one for all.
  learning_rate <- rep(if (learned) lambda_start else lambda, columns)
  p <- matrix(0, k, columns)                 # Weights; zero before an expert's birth.
  log_weights <- matrix(-Inf, k, columns)    # Kept across updates, even if exp() underflows.
  moments <- matrix(NA_real_, length(unlist(feedback_by_week)), columns)
  n_updates <- 0L
  weights <- array(0, c(n, k, columns), dimnames = list(dimnames(forecasts)[[1]],
    dimnames(forecasts)[[2]], if (per_quantile) dimnames(forecasts)[[3]] else "all levels"))
  lambda_used <- matrix(NA_real_, n, columns)
  log_mass <- function(x) apply(x, 2, function(v) max(v) + log(sum(exp(v - max(v)))))
  retired <- rep(FALSE, k)
  q <- base
  for (week in seq_len(n)) {
    admitted <- which(birth < week & !retired)
    if (learn && week > warmup_weeks) for (past in feedback_by_week[[week]]) {
      eligible <- admitted[birth[admitted] <= past]
      if (length(eligible) < 2L) next          # A single expert keeps all of its mass.
      sleeping <- setdiff(admitted, eligible)
      if (learned && n_updates >= H) {
        # Use the PREVIOUS H loss moments, excluding the feedback being processed now.
        recent <- seq.int(n_updates - H + 1L, n_updates)
        denominator <- colSums(moments[recent, , drop = FALSE])
        positive <- denominator > 0
        constant <- log(2 * length(eligible) * H) + 1
        learning_rate[positive] <- sqrt(constant / denominator[positive])
        # Tiny positive moments can overflow the ratio before its square root.
        overflow <- positive & !is.finite(learning_rate)
        learning_rate[overflow] <- sqrt(constant) / sqrt(denominator[overflow])
        # If recent losses are all zero, keep the previous finite learning rate.
      }
      issued <- matrix(forecasts[past, eligible, ], length(eligible), a)
      pinball <- function(y) {
        error <- y - issued
        (pmax(error, 0) * rep(alpha, each = length(eligible)) +
          pmax(-error, 0) * rep(1 - alpha, each = length(eligible))) / loss_scale
      }
      # Corrected labels: the mean quantile loss over the plausible counts.
      plausible <- feedback_samples[[past]]
      loss <- if (is.null(plausible)) pinball(feedback[past]) else
        Reduce(`+`, lapply(plausible, pinball)) / length(plausible)
      if (!per_quantile) loss <- matrix(rowMeans(loss), ncol = 1L)
      stopifnot(all(is.finite(loss)))
      # Experts without a forecast for this row sleep and keep their mass.
      mass <- 1 - colSums(p[sleeping, , drop = FALSE])
      # Second moment of the expert losses under the pre-update weights (the paper's), or
      # their weighted variance (centered). For the variance, losses that differ only by
      # rounding (relative spread below 1e-6) carry no information: no moment, as for equal
      # losses, instead of a lambda that turns rounding into weight changes.
      n_updates <- n_updates + 1L
      prior <- sweep(p[eligible, , drop = FALSE], 2, mass, `/`)
      second <- colSums(prior * loss^2)
      moments[n_updates, ] <- if (!centered) second else {
        variance <- colSums(prior * sweep(loss, 2, colSums(prior * loss))^2)
        ifelse(variance > 1e-12 * second, variance, 0)
      }
      updated <- log_weights[eligible, , drop = FALSE] - sweep(loss, 2, learning_rate, `*`)
      if (length(sleeping)) {
        # Shift the eligible experts' log weights so that their total mass is unchanged.
        updated <- sweep(updated, 2, log_mass(log_weights[eligible, , drop = FALSE]) -
          log_mass(updated), `+`)
      }
      log_weights[eligible, ] <- updated
      # Exponential weights, stabilized and normalized separately for each column.
      log_weights[admitted, ] <- sweep(log_weights[admitted, , drop = FALSE], 2,
        apply(log_weights[admitted, , drop = FALSE], 2, max), `-`)
      p[admitted, ] <- exp(log_weights[admitted, , drop = FALSE])
      p[admitted, ] <- sweep(p[admitted, , drop = FALSE], 2,
        colSums(p[admitted, , drop = FALSE]), `/`)
      if (sharing) {
        p[eligible, ] <- (1 - share) * p[eligible, , drop = FALSE] +
          rep(share * mass / length(eligible), each = length(eligible))
        log_weights[admitted, ] <- log(p[admitted, , drop = FALSE])
      }
    }
    newborn <- which(birth == week)
    active <- c(admitted, newborn)
    if (length(active) > keep) {
      # The oldest retire; the incumbents that stay share their mass (uniformly if none had any).
      retiring <- head(active[order(birth[active], active)], length(active) - keep)
      retired[retiring] <- TRUE
      p[retiring, ] <- 0
      log_weights[retiring, ] <- -Inf
      active <- setdiff(active, retiring)
      incumbents <- setdiff(active, newborn)
      if (length(incumbents)) {
        total <- colSums(p[incumbents, , drop = FALSE])
        p[incumbents, ] <- sweep(p[incumbents, , drop = FALSE], 2, total, `/`)
        p[incumbents, total == 0] <- 1 / length(incumbents)
        log_weights[incumbents, ] <- log(p[incumbents, , drop = FALSE])
      }
    }
    if (!learn) {
      p[active, ] <- 1 / length(active)
    } else if (length(newborn)) {
      # Each newborn gets 1/K among K active experts, before issuing its first forecast.
      # Preserve the incumbents' relative weights after processing released losses.
      new_mass <- length(newborn) / length(active)
      p <- p * (1 - new_mass)
      p[newborn, ] <- new_mass / length(newborn)
      log_weights[active, ] <- log(p[active, , drop = FALSE])
    }
    if (length(active) && week > warmup_weeks) {
      stopifnot(all(is.finite(p)), max(abs(colSums(p) - 1)) < 1e-10)
      current <- matrix(forecasts[week, active, ], length(active), a)
      # Weight first, then repair crossings of per-level mixtures. Keep an all-base mixture exact.
      if (all(current == rep(base[week, ], each = length(active)))) {
        q[week, ] <- base[week, ]
      } else {
        mixture <- colSums(current * p[active, rep(seq_len(columns), length.out = a), drop = FALSE])
        q[week, ] <- if (per_quantile) isoreg(mixture)$yf else mixture
      }
    }
    weights[week, , ] <- p
    lambda_used[week, ] <- learning_rate
  }
  list(q = q, weights = weights, lambda = lambda_used, experts = forecasts)
}

# 7. THE RECALIBRATION METHODS OF SCRIPT 04 --------------------------------------
# Settings of the methods, none tuned (no tuning since 13 November 2026). Every tracker
# learns at Delphi's adaptive rate on its own scores, adaptive_multiplier times the 90%
# quantile of the absolute scores of the last adaptive_rounds learned rows, the first rows
# replayed (rolling_rate()); the rate grid and the weekly starts use multiples of the same
# rolling scale. The leak is a fraction of the offsets per week. by_method overrides a
# setting for one method, e.g. list(shared = list(H = 4L)).
method_parameters <- function() list(
  adaptive_multiplier = .03, adaptive_rounds = 20L,
  # Delphi's REF-op as in covid prod, cold (no burn-in season), with Delphi's constants,
  # not tuned (pipelines/covid_hosp_prod.R, cal_ref_args("op") in R/calibration/views.R of
  # cmu-delphi/delphi-hub-forecasting): MultiQT on sqrt(counts) for either tracker, with
  # the fast rate ref_op_multiplier * Q90 of the last ref_op_window learned rows (at least
  # ref_op_floor) and a slow term with ref_op_slow_multiplier * Q90 of all learned rows;
  # the fast offsets shrink by ref_op_fast_decay every week, the slow ones never.
  ref_op_transform = "sqrt", ref_op_multiplier = .03, ref_op_window = 20L,
  ref_op_floor = 1e-3, ref_op_slow_multiplier = .003, ref_op_fast_decay = .1,
  # Delphi's flu prod since 2026-10-07 (pipelines/flu_hosp_prod.R, "sqrt constant 0.018 per
  # 100k", cold), not tuned: Fixed MultiQT on sqrt(counts per 100k residents) for either
  # tracker, with this constant rate in those units (Delphi chose it on flu 2024/25-2025/26).
  sqrt_constant_rate = .018,
  # Learned mixtures (exp_weights()): adaptive lambda from lambda_start, look-back H (fixed
  # for now, not tuned); sharing, where used, 1/(2H).
  H = 8L, lambda_start = 2.72,
  outer_multiplier = .001,               # Outer MultiQT on a mixture: eta = m * S0.
  # Weekly restarts: a fresh copy at every issue week learning at restart_multiplier times
  # the rolling scale (the adaptive Fixed tracker, started anew every week);
  # restart_recent: learned weights (H, sharing 1/(2H)) on the restart_recent_k most recent
  # copies. restart_recent_dual: two copies per start, a slow and a fast one
  # (restart_dual_multipliers), with learned weights and sharing on the copies of the
  # restart_recent_k most recent starts. None tuned.
  restart_multiplier = .03, restart_recent_k = 8L, restart_dual_multipliers = c(.03, .3),
  # Leaky MultiQT (Delphi's leak): Fixed MultiQT whose offsets shrink by the fraction leak
  # at every issue week (Delphi's fast_decay, .1 as in REF-op). leak_grid: the leak bank of
  # leaky_select (the leak with the smallest quantile loss over the last leak_window learned
  # rows) and leaky_weighted (learned weights).
  # Polar MultiQT (polar_multiqt()): the rate of the shared angle (unit-free) and its
  # distance eps from 0 and 180 degrees; not tuned.
  polar_angle_rate = .3, polar_eps = .1,
  leak = .1, leak_grid = c(0, .02, .05, .1, .2, .3, .5), leak_window = 20L,
  # phase_power: the decline probability from the last phase_rows learned rows
  # (decline_probability()).
  phase_rows = 4L,
  # phase_angle: the rate factor of its slow per-level offset (as REF-op's slow term, .003
  # against .03); not tuned.
  phase_slow = .1,
  # The rate multipliers of phase_angle's one-offset controls (angle_controls); fixed before
  # scoring (20 Nov 2026), a diagnostic grid, not tuned.
  angle_control_rates = c(angle_slow_only = .1, angle_rate_025 = .25, angle_rate_05 = .5,
    angle_rate_1 = 1),
  # Score-based experts (score_experts, score_methods): learned weights per level with H and
  # sharing 1/(2H) as the other learned mixtures; not tuned.
  # Studentized MultiQT: sigma_t from the RMS of the last median_window released
  # residuals (the base until median_min_rows exist) or from the base width between
  # median_base_levels; median_start: "normal" or "t" start of the residual median methods
  # (their *_zero variants start all quantiles at the median).
  median_window = 5L, median_min_rows = 5L, median_base_levels = c(.1, .9),
  median_start = "normal",
  # Residual scales (studentized median, constant sigma, level RMS): "rms" or "median"
  # (robust to one outlying residual; see residual_scale()).
  residual_estimator = "rms",
  # Every studentized scale is clipped to these multiples of S0 (tracker units). Near 0
  # (e.g. a base of zero width) the offsets grow without effect for weeks; a jump (e.g.
  # a broken base forecast or report) then multiplies all of them. c(0, Inf): no limits.
  sigma_bounds = c(.05, 1),
  # Scales of the studentized methods with the most stationary scores (scaled_methods):
  # power_sigma() exponent, spread_sigma() lambda; each normalized to S0 and guarded
  # (normalized_scale()).
  power_p = .5, spread_lambda = .4, scale_guard = 20,
  by_method = list())
method_parameter <- function(params, method, name) {
  value <- params$by_method[[method]][[name]]
  # Variants use the values set for their method (parameter_source).
  if (is.null(value) && method %in% names(parameter_source))
    value <- params$by_method[[parameter_source[[method]]]][[name]]
  if (is.null(value)) params[[name]] else value
}
# Score-based experts: one bank of the forecasts of these five methods (score_experts: the
# base, the additive trackers, the best studentized method and Delphi's square-root tracker;
# the base-width studentized methods left out), weighted by equal weights (score_equal) or
# by learned weights per level (score_learned). The members are the same in every week and
# re-run; members that issue identical forecasts in some weeks stay, so later forecasts never
# change earlier weights.
score_experts <- c("original", "fixed", "leaky", "studentized_level_rms", "sqrt_constant")
score_methods <- c("score_equal", "score_learned")
# Studentized MultiQT with two scales of the most stationary scores (power_sigma(),
# spread_sigma()): per scale the tracker at the adaptive rate, and the rate grid and the
# most recent weekly starts, each with equal and with learned weights.
scaled_methods <- unlist(lapply(c("power", "spread"), function(scale)
  c(paste0("studentized_", scale), paste0(scale, c("_uniform", "_weighted_common",
    "_restart_equal", "_restart_recent")))))
# Methods that use the settings of another (method_parameter()): the mixtures of the scaled
# trackers those of their unscaled counterpart, the outer MultiQT H of the mixture it
# corrects.
parameter_source <- c(studentized_zero = "studentized", median_constant_zero = "median_constant",
  power_uniform = "uniform", power_weighted_common = "weighted_common",
  power_restart_equal = "restart_recent", power_restart_recent = "restart_recent",
  spread_uniform = "uniform", spread_weighted_common = "weighted_common",
  spread_restart_equal = "restart_recent", spread_restart_recent = "restart_recent",
  shared_outer = "shared", shared_common_outer = "shared_common",
  restart_uniform_outer = "restart_uniform", restart_weighted_outer = "restart_weighted")
# Sign-flip experts and the two-phase trackers (notes, "Sign-flip experts" and "A phase
# probability in the score").
flip_methods <- c("power_flip", "level_flip", "width_flip")
phase_methods <- c("phase_power", "phase_oracle")
# Bias and phase probabilities (R/05_phase_probabilities.R; the models of 23 November 2026).
outside_methods <- c("phase_angle", "phase_hidden")
# Controls of phase_angle (20 November 2026; diagnostic, not in method_groups, so neither in
# 06's tree nor in 04's trade-off pages; their page: diagnostic_groups): one power-scale
# offset at a multiple of phase_angle's rate (angle_control_rates; .1 = p fixed at 1/2: the
# slow offset alone), the angle with |1 - 2p| in forecast and update (no reversals), and one
# offset learning each released error at {.1 + (1 - 2 p_t)^2} times the rate of the learning
# week t, issued without a probability multiplier. Same power scale, adaptive rate, warmup and
# phase inputs as phase_angle; each learns from its own issued forecasts.
angle_controls <- c("angle_slow_only", "angle_rate_025", "angle_rate_05", "angle_rate_1",
  "angle_abs", "angle_rate_only")
# Controls of phase_hidden (22 November 2026; diagnostic like angle_controls): equal
# probabilities (1/3, 1/3, 1/3) at issuance and learning (no phase information; from zero
# states it is one tracker at a third of the rate), and phase_hidden learning each report with
# the probabilities P_s its row was issued with instead of W_{s,t}, those of its target week
# as known at the learning cutoff (no revised assignments; also in re-runs). Same base, power
# scale, adaptive rate, warmup, projection and feedback as phase_hidden; each learns from its
# own issued forecasts. The phases of the hidden-phase model, in the order of R/05's columns
# rise, over, tail.
hidden_controls <- c("hidden_uniform", "hidden_issued_weights")
hidden_phases <- c("rise", "overshoot", "tail")
# The true three phases (23 November 2026; a hindsight benchmark like phase_oracle, drawn as a
# reference and never ranked in 06): phase_hidden issued and learned with the phase of each
# row's target week relative to the nearest peak of the latest counts instead of the
# hidden-phase model's probabilities. Needs no phase inputs (every location).
hindsight_hidden <- "hidden_oracle"
method_names <- c("original", "fixed", "leaky", "fixed_log", "leaky_log", "ref_op",
  "sqrt_constant", "sqrt_adaptive", "leaky_select", "leaky_weighted", "studentized",
  "studentized_base", "median_constant",
  "studentized_zero", "median_constant_zero", "studentized_full", "studentized_level_rms",
  "studentized_pit", scaled_methods, "oracle_target", "uniform",
  "weighted", "weighted_common", "shared",
  "shared_common", "restart_uniform", "restart_weighted", "restart_recent", "restart_recent_dual",
  score_methods, "restart_persistent", "shared_outer", "shared_common_outer",
  "restart_uniform_outer", "restart_weighted_outer", flip_methods, phase_methods, "polar",
  outside_methods, hindsight_hidden, angle_controls, hidden_controls)
# The mixture each outer MultiQT corrects.
outer_inner <- c(shared_outer = "shared", shared_common_outer = "shared_common",
  restart_uniform_outer = "restart_uniform", restart_weighted_outer = "restart_weighted")
# Methods kept in the code but not run (13 November 2026, from the run of 12 November; notes,
# "Methods not run"): apply_methods() leaves them out unless asked for, and so do the tuning
# (R/05), the plots (method_groups) and 06's tree. Delete a name to run that method again.
paused_methods <- c(
  "oracle_target", "sqrt_adaptive",                                      # Single trackers.
  "studentized", "studentized_zero", "studentized_base", "median_constant",
  "median_constant_zero", "studentized_full",                            # Studentized.
  "power_uniform", "power_restart_equal", "spread_uniform", "spread_restart_equal",
  "shared", "shared_common", "shared_outer", "shared_common_outer",      # Rate grid.
  "restart_persistent", "restart_weighted", "restart_uniform_outer",
  "restart_weighted_outer",                                              # Weekly starts.
  "score_equal",                                                         # Score bank.
  # 17 November 2026: the rate-grid mixtures and the PIT tracker.
  "uniform", "weighted", "weighted_common", "power_weighted_common", "spread_weighted_common",
  "studentized_pit",
  # 19 November 2026: the score bank and the slow + fast weekly starts (close to the base).
  "score_learned", "restart_recent_dual",
  # 20 November 2026: the level RMS of recent errors, its tracker and its sign-flip experts.
  "studentized_level_rms", "level_flip")
active_methods <- setdiff(method_names, paused_methods)
# The method families of 04's plots (forecast pages, trade-off colours); their methods are
# also the candidates of 06's tree.
method_groups <- list(
  "Single trackers" = c("original", "fixed", "leaky", "ref_op", "sqrt_constant",
    "oracle_target", "sqrt_adaptive"),
  "Log1p tracker" = c("original", "fixed", "leaky", "fixed_log", "leaky_log"),
  "Leak chosen online" = c("original", "leaky", "leaky_select", "leaky_weighted"),
  "Rate-grid experts" = c("original", "fixed", "uniform", "weighted", "weighted_common",
    "shared", "shared_common", "shared_outer", "shared_common_outer"),
  "Weekly starts: immediate admission" = c("original", "restart_persistent", "restart_uniform",
    "restart_weighted", "restart_recent", "restart_recent_dual", "restart_uniform_outer",
    "restart_weighted_outer"),
  "Score-based experts" = c("original", score_methods),
  "Studentized median" = c("original", "fixed", "studentized", "studentized_zero",
    "studentized_base", "median_constant", "median_constant_zero"),
  "Studentized full base" = c("original", "fixed", "studentized_full",
    "studentized_level_rms", "studentized_pit"),
  "Power scale" = c("original", scaled_methods[startsWith(scaled_methods, "studentized_power") |
    startsWith(scaled_methods, "power_")]),
  "Spread scale" = c("original", scaled_methods[startsWith(scaled_methods, "studentized_spread") |
    startsWith(scaled_methods, "spread_")]),
  "Sign-flip experts" = c("original", "studentized_power", "studentized_level_rms",
    flip_methods, "polar"),
  "Two-phase trackers" = c("original", "studentized_power", phase_methods),
  "Bias and phase probabilities" = c("original", "phase_oracle", outside_methods,
    hindsight_hidden))
# Only the methods that run: a group keeps its active methods and is dropped when only the
# base and Fixed MultiQT are left.
method_groups <- Filter(function(group) length(setdiff(group, c("original", "fixed"))) > 0L,
  lapply(method_groups, intersect, active_methods))
# Diagnostic pages of 04's forecast plots (not candidates of 06's tree, not on the trade-off
# pages), with their own colours and line types: the rate multipliers dashed in blues from
# light (slow) to dark (1x), no reversals solid vermillion, rate only dot-dashed green.
diagnostic_groups <- list("Angle + slow offset: controls" = intersect(c("original",
  "phase_angle", angle_controls), active_methods),
  "Hidden-phase model: controls" = intersect(c("original", "phase_hidden", hidden_controls),
    active_methods))
# Hidden-phase controls: equal probabilities dashed turquoise, issued weights dot-dashed
# magenta (both far from every other method's colour).
diagnostic_colours <- c(phase_angle = "#E69F00", angle_slow_only = "#9ECAE1",
  angle_rate_025 = "#4292C6", angle_rate_05 = "#2171B5", angle_rate_1 = "#08306B",
  angle_abs = "#D55E00", angle_rate_only = "#009E73", phase_hidden = "#7815C4",
  hidden_uniform = "#00CED1", hidden_issued_weights = "#E7298A")
diagnostic_line_types <- c(phase_angle = 1, angle_slow_only = 2, angle_rate_025 = 2,
  angle_rate_05 = 2, angle_rate_1 = 2, angle_abs = 1, angle_rate_only = 4, phase_hidden = 1,
  hidden_uniform = 2, hidden_issued_weights = 4)
# Short method names without settings (06's tree; 04's labels add the tuned values).
method_short_labels <- c(original = "Base", fixed = "Fixed MultiQT", leaky = "Leaky MultiQT",
  fixed_log = "Fixed MultiQT (log1p tracker)", leaky_log = "Leaky MultiQT (log1p tracker)",
  ref_op = "Delphi REF-op (covid prod)", sqrt_constant = "Delphi sqrt 0.018 (flu prod)",
  oracle_target = "Target oracle (hindsight)",
  sqrt_adaptive = "Sqrt tracker", leaky_select = "Leaky MultiQT, leak chosen online",
  leaky_weighted = "Leaky MultiQT, leak bank with learned weights",
  polar = "Polar MultiQT (shared angle)",
  power_flip = "Sign-flip experts: power scale", level_flip = "Sign-flip experts: level RMS",
  width_flip = "Sign-flip experts: half the base IQR",
  phase_power = "Two-phase power scale (decline probability)",
  phase_angle = "Angle + slow offset (bias probability)",
  phase_hidden = "Hidden-phase model (three phases, outside data)",
  phase_oracle = "Two-phase power scale (phase in hindsight)",
  hidden_oracle = "Three phases (phase in hindsight)",
  angle_slow_only = "Power MultiQT (0.1x rate; slow only)",
  angle_rate_025 = "Power MultiQT (0.25x rate)", angle_rate_05 = "Power MultiQT (0.5x rate)",
  angle_rate_1 = "Power MultiQT (1x rate)",
  angle_abs = "Angle magnitude + slow offset (no reversals)",
  angle_rate_only = "Power MultiQT (probability-dependent rate only)",
  hidden_uniform = "Hidden phases, equal probabilities (one tracker at rate/3)",
  hidden_issued_weights = "Hidden-phase model, learned with issued probabilities",
  uniform = "Rate grid: equal", weighted = "Rate grid: learned, no sharing",
  weighted_common = "Rate grid: learned, common weights",
  shared = "Rate grid: learned + sharing",
  shared_common = "Rate grid: learned + sharing, common weights",
  shared_outer = "Rate grid: sharing + outer",
  shared_common_outer = "Rate grid: sharing, common weights + outer",
  restart_persistent = "First start only", restart_uniform = "Weekly starts: equal",
  restart_weighted = "Weekly starts: learned + sharing",
  restart_uniform_outer = "Starts equal + outer",
  restart_weighted_outer = "Starts learned + sharing + outer",
  restart_recent = "Weekly starts: most recent, learned + sharing",
  restart_recent_dual = "Weekly starts: most recent, slow + fast, learned + sharing",
  score_equal = "Score experts: equal weights",
  score_learned = "Score experts: learned weights",
  studentized = "Studentized median (residual RMS)",
  studentized_base = "Studentized median (base width)",
  median_constant = "Median, constant sigma",
  studentized_zero = "Studentized median (residual RMS), median start",
  median_constant_zero = "Median, constant sigma, median start",
  studentized_full = "Studentized full base (base width)",
  studentized_level_rms = "Studentized full base (level RMS)",
  studentized_pit = "Studentized full base (PIT)")
scale_titles <- c(power = "Power scale", spread = "Spread scale")
method_short_labels <- c(method_short_labels, unlist(lapply(names(scale_titles), function(scale)
  setNames(c(paste0("Studentized full base (", sub("^(Power|Spread)", "\\L\\1",
    scale_titles[[scale]], perl = TRUE), ")"),
    paste0(scale_titles[[scale]], ": ", c("rate grid, equal", "rate grid, learned, common weights",
      "most recent starts, equal", "most recent starts, learned + sharing"))),
    c(paste0("studentized_", scale), paste0(scale, c("_uniform", "_weighted_common",
      "_restart_equal", "_restart_recent")))))))
# Colours of 04's trade-off plots and 06's tree: a method takes the colour family of the
# first group of method_groups that lists it, one shade per method; common-weight and outer
# variants share their method's shade. Every method group needs a palette (an hcl.colors()
# palette or one colour per method): methods of a group without one get no colour and are
# not drawn. The bias and phase probability methods in amber and violet, hues no other method
# has: in shades of one palette they could not be told apart from each other and from the
# phase oracle drawn with them.
group_palettes <- list("Single trackers" = "Blues 3", "Log1p tracker" = "Earth",
  "Rate-grid experts" = "Oranges",
  "Weekly starts: immediate admission" = "Greens 3", "Score-based experts" = "Purples 3",
  "Studentized median" = "Burg",
  "Studentized full base" = "BluGrn", "Power scale" = "Peach",
  "Spread scale" = "Purp", "Leak chosen online" = "Reds 3", "Sign-flip experts" = "Magenta",
  "Two-phase trackers" = "Dark Mint",
  "Bias and phase probabilities" = c(phase_angle = "#E69F00", phase_hidden = "#7815C4",
    hidden_oracle = "#C9A0F0"))
home_group <- local({
  home <- character()
  for (group in names(method_groups)) for (method in method_groups[[group]])
    if (!method %in% c("original", names(home))) home[method] <- group
  home
})
tradeoff_colours <- local({
  base_method <- sub("_outer$", "", sub("_common", "", names(home_group)))
  colours <- c(original = "black")
  for (group in names(group_palettes)) {
    members <- home_group == group
    bases <- unique(base_method[members])
    # Sequential palettes run dark to light; drop the two lightest shades.
    palette <- group_palettes[[group]]
    shades <- if (is.null(names(palette)))
      hcl.colors(length(bases) + 2L, palette)[seq_along(bases)] else palette[bases]
    colours[names(home_group)[members]] <- shades[match(base_method[members], bases)]
  }
  colours
})

# Issue the forecasts of the given methods (default: those that run, active_methods) for one
# case, with every location fitted independently: offsets, expert weights, lambda and loss
# history start fresh. Returns
# the predictions in counts / count_scale (original first) and what 04 plots and saves.
# cache (R/05_tuning.R): an environment for one case and fixed settings apart from the
# transform. Scales, expert banks and mixtures computed for one grid value are reused for
# the next; each key holds every setting the stored object depends on. NULL: no reuse.
# Re-run design (get_case(revisions = TRUE)): rerun_methods() below; weeks: the issue
# weeks it re-runs (default: all), e.g. only the scored ones.
apply_methods <- function(case, settings, params = method_parameters(), methods = active_methods,
                          cache = NULL, weeks = NULL) {
  if (isTRUE(case$feedback_revised))
    return(rerun_methods(case, settings, params, methods, cache, weeks))
  par <- function(method, name) method_parameter(params, method, name)
  methods <- union(methods, outer_inner[intersect(names(outer_inner), methods)])
  want <- function(...) any(c(...) %in% methods)
  scored <- want(score_methods)                    # Needs the six score_experts.
  transform <- settings$transform
  memo <- function(key, value) {
    if (is.null(cache)) return(value)
    key <- paste(transform, key)
    if (is.null(cache[[key]])) cache[[key]] <- value
    cache[[key]]
  }
  scaling <- choose_count_scale(case, settings$scale_method)
  count_scale <- scaling$scale
  s0 <- if (!is.null(case$frozen_s0)) case$frozen_s0 else memo("s0", error_scale(case,
    settings$warmup_rounds, settings$scale_floor, transform, settings$replay_warmup))
  # Keep stored predictions and the scoring/weighting inputs in counts / count_scale.
  base <- case$base / count_scale
  feedback <- case$feedback_observed / count_scale
  adjust <- function(counts, eta, warmup_weeks) multiqt(counts, case$feedback_observed,
    case$feedback_by_week, case$alpha, eta, warmup_weeks, transform,
    feedback_samples = case$feedback_samples) / count_scale

  # Every tracker learns at Delphi's adaptive rate (REF-op's fast term), untuned:
  # adaptive_multiplier times the 90% quantile of its own absolute scores over the last
  # adaptive_rounds learned rows, the first rows replayed in hindsight (rolling_rate(),
  # tracker_score()). The rate grid and the weekly starts use multiples of the same rolling
  # scale (unit()); there is no tuning.
  unit <- function(score) rolling_rate(case, score, 1, params$adaptive_rounds,
    settings$replay_warmup)
  rate_of <- function(unit) params$adaptive_multiplier * unit
  count_unit <- memo("count unit", unit(tracker_score(case, transform)))
  predictions <- list(original = base)
  if (want("fixed", "oracle_target") || scored)
    predictions$fixed <- adjust(case$base, rate_of(count_unit), 0L)
  # Leaky MultiQT (Delphi's leak; REF-op's fast term without its slow term): Fixed MultiQT
  # whose offsets shrink by the fraction leak at every issue week before its update
  # (multiqt(): fast_decay), so a correction that new reports do not refresh fades back to
  # the base: a ridge penalty on the offsets (notes, "The leak as a ridge penalty").
  if (want("leaky") || scored)
    predictions$leaky <- multiqt(case$base, case$feedback_observed, case$feedback_by_week,
      case$alpha, rate_of(count_unit), 0L, transform, fast_decay = params$leak,
      feedback_samples = case$feedback_samples) / count_scale
  # Polar MultiQT (polar_multiqt()): Fixed MultiQT's lengths at its adaptive rate and one
  # angle for all levels that can turn the sign of the whole correction.
  if (want("polar"))
    predictions$polar <- polar_multiqt(case$base, case$feedback_observed,
      case$feedback_by_week, case$alpha, rate_of(count_unit), count_unit,
      params$polar_angle_rate, params$polar_eps, case$feedback_samples) / count_scale
  # Fixed and Leaky MultiQT on the log1p tracker, whatever the tracker of the other methods,
  # at the adaptive rate of their log1p scores. The log1p tracker issues (1 + b) * exp(offset)
  # - 1: its correction is multiplicative.
  if (want("fixed_log", "leaky_log")) {
    log_rate <- rate_of(unit(tracker_score(case, "log1p")))
    log_tracker <- function(...) multiqt(case$base, case$feedback_observed,
      case$feedback_by_week, case$alpha, log_rate, 0L, "log1p", ...,
      feedback_samples = case$feedback_samples) / count_scale
    if (want("fixed_log")) predictions$fixed_log <- log_tracker()
    if (want("leaky_log")) predictions$leaky_log <- log_tracker(fast_decay = params$leak)
  }
  # Delphi's REF-op (method_parameters()): Fixed MultiQT's update on its own sqrt tracker,
  # with Delphi's rolling rates from the rows learned so far (adaptive_rate(), no replay;
  # re-run: as known at the issue week, as Delphi's prod re-runs), the slow term and the
  # leak; Delphi's constants.
  if (want("ref_op")) {
    tracker <- par("ref_op", "ref_op_transform")
    rate <- function(multiplier, rounds, floor) adaptive_rate(case,
      par("ref_op", multiplier), rounds, floor, tracker)
    predictions$ref_op <- multiqt(case$base, case$feedback_observed, case$feedback_by_week,
      case$alpha, rate("ref_op_multiplier", par("ref_op", "ref_op_window"),
        par("ref_op", "ref_op_floor")), 0L, tracker,
      eta_slow = rate("ref_op_slow_multiplier", Inf, 0),
      fast_decay = par("ref_op", "ref_op_fast_decay"),
      feedback_samples = case$feedback_samples) / count_scale
  }
  # Fixed MultiQT on the sqrt tracker at the adaptive rate (REF-op without its slow term and
  # leak).
  if (want("sqrt_adaptive")) predictions$sqrt_adaptive <- multiqt(case$base,
    case$feedback_observed, case$feedback_by_week, case$alpha,
    rate_of(unit(tracker_score(case, "sqrt"))), 0L, "sqrt",
    feedback_samples = case$feedback_samples) / count_scale
  # Delphi's flu prod (method_parameters()): Fixed MultiQT on sqrt(counts per 100k) with a
  # constant rate; no S0 and no warmup. Unchanged forecasts keep the exact base.
  if (want("sqrt_constant") || scored) {
    per_100k <- choose_count_scale(case, "population")$scale
    tracked <- multiqt(case$base / per_100k, case$feedback_observed / per_100k,
      case$feedback_by_week, case$alpha, par("sqrt_constant", "sqrt_constant_rate"), 0L, "sqrt",
      feedback_samples = scale_samples(case$feedback_samples, per_100k))
    q <- tracked * per_100k
    unchanged <- tracked == case$base / per_100k
    q[unchanged] <- case$base[unchanged]
    predictions$sqrt_constant <- q / count_scale
  }
  # Studentized median MultiQT (median_multiqt()): Fixed MultiQT's coverage indicator
  # and update on the base median b_t, issued as b_t + sigma_t * offsets.
  # Studentized scales stay within sigma_bounds * S0 (missing scales stay missing).
  bound <- function(scale) if (!s0$initialized) scale else
    pmin(pmax(scale, params$sigma_bounds[1] * s0$scale), params$sigma_bounds[2] * s0$scale)
  bounds_key <- paste(params$sigma_bounds, collapse = ",")
  sigma <- memo(paste("sigma", params$median_window, params$median_min_rows,
    params$residual_estimator, bounds_key), {
    sigma <- median_sigma(case, "residuals", params$median_window, params$median_min_rows,
      settings$scale_floor, transform, estimator = params$residual_estimator)
    sigma$scale <- bound(sigma$scale)
    sigma
  })
  base_sigma <- memo(paste("base_sigma", paste(params$median_base_levels, collapse = ","),
    bounds_key), {
    base_sigma <- median_sigma(case, "base", floor = settings$scale_floor,
      transform = transform, levels = params$median_base_levels)
    base_sigma$scale <- bound(base_sigma$scale)
    base_sigma
  })
  # The original base is issued until sigma_t exists (and during the S0 warmup without replay).
  median_warmup_weeks <- max(s0$warmup_weeks, sigma$warmup_weeks)
  # Offsets before the first update: normal quantiles, or t quantiles of a residual divided
  # by the RMS of median_window other normal residuals (wider tails for small windows).
  start <- switch(params$median_start, normal = qnorm(case$alpha),
    t = qt(case$alpha, params$median_window), stop("Unknown median start ", params$median_start))
  # The median methods learn at the adaptive rate of their median scores (y - b^.5) / sigma_t.
  median_matrix <- matrix(case$base[, match(.5, case$alpha)], nrow(case$base),
    length(case$alpha))
  median_rate <- function(scale) rate_of(unit(tracker_score(case, transform, scale,
    median_matrix)))
  if (want("studentized"))
    predictions$studentized <- median_multiqt(case, median_rate(sigma$scale),
      sigma$scale, start, median_warmup_weeks, transform) / count_scale
  # sigma_t from the width of the base: no feedback lag, so no warmup beyond S0's. This
  # sigma_t is a normal sd, so it starts from normal quantiles: the base's own interval.
  if (want("studentized_base"))
    predictions$studentized_base <- median_multiqt(case, median_rate(base_sigma$scale),
      base_sigma$scale, qnorm(case$alpha), max(s0$warmup_weeks, base_sigma$warmup_weeks),
      transform) / count_scale
  # Residual sigma frozen at its first issued value: additive MultiQT on the median.
  # It separates reducing the base to its median from studentizing.
  constant_sigma <- sigma$scale[median_warmup_weeks + 1L]
  if (want("median_constant"))
    predictions$median_constant <- median_multiqt(case, median_rate(constant_sigma),
      constant_sigma, start, median_warmup_weeks, transform) / count_scale
  # Median start: the same two methods (rate, scale, feedback and warmup) with all quantiles
  # starting at the median (offsets 0) instead of the normal (or t) quantiles.
  if (want("studentized_zero"))
    predictions$studentized_zero <- median_multiqt(case, median_rate(sigma$scale),
      sigma$scale, 0, median_warmup_weeks, transform) / count_scale
  if (want("median_constant_zero"))
    predictions$median_constant_zero <- median_multiqt(case, median_rate(constant_sigma),
      constant_sigma, 0, median_warmup_weeks, transform) / count_scale
  # Studentized MultiQT on the full base: offsets start at 0 (the base) and are issued as
  # b_t^alpha + sigma_t * offsets, with one scale per week (the base width; isoreg) or one
  # per level (weighted isoreg): the RMS of that level's recent residuals y - b^alpha.
  # Each learns at the adaptive rate of its own studentized scores (y - b^alpha) / sigma_t.
  full <- function(scales, warmup_weeks = s0$warmup_weeks) multiqt(case$base,
    case$feedback_observed, case$feedback_by_week, case$alpha,
    rate_of(unit(tracker_score(case, transform, scales))), warmup_weeks, transform,
    sigma = scales, feedback_samples = case$feedback_samples) / count_scale
  if (want("studentized_full"))
    predictions$studentized_full <- full(base_sigma$scale)
  # The tracker of width_flip: half the base IQR (half_iqr_sigma()), bounded like the base
  # width; a helper, not issued.
  if (want("width_flip"))
    predictions$studentized_iqr <- full(bound(half_iqr_sigma(case)))
  # The residual scales exist once median_min_rows residuals are released, as for the median.
  if (want("studentized_level_rms", "level_flip") || scored)
    predictions$studentized_level_rms <- full(bound(level_residual_sigma(case,
      params$median_window, params$median_min_rows, settings$scale_floor, transform,
      params$residual_estimator)), median_warmup_weeks)
  # The PIT of the base: offsets on the probability scale (pit_multiqt()) at the adaptive
  # rate of |PIT - alpha|; the same for either tracker, the same warmup as Fixed MultiQT.
  if (want("studentized_pit"))
    predictions$studentized_pit <- pit_multiqt(case$base, case$feedback_observed,
      case$feedback_by_week, case$alpha, rate_of(unit(pit_score(case))), s0$warmup_weeks,
      feedback_samples = case$feedback_samples) / count_scale

  # What if we knew the true peak? Hindsight benchmark: bypass correction for postpeak
  # targets, even if issued prepeak (each week's nearest peak). Setting eta=0 instead
  # would only freeze the accumulated offset.
  if (want("oracle_target")) {
    predictions$oracle_target <- predictions$fixed
    after_peak <- which(case$target_date > case$peak_date)
    predictions$oracle_target[after_peak, ] <- base[after_peak, ]
  }

  # Several trackers: weight their forecasts. Rate grid: one persistent expert per
  # multiplier of the rolling scale. Restarts: a fresh copy per multiplier at every issue
  # week; no permanent base expert. Leak bank: Leaky MultiQT for every leak of leak_grid.
  # Score-based experts: the forecasts of five methods (score_experts). Sign-flip experts:
  # the base, a studentized tracker and its mirror. Each learned weighting learns from
  # released quantile losses of issued forecasts, divided by the base's mean loss over the S0
  # rows, so the count scale does not change the weights; adaptive lambda and look-back H
  # (par(method, "H")) for every one. Equal averaging is a control.
  loss_scale <- if (!is.null(case$frozen_loss_scale)) case$frozen_loss_scale else
    weighting_loss_scale(case, s0, count_scale, settings$scale_floor)
  weigh <- function(method, experts, ...) exp_weights(experts, base, feedback,
    case$feedback_by_week, case$alpha, H = par(method, "H"),
    lambda_start = par(method, "lambda_start"), warmup_weeks = s0$warmup_weeks,
    loss_scale = loss_scale, feedback_samples = scale_samples(case$feedback_samples,
      count_scale), ...)
  # Without a cache (R/04) memo() only returns the value.
  mix <- function(method, value) memo(paste("mixture", method), value)
  mixtures <- list()
  rate_experts <- restart_experts <- fixed_restarts <- NULL
  if (want("uniform", "weighted", "weighted_common", "shared", "shared_common")) {
    rate_experts <- make_rate_experts(case, settings$expert_multipliers, count_unit,
      count_scale, transform)
    if (want("uniform")) mixtures$uniform <- mix("uniform", weigh("uniform", rate_experts,
      lambda = 0))
    if (want("weighted")) mixtures$weighted <- mix("weighted", weigh("weighted", rate_experts,
      sharing = FALSE))
    if (want("weighted_common")) mixtures$weighted_common <- mix("weighted_common",
      weigh("weighted_common", rate_experts, sharing = FALSE, per_quantile = FALSE))
    if (want("shared")) mixtures$shared <- mix("shared", weigh("shared", rate_experts))
    if (want("shared_common")) mixtures$shared_common <- mix("shared_common",
      weigh("shared_common", rate_experts, per_quantile = FALSE))
  }
  if (want("restart_uniform", "restart_weighted", "restart_persistent", "restart_recent")) {
    restart_experts <- make_restart_experts(case, params$restart_multiplier, count_unit,
      count_scale, transform)
    fixed_restarts <- restart_experts
    if (want("restart_uniform")) mixtures$restart_uniform <- mix("restart_uniform",
      weigh("restart_uniform", fixed_restarts, lambda = 0))
    if (want("restart_weighted")) mixtures$restart_weighted <- mix("restart_weighted",
      weigh("restart_weighted", fixed_restarts))
    if (want("restart_recent")) mixtures$restart_recent <- mix("restart_recent",
      weigh("restart_recent", fixed_restarts, keep = par("restart_recent", "restart_recent_k")))
  }
  # Two copies per weekly start, a slow and a fast one: learned weights with sharing on the
  # copies of the restart_recent_k most recent starts; both copies of the oldest start
  # retire together.
  if (want("restart_recent_dual")) {
    dual_multipliers <- par("restart_recent_dual", "restart_dual_multipliers")
    dual_restarts <- make_restart_experts(case, dual_multipliers, count_unit, count_scale,
      transform)
    mixtures$restart_recent_dual <- mix("restart_recent_dual", weigh("restart_recent_dual",
      dual_restarts, keep = length(dual_multipliers) *
        par("restart_recent_dual", "restart_recent_k")))
  }
  # The leak chosen online (notes, "The leak as a ridge penalty"): Leaky MultiQT at the
  # adaptive rate for every leak of leak_grid. leaky_select issues each week the forecast of
  # the leak whose forecasts had the smallest quantile loss over the last leak_window
  # learned rows (prequential choice, select_expert()); leaky_weighted a learned mixture with
  # one weight vector for all levels and sharing.
  if (want("leaky_select", "leaky_weighted")) {
    leak_bank <- make_leak_experts(case, params$leak_grid, rate_of(count_unit), count_scale,
      transform)
    if (want("leaky_select")) predictions$leaky_select <- select_expert(leak_bank, case,
      params$leak_window, count_scale, match(params$leak, params$leak_grid))$q
    if (want("leaky_weighted")) mixtures$leaky_weighted <- mix("leaky_weighted",
      weigh("leaky_weighted", leak_bank, per_quantile = FALSE))
  }
  # Studentized MultiQT on the full base with two scales of the most stationary scores
  # (scaled_methods): each alone at the adaptive rate of its scores, and in two expert banks,
  # the rate grid (expert_multipliers) and weekly starts (restart_multiplier) in multiples of
  # the rolling scale of its scores, each with equal and with learned weights as
  # weighted_common and restart_recent (parameter_source); the equal-weight starts also keep
  # the restart_recent_k most recent.
  scaled <- list()
  for (scale_name in c("power", "spread")) {
    names_of <- c(paste0("studentized_", scale_name), paste0(scale_name, c("_uniform",
      "_weighted_common", "_restart_equal", "_restart_recent")))
    needed <- want(names_of) || (scale_name == "power" &&
      want("power_flip", "phase_power", "phase_oracle", outside_methods, angle_controls,
        hidden_controls, hindsight_hidden))
    if (!needed) next
    scale_key <- paste(scale_name, params$power_p, params$spread_lambda, params$scale_guard)
    studentized <- memo(paste("scale", scale_key), normalized_scale(studentized_scale(case,
      scale_name, params), s0, params$scale_guard, case$frozen_scale_reference[[scale_name]]))
    warmup <- max(s0$warmup_weeks, studentized$warmup_weeks)
    scaled_unit <- unit(tracker_score(case, transform, studentized$scale))
    scaled[[scale_name]] <- list(scale = studentized$scale, warmup = warmup, unit = scaled_unit)
    if (want(names_of[1]) || (scale_name == "power" && want("power_flip")))
      predictions[[names_of[1]]] <- full(studentized$scale, warmup)
    if (want(names_of[2:3])) {
      grid <- make_rate_experts(case, settings$expert_multipliers, scaled_unit, count_scale,
        transform, studentized$scale, warmup)
      if (want(names_of[2])) mixtures[[names_of[2]]] <- mix(names_of[2], weigh(names_of[2],
        grid, lambda = 0))
      if (want(names_of[3])) mixtures[[names_of[3]]] <- mix(names_of[3], weigh(names_of[3],
        grid, sharing = FALSE, per_quantile = FALSE))
    }
    if (want(names_of[4:5])) {
      starts <- make_restart_experts(case, params$restart_multiplier, scaled_unit,
        count_scale, transform, studentized$scale, warmup)
      keep <- par(names_of[5], "restart_recent_k")
      if (want(names_of[4])) mixtures[[names_of[4]]] <- mix(names_of[4], weigh(names_of[4],
        starts, lambda = 0, per_quantile = FALSE, keep = keep))
      if (want(names_of[5])) mixtures[[names_of[5]]] <- mix(names_of[5], weigh(names_of[5],
        starts, keep = keep))
    }
  }
  # Score-based experts: one bank of the forecasts of the score_experts. score_equal: equal
  # weights. score_learned: learned weights per level as the other learned mixtures (H,
  # sharing 1/(2H)).
  score_bank <- NULL
  if (scored) {
    experts <- score_experts
    forecasts <- aperm(vapply(predictions[experts], identity, base), c(1, 3, 2))
    dimnames(forecasts) <- list(as.character(case$reference_date), experts,
      as.character(case$alpha))
    score_bank <- list(forecasts = forecasts, birth_week = rep(1L, length(experts)),
      experts = data.frame(expert = experts, start_date = case$reference_date[1]))
    if (want("score_equal")) mixtures$score_equal <- mix("score_equal",
      weigh("score_equal", score_bank, lambda = 0, per_quantile = FALSE))
    if (want("score_learned")) mixtures$score_learned <- mix("score_learned",
      weigh("score_learned", score_bank, centered = TRUE))
  }
  # Sign-flip experts (notes, "Sign-flip experts"): the base, a studentized tracker and its
  # mirror, the base minus the tracker's correction (flip_bank()), with learned weights per
  # level and sharing. The bias of the base changes sign at a peak; the mirror issues the
  # learned correction with the other sign, and the weights move to it when it starts to
  # lose less. Trackers: the power scale, the level RMS (paused) and half the base IQR.
  flips <- c(power_flip = "studentized_power", level_flip = "studentized_level_rms",
    width_flip = "studentized_iqr")
  for (method in intersect(names(flips), methods))
    mixtures[[method]] <- mix(method, weigh(method, flip_bank(base,
      predictions[[flips[[method]]]], case)))
  # Two-phase MultiQT on the power scale (phase_multiqt()): one state for the rise and one
  # for the decline of a wave, mixed by the probability of a decline. phase_power: that
  # probability from the recent errors of the base (decline_probability()); phase_oracle:
  # the phase in hindsight (after the nearest peak of the latest counts), a benchmark.
  if (want("phase_power", "phase_oracle")) {
    power <- scaled$power
    phase <- function(p) phase_multiqt(case$base, case$feedback_observed,
      case$feedback_by_week, case$alpha, rate_of(power$unit), p, power$scale,
      power$warmup, feedback_samples = case$feedback_samples) / count_scale
    if (want("phase_power"))
      predictions$phase_power <- phase(decline_probability(case, params$phase_rows))
    if (want("phase_oracle")) predictions$phase_oracle <- phase(ifelse(
      is.finite(case$relative_week), as.numeric(case$relative_week > 0), .5))
  }
  # Bias and phase probabilities (R/05_phase_probabilities.R; notes, "The angle from a phase
  # probability", "Bias and phase probabilities"), on the power scale at the
  # adaptive rate of phase_power: phase_angle, the angle with the probability that the base is
  # too high plus a slow per-level offset (phase_slow); phase_hidden, three free phases with the
  # hidden-phase model's probabilities, learned with the probabilities as known when each
  # report arrives. Locations without phase inputs (outside the pooled locations): NA.
  # phase_angle, phase_hidden and their controls (angle_controls, hidden_controls) also return
  # their issued components (diagnostics, in counts / count_scale like the forecasts, with the
  # probabilities, p_t or P_t, and the probability source); without phase inputs: NA forecasts
  # and diagnostics marked unavailable.
  diagnostics <- list()
  if (want(outside_methods, angle_controls, hidden_controls, hindsight_hidden)) {
    power <- scaled$power
    outside <- function(P, weights, V, rate, diagnostics = FALSE, with_p = TRUE) {
      x <- phase_weighted_multiqt(case$base, case$feedback_observed, case$feedback_by_week,
        case$alpha, rate_of(power$unit), power$scale, P, weights, V, rate, power$warmup,
        case$feedback_samples, diagnostics)
      if (!diagnostics) return(x / count_scale)
      d <- x$diagnostics
      for (part in c("components", "raw", "issued", "adjustment", "sigma"))
        d[[part]] <- d[[part]] / count_scale
      list(q = x$q / count_scale, diagnostics = c(list(available = TRUE),
        if (with_p) list(p = case$phase$p_bias), list(source = case$phase$source,
          md5 = case$phase$md5, base = base, count_scale = count_scale), d))
    }
    ph <- case$phase
    unavailable <- list(available = FALSE, reason = "no phase inputs for this location")
    angle <- cbind(shared = c(1, 1), flip = c(1, -1))
    one <- matrix(1, 1, 1, dimnames = list(NULL, "offset"))
    run_angle <- function(method, x) {
      if (is.null(ph)) { predictions[[method]] <<- base * NA; diagnostics[[method]] <<- unavailable }
      else { r <- x(); predictions[[method]] <<- r$q; diagnostics[[method]] <<- r$diagnostics }
    }
    angle_with <- function(p) function() {
      P2 <- cbind(1 - p, p)
      outside(P2, function(t, s) P2[s, ], angle, c(params$phase_slow, 1), TRUE)
    }
    if (want("phase_angle")) run_angle("phase_angle", angle_with(ph$p_bias))
    # The hidden-phase model and its controls (hidden_controls): three free phases (V = I,
    # named hidden_phases), issued with P and learning row s at week t with weights(t, s).
    n_rows <- nrow(case$base)
    phases_V <- diag(3)
    dimnames(phases_V) <- list(hidden_phases, hidden_phases)
    hidden_with <- function(P, weights) function() {
      colnames(P) <- hidden_phases
      outside(P, weights, phases_V, rep(1, 3), TRUE, with_p = FALSE)
    }
    # phase_hidden: issued with P_t, learned with W_{s,t} (the probabilities of row s's target
    # week as known at the learning cutoff t).
    if (want("phase_hidden")) run_angle("phase_hidden",
      hidden_with(ph$P, function(t, s) ph$W[s, t, ]))
    # Equal probabilities at issuance and learning.
    if (want("hidden_uniform")) run_angle("hidden_uniform",
      hidden_with(matrix(1 / 3, n_rows, 3), function(t, s) rep(1 / 3, 3)))
    # Learned with the issued probabilities P_s of the learned row.
    if (want("hidden_issued_weights")) run_angle("hidden_issued_weights",
      hidden_with(ph$P, function(t, s) ph$P[s, ]))
    # The true three phases (hidden_oracle, hindsight): rise up to and including the nearest
    # peak of the latest counts (relative week <= 0), overshoot weeks 1-4 after it, tail from
    # week 5; no nearest peak: thirds. Floored at .005 as phase_inputs() floors P and W; the
    # same weights at issuance and learning (the true phase of a week does not change).
    if (want(hindsight_hidden)) {
      rw <- case$relative_week
      truth <- cbind(rw <= 0, rw > 0 & rw <= 4, rw > 4) * 1
      truth[!is.finite(rw), ] <- 1 / 3
      truth <- pmin(pmax(truth, .005), 1)
      r <- hidden_with(truth, function(t, s) truth[s, ])()
      predictions$hidden_oracle <- r$q
      diagnostics$hidden_oracle <- r$diagnostics
    }
    for (method in intersect(names(params$angle_control_rates), methods))
      run_angle(method, local({ r <- params$angle_control_rates[[method]]; function()
        outside(matrix(1, n_rows, 1), function(t, s) 1, one, r, TRUE) }))
    # No reversals: p replaced by min(p, 1 - p), so 1 - 2p becomes |1 - 2p| in the forecast
    # (P_t V) and in the update (the saved P_s of the learned row).
    if (want("angle_abs")) run_angle("angle_abs", angle_with(pmin(ph$p_bias, 1 - ph$p_bias)))
    # Rate only: the weight of a report learned at week t is .1 + (1 - 2 p_t)^2 (p_t of the
    # learning week, not of the learned row); issued without a multiplier.
    if (want("angle_rate_only")) run_angle("angle_rate_only", function()
      outside(matrix(1, n_rows, 1), function(t, s) params$phase_slow + (1 - 2 * ph$p_bias[t])^2,
        one, 1, TRUE))
  }
  for (method in names(mixtures)) predictions[[method]] <- mixtures[[method]]$q
  # The first-start control has exactly the same birth rule as the restart copies.
  if (want("restart_persistent")) {
    predictions$restart_persistent <- base
    if (length(fixed_restarts$birth_week)) {
      issued <- seq.int(fixed_restarts$birth_week[1], nrow(base))
      predictions$restart_persistent[issued, ] <- fixed_restarts$forecasts[issued, 1, ]
    }
  }
  # Each outer tracker has its own offsets and learns from its OWN final forecasts.
  # Its state stays continuous when new restart copies join the underlying mixture.
  for (method in intersect(names(outer_inner), methods))
    predictions[[method]] <- adjust(predictions[[outer_inner[[method]]]] * count_scale,
      par(method, "outer_multiplier") * s0$scale, s0$warmup_weeks)
  # Only the methods asked for (and the score bank's members); helpers computed for another
  # method (e.g. the half-IQR tracker of width_flip) are not issued.
  predictions <- predictions[intersect(method_names, union(c("original", methods),
    if (scored) score_experts))]
  # expert_forecasts: the experts each mixture weighs [week x expert x level] in counts /
  # count_scale (04's weight pages).
  list(case = case, scaling = scaling, s0 = s0, sigma = sigma, base_sigma = base_sigma,
    median_warmup_weeks = median_warmup_weeks,
    constant_sigma = constant_sigma, predictions = predictions,
    diagnostics = diagnostics[intersect(names(diagnostics), names(predictions))],
    weights = lapply(mixtures, `[[`, "weights"), lambda = lapply(mixtures, `[[`, "lambda"),
    expert_forecasts = lapply(mixtures, `[[`, "experts"),
    rate_experts = rate_experts, restart_experts = restart_experts, score_experts = score_bank)
}

# Re-run design, as in Delphi's calibration reports: at every issue week t of `weeks`
# (default: all), every method is re-run from the first week with the reports current at
# t's cutoff (case_as_of()): the same rows are consumed in the same rounds as in the
# stateful design with this delay, but each with its report as known at t. Only the week-t
# forecasts, weights and expert forecasts are issued; the earlier forecasts
# of a re-run are reconstructions and generally differ from those issued then. S0 and the
# weighting loss scale are those of the stateful design with this delay (from its stateful
# feedback), frozen for every re-run. Rows outside `weeks` are NA; experts not yet started
# have weight 0. The other results (scales, experts) are those of the run at the last week;
# rerun_omitted: the rows left out of the re-run at each week because their current report
# was invalid.
# cache: one environment per week (see apply_methods()).
rerun_methods <- function(case, settings, params, methods, cache = NULL, weeks = NULL) {
  n <- nrow(case$base)
  weeks <- sort(if (is.null(weeks)) seq_len(n) else weeks)
  stateful <- case
  stateful$feedback_revised <- FALSE
  s0 <- error_scale(stateful, settings$warmup_rounds, settings$scale_floor,
    settings$transform, settings$replay_warmup)
  loss_scale <- weighting_loss_scale(stateful, s0,
    choose_count_scale(case, settings$scale_method)$scale, settings$scale_floor)
  # The scaled methods (scaled_methods) freeze the reference of normalized_scale() likewise.
  scale_reference <- lapply(setNames(nm = c("power", "spread")), function(name)
    if (any(grepl(name, intersect(methods, scaled_methods))))
      normalized_scale(studentized_scale(stateful, name, params), s0, params$scale_guard)$reference)
  run_at <- function(t) {
    week_cache <- NULL
    if (!is.null(cache)) {
      key <- paste("re-run week", t)
      if (is.null(cache[[key]])) cache[[key]] <- new.env()
      week_cache <- cache[[key]]
    }
    as_of <- case_as_of(case, t, s0, loss_scale)
    as_of$frozen_scale_reference <- scale_reference
    apply_methods(as_of, settings, params, methods, week_cache)
  }
  # n rows of NA with the other dimensions of x (one row per issue week).
  blank <- function(x) {
    if (is.data.frame(x)) {
      x <- x[rep(NA_integer_, n), , drop = FALSE]
      rownames(x) <- NULL
      return(x)
    }
    d <- if (is.null(dim(x))) c(length(x), 1L) else dim(x)
    array(NA_real_, c(n, d[-1]), c(list(NULL), dimnames(x)[-1]))
  }
  # Row t of x from the run at t (experts matched by name: weekly restarts join over time);
  # empty: the value of the experts missing at t (weight 0, forecast NA). Every expert of
  # the run at t must have a column (allocated from the run at the last week).
  fill <- function(target, x, t, empty = 0) {
    if (is.data.frame(target)) {
      target[t, ] <- x[t, ]
    } else if (length(dim(target)) == 3L) {
      target[t, , ] <- empty
      columns <- match(dimnames(x)[[2]], dimnames(target)[[2]])
      stopifnot(!anyNA(columns))
      target[t, columns, ] <- x[t, , ]
    } else {
      target[t, ] <- if (is.null(dim(x))) x[t] else x[t, ]
    }
    target
  }
  last <- run_at(max(weeks))
  parts <- c("predictions", "weights", "lambda", "expert_forecasts")
  issued <- lapply(last[parts], function(part) lapply(part, blank))
  # Diagnostics (phase_angle, phase_hidden and their controls): every per-week part from row t
  # of the replay at t, as the forecasts, including the scale sigma that replay used; the
  # learning summary (replay_parts) at t: its totals over the replay at t, which learns every
  # row reported by t with the weights known at t (case_as_of()).
  row_parts <- c("p", "P", "base", "states", "multiplier", "sigma", "eta", "components", "raw",
    "issued", "adjustment", "warmup")
  replay_parts <- c("learn_rows", "learn_issued", "learn_weights")
  issued_diagnostics <- lapply(last$diagnostics, function(d) if (!isTRUE(d$available)) d else
    c(d[setdiff(names(d), c(row_parts, replay_parts))],
      lapply(d[intersect(c(row_parts, replay_parts), names(d))], blank)))
  omitted <- list()
  for (t in weeks) {
    run <- if (t == max(weeks)) last else run_at(t)
    for (method in names(issued_diagnostics)) if (isTRUE(issued_diagnostics[[method]]$available)) {
      d <- issued_diagnostics[[method]]
      for (part in intersect(row_parts, names(d)))
        d[[part]] <- fill(d[[part]], run$diagnostics[[method]][[part]], t, NA)
      for (part in intersect(replay_parts, names(d)))
        d[[part]][t, ] <- colSums(as.matrix(run$diagnostics[[method]][[part]]))
      issued_diagnostics[[method]] <- d
    }
    # The score bank keeps its members in every re-run; an invalid report only removes its
    # feedback (rerun_omitted).
    for (method in intersect(names(run$weights), score_methods))
      stopifnot(identical(dimnames(run$weights[[method]])[[2]], score_experts),
        identical(dimnames(run$expert_forecasts[[method]])[[2]], score_experts))
    for (part in parts) for (method in names(issued[[part]]))
      issued[[part]][[method]] <- fill(issued[[part]][[method]], run[[part]][[method]], t,
        if (part == "expert_forecasts") NA else 0)
    if (length(run$case$rerun_omitted))
      omitted[[length(omitted) + 1L]] <- data.frame(week = t, row = run$case$rerun_omitted)
  }
  last[parts] <- issued
  last$diagnostics <- issued_diagnostics
  last$case <- case
  last$rerun_omitted <- do.call(rbind, c(list(data.frame(week = integer(), row = integer())),
    omitted))
  last$frozen <- list(s0 = s0, loss_scale = loss_scale)
  last
}

# 8. SAVED RESULTS OF LONG RUNS --------------------------------------------
# R/04 and R/05 save every finished piece of work (one location's run, one tuning stage) in
# <output_dir>/cache/<name>/, under a key of its arguments and of everything the function reads:
# the code of every function it calls (without comments) and the values of the global
# objects these use (codetools::findGlobals()). A restarted run, after a crash or a change,
# loads every piece whose key is unchanged and recomputes the others; changes to plots,
# colours or labels recompute nothing. Not in the key: the R and package versions and files
# read inside the function. Delete the cache folder to recompute everything.
cache_dir <- file.path(multiqt_settings()$output_dir, "cache")
cache_ignore <- "cores"                   # Global objects that never change a result.
# The plotted locations (settings$location_ids) only choose pages, and the base model and its
# output folder only choose the folders (each base has its own cache): left out of keys.
key_value <- function(x) {
  if (is.list(x) && all(c("location_ids", "pooled_locations") %in% names(x)))
    x[c("location_ids", "base_model", "output_dir")] <- NULL
  x
}
# The code and the global values a function reads, followed through the functions it calls
# (also functions inside lists).
code_key <- function(fun) {
  key <- list(deparse(fun))
  queue <- list(fun)
  while (length(queue)) {
    names_used <- setdiff(codetools::findGlobals(queue[[1]]), c(names(key), cache_ignore))
    queue <- queue[-1]
    for (name in names_used) if (exists(name, globalenv(), inherits = FALSE)) {
      value <- get(name, globalenv())
      functions <- if (is.function(value)) list(value) else if (is.list(value))
        Filter(is.function, value) else list()
      key[[name]] <- if (is.function(value)) deparse(value) else
        rapply(list(key_value(value)), deparse, classes = "function", how = "replace")
      queue <- c(queue, functions)
    }
  }
  key
}
# fun(...), loaded from <output_dir>/cache/<name>/ if saved under the same key, else computed and
# saved (written to a temporary file first, so an interrupted write leaves no result).
cached <- function(name, fun, ...) {
  key <- tools::md5sum(bytes = serialize(list(code_key(fun), lapply(list(...), key_value)),
    NULL))
  file <- file.path(cache_dir, name, paste0(key, ".rds"))
  if (file.exists(file)) return(readRDS(file))
  value <- fun(...)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  temporary <- file.path(dirname(file), paste0("tmp_", Sys.getpid(), "_", basename(file)))
  saveRDS(value, temporary)
  file.rename(temporary, file)
  value
}

# PHASE AND BIAS PROBABILITIES (R/05_phase_probabilities.R; notes, "Bias and phase
# probabilities") -------------------------------------------------------------------------
# The inputs of phase_angle (p_t, the probability that the report falls below the base median)
# and phase_hidden (P_t^k and W_{s,t}^k, the probabilities of the phases rise, overshoot, tail),
# computed by R/05 at every issuance cutoff of the base from what is known then, for the pooled
# locations. Every function takes the data it uses as arguments (`pp`, the list built by R/05:
# states, incomplete, reports, nssp, nrevss, erviss); none reads files or global data.
#   bias   one binomial GAM per disease and base (G0): logit p = Xlin beta + f(year_phase) +
#          f(g1) + a_state, refitted at every cutoff on the labels known then; the label is the
#          report as of target-week end + 11 days (COVID) or + 4 days (flu);
#   phase  a hidden semi-Markov chain (40 rise ages, 4 overshoot weeks, 40 tail ages), reset at
#          reporting gaps, with disease-specific transitions and observations: COVID, exit hazards
#          from the durations of the European waves (ERVISS) and conditional-growth GAMs of the
#          hospital and NSSP growth (O2) at the NSSP lead chosen per 13-week block from earlier
#          blocks; flu, exit-hazard GAMs of duration and calendar on the completed episodes of the
#          hospital reports and NREVSS, and phase-mean Gaussians of both growths at lead 0 (O1).
pp_classes <- c("rise", "over", "tail")
pp_settings <- list(
  label_days = c(covid = 11, flu = 4),          # Label: report as of target-week end + days.
  clip = .02,                                   # Operational p_t in [.02, .98].
  bias_terms = c("base_growth", "horizon", "hol_break", "hol_post", "hol_thanks"),
  bias_rule = list(min_n = 200, min_class = 20, min_weeks = 8, min_states = 3,
    max_coef_ratio = .2),
  bias_distinct = list(year_phase = 12, g1 = 8),
  peak = list(halfwidth = 8L, drop = .3, min_prominence = .1, smooth = 3L, min_age = 56,
    pair_age = 84),
  D = 40L, ahead = 6L, window_weeks = 80L,
  hazard_rule = list(min_events = 20L, min_nonevents = 20L, min_states = 3L,
    min_span_days = 730, min_ages = 4L, min_year_phase = 12L),
  hazard_clip = 1e-6,
  sd_floor = .05, min_rows = 20L, min_lag_distinct = 8L, rho_prior = 50, rho_cap = .9,
  nrevss_lag = 14, erviss_known = 182, erviss_min_weeks = 100, erviss_min_median = 5,
  covid_variants = paste0("O2_L", 0:2), flu_variant = "O1_L0", block_weeks = 13L)

# Peaks of a weekly series x: local maxima of its centred `smooth`-week mean s that are
# the largest value of s within `halfwidth` weeks on either side (ties: the earliest) and
# whose prominence is at least `drop` times their height and `min_prominence` times
# `scale`. Prominence: the height above the higher of the two lowest points between the
# week and the nearest higher week on each side (or the end of the series), so a
# maximum at the end of the series has none. Each peak is placed at the largest x within
# one week. Returns the rows of x and the relative prominence. (R/01's peaks and the episodes
# of the phase model.)
find_peaks <- function(x, halfwidth = 8L, drop = .3, min_prominence = .1, smooth = 3L,
                       scale = max(x)) {
  n <- length(x)
  s <- if (smooth > 1L) as.numeric(stats::filter(x, rep(1 / smooth, smooth), sides = 2)) else x
  s[is.na(s)] <- x[is.na(s)]
  prominence <- vapply(seq_len(n), function(i) {
    higher_left <- which(s[seq_len(i - 1L)] > s[i])
    higher_right <- i + which(s[i + seq_len(n - i)] > s[i])
    left <- min(s[(if (length(higher_left)) max(higher_left) else 1L):i])
    right <- min(s[i:(if (length(higher_right)) min(higher_right) else n)])
    s[i] - max(left, right)
  }, numeric(1))
  local_max <- vapply(seq_len(n), function(i) {
    around <- max(1L, i - halfwidth):min(n, i + halfwidth)
    s[i] == max(s[around]) && !any(s[around[around < i]] == s[i])
  }, logical(1))
  top <- which(local_max & prominence > 0 & prominence >= drop * s &
    prominence >= min_prominence * scale)
  row <- vapply(top, function(i) {
    near <- max(1L, i - 1L):min(n, i + 1L)
    near[which.max(x[near])]
  }, integer(1))
  data.frame(row, relative_prominence = prominence[top] / s[top])
}

# REPORTS AS KNOWN AT A DATE --------------------------------------------------------------------
# r: the versions of one disease and location (target_end_date, version, value; a masked
# version has value NA and stays, so it is never replaced by an older one).
# The reports known at `date`: the latest version per target week.
report_known_at <- function(r, date) {
  r <- r[r$version <= date, ]
  r <- r[!duplicated(r$target_end_date, fromLast = TRUE), ]
  r[order(r$target_end_date), c("target_end_date", "value")]
}
# Value of the targets `target` as of the dates `date` (one each): the latest version by then.
report_value_at <- function(r, target, date) {
  idx <- split(seq_len(nrow(r)), as.numeric(r$target_end_date))
  vapply(seq_along(target), function(i) {
    j <- idx[[as.character(as.numeric(target[i]))]]
    if (is.null(j)) return(NA_real_)
    k <- findInterval(as.numeric(date[i]), as.numeric(r$version[j]))
    if (k == 0L) NA_real_ else r$value[j][k]
  }, numeric(1))
}
# Known at `date` with a value (a masked latest version leaves its week out); _complete also
# without the weeks of incomplete reporting (pp$incomplete, inclusive).
pp_known_valid <- function(r, date) {
  x <- report_known_at(r, date)
  x[is.finite(x$value), ]
}
pp_known_complete_valid <- function(pp, r, date) {
  x <- report_known_at(r, date)
  x <- x[x$target_end_date < pp$incomplete[1] | x$target_end_date > pp$incomplete[2], ]
  x[is.finite(x$value), ]
}
# Revision shift at a cutoff: for the ages a = 0, 1, 2 weeks after the newest report (lag_days
# after its target week), the median of log1p(now) - log1p(at that age) over the reports of the
# pooled states whose target week lies 9 to 61 weeks before the cutoff and after the weeks of
# incomplete reporting; 0 with fewer than 8 such pairs. Added to log1p of the newest three
# reports (a = 0 the newest).
pp_revision_shift <- function(pp, disease, cutoff, lag_days) {
  vapply(0:2, function(a) {
    d <- unlist(lapply(pp$states, function(l) {
      r <- pp$reports[[paste(disease, l)]]
      targets <- unique(r$target_end_date[r$target_end_date <= cutoff - 63 &
        r$target_end_date > cutoff - 63 - 364 & r$target_end_date > pp$incomplete[2]])
      if (!length(targets)) return(NULL)
      now <- report_value_at(r, targets, rep(cutoff, length(targets)))
      young <- report_value_at(r, targets, targets + lag_days + 7 * a)
      log1p(now) - log1p(young)
    }))
    d <- d[is.finite(d)]
    if (length(d) < 8L) 0 else median(d)
  }, numeric(1))
}
# Report age of the newest report at a cutoff (days from its target week end), from the first
# pooled state: with all versions (bias inputs) or without incomplete weeks (phase series).
pp_lag <- function(pp, disease, cutoff) {
  r <- pp_known_valid(pp$reports[[paste(disease, pp$states[1])]], cutoff)
  as.numeric(cutoff - max(r$target_end_date))
}
pp_hidden_shift <- function(pp, disease, cutoff) {
  lag <- as.numeric(cutoff - max(pp_known_complete_valid(pp,
    pp$reports[[paste(disease, pp$states[1])]], cutoff)$target_end_date))
  pp_revision_shift(pp, disease, cutoff, lag)
}
# Holidays of target weeks ending t: hol_break, the share of the week's 7 days in 22 December -
# 3 January; hol_post, the week ends 1-21 days after 3 January; hol_thanks, it ends 0-6 days
# after the fourth Thursday of November.
holiday_covariates <- function(t) {
  in_break <- function(d) { md <- format(d, "%m-%d"); md >= "12-22" | md <= "01-03" }
  jan3 <- as.Date(paste0(format(t, "%Y"), "-01-03"))
  thanksgiving <- vapply(as.integer(format(t, "%Y")), function(y) {
    nov <- as.Date(paste0(y, "-11-01")) + 0:29
    as.numeric(nov[format(nov, "%u") == "4"][4])
  }, numeric(1))
  days <- matrix(vapply(0:6, function(k) in_break(t - k), logical(length(t))), length(t), 7)
  data.frame(hol_break = rowMeans(days),
    hol_post = as.numeric(as.numeric(t - jan3) %in% 1:21),
    hol_thanks = as.numeric((as.numeric(t) - thanksgiving) %in% 0:6))
}
# Fraction of the calendar year elapsed at date d (days since 1 January over the days of the year).
year_phase <- function(d) {
  y <- as.integer(format(d, "%Y"))
  start <- as.Date(paste0(y, "-01-01")); end <- as.Date(paste0(y + 1L, "-01-01"))
  as.numeric(d - start) / as.numeric(end - start)
}

# BIAS PROBABILITY ------------------------------------------------------------------------------
# The rows of every case (pooled states, every horizon of the base): base growth (horizon-3 minus
# horizon-0 base median of the reference week, log1p), the label B = 1{report as of target +
# label_days < base median} (missing when no report or a masked one), the mature label (latest
# report; evaluation only), holidays, year_phase, and g1 at R/05's feature cutoffs (the cutoffs of
# the horizon-0 cases): the newest weekly growth of the revision-shifted log1p reports known at
# the cutoff, with at least four reports and the newest two 7 days apart.
pp_bias_rows <- function(pp, prepared, horizons, label_days = pp_settings$label_days,
                         cores = 1L) {
  med <- 12L
  grid <- expand.grid(disease = c("covid", "flu"), location = pp$states, horizon = horizons,
    stringsAsFactors = FALSE)
  rows <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), function(i) {
    disease <- grid$disease[i]; location <- grid$location[i]; horizon <- grid$horizon[i]
    case <- get_case(prepared, disease, location, horizon, 0L, FALSE)
    stopifnot(case$alpha[med] == .5)
    h <- case$context_forecasts
    h <- h[!is.na(h$as_of) & is.finite(h$median) & h$median >= 0, ]
    m <- function(ref, hz) log1p(h$median[match(paste(ref, hz), paste(h$reference_date,
      h$horizon))])
    r <- pp$reports[[paste(disease, location)]]
    label_date <- case$target_date + label_days[[disease]]
    now <- report_value_at(r, case$target_date, label_date)
    latest <- report_value_at(r, case$target_date, rep(as.Date("2100-01-01"),
      length(case$target_date)))
    data.frame(disease, location, horizon, row = seq_len(nrow(case$base)),
      reference_date = case$reference_date, target_date = case$target_date,
      cutoff = case$as_of, period = case$period,
      base_growth = m(case$reference_date, 3) - m(case$reference_date, 0), label_date,
      label_report = now, base_median = case$base[, med],
      bias = as.numeric(now < case$base[, med]), bias_mature = as.numeric(latest <
        case$base[, med]), holiday_covariates(case$target_date),
      year_phase = year_phase(case$target_date))
  }, mc.cores = cores))
  jobs <- do.call(rbind, lapply(c("covid", "flu"), function(d) data.frame(disease = d,
    cutoff = sort(unique(do.call(c, lapply(pp$states, function(l)
      prepared$cases[[paste(d, l, 0, sep = "_")]]$as_of)))))))
  shifts <- parallel::mclapply(seq_len(nrow(jobs)), function(i) pp_revision_shift(pp,
    jobs$disease[i], jobs$cutoff[i], pp_lag(pp, jobs$disease[i], jobs$cutoff[i])),
    mc.cores = cores)
  shift_of <- setNames(shifts, paste(jobs$disease, jobs$cutoff))
  feats <- merge(jobs, data.frame(location = pp$states), by = NULL)
  feats$g1 <- unlist(parallel::mclapply(seq_len(nrow(feats)), function(i) {
    r <- pp_known_valid(pp$reports[[paste(feats$disease[i], feats$location[i])]],
      feats$cutoff[i])
    last <- nrow(r)
    if (last < 4L) return(NA_real_)
    shift <- shift_of[[paste(feats$disease[i], feats$cutoff[i])]]
    y <- log1p(r$value[last - 0:3]) + c(shift, shift[3])
    if (as.numeric(r$target_end_date[last] - r$target_end_date[last - 1]) != 7) NA_real_
    else y[1] - y[2]
  }, mc.cores = cores))
  rows <- merge(rows, feats, by = c("disease", "location", "cutoff"), all.x = TRUE, sort = FALSE)
  rows <- rows[order(rows$disease, rows$location, rows$horizon, rows$row), ]
  rownames(rows) <- NULL
  list(rows = rows, shifts = shift_of)
}
# Runs expr, recording every warning (text) and an error instead of stopping.
capture <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(tryCatch(expr, error = function(e) structure(conditionMessage(e),
    class = "capture_error")), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  error <- if (inherits(value, "capture_error")) as.character(value) else NA_character_
  if (!is.na(error)) value <- NULL
  list(value = value, warnings = unique(warnings), error = error)
}
# Linear design columns and their training-only standardization (columns with sd <= 1e-10 of the
# training rows are dropped and recorded).
numeric_design <- function(data, terms) {
  X <- vapply(terms, function(term) {
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    Reduce(`*`, lapply(parts, function(v) as.numeric(data[[v]])))
  }, numeric(nrow(data)))
  if (is.null(dim(X))) X <- matrix(X, nrow(data), length(terms), dimnames = list(NULL, terms))
  X
}
fit_scaler <- function(X) {
  if (!ncol(X)) return(list(centre = numeric(), scale = numeric(), columns = character(),
    dropped = character()))
  centre <- colMeans(X)
  scale <- apply(X, 2, sd)
  keep <- is.finite(scale) & scale > 1e-10
  list(centre = centre[keep], scale = scale[keep], columns = colnames(X)[keep],
    dropped = colnames(X)[!keep])
}
apply_scaler <- function(s, X) {
  X <- X[, s$columns, drop = FALSE]
  sweep(sweep(X, 2, s$centre, "-"), 2, s$scale, "/")
}
# Eligibility of the training rows: both labels with min_class rows each, min_n rows, min_weeks
# target weeks, min_states states and the distinct values of the listed smooth variables.
train_check <- function(train, rule, y = "bias", week = "target_date", state = "location",
                        distinct = NULL) {
  n <- nrow(train)
  ones <- sum(train[[y]] == 1); zeros <- sum(train[[y]] == 0)
  why <- character()
  if (n < rule$min_n) why <- c(why, sprintf("rows %d < %d", n, rule$min_n))
  if (min(ones, zeros) < rule$min_class) why <- c(why, sprintf("smaller class %d < %d",
    min(ones, zeros), rule$min_class))
  weeks <- length(unique(train[[week]]))
  if (weeks < rule$min_weeks) why <- c(why, sprintf("target weeks %d < %d", weeks, rule$min_weeks))
  states <- length(unique(train[[state]]))
  if (states < rule$min_states) why <- c(why, sprintf("states %d < %d", states, rule$min_states))
  for (v in names(distinct)) {
    k <- length(unique(train[[v]]))
    if (k < distinct[[v]]) why <- c(why, sprintf("distinct %s %d < %d", v, k, distinct[[v]]))
  }
  list(ok = !length(why), reason = paste(why, collapse = "; "), n = n, ones = ones,
    zeros = zeros, weeks = weeks, states = states)
}
# Coefficients, smoothing parameters, edf per term (Xlin summed), penalty ranks of a GAM.
gam_summary <- function(m) {
  labels <- vapply(m$smooth, function(s) s$label, "")
  edf <- setNames(vapply(m$smooth, function(s) sum(m$edf[s$first.para:s$last.para]), 1), labels)
  lin <- grep("^Xlin", names(coef(m)))
  list(coef = coef(m), sp = m$sp, edf = c(Xlin = sum(m$edf[lin]), edf), edf_total = sum(m$edf),
    penalty_rank = setNames(vapply(m$smooth, function(s) sum(vapply(s$S, function(S)
      qr(S)$rank, 1)), 1), labels),
    n_coef_smooth = setNames(vapply(m$smooth, function(s) s$last.para - s$first.para + 1, 1),
      labels),
    design_rank = qr(model.matrix(m))$rank, n_coef = length(coef(m)),
    outer = if (is.null(m$outer.info$conv)) NA_character_ else m$outer.info$conv)
}
# The bias GAM: B ~ Xlin + s(year_phase, cc, 6) + s(g1, cr, 4) + s(state, re), binomial, logit,
# REML, select = TRUE, cyclic year_phase on [0, 1], Xlin the standardized linear terms with an
# identity (ridge) penalty whose strength REML estimates; an unpenalized intercept.
pp_bias_formula <- B ~ Xlin + s(year_phase, bs = "cc", k = 6) + s(g1, bs = "cr", k = 4) +
  s(state, bs = "re")
pp_bias_data <- function(data, states, scaler = NULL) {
  X <- numeric_design(data, pp_settings$bias_terms)
  if (is.null(scaler)) scaler <- fit_scaler(X)
  d <- data.frame(B = if (is.null(data$bias)) NA_real_ else data$bias,
    year_phase = data$year_phase, g1 = data$g1, state = factor(data$location, levels = states))
  d$Xlin <- apply_scaler(scaler, X)
  list(data = d, scaler = scaler)
}
pp_bias_gam <- function(train, states, max_coef_ratio = pp_settings$bias_rule$max_coef_ratio) {
  g <- pp_bias_data(train, states)
  p <- ncol(g$data$Xlin)
  if (!p) return(list(error = "no non-constant linear column", warnings = character()))
  r <- capture(mgcv::gam(pp_bias_formula, family = binomial(link = "logit"), data = g$data,
    method = "REML", select = TRUE, knots = list(year_phase = c(0, 1)),
    paraPen = list(Xlin = list(diag(p))), drop.unused.levels = FALSE))
  m <- r$value
  if (!is.null(m)) {
    conv <- m$outer.info$conv
    if (!isTRUE(m$converged)) r$error <- "inner iteration did not converge"
    else if (!is.null(conv) && !conv %in% c("full convergence", "step failed"))
      r$error <- paste("outer:", conv)
    else if (length(coef(m)) > max_coef_ratio * nrow(train))
      r$error <- sprintf("%d coefficients > %.2f x %d rows", length(coef(m)), max_coef_ratio,
        nrow(train))
    else if (!all(is.finite(coef(m)))) r$error <- "non-finite coefficients"
    if (!is.null(conv) && conv == "step failed")
      r$warnings <- c(r$warnings, "outer REML step failed (smoothing parameter at a boundary)")
  }
  list(fit = if (is.na(r$error)) m else NULL, scaler = g$scaler, states = states,
    train_states = sort(unique(train$location)), warnings = r$warnings, error = r$error)
}
# Prediction; a state without training rows gets zero state deviation.
pp_bias_predict <- function(obj, newdata) {
  g <- pp_bias_data(newdata, obj$states, obj$scaler)
  labels <- vapply(obj$fit$smooth, function(s) s$label, "")
  st <- labels[grepl("state", labels)]
  eta <- predict(obj$fit, g$data, type = "link")
  unseen <- !newdata$location %in% obj$train_states
  if (any(unseen)) eta[unseen] <- predict(obj$fit, g$data[unseen, , drop = FALSE],
    type = "link", exclude = st)
  as.numeric(plogis(eta))
}
# One disease at one issuance cutoff: trained on the rows issued before the cutoff whose label has
# arrived (label_date <= cutoff) and whose inputs are complete; predicted for the rows issued at
# the cutoff. Ineligible or failed fit, or missing prediction inputs: the smoothed label frequency
# (sum B + 1) / (n + 2) over every label known by then (1/2 before any). p_raw and the operational
# p (clipped to [.02, .98]), source ("gam" / "frequency"), reason, record, object (keep).
pp_bias_fit_at <- function(rows, disease, cutoff, states, keep = FALSE) {
  x <- rows$disease == disease
  inputs <- c("year_phase", pp_settings$bias_terms, "g1")
  complete <- Reduce(`&`, lapply(inputs, function(v) is.finite(rows[[v]])))
  labels <- which(x & is.finite(rows$bias) & rows$label_date <= cutoff & rows$cutoff < cutoff)
  train <- rows[labels[complete[labels]], ]
  now <- which(x & rows$cutoff == cutoff)
  frequency <- (sum(rows$bias[labels]) + 1) / (length(labels) + 2)
  chk <- train_check(train, pp_settings$bias_rule, distinct = pp_settings$bias_distinct)
  record <- data.frame(disease, cutoff, labels_available = length(labels), frequency,
    n = chk$n, ones = chk$ones, zeros = chk$zeros, target_weeks = chk$weeks,
    states = chk$states, eligible = chk$ok, reason = chk$reason, fitted = FALSE,
    error = NA_character_, warnings = "", dropped_columns = "", n_coef = NA_integer_, edf = "",
    sp = "", outer = NA_character_, issued_rows = length(now),
    complete_issued_rows = sum(complete[now]))
  stopifnot(!nrow(train) || (max(train$label_date) <= cutoff && max(train$cutoff) < cutoff))
  p <- rep(NA_real_, length(now))
  why <- rep(if (!chk$ok) paste("training rows:", chk$reason) else "", length(now))
  obj <- NULL
  if (chk$ok) {
    obj <- pp_bias_gam(train, states)
    record$warnings <- paste(obj$warnings, collapse = "; ")
    record$dropped_columns <- paste(obj$scaler$dropped, collapse = " ")
    if (!is.null(obj$fit)) {
      use <- complete[now]
      pg <- if (any(use)) pp_bias_predict(obj, rows[now[use], ]) else numeric()
      if (!all(is.finite(pg))) obj$error <- "non-finite predictions" else {
        p[use] <- pg
        why[!use] <- "prediction inputs missing"
        s <- gam_summary(obj$fit)
        record$fitted <- TRUE
        record$n_coef <- s$n_coef
        record$edf <- paste(sprintf("%s=%.3g", names(s$edf), s$edf), collapse = " ")
        record$sp <- paste(sprintf("%s=%.4g", names(obj$fit$sp), obj$fit$sp), collapse = " ")
        record$outer <- s$outer
      }
    }
    if (is.null(obj$fit) || !is.na(obj$error)) {
      record$error <- obj$error
      why[] <- paste("fit failed:", obj$error)
      p[] <- NA_real_
    }
  }
  source <- ifelse(is.finite(p), "gam", "frequency")
  p_raw <- ifelse(is.finite(p), p, frequency)
  list(record = record, rows = now, p_raw = p_raw,
    p = pmin(pmax(p_raw, pp_settings$clip), 1 - pp_settings$clip), source = source,
    reason = why, object = if (keep && record$fitted) obj)
}

# PHASE PROBABILITIES: SERIES AND EPISODES ---------------------------------------------------------
# One disease and state as known at the cutoff: the complete weeks (incomplete reporting and
# masked versions left out), x = log1p(count) with the newest three shifted by the revision
# shift, contiguous weekly segments, growth G = diff(x) within a segment, holidays, year_phase.
pp_series <- function(pp, disease, location, cutoff, shift) {
  s <- pp_known_complete_valid(pp, pp$reports[[paste(disease, location)]], cutoff)
  s$x <- log1p(s$value)
  n <- nrow(s)
  if (n >= 3L) s$x[n - 0:2] <- s$x[n - 0:2] + shift
  gap <- c(TRUE, diff(as.numeric(s$target_end_date)) != 7)[seq_len(n)]
  s$segment <- cumsum(gap)
  s$growth <- c(NA, diff(s$x))[seq_len(n)]
  s$growth[gap] <- NA
  h <- holiday_covariates(s$target_end_date)
  s$hol_break <- h$hol_break
  s$hol_post <- h$hol_post
  s$year_phase <- year_phase(s$target_end_date)
  s
}
# Peaks within each contiguous segment (find_peaks(value, 8, .3, .1, 3, scale = the segment's
# maximum); a segment shorter than 3 weeks has none), at least 56 days old at the cutoff.
pp_peaks <- function(s, cutoff, peak = pp_settings$peak) {
  out <- lapply(unique(s$segment), function(g) {
    k <- which(s$segment == g)
    v <- s$value[k]
    if (length(v) < 3L || !any(v > 0)) return(NULL)
    rows <- find_peaks(v, peak$halfwidth, peak$drop, peak$min_prominence, peak$smooth,
      scale = max(v))$row
    if (!length(rows)) return(NULL)
    data.frame(week = s$target_end_date[k[rows]], segment = g)
  })
  p <- do.call(rbind, out)
  if (is.null(p)) return(data.frame(week = as.Date(character()), segment = integer()))
  p[p$week <= cutoff - peak$min_age, ]
}
# Adjacent retained peaks a < b give one completed tail episode (a + 35 .. m, exit at m) and one
# completed rise episode (m + 7 .. b, exit at b), m = a + 7 floor((b - a) / 14) days; overshoot a + 7
# .. a + 28. Admitted only if a and b lie in one segment, b <= t - 84 and m >= a + 35.
pp_pairs <- function(peaks, cutoff, pair_age = pp_settings$peak$pair_age) {
  if (nrow(peaks) < 2L) return(data.frame(a = as.Date(character()), b = as.Date(character()),
    m = as.Date(character()), admitted = logical(), reason = character()))
  peaks <- peaks[order(peaks$week), ]
  a <- peaks$week[-nrow(peaks)]; b <- peaks$week[-1]
  m <- a + 7 * floor(as.numeric(b - a) / 14)
  reason <- ifelse(peaks$segment[-nrow(peaks)] != peaks$segment[-1], "gap between the peaks",
    ifelse(b > cutoff - pair_age, "second peak less than 84 days old",
      ifelse(m < a + 35, "peaks too close (midpoint before week 5 after the first)", "")))
  data.frame(a, b, m, admitted = reason == "", reason)
}
# Labelled weeks and risk rows of the admitted pairs of one series s. Risk rows: every week of an
# episode with its age d (1 = first week), the capped age min(d, D), the year fraction of the
# risk week and event = 1 only in the episode's last week (the exit at the end of that week).
pp_episode_rows <- function(s, pairs, location, D = pp_settings$D) {
  pairs <- pairs[pairs$admitted, , drop = FALSE]
  if (!nrow(pairs)) return(list(labels = NULL, risk = NULL))
  labels <- list(); risk <- list()
  for (i in seq_len(nrow(pairs))) {
    a <- pairs$a[i]; b <- pairs$b[i]; m <- pairs$m[i]
    over <- a + 7 * (1:4)
    tail <- seq(a + 35, m, by = 7)
    rise <- seq(m + 7, b, by = 7)
    weeks <- c(over, tail, rise)
    stopifnot(all(weeks %in% s$target_end_date))
    labels[[i]] <- data.frame(location, pair = i, week = weeks,
      phase = factor(rep(c("over", "tail", "rise"), c(4L, length(tail), length(rise))),
        pp_classes))
    episode <- function(w, phase) data.frame(location, pair = i, phase, week = w,
      age = seq_along(w), age_cap = pmin(seq_along(w), D),
      year_phase = year_phase(w), event = as.numeric(seq_along(w) == length(w)))
    risk[[i]] <- rbind(episode(tail, "tail"), episode(rise, "rise"))
  }
  labels <- do.call(rbind, labels)
  k <- match(labels$week, s$target_end_date)
  labels <- cbind(labels, s[k, c("growth", "hol_break", "hol_post")])
  list(labels = labels, risk = do.call(rbind, risk))
}
# An outside series (NREVSS) as known at the cutoff: the weeks at least `lag` days old with a
# value, contiguous weekly segments.
pp_long_series <- function(x, cutoff, lag) {
  x <- x[x$target_end_date <= cutoff - lag & is.finite(x$value), ]
  gap <- c(TRUE, diff(as.numeric(x$target_end_date)) != 7)[seq_len(nrow(x))]
  data.frame(target_end_date = x$target_end_date, value = x$value, segment = cumsum(gap),
    growth = rep(NA_real_, nrow(x)), hol_break = rep(NA_real_, nrow(x)),
    hol_post = rep(NA_real_, nrow(x)))
}
# The NREVSS episode risk rows at a cutoff (same episode rules), series ids "nrevss:<state>".
pp_nrevss_risk <- function(pp, cutoff, D = pp_settings$D) {
  risk <- list(); log <- list()
  for (id in names(pp$nrevss)) {
    s <- pp_long_series(pp$nrevss[[id]], cutoff, pp_settings$nrevss_lag)
    if (nrow(s) < 3L) next
    pr <- pp_pairs(pp_peaks(s, cutoff), cutoff)
    sid <- paste0("nrevss:", id)
    if (nrow(pr)) log[[id]] <- cbind(location = sid, pr)
    risk[[id]] <- pp_episode_rows(s, pr, sid, D)$risk
  }
  list(risk = do.call(rbind, risk), pairs = do.call(rbind, log))
}

# PHASE PROBABILITIES: TRANSITIONS --------------------------------------------------------------
# Exit hazards of one phase from its risk rows: the GAM event ~ s(age_cap, cr, 4) +
# s(year_phase, cc, 6), binomial, REML, select = TRUE, cyclic on [0, 1], when eligible (20
# events, 20 non-events, 3 series, a 730-day span, 4 capped ages, 12 year fractions, a full-rank
# design, a converged fit); else the constant (E + 1) / (N + 2) when both classes exist; else
# unavailable.
pp_hazard_fit <- function(risk, rule = pp_settings$hazard_rule) {
  E <- sum(risk$event); N <- nrow(risk)
  rec <- data.frame(events = E, risk_rows = N, nonevents = N - E,
    states = length(unique(risk$location)),
    span_days = if (N) as.numeric(max(risk$week) - min(risk$week)) else 0,
    first_week = if (N) min(risk$week) else as.Date(NA),
    last_week = if (N) max(risk$week) else as.Date(NA),
    distinct_ages = length(unique(risk$age_cap)),
    distinct_year_phase = length(unique(risk$year_phase)),
    max_age = if (N) max(risk$age) else NA_integer_,
    episodes = length(unique(paste(risk$location, risk$pair))),
    kind = NA_character_, reason = "", error = NA_character_, warnings = "",
    edf_age = NA_real_, edf_year_phase = NA_real_, sp = "", constant = NA_real_)
  if (E == 0L || N - E == 0L) {
    rec$kind <- "unavailable"
    rec$reason <- "no events or no non-events"
    return(list(record = rec))
  }
  rec$constant <- (E + 1) / (N + 2)
  why <- character()
  if (E < rule$min_events) why <- c(why, sprintf("events %d < 20", E))
  if (N - E < rule$min_nonevents) why <- c(why, sprintf("non-events %d < 20", N - E))
  if (rec$states < rule$min_states) why <- c(why, sprintf("states %d < 3", rec$states))
  if (rec$span_days < rule$min_span_days)
    why <- c(why, sprintf("span %.0f days < 730", rec$span_days))
  if (rec$distinct_ages < rule$min_ages)
    why <- c(why, sprintf("distinct ages %d < 4", rec$distinct_ages))
  if (rec$distinct_year_phase < rule$min_year_phase)
    why <- c(why, sprintf("distinct year fractions %d < 12", rec$distinct_year_phase))
  fit <- NULL
  if (!length(why)) {
    d <- data.frame(event = risk$event, age_cap = risk$age_cap, year_phase = risk$year_phase)
    r <- capture(mgcv::gam(event ~ s(age_cap, bs = "cr", k = 4) +
      s(year_phase, bs = "cc", k = 6), family = binomial(link = "logit"), data = d,
      method = "REML", select = TRUE, knots = list(year_phase = c(0, 1))))
    m <- r$value
    if (!is.null(m)) {
      X <- model.matrix(m)
      if (qr(X)$rank < ncol(X)) r$error <- sprintf("design rank %d < %d", qr(X)$rank, ncol(X))
      else if (!isTRUE(m$converged)) r$error <- "inner iteration did not converge"
      else if (!all(is.finite(coef(m)))) r$error <- "non-finite coefficients"
    }
    rec$warnings <- paste(r$warnings, collapse = "; ")
    if (is.na(r$error)) {
      fit <- m
      labels <- vapply(m$smooth, function(s) s$label, "")
      edf <- vapply(m$smooth, function(s) sum(m$edf[s$first.para:s$last.para]), 1)
      rec$edf_age <- edf[labels == "s(age_cap)"]
      rec$edf_year_phase <- edf[labels == "s(year_phase)"]
      rec$sp <- paste(sprintf("%s=%.4g", names(m$sp), m$sp), collapse = " ")
    } else {
      rec$error <- r$error
      why <- paste("fit failed:", r$error)
    }
  }
  rec$kind <- if (!is.null(fit)) "gam" else "constant"
  rec$reason <- paste(why, collapse = "; ")
  list(record = rec, fit = fit)
}
# Hazards h[d, week] (rows: ages 1..D; columns: weeks), raw and clipped to [1e-6, 1 - 1e-6].
pp_hazard_grid <- function(h, weeks, D = pp_settings$D) {
  if (h$record$kind == "gam") {
    nd <- expand.grid(age_cap = pmin(seq_len(D), D), week = seq_along(weeks))
    nd$year_phase <- year_phase(weeks[nd$week])
    raw <- matrix(as.numeric(predict(h$fit, nd, type = "response")), D, length(weeks))
  } else raw <- matrix(h$record$constant, D, length(weeks))
  c0 <- pp_settings$hazard_clip
  list(raw = raw, h = pmin(pmax(raw, c0), 1 - c0))
}
# COVID's duration prior (no calendar term): rise and tail durations of the waves of the European
# countries (ERVISS; series of at least 100 weeks with a median of at least 5 admissions, judged on
# the whole snapshot), each series truncated at cutoff - 182 days before its peaks are found
# (find_peaks, the segment-free series); a rise from the midpoint between two peaks to the peak, a
# tail from week 5 after a peak to the next midpoint. Hazard of duration d: events / at risk pooled
# over d - 1 .. d + 1, clipped to [.01, .6], .3 with fewer than 3 at risk.
pp_duration_hazard <- function(durations, D = pp_settings$D) {
  durations <- pmin(pmax(durations[is.finite(durations)], 1), D)
  e <- tabulate(durations, D); r <- rev(cumsum(rev(e)))
  vapply(seq_len(D), function(d) {
    k <- max(1, d - 1):min(D, d + 1)
    if (sum(r[k]) < 3) .3 else min(max(sum(e[k]) / sum(r[k]), .01), .6)
  }, numeric(1))
}
pp_covid_prior <- function(pp, cutoff, D = pp_settings$D) {
  known <- cutoff - pp_settings$erviss_known
  pk_ <- pp_settings$peak
  w <- do.call(rbind, lapply(pp$erviss, function(s) {
    s <- s[s$target_end_date <= known & is.finite(s$value), ]
    if (nrow(s) < 3L) return(NULL)
    pk <- s$target_end_date[find_peaks(s$value, pk_$halfwidth, pk_$drop, pk_$min_prominence,
      pk_$smooth, scale = max(s$value))$row]
    if (length(pk) < 2L) return(NULL)
    mid <- as.Date(floor((as.numeric(pk[-1]) + as.numeric(pk[-length(pk)])) / 2),
      origin = "1970-01-01")
    data.frame(rise = round(as.numeric(pk[-1] - mid) / 7),
      tail = round(as.numeric(mid - pk[-length(pk)]) / 7) - 4)
  }))
  list(hr = pp_duration_hazard(w$rise, D), ht = pp_duration_hazard(w$tail, D),
    record = data.frame(kind = "duration prior", waves = if (is.null(w)) 0L else nrow(w),
      known_until = known))
}
# The transitions of one disease at a cutoff: a function of grid weeks giving HR and HT (D x
# weeks), a record, and whether they are available. own_risk: the pooled states' hospital episode
# rows. Flu: the hazard GAMs on the hospital and NREVSS episodes; COVID: the duration prior.
pp_transitions <- function(pp, disease, cutoff, own_risk, D = pp_settings$D) {
  if (disease == "flu") {
    ex <- pp_nrevss_risk(pp, cutoff, D)
    risk <- if (is.null(ex$risk)) own_risk else rbind(own_risk, ex$risk[names(own_risk)])
    haz <- lapply(c(rise = "rise", tail = "tail"), function(k)
      pp_hazard_fit(risk[risk$phase == k, ]))
    rec <- data.frame(kind = paste(haz$rise$record$kind, haz$tail$record$kind, sep = "/"),
      events_rise = haz$rise$record$events, events_tail = haz$tail$record$events,
      extra_pairs = if (is.null(ex$pairs)) 0L else sum(ex$pairs$admitted), waves = NA_integer_,
      edf_rise = paste(round(c(haz$rise$record$edf_age, haz$rise$record$edf_year_phase), 3),
        collapse = "/"), edf_tail = paste(round(c(haz$tail$record$edf_age,
        haz$tail$record$edf_year_phase), 3), collapse = "/"),
      sp = paste(haz$rise$record$sp, haz$tail$record$sp, sep = " | "),
      warnings = paste(haz$rise$record$warnings, haz$tail$record$warnings, sep = " | "),
      reason = paste(haz$rise$record$reason, haz$tail$record$reason, sep = " | "))
    available <- all(vapply(haz, function(h) h$record$kind != "unavailable", TRUE))
    grid <- function(weeks) list(HR = pp_hazard_grid(haz$rise, weeks, D)$h,
      HT = pp_hazard_grid(haz$tail, weeks, D)$h)
  } else {
    pr <- pp_covid_prior(pp, cutoff, D)
    rec <- data.frame(kind = pr$record$kind, events_rise = NA_integer_, events_tail = NA_integer_,
      extra_pairs = NA_integer_, waves = pr$record$waves, edf_rise = "", edf_tail = "", sp = "",
      warnings = "", reason = "")
    available <- pr$record$waves > 0
    grid <- function(weeks) list(HR = matrix(pr$hr, D, length(weeks)),
      HT = matrix(pr$ht, D, length(weeks)))
  }
  list(grid = grid, record = rec, available = available)
}

# PHASE PROBABILITIES: THE CHAIN ---------------------------------------------------------------
# States 1..D rise ages, D + 1..D + 4 overshoot weeks, D + 5..2D + 4 tail ages. step(a, hr, ht):
# the distribution one week later given a and the exit hazards of the current week (an exit at
# week w changes the phase at w + 1); age D stays at D and can still exit; back(): its transpose.
pp_index <- function(D) list(R = seq_len(D), O = D + 1:4, T = D + 4L + seq_len(D), S = 2L * D + 4L)
pp_step <- function(a, hr, ht, D) {
  i <- pp_index(D); out <- numeric(i$S)
  out[i$R[1]] <- sum(a[i$T] * ht)
  if (D > 1L) out[i$R[-1]] <- a[i$R[-D]] * (1 - hr[-D])
  out[i$R[D]] <- out[i$R[D]] + a[i$R[D]] * (1 - hr[D])
  out[i$O[1]] <- sum(a[i$R] * hr)
  out[i$O[2:4]] <- a[i$O[1:3]]
  out[i$T[1]] <- out[i$T[1]] + a[i$O[4]]
  if (D > 1L) out[i$T[-1]] <- out[i$T[-1]] + a[i$T[-D]] * (1 - ht[-D])
  out[i$T[D]] <- out[i$T[D]] + a[i$T[D]] * (1 - ht[D])
  out
}
pp_back <- function(v, hr, ht, D) {
  i <- pp_index(D); out <- numeric(i$S)
  nextR <- if (D > 1L) c(i$R[-1], i$R[D]) else i$R[D]
  nextT <- if (D > 1L) c(i$T[-1], i$T[D]) else i$T[D]
  out[i$R] <- (1 - hr) * v[nextR] + hr * v[i$O[1]]
  out[i$O[1:3]] <- v[i$O[2:4]]
  out[i$O[4]] <- v[i$T[1]]
  out[i$T] <- (1 - ht) * v[nextT] + ht * v[i$R[1]]
  out
}
pp_transition_matrix <- function(hr, ht, D) {
  S <- pp_index(D)$S
  t(vapply(seq_len(S), function(z) pp_step(replace(numeric(S), z, 1), hr, ht, D), numeric(S)))
}
# Initial distribution of a segment: 1/3 per coarse phase, uniform over its internal states.
pp_initial <- function(D) {
  i <- pp_index(D); a <- numeric(i$S)
  a[i$R] <- 1 / (3 * D); a[i$O] <- 1 / 12; a[i$T] <- 1 / (3 * D)
  a
}
pp_coarse <- function(A, D) {
  i <- pp_index(D)
  if (is.null(dim(A))) A <- matrix(A, 1)
  cbind(rise = rowSums(A[, i$R, drop = FALSE]), over = rowSums(A[, i$O, drop = FALSE]),
    tail = rowSums(A[, i$T, drop = FALSE]))
}
# Forward filter and backward smoother within each contiguous segment of s (reset at every gap),
# with the weekly likelihoods L3 (weeks x 3 coarse phases) and hazards HR, HT (D x weeks, the
# hazards of the week the transition leaves); then `ahead` predictive weeks after the last row.
pp_filter <- function(s, L3, HR, HT, D = pp_settings$D, ahead = 6L, HR_ahead = NULL,
                      HT_ahead = NULL) {
  i <- pp_index(D); n <- nrow(s)
  Lf <- cbind(matrix(L3[, 1], n, D), matrix(L3[, 2], n, 4), matrix(L3[, 3], n, D))
  alpha <- matrix(NA_real_, n, i$S); smooth <- alpha
  for (g in unique(s$segment)) {
    k <- which(s$segment == g)
    a <- pp_initial(D) * Lf[k[1], ]
    alpha[k[1], ] <- a / sum(a)
    if (length(k) > 1L) for (j in 2:length(k)) {
      a <- pp_step(alpha[k[j - 1], ], HR[, k[j - 1]], HT[, k[j - 1]], D) * Lf[k[j], ]
      alpha[k[j], ] <- a / sum(a)
    }
    beta <- rep(1, i$S)
    smooth[k[length(k)], ] <- alpha[k[length(k)], ]
    if (length(k) > 1L) for (j in (length(k) - 1L):1) {
      b <- pp_back(Lf[k[j + 1], ] * beta, HR[, k[j]], HT[, k[j]], D)
      beta <- b / sum(b)
      sm <- alpha[k[j], ] * beta
      smooth[k[j], ] <- sm / sum(sm)
    }
  }
  pred <- matrix(NA_real_, ahead, i$S)
  a <- alpha[n, ]
  hr <- HR[, n]; ht <- HT[, n]
  for (m in seq_len(ahead)) {
    a <- pp_step(a, hr, ht, D)
    pred[m, ] <- a
    if (m < ahead) { hr <- HR_ahead[, m]; ht <- HT_ahead[, m] }
  }
  list(filtered = pp_coarse(alpha, D), smoothed = pp_coarse(smooth, D),
    ahead = pp_coarse(pred, D), last_state = alpha[n, ])
}

# PHASE PROBABILITIES: OBSERVATIONS ----------------------------------------------------------------
# NSSP as known at a cutoff (published with issue + 6 <= cutoff; the latest issue of each week)
# and its weekly log growth (consecutive weeks), lag growth and holidays of the NSSP week.
pp_nssp_known <- function(pp, disease, location, cutoff) {
  n <- pp$nssp[[paste(disease, location)]]
  n <- n[n$issue + 6 <= cutoff, ]
  n <- n[order(n$week_end, n$issue), ]
  n <- n[!duplicated(n$week_end, fromLast = TRUE), ]
  n[order(n$week_end), ]
}
pp_nssp_growth <- function(pp, disease, location, cutoff) {
  n <- pp_nssp_known(pp, disease, location, cutoff)
  if (nrow(n) < 2L) return(data.frame(week = as.Date(character()), ng = numeric(),
    ng_lag = numeric(), n_hol_break = numeric(), n_hol_post = numeric()))
  lv <- log(n$value)
  consecutive <- c(FALSE, diff(as.numeric(n$week_end)) == 7)
  ng <- ifelse(consecutive, c(NA, diff(lv)), NA)
  ng_lag <- c(NA, ng[-length(ng)])
  ng_lag[!consecutive] <- NA
  h <- holiday_covariates(n$week_end)
  data.frame(week = n$week_end, ng, ng_lag, n_hol_break = h$hol_break, n_hol_post = h$hol_post)
}
# The observation rows of one series s at NSSP lead L: hospital growth and lag growth
# (consecutive), holidays; the NSSP growth of week w - L, its lag and holidays; plus the weeks
# after the newest report up to the newest NSSP week + L (at most 6), informed by NSSP alone.
pp_obs_frame <- function(s, nssp, L, extend = TRUE) {
  hosp_lag <- c(NA, s$growth[-nrow(s)])
  hosp_lag[c(TRUE, diff(s$segment) != 0)] <- NA
  x <- data.frame(week = s$target_end_date, segment = s$segment, growth = s$growth,
    lag_growth = hosp_lag, hol_break = s$hol_break, hol_post = s$hol_post, extension = FALSE)
  if (extend && nrow(nssp)) {
    last <- max(s$target_end_date)
    ext <- last + 7 * seq_len(6)
    ext <- ext[ext - 7 * L <= max(nssp$week)]
    if (length(ext)) {
      h <- holiday_covariates(ext)
      x <- rbind(x, data.frame(week = ext, segment = s$segment[nrow(s)], growth = NA,
        lag_growth = c(s$growth[nrow(s)], rep(NA, length(ext) - 1L)), hol_break = h$hol_break,
        hol_post = h$hol_post, extension = TRUE))
    }
  }
  k <- match(x$week - 7 * L, nssp$week)
  x$ng <- nssp$ng[k]; x$ng_lag <- nssp$ng_lag[k]
  x$n_hol_break <- nssp$n_hol_break[k]; x$n_hol_post <- nssp$n_hol_post[k]
  x
}
# One channel's phase-mean model: growth ~ 0 + phase + holidays (least squares), phase-specific
# residual sds (floor .05); at least 20 rows per phase; an all-zero holiday column is dropped.
pp_channel_lm <- function(x, y, hol) {
  x <- x[is.finite(x[[y]]), ]
  counts <- table(factor(x$phase, pp_classes))
  if (any(counts < pp_settings$min_rows)) return(list(available = FALSE,
    reason = sprintf("rows per phase %s < %d", paste(counts, collapse = "/"),
      pp_settings$min_rows)))
  dropped <- hol[colSums(abs(x[hol])) == 0]
  fit <- lm(reformulate(c("0", "phase", setdiff(hol, dropped)), y), x)
  if (fit$rank < length(coef(fit)) || anyNA(coef(fit)))
    return(list(available = FALSE, reason = "phase-mean regression not estimable"))
  cf <- coef(fit)
  delta <- setNames(vapply(hol, function(h) if (h %in% names(cf)) cf[[h]] else 0, 1), hol)
  list(available = TRUE, kind = "lm", mu = cf[paste0("phase", pp_classes)], delta = delta,
    sd = pmax(tapply(residuals(fit), x$phase, sd), pp_settings$sd_floor)[pp_classes], hol = hol,
    rows = nrow(x), dropped = paste(dropped, collapse = " "))
}
# One channel's conditional-growth GAM: growth ~ 0 + phase + s(lag growth, cr, 4) + holidays
# (Gaussian, REML, select = TRUE); at least 20 rows with a lag per phase and 8 distinct lags.
pp_channel_gam <- function(x, y, lag, hol) {
  x <- x[is.finite(x[[y]]) & is.finite(x[[lag]]), ]
  counts <- table(factor(x$phase, pp_classes))
  if (any(counts < pp_settings$min_rows)) return(list(available = FALSE,
    reason = sprintf("rows with lag per phase %s < %d", paste(counts, collapse = "/"),
      pp_settings$min_rows)))
  if (length(unique(x[[lag]])) < pp_settings$min_lag_distinct) return(list(available = FALSE,
    reason = "too few distinct lag growths"))
  dropped <- hol[colSums(abs(x[hol])) == 0]
  d <- data.frame(g = x[[y]], phase = x$phase, lag = x[[lag]], x[setdiff(hol, dropped)])
  form <- as.formula(paste("g ~ 0 + phase + s(lag, bs = \"cr\", k = 4)",
    paste(c("", setdiff(hol, dropped)), collapse = " + ")))
  r <- capture(mgcv::gam(form, data = d, method = "REML", select = TRUE))
  m <- r$value
  if (is.null(m) || !isTRUE(m$converged) || !all(is.finite(coef(m))))
    return(list(available = FALSE, reason = paste("GAM failed:", r$error)))
  list(available = TRUE, kind = "gam", fit = m, hol = setdiff(hol, dropped),
    sd = pmax(tapply(residuals(m), x$phase, sd), pp_settings$sd_floor)[pp_classes],
    rows = nrow(x), edf_lag = sum(m$edf[grep("s\\(lag\\)", names(coef(m)))]),
    sp = paste(sprintf("%.4g", m$sp), collapse = " "), warnings = paste(r$warnings,
      collapse = "; "))
}
# Channel means (n x 3) of a channel model for the rows of x (NA where its inputs are missing).
pp_channel_means <- function(ch, x, y, lag, hol_names) {
  n <- nrow(x)
  if (ch$kind == "lm") {
    M <- matrix(ch$mu, n, 3, byrow = TRUE)
    for (h in names(ch$delta)) M <- M + ch$delta[[h]] * x[[hol_names[[h]]]]
    return(M)
  }
  M <- matrix(NA_real_, n, 3)
  ok <- is.finite(x[[lag]]) & Reduce(`&`, lapply(ch$hol, function(h)
    is.finite(x[[hol_names[[h]]]])), rep(TRUE, n))
  if (!any(ok)) return(M)
  for (k in seq_along(pp_classes)) {
    nd <- data.frame(phase = factor(pp_classes[k], pp_classes), lag = x[[lag]][ok])
    for (h in ch$hol) nd[[h]] <- x[[hol_names[[h]]]][ok]
    M[ok, k] <- as.numeric(predict(ch$fit, nd))
  }
  M
}
pp_hosp_names <- c(hol_break = "hol_break", hol_post = "hol_post")
pp_nssp_names <- c(n_hol_break = "n_hol_break", n_hol_post = "n_hol_post")
# Means and sds of both channels (n x 3 each): the GAM where its inputs exist (O2), else the
# channel's phase-mean model.
pp_channel_moments <- function(model, x) {
  mh <- pp_channel_means(model$hosp, x, "growth", "lag_growth", pp_hosp_names)
  sh <- matrix(model$hosp$sd, nrow(x), 3, byrow = TRUE)
  if (!is.null(model$hosp_gam)) {
    g <- pp_channel_means(model$hosp_gam, x, "growth", "lag_growth", pp_hosp_names)
    use <- is.finite(g[, 1])
    mh[use, ] <- g[use, ]
    sh[use, ] <- matrix(model$hosp_gam$sd, sum(use), 3, byrow = TRUE)
  }
  out <- list(mh = mh, sh = sh, mn = NULL, sn = NULL)
  if (!is.null(model$nssp)) {
    xn <- x
    xn$n_hol_break[!is.finite(xn$n_hol_break)] <- 0
    xn$n_hol_post[!is.finite(xn$n_hol_post)] <- 0
    mn <- pp_channel_means(model$nssp, xn, "ng", "ng_lag", pp_nssp_names)
    sn <- matrix(model$nssp$sd, nrow(x), 3, byrow = TRUE)
    if (!is.null(model$nssp_gam)) {
      g <- pp_channel_means(model$nssp_gam, xn, "ng", "ng_lag", pp_nssp_names)
      use <- is.finite(g[, 1])
      mn[use, ] <- g[use, ]
      sn[use, ] <- matrix(model$nssp_gam$sd, sum(use), 3, byrow = TRUE)
    }
    out$mn <- mn; out$sn <- sn
  }
  out
}
# Standardized residuals of the labelled rows at their own phase (for the correlation).
pp_standardized <- function(model, lab) {
  m <- pp_channel_moments(model, lab)
  k <- match(as.character(lab$phase), pp_classes)
  idx <- cbind(seq_len(nrow(lab)), k)
  list(h = (lab$growth - m$mh[idx]) / m$sh[idx],
    n = if (is.null(m$mn)) rep(NA_real_, nrow(lab)) else (lab$ng - m$mn[idx]) / m$sn[idx])
}
# The record of an observation model (the same columns whether or not it could be fitted).
pp_obs_record <- function(variant, reason = "") data.frame(variant, hospital = FALSE,
  hospital_reason = reason, nssp = FALSE, nssp_reason = "", hospital_gam = FALSE,
  nssp_gam = FALSE, rho_raw = NA_real_, rho = 0, rho_rows = 0L, gam_note = "")
# Fit an observation model on the labelled weeks: "O1" hospital and NSSP phase-mean models; "O2"
# their conditional-growth GAMs on top (each falling back to its phase-mean model). Correlation of
# the paired standardized residuals, shrunk by n / (n + 50), capped at |rho| <= .9 (0 with fewer
# than 3 pairs).
pp_obs_fit <- function(lab, variant) {
  hosp <- pp_channel_lm(lab, "growth", c("hol_break", "hol_post"))
  rec <- pp_obs_record(variant, if (hosp$available) "" else hosp$reason)
  rec$hospital <- hosp$available
  if (!hosp$available) return(list(available = FALSE, record = rec))
  model <- list(variant = variant, hosp = hosp, nssp = NULL, rho = 0)
  zero_na <- function(x) transform(x, n_hol_break = ifelse(is.finite(n_hol_break), n_hol_break, 0),
    n_hol_post = ifelse(is.finite(n_hol_post), n_hol_post, 0))
  if (variant %in% c("O1", "O2")) {
    nssp <- pp_channel_lm(zero_na(lab), "ng", c("n_hol_break", "n_hol_post"))
    rec$nssp <- nssp$available
    rec$nssp_reason <- if (nssp$available) "" else nssp$reason
    if (nssp$available) model$nssp <- nssp
  }
  if (variant == "O2") {
    hg <- pp_channel_gam(lab, "growth", "lag_growth", c("hol_break", "hol_post"))
    rec$hospital_gam <- hg$available
    if (hg$available) model$hosp_gam <- hg else rec$gam_note <- paste("hospital:", hg$reason)
    if (!is.null(model$nssp)) {
      ng <- pp_channel_gam(zero_na(lab), "ng", "ng_lag", c("n_hol_break", "n_hol_post"))
      rec$nssp_gam <- ng$available
      if (ng$available) model$nssp_gam <- ng else rec$gam_note <- paste(rec$gam_note,
        "NSSP:", ng$reason)
    }
  }
  if (!is.null(model$nssp)) {
    z <- pp_standardized(model, lab)
    both <- is.finite(z$h) & is.finite(z$n)
    if (sum(both) >= 3L) {
      raw <- cor(z$h[both], z$n[both])
      rec$rho_raw <- raw
      rec$rho_rows <- sum(both)
      model$rho <- max(min(raw * sum(both) / (sum(both) + pp_settings$rho_prior),
        pp_settings$rho_cap), -pp_settings$rho_cap)
      rec$rho <- model$rho
    }
  }
  list(available = TRUE, model = model, record = rec)
}
# Log-likelihood (n x 3) of the observations of x under each phase: the joint Gaussian when both
# channels are observed, the marginal of the observed one, 0 (likelihood 1) when none is.
pp_loglik <- function(model, x) {
  m <- pp_channel_moments(model, x)
  n <- nrow(x)
  gh <- matrix(x$growth, n, 3); oh <- is.finite(x$growth)
  ll <- matrix(0, n, 3)
  zh <- (gh - m$mh) / m$sh
  if (is.null(m$mn)) {
    ll[oh, ] <- dnorm(gh[oh, , drop = FALSE], m$mh[oh, , drop = FALSE], m$sh[oh, , drop = FALSE],
      log = TRUE)
    return(ll)
  }
  gn <- matrix(x$ng, n, 3); on <- is.finite(x$ng) & is.finite(m$mn[, 1])
  zn <- (gn - m$mn) / m$sn
  rho <- model$rho
  both <- oh & on; honly <- oh & !on; nonly <- on & !oh
  if (any(both)) ll[both, ] <- -log(2 * pi) - log(m$sh[both, , drop = FALSE]) -
    log(m$sn[both, , drop = FALSE]) - .5 * log(1 - rho^2) -
    (zh[both, , drop = FALSE]^2 - 2 * rho * zh[both, , drop = FALSE] * zn[both, , drop = FALSE] +
      zn[both, , drop = FALSE]^2) / (2 * (1 - rho^2))
  if (any(honly)) ll[honly, ] <- dnorm(gh[honly, , drop = FALSE], m$mh[honly, , drop = FALSE],
    m$sh[honly, , drop = FALSE], log = TRUE)
  if (any(nonly)) ll[nonly, ] <- dnorm(gn[nonly, , drop = FALSE], m$mn[nonly, , drop = FALSE],
    m$sn[nonly, , drop = FALSE], log = TRUE)
  ll
}
# The hospital growth's moments per phase given, where observed, the same week's NSSP growth
# (the conditional Gaussian): the predictive density of the lead-selection score.
pp_hosp_conditional <- function(model, x) {
  m <- pp_channel_moments(model, x)
  mu <- m$mh; sd <- m$sh
  if (!is.null(m$mn)) {
    on <- is.finite(x$ng) & is.finite(m$mn[, 1])
    if (any(on)) {
      zn <- (matrix(x$ng, nrow(x), 3)[on, , drop = FALSE] - m$mn[on, , drop = FALSE]) /
        m$sn[on, , drop = FALSE]
      mu[on, ] <- mu[on, , drop = FALSE] + model$rho * sd[on, , drop = FALSE] * zn
      sd[on, ] <- sd[on, , drop = FALSE] * sqrt(1 - model$rho^2)
    }
  }
  list(mu = mu, sd = sd)
}

# PHASE PROBABILITIES: ONE DISEASE AND CUTOFF ----------------------------------------------------
# Episodes, transitions and each observation variant ("O1_L<lead>", "O2_L<lead>") at a cutoff.
# Per variant the probabilities of the 80-week grid (the pooled states' weeks up to 80 weeks back,
# "smoothed"; the 6 weeks after the newest hospital report: "nssp" where shifted NSSP informs
# them, else "ahead"), the predictive inputs of the next hospital week, and the records.
# Unavailable transitions or observation model: (1/3, 1/3, 1/3).
pp_hidden_at <- function(pp, disease, cutoff, variants, D = pp_settings$D,
                         ahead = pp_settings$ahead, window_weeks = pp_settings$window_weeks,
                         keep = FALSE) {
  shift <- pp_hidden_shift(pp, disease, cutoff)
  ser <- lapply(setNames(nm = pp$states), function(l) pp_series(pp, disease, l, cutoff, shift))
  pairs <- lapply(ser, function(s) pp_pairs(pp_peaks(s, cutoff), cutoff))
  ep <- lapply(pp$states, function(l) pp_episode_rows(ser[[l]], pairs[[l]], l, D))
  labels <- do.call(rbind, lapply(ep, `[[`, "labels"))
  risk <- do.call(rbind, lapply(ep, `[[`, "risk"))
  if (is.null(risk)) risk <- data.frame(location = character(), pair = integer(),
    phase = character(), week = as.Date(character()), age = integer(), age_cap = integer(),
    year_phase = numeric(), event = numeric())
  tr <- pp_transitions(pp, disease, cutoff, risk, D)
  nssp <- lapply(setNames(nm = pp$states), function(l) pp_nssp_growth(pp, disease, l, cutoff))
  pair_log <- do.call(rbind, lapply(pp$states, function(l) if (nrow(pairs[[l]]))
    cbind(location = l, pairs[[l]])))
  out <- list(); scores <- list(); records <- list(); models <- list()
  for (v in variants) {
    L <- as.integer(sub(".*_L", "", v))
    base_variant <- sub("_L.*", "", v)
    frames <- lapply(setNames(nm = pp$states), function(l) pp_obs_frame(ser[[l]], nssp[[l]], L))
    lab <- NULL
    if (!is.null(labels)) {
      lab <- do.call(rbind, lapply(pp$states, function(l) {
        x <- labels[labels$location == l, c("location", "week", "phase")]
        f <- frames[[l]]
        cbind(x, f[match(x$week, f$week), c("growth", "lag_growth", "hol_break", "hol_post",
          "ng", "ng_lag", "n_hol_break", "n_hol_post")])
      }))
    }
    om <- if (is.null(lab)) list(available = FALSE, record = pp_obs_record(base_variant,
      "no admitted pair")) else pp_obs_fit(lab, base_variant)
    available <- isTRUE(om$available) && tr$available
    rec <- data.frame(disease, cutoff, variant = v, lead = L, transition_kind = tr$record$kind[1],
      available, labelled_weeks = if (is.null(lab)) 0L else nrow(lab),
      pairs_found = if (is.null(pair_log)) 0L else nrow(pair_log),
      pairs_admitted = if (is.null(pair_log)) 0L else sum(pair_log$admitted))
    records[[v]] <- cbind(rec, om$record[setdiff(names(om$record), "variant")])
    if (keep && isTRUE(om$available)) models[[v]] <- om$model
    probs <- list(); sc <- list()
    for (l in pp$states) {
      s <- ser[[l]]; f <- frames[[l]]
      n_h <- nrow(s); last <- s$target_end_date[n_h]
      grid_after <- last + 7 * seq_len(ahead)
      if (!available) {
        keep_w <- s$target_end_date[s$target_end_date >= cutoff - 7 * window_weeks]
        third <- matrix(1 / 3, length(keep_w) + ahead, 3, dimnames = list(NULL, pp_classes))
        probs[[l]] <- data.frame(disease, location = l, cutoff, week = c(keep_w, grid_after),
          kind = rep(c("smoothed", "ahead"), c(length(keep_w), ahead)), third)
        next
      }
      ll <- pp_loglik(om$model, f)
      L3 <- exp(ll - apply(ll, 1, max))
      n_f <- nrow(f); n_ext <- sum(f$extension); n_ahead <- ahead - n_ext
      weeks_all <- c(f$week, if (n_ahead > 0) f$week[n_f] + 7 * seq_len(n_ahead))
      g <- tr$grid(weeks_all)
      fs <- data.frame(target_end_date = f$week, segment = f$segment)
      res <- suppressWarnings(pp_filter(fs, L3, g$HR[, seq_len(n_f), drop = FALSE],
        g$HT[, seq_len(n_f), drop = FALSE], D, n_ahead,
        g$HR[, n_f + seq_len(max(n_ahead, 0)), drop = FALSE],
        g$HT[, n_f + seq_len(max(n_ahead, 0)), drop = FALSE]))
      k <- which(!f$extension & f$week >= cutoff - 7 * window_weeks)
      e <- which(f$extension)
      probs[[l]] <- rbind(
        data.frame(disease, location = l, cutoff, week = f$week[k], kind = "smoothed",
          res$smoothed[k, , drop = FALSE]),
        if (length(e)) data.frame(disease, location = l, cutoff, week = f$week[e], kind = "nssp",
          res$smoothed[e, , drop = FALSE]),
        if (n_ahead > 0) data.frame(disease, location = l, cutoff,
          week = f$week[n_f] + 7 * seq_len(n_ahead), kind = "ahead", res$ahead))
      nxt <- last + 7
      j <- match(nxt, f$week)
      pi_next <- if (!is.na(j)) res$filtered[j, ] else res$ahead[1, ]
      xr <- if (!is.na(j)) f[j, ] else {
        h <- holiday_covariates(nxt)
        data.frame(week = nxt, segment = NA, growth = NA, lag_growth = s$growth[n_h],
          hol_break = h$hol_break, hol_post = h$hol_post, extension = TRUE, ng = NA, ng_lag = NA,
          n_hol_break = NA, n_hol_post = NA)
      }
      hm <- pp_hosp_conditional(om$model, xr)
      sc[[l]] <- data.frame(disease, location = l, cutoff, variant = v, last_week = last,
        last_growth = s$growth[n_h], next_week = nxt, nssp_informed = !is.na(j),
        p_rise = pi_next[1], p_over = pi_next[2], p_tail = pi_next[3], mu_rise = hm$mu[1, 1],
        mu_over = hm$mu[1, 2], mu_tail = hm$mu[1, 3], sd_rise = hm$sd[1, 1],
        sd_over = hm$sd[1, 2], sd_tail = hm$sd[1, 3])
    }
    out[[v]] <- do.call(rbind, probs)
    rownames(out[[v]]) <- NULL
    scores[[v]] <- do.call(rbind, sc)
  }
  list(probabilities = out, scores = do.call(rbind, scores), records = do.call(rbind, records),
    transition = cbind(disease, cutoff, tr$record), pairs = pair_log,
    models = if (keep) models)
}

# PHASE PROBABILITIES: LEAD SELECTION (COVID) -------------------------------------------------------
# Blocks of 13 weeks of cutoffs from the first cutoff of the base. The score of a cutoff's next
# hospital week, log sum_k pi_k N(G; mu_k, sd_k), counts once that week's growth is reported (its
# growth as the newest report of a later cutoff, revision-shifted). The lead of block j: the largest
# summed score over the weeks reported before block j begins (ties: the smallest lead); lead 0
# while a lead has no score.
pp_block_of <- function(cutoff, first, weeks = pp_settings$block_weeks)
  as.integer(as.numeric(cutoff - first) %/% (7 * weeks))
pp_predictive_scores <- function(sc) {
  realized <- unique(sc[c("disease", "location", "last_week", "last_growth", "cutoff")])
  realized <- realized[order(realized$cutoff), ]
  realized <- realized[!duplicated(realized[c("disease", "location", "last_week")]), ]
  names(realized)[names(realized) == "cutoff"] <- "reported_at"
  k <- match(paste(sc$disease, sc$location, sc$next_week), paste(realized$disease,
    realized$location, realized$last_week))
  sc$realized <- realized$last_growth[k]
  sc$reported_at <- realized$reported_at[k]
  P <- as.matrix(sc[c("p_rise", "p_over", "p_tail")])
  dens <- P[, 1] * dnorm(sc$realized, sc$mu_rise, sc$sd_rise) + P[, 2] * dnorm(sc$realized,
    sc$mu_over, sc$sd_over) + P[, 3] * dnorm(sc$realized, sc$mu_tail, sc$sd_tail)
  sc$log_score <- log(pmax(dens, 1e-300))
  sc
}
pp_select_leads <- function(sc, family, first, cutoffs, leads = 0:2) {
  x <- sc[startsWith(sc$variant, paste0(family, "_L")) & is.finite(sc$log_score), ]
  blocks <- sort(unique(pp_block_of(cutoffs, first)))
  do.call(rbind, lapply(unique(sc$disease), function(d) do.call(rbind, lapply(blocks,
    function(b) {
      start <- first + 7 * pp_settings$block_weeks * b
      y <- x[x$disease == d & x$reported_at < start, ]
      tot <- vapply(leads, function(L) sum(y$log_score[y$variant == paste0(family, "_L", L)]), 1)
      n <- vapply(leads, function(L) sum(y$variant == paste0(family, "_L", L)), 1)
      lead <- if (!nrow(y) || any(n == 0)) 0L else leads[which.max(tot)]
      data.frame(disease = d, family, block = b, block_start = start, lead, scored_weeks = min(n),
        score_L0 = tot[1], score_L1 = tot[2], score_L2 = tot[3])
    }))))
}

# PHASE PROBABILITIES: OUTSIDE SERIES ----------------------------------------------------------------
# NREVSS (WHO/NREVSS clinical-laboratory files of Delphi's flu auxiliary data, external/delphi/
# aux_data/flusion_data): positives / specimens per pooled state and MMWR week (Saturday week ends),
# New York City added to New York by counts (New Jersey has no series before October 2023).
pp_nrevss_series <- function(dir, states) {
  num <- function(x) suppressWarnings(as.numeric(x))
  pooled_names <- c("06" = "California", "48" = "Texas", "36" = "New York",
    "42" = "Pennsylvania", "17" = "Illinois", "39" = "Ohio", "13" = "Georgia",
    "37" = "North Carolina", "26" = "Michigan", "34" = "New Jersey")
  stopifnot(setequal(names(pooled_names), states))
  old <- read.csv(file.path(dir, "WHO_NREVSS_Combined_prior_to_2015_16_State.csv"), skip = 1,
    check.names = FALSE)
  old$positives <- num(old[["TOTAL SPECIMENS"]]) * num(old[["PERCENT POSITIVE"]]) / 100
  new <- read.csv(file.path(dir, "WHO_NREVSS_Clinical_Labs_State.csv"), skip = 1,
    check.names = FALSE)
  new$positives <- num(new[["TOTAL A"]]) + num(new[["TOTAL B"]])
  lab <- rbind(old[c("REGION", "YEAR", "WEEK", "TOTAL SPECIMENS", "positives")],
    new[c("REGION", "YEAR", "WEEK", "TOTAL SPECIMENS", "positives")])
  names(lab) <- c("region", "year", "week", "specimens", "positives")
  lab$specimens <- num(lab$specimens)
  lab$region[lab$region == "New York City"] <- "New York"
  lab <- lab[lab$region %in% pooled_names & is.finite(lab$specimens) &
    is.finite(lab$positives), ]
  lab <- aggregate(cbind(specimens, positives) ~ region + year + week, lab, sum)
  lab <- lab[lab$specimens > 0, ]
  jan4 <- as.Date(paste0(lab$year, "-01-04"))
  lab$target_end_date <- jan4 - as.integer(format(jan4, "%w")) + 7 * (lab$week - 1) + 6
  lab$value <- lab$positives / lab$specimens
  lapply(split(lab[c("target_end_date", "value")], lab$region), function(x)
    x[order(x$target_end_date), ])
}
# ERVISS (ECDC, weekly SARS-CoV-2 hospital admissions of the EU/EEA countries; a snapshot): one
# series per country with at least 100 weeks and a median of at least 5 admissions (ISO week ends).
pp_erviss_series <- function(file) {
  eu <- as.data.frame(arrow::read_parquet(file))
  y <- as.integer(substr(eu$yearweek, 1, 4)); w <- as.integer(substr(eu$yearweek, 7, 8))
  jan4 <- as.Date(paste0(y, "-01-04"))
  eu$target_end_date <- jan4 - (as.integer(format(jan4, "%u")) - 1) + 7 * (w - 1) + 6
  eu <- eu[order(eu$country, eu$target_end_date), ]
  eu$value <- as.numeric(eu$value)
  Filter(function(s) nrow(s) >= pp_settings$erviss_min_weeks &&
    median(s$value) >= pp_settings$erviss_min_median,
    split(eu[c("target_end_date", "value")], eu$country))
}
