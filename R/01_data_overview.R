# Purpose: understand the two data files and how reported outcomes are revised,
# and prepare the data for scripts 02 and 04.
# Run from the project root after R/00_base_forecasts.R: Rscript R/01_data_overview.R
# install.packages(c("nanoparquet", "pak", "epidatr", "epidatasets"))
# pak::pkg_install("cmu-delphi/epiprocess@main")
.libPaths(c("agent/tmp/multiqt-lib", .libPaths()))  # Same local package folder as script 02.
library(epiprocess)
library(dplyr)  # glimpse() used in this overview.

# 0. THE DATA FILES ------------------------------------------------------------
# R/00_base_forecasts.R downloads both tables once and stores them as Parquet in
# the project root; this script only reads them.
forecasts_path <- "cmu_timeseries_hub_forecasts.parquet"
vintages_path <- "nhsn_vintages_epidata.parquet"
if (!all(file.exists(c(forecasts_path, vintages_path))))
  stop("Run Rscript R/00_base_forecasts.R first: it downloads the data files.")

# 1. FORECAST TABLE: WHAT WAS PREDICTED? ---------------------------------------
# Each row contains ONE quantile of a submitted forecast. For a fixed disease,
# place, reference date and horizon, the 23 rows describe 23 quantile thresholds
# for the SAME weekly outcome. They are not 23 different observed counts.
#
# Column            Meaning and examples
# ---------------------------------------------------------------------------
# disease           Which illness is being predicted: "covid" or "flu".
#
# reference_date    Saturday ending the forecast's reference week. This is the
#                   date from which horizon is measured. It labels the forecast
#                   week.
#
# target            The quantity being predicted:
#                     "wk inc covid hosp" = weekly new COVID hospital admissions;
#                     "wk inc flu hosp"   = weekly new flu hospital admissions.
#
# horizon           How far the TARGET WEEK is from the REFERENCE WEEK:
#                     -1 = the previous week;
#                      0 = the reference week itself;
#                      1 = the following week (and similarly for 2 and 3).
#                   The relationship is:
#                     target_end_date = reference_date + 7 * horizon days.
#
#                   A negative horizon is useful because last week's admissions
#                   may still be unreported or only provisionally reported.
#                   Predicting the eventual count for that already-ended week
#                   is a nowcast: a past event can still have an uncertain total.
#
# target_end_date   Saturday ending the week whose admissions are predicted.
#                   The weekly count covers Sunday through that Saturday.
#
# location          Geographic identifier, stored as text:
#                     "01" = Alabama; "06" = California; "11" = Washington DC;
#                     "72" = Puerto Rico; "US" = the national total.
#                   The two-digit codes are state/territory FIPS-style codes.
#                   Retain leading zeros: converting "01" to 1 can break joins
#                   to another table that still identifies Alabama as "01".
#
# quantile_level    Which percentile of the forecast distribution this row gives.
#
# value             The predicted hospitalization threshold at this quantile level.
#
# truth_finalized   The latest stored report of the target week (section 2),
#                   added at download for evaluation. It is the SAME outcome
#                   across the 23 quantile rows for the same disease, location
#                   and target week, and missing while the week is unreported.
#
forecasts <- nanoparquet::read_parquet(forecasts_path)
glimpse(forecasts)
forecast_unique <- lapply(forecasts, unique)
for (column in names(forecast_unique)) {
  cat(
    "\nColumn:", column, "| Unique values:",
    length(forecast_unique[[column]]), "\n"
  )
  print(forecast_unique[[column]], max = 100)
}

# 2. REPORT TABLE: WHAT WAS REPORTED, AND WHEN? --------------------------------
# A "vintage" ("=" Jahrgang/ Kohorte) is a stored version of an observation. Each row contains ONE
# reported count for ONE disease/place/target week, as recorded in ONE version.
# Multiple rows for that week describe the history of reporting and revisions.
#
# Column            Meaning and examples
# ---------------------------------------------------------------------------
# disease           "covid" or "flu", with the same meaning as in forecasts.
#
# geo_value         Lowercase geographic abbreviation used by this report table:
#                   "al" = Alabama, "ca" = California, "pr" = Puerto Rico,
#                   "us" = national total. 
#
# location          The matching forecast-table code: "al" corresponds to "01",
#                   "ca" to "06", "pr" to "72", and "us" to "US".
#                   geo_value and location are two ways of identifying a place.
#                   Some territories have geo_value but a missing location in this
#                   file; in the forecast data these territories are not included.
#
# target_end_date   End of the week whose admissions are being counted. It has
#                   the same meaning as target_end_date in the forecast table.
#
# version           Date assigned to this archived report. This tells us which
#                   reporting version a row represents. Two rows with the same
#                   target week but different version dates are different reports
#                   of that week's count, not observations from different weeks.
#
# value             Hospital admissions reported in THIS version for THIS week.
#                   It can be revised up or down in a later version. 
#
# Why horizon -1 can still require a prediction: Alabama COVID admissions for
# the week ending 2024-11-16 were reported as 80 on 2024-11-20, 81 on 2024-11-22,
# 83 on 2024-12-06, and 82 on 2025-11-26. A forecast with reference_date
# 2024-11-23 and horizon -1 concerns that already-ended week, but its eventual
# count was still subject to revision. These dates describe the archive; without
# a submission timestamp we cannot say which reports the model had actually seen.

# The report archive (downloaded by R/00_base_forecasts.R) covers target weeks from Aug 2020
# on. It begins on Nov 19, 2024 with a bulk load of older target weeks; it does
# not show when those older counts were originally published.

vintages <- nanoparquet::read_parquet(vintages_path)
glimpse(vintages)
vintage_unique <- lapply(vintages, unique)
for (column in names(vintage_unique)) {
  cat(
    "\nColumn:", column, "| Unique values:",
    length(vintage_unique[[column]]), "\n"
  )
  print(vintage_unique[[column]], max = 100)
}

# 3. Put the reports in an epiprocess archive so we can ask "what was known when?"
# epiprocess calls the observation's target date time_value; rename it here.
# It already recognizes geo_value and version. other_keys keeps the diseases
# separate, so COVID and flu reports for the same place/week are distinct series.
# In the printed archive, time_value therefore means exactly the same week as
# target_end_date did above. Creating the archive does not change any counts.
archive <- as_epi_archive(
  vintages,
  time_value = "target_end_date",
  other_keys = "disease"
)
print(archive)

# 4. Compare two snapshots. Each has at most one row per disease/place/target week.
# "Latest" selects the last stored version of each count across the whole archive.
# These are the latest values in this file, not a promise that revisions are over.
# Snapshot rows use time_value (the renamed target_end_date) and omit version
# after choosing a report. The printed as_of date is the snapshot cutoff, not
# the individual release date of every observation in that snapshot.
# archive$versions_end is the last version date stored in this archive; using it
# as the cutoff lets every recorded revision participate in choosing the latest.
latest <- epix_as_of(archive, archive$versions_end)
print(latest)


# For a historical snapshot, choose the last version ON OR BEFORE the cutoff.
# The cutoff applies to version (report date), not to the target week's date.
# Later revisions are excluded; counts with no report by then are absent.
# Change this date to inspect another point within the stored version history.
# For the horizon -1 example above (target week ending 2024-11-16), this cutoff
# selects 83 from the 2024-12-06 report. The latest snapshot selects 82 instead,
# because it can also see the later revision dated 2025-11-26.
cutoff <- as.Date("2025-01-11")
available_then <- epix_as_of(archive, cutoff)
print(available_then)

# 5. PREPARE THE DATA FOR SCRIPTS 02 AND 04 ------------------------------------
# Keep every location's full forecast chronology in raw count units, from the first
# 2023/24 round on. Seasons run from October to September and are set by the target
# week: 2023/24 tunes the methods' settings, 2024/25 and 2025/26 evaluate them.
# Peaks only label weeks in hindsight; they do not affect online feedback. A peak is
# the largest latest-reported weekly count of a season (flu) or a half season (covid:
# October to March and April to September, one for the winter and one for the summer
# wave). Only peaks whose week was forecast are kept. Boundary maxima are retained;
# ties use the earliest date. Each target week belongs to its nearest peak.
#
# The output contains:
#   settings : the settings below;
#   peaks    : one row per disease, location and (half) season with a forecast peak;
#   cases    : named disease_location_horizon lists, e.g. covid_18_0.
# Each case contains its raw forecast matrix, quantiles, online feedback per delay,
# latest-vintage truth for scoring, dates, seasons and periods, peak-relative week
# labels, the population (counts per 100k residents) and the historical counts for
# other scalings. Scripts 02/04 choose the scaling without reopening the data files.
# US remains a separate case and must not be counted as an additional state.
#
# base: the base forecasts per disease. "windowed_seasonal": R/00_base_forecasts.R,
#   Delphi's forecaster, made on Wednesdays (the generation date) from the reports as
#   of that day. "cmu_hub": the CMU-TimeSeries submissions (covid from 2024-11-23 only).
#   The online feedback of a round uses the reports as of the base's generation date
#   (windowed_seasonal) or the day before the reference date (cmu_hub).
# final_data_through: there are no report versions before 2024-11-19. Target weeks
#   ending by this date get a first report with their latest value on the Wednesday
#   after the week (target_end_date + 4), so in 2023/24 MultiQT learns from final
#   counts, as the windowed_seasonal forecasts were made from final counts.
# incomplete_weeks: mandatory hospital reporting paused from 2024-05-01 to 2024-10-31,
#   so the counts of these target weeks are incomplete (Texas covid: 643, 274, 8, 2 in
#   the weeks ending 2024-04-27 to 2024-05-18). Forecasts of them are dropped: never
#   learned from, never scored.
# feedback_delays: the delays (weeks after the first report) prepared for the
#   online feedback; feedback_delay in R/03_functions.R chooses one of them.
# fill_gaps_up_to: fill gaps of at most this many missing rounds (0: no filling),
#   see fill_missing_rounds() below; the gap from May to November 2024 stays.
prepare_settings <- list(base = c(covid = "windowed_seasonal", flu = "cmu_hub"),
                         start_date = as.Date("2023-10-14"), season_start_month = 10L,
                         tuning_seasons = "2023/24",
                         evaluation_seasons = c("2024/25", "2025/26"),
                         half_season_peaks = c(covid = TRUE, flu = FALSE),
                         final_data_through = as.Date("2024-04-27"),
                         incomplete_weeks = as.Date(c("2024-05-04", "2024-11-02")),
                         horizons = -1:3, feedback_delays = 0:4, fill_gaps_up_to = 6L)

# Season of a date (October to September, e.g. "2024/25") and its half season.
season_of <- function(date, start_month = 10L) {
  year <- as.integer(format(date, "%Y")) - (as.integer(format(date, "%m")) < start_month)
  paste0(year, "/", substr(year + 1L, 3L, 4L))
}
half_season_of <- function(date, start_month = 10L)
  ifelse((as.integer(format(date, "%m")) - start_month) %% 12L < 6L, "winter", "summer")
# A missing horizon-h forecast of round r is filled with the latest earlier forecast
# of the same target week: the horizon-(h + k) forecast issued k weeks earlier, with
# the smallest k. It uses no information after r. Only gaps of at most max_gap
# consecutive missing rounds are filled (not the flu off-season), and the longest
# horizon has no earlier forecast to fill from. filled_from is the source round
# (NA for issued forecasts).
fill_missing_rounds <- function(forecasts, max_gap) {
  forecasts$filled_from <- as.Date(NA)
  if (max_gap < 1L) return(forecasts)
  key <- paste(forecasts$disease, forecasts$location)
  filled <- list()
  for (k in unique(key)) {
    f <- forecasts[key == k, , drop = FALSE]
    for (h in sort(unique(f$horizon))) {
      rounds <- sort(unique(f$reference_date[f$horizon == h]))
      weeks <- seq(min(rounds), max(rounds), by = 7)
      missing <- weeks[!weeks %in% rounds]
      if (!length(missing)) next
      gap <- cumsum(c(1L, diff(as.numeric(missing)) != 7))
      missing <- missing[ave(seq_along(missing), gap, FUN = length) <= max_gap]
      for (i in seq_along(missing)) for (later in seq_len(max(f$horizon) - h)) {
        source <- f[f$reference_date == missing[i] - 7 * later & f$horizon == h + later, , drop = FALSE]
        if (!nrow(source)) next
        stopifnot(all(source$target_end_date == missing[i] + 7 * h))
        source$filled_from <- source$reference_date
        source$reference_date <- missing[i]
        source$horizon <- h
        filled[[length(filled) + 1L]] <- source
        break
      }
    }
  }
  rbind(forecasts, do.call(rbind, filled))
}

prepare_data <- function(forecasts, vintages, settings) {
  forecast_columns <- c("disease", "location", "horizon", "reference_date",
                        "target_end_date", "quantile_level", "value", "truth_finalized")
  selected <- as.data.frame(forecasts[forecasts$horizon %in% settings$horizons &
    forecasts$reference_date >= settings$start_date, forecast_columns])
  # Invalid submissions: a round and horizon whose forecasts are 0 at every quantile
  # and location are dropped and filled like a missing round. CMU's flu submission of
  # 2025-11-22 has such zeros at horizon -1 (its horizons 0-3 are fine).
  block <- paste(selected$disease, selected$reference_date, selected$horizon)
  all_zero <- tapply(selected$value == 0, block, all)
  invalid <- names(all_zero)[all_zero]
  if (length(invalid)) message("All-zero forecasts dropped and filled (disease, round, ",
    "horizon): ", paste(invalid, collapse = "; "))
  selected <- selected[!block %in% invalid, , drop = FALSE]
  selected <- fill_missing_rounds(selected, settings$fill_gaps_up_to)
  selected <- selected[selected$target_end_date < settings$incomplete_weeks[1] |
    selected$target_end_date > settings$incomplete_weeks[2], , drop = FALSE]
  # Reports as of each round: the base's generation date, or the day before the
  # reference date. Filled rounds: the Wednesday before the reference date.
  generation <- unique(as.data.frame(forecasts[!is.na(forecasts$generation_date),
    c("disease", "reference_date", "generation_date")]))
  stopifnot(!anyDuplicated(generation[c("disease", "reference_date")]))
  round_as_of <- function(disease, reference_date) {
    if (settings$base[[disease]] == "cmu_hub") return(reference_date - 1)
    as_of <- generation$generation_date[match(paste(disease, reference_date),
      paste(generation$disease, generation$reference_date))]
    as.Date(ifelse(is.na(as_of), reference_date - 3, as_of), origin = "1970-01-01")
  }
  reports <- as.data.frame(vintages[!is.na(vintages$location),
    c("disease", "location", "geo_value", "target_end_date", "version", "value")])
  # Final counts before the report versions: a first report of the latest value on the
  # Wednesday after the week, before every stored version.
  reports <- reports[order(reports$disease, reports$location, reports$target_end_date,
    reports$version), , drop = FALSE]
  final <- reports[!duplicated(reports[c("disease", "location", "target_end_date")],
    fromLast = TRUE) & reports$target_end_date <= settings$final_data_through, , drop = FALSE]
  final$version <- final$target_end_date + 4
  stopifnot(all(final$version < min(reports$version)))
  reports <- rbind(final, reports)
  group_key <- function(x) paste(x$disease, x$location, sep = "_")
  forecast_groups <- split(seq_len(nrow(selected)), group_key(selected))
  report_groups <- split(seq_len(nrow(reports)), group_key(reports))
  periods <- function(season) ifelse(season %in% settings$tuning_seasons, "tuning",
    ifelse(season %in% settings$evaluation_seasons, "evaluation", "other"))
  cases <- peaks <- list()
  for (key in names(forecast_groups)) {
    f <- selected[forecast_groups[[key]], , drop = FALSE]
    report_rows <- report_groups[[key]]
    if (!length(report_rows)) stop("No reports for ", key)
    r <- reports[report_rows, , drop = FALSE]
    r <- r[order(r$target_end_date, r$version), , drop = FALSE]
    stopifnot(!anyDuplicated(r[c("target_end_date", "version")]),
              length(unique(r$geo_value)) == 1L)
    disease <- r$disease[1]
    location_name <- state.name[match(toupper(r$geo_value[1]), state.abb)]
    if (r$location[1] == "11") location_name <- "District of Columbia"
    if (r$location[1] == "72") location_name <- "Puerto Rico"
    if (r$location[1] == "US") location_name <- "United States"
    stopifnot(!is.na(location_name))
    # 2019 census population (FIPS "00" is the US), for counts per 100k residents.
    census <- epidatasets::state_census
    population <- census$pop[match(if (r$location[1] == "US") "00" else r$location[1], census$fips)]
    stopifnot(is.finite(population), population > 0)
    # Choose latest vintages first; missing newest values never fall back to an
    # older finite vintage. This ordering also applies to scaling and feedback.
    latest <- r[!duplicated(r$target_end_date, fromLast = TRUE), , drop = FALSE]

    # Peaks: the largest latest report of each tuning or evaluation (half) season,
    # among the reported weeks of that (half) season; kept if its week was forecast.
    half <- settings$half_season_peaks[[disease]]
    part_of <- function(date) {
      season <- season_of(date, settings$season_start_month)
      if (half) paste(season, half_season_of(date, settings$season_start_month)) else season
    }
    latest_part <- part_of(latest$target_end_date)
    location_peaks <- list()
    for (part in unique(part_of(f$target_end_date))) {
      window <- latest[latest_part == part & is.finite(latest$value), , drop = FALSE]
      season <- season_of(window$target_end_date[1], settings$season_start_month)
      if (!nrow(window) || periods(season) == "other") next
      best <- which.max(window$value)
      if (!window$target_end_date[best] %in% f$target_end_date) next
      tied <- which(window$value == window$value[best])
      location_peaks[[part]] <- data.frame(disease, location = r$location[1],
        geo_value = r$geo_value[1], location_name, season,
        half = if (half) half_season_of(window$target_end_date[best],
          settings$season_start_month) else "season",
        period = periods(season), peak_date = window$target_end_date[best],
        peak_value = window$value[best], boundary_peak = best %in% c(1L, nrow(window)),
        tied_peak_n = length(tied), last_tied_date = max(window$target_end_date[tied]),
        window_start = min(window$target_end_date), window_end = max(window$target_end_date),
        n_window_weeks = nrow(window))
    }
    location_peaks <- do.call(rbind, location_peaks)
    if (is.null(location_peaks)) {
      warning("No forecast peak for ", key, ": no before/after windows.", call. = FALSE)
    } else {
      location_peaks <- location_peaks[order(location_peaks$peak_date), , drop = FALSE]
      rownames(location_peaks) <- NULL
      peaks[[key]] <- location_peaks
    }

    # One raw report lookup table is shared by this location's horizons.
    # The median identifies the issue/target rows; all quantiles enter base below.
    medians <- f[is.finite(f$truth_finalized) & f$quantile_level == .5, , drop = FALSE]
    issue_dates <- sort(unique(medians$reference_date))
    issue_as_of <- round_as_of(disease, issue_dates)
    target_dates <- sort(unique(medians$target_end_date))
    available <- matrix(NA_real_, length(issue_dates), length(target_dates))
    report_targets <- split(seq_len(nrow(r)), as.character(r$target_end_date))
    for (j in seq_along(target_dates)) {
      rows <- report_targets[[as.character(target_dates[j])]]
      if (!length(rows)) next
      # The latest report with version <= the round's as-of date.
      latest_index <- findInterval(as.numeric(issue_as_of), as.numeric(r$version[rows]))
      found <- latest_index > 0L
      available[found, j] <- r$value[rows[latest_index[found]]]
    }
    for (horizon in sort(unique(f$horizon))) {
      fh <- f[f$horizon == horizon, , drop = FALSE]
      z <- medians[medians$horizon == horizon, , drop = FALSE]
      z <- z[order(z$reference_date), , drop = FALSE]
      if (!nrow(z)) stop("No finite finalized truth for ", key, " horizon ", horizon)
      stopifnot(!anyDuplicated(z$reference_date), !anyDuplicated(z$target_end_date))
      reference_date <- as.Date(as.numeric(z$reference_date), origin = "1970-01-01")
      target_date <- as.Date(as.numeric(z$target_end_date), origin = "1970-01-01")
      alpha <- sort(unique(fh$quantile_level))
      finite <- fh[is.finite(fh$truth_finalized), , drop = FALSE]
      row <- match(finite$reference_date, reference_date)
      column <- match(finite$quantile_level, alpha)
      stopifnot(!anyNA(row), !anyNA(column),
                !anyDuplicated(data.frame(row, column)),
                nrow(finite) == nrow(z) * length(alpha),
                all(finite$target_end_date == target_date[row]),
                all(finite$truth_finalized == z$truth_finalized[row]))
      base <- matrix(NA_real_, nrow(z), length(alpha), dimnames = list(NULL, as.character(alpha)))
      base[cbind(row, column)] <- finite$value
      stopifnot(all(is.finite(base)))

      # Each earlier forecast is consumed once. With delay d it waits d weeks after
      # the first later issue with a finite report, and is consumed at the first
      # issue from then on, with the report available at that issue (d = 0: the
      # first report). The saved value never changes when later revisions arrive.
      issues <- match(reference_date, issue_dates)
      targets <- match(target_date, target_dates)
      feedback_by_delay <- lapply(setNames(settings$feedback_delays, settings$feedback_delays),
        function(delay) {
          observed <- rep(NA_real_, nrow(z))
          date <- rep(as.Date(NA), nrow(z))
          for (i in seq_len(nrow(z))) {
            if (i == nrow(z)) next
            later <- seq.int(i + 1L, nrow(z))
            counts <- available[issues[later], targets[i]]
            first <- which(is.finite(counts))[1L]
            if (is.na(first)) next
            due <- reference_date[later[first]] + 7 * delay
            use <- which(is.finite(counts) & reference_date[later] >= due)[1L]
            if (is.na(use)) next
            observed[i] <- counts[use]
            date[i] <- reference_date[later[use]]
          }
          used_at <- match(date, reference_date)
          list(observed = observed, date = date,
            by_week = lapply(seq_len(nrow(z)), function(i) which(!is.na(used_at) & used_at == i)))
        })

      season <- season_of(target_date, settings$season_start_month)
      nearest <- if (is.null(location_peaks)) rep(NA_integer_, length(target_date)) else
        apply(abs(outer(as.numeric(target_date), as.numeric(location_peaks$peak_date), "-")),
          1, which.min)
      relative_week <- as.numeric(target_date - location_peaks$peak_date[nearest]) / 7

      # Save the exact pre-first-forecast history, before filtering finalized
      # truth. Scripts 02/03 can use its median, maximum, or no scaling.
      scale_as_of <- min(fh$reference_date) - 1
      history <- r[r$version <= scale_as_of, , drop = FALSE]
      history <- history[!duplicated(history$target_end_date, fromLast = TRUE), , drop = FALSE]
      history <- history[history$target_end_date <= scale_as_of & is.finite(history$value), , drop = FALSE]
      case_id <- paste(key, horizon, sep = "_")
      cases[[case_id]] <- list(disease = z$disease[1], location = z$location[1],
        geo_value = r$geo_value[1], location_name = location_name, population = population,
        horizon = as.integer(horizon),
        alpha = alpha, base = base, reference_date = reference_date, target_date = target_date,
        truth_finalized = z$truth_finalized,
        # Filled rounds: the source round of the earlier forecast (see fill_missing_rounds()).
        filled = !is.na(z$filled_from), filled_from = z$filled_from,
        # Latest stored report of each target week: scoring target only, never feedback.
        truth_latest = latest$value[match(target_date, latest$target_end_date)],
        # Online feedback per delay; get_case() in R/03_functions.R selects one.
        feedback_by_delay = feedback_by_delay,
        # The base, and the reports' as-of date of each round (also for the feedback).
        base_name = settings$base[[disease]], as_of = issue_as_of[issues],
        # Season of the target week and its period: "tuning", "evaluation" or "other".
        season = season, period = periods(season),
        # This location's peaks, and per week the nearest one (hindsight labels only).
        peaks = location_peaks, peak_index = nearest,
        peak_date = location_peaks$peak_date[nearest], relative_week = relative_week,
        scale_history = history$value, scale_as_of = scale_as_of,
        scale_reference_start = if (nrow(history)) min(history$target_end_date) else as.Date(NA),
        scale_reference_end = if (nrow(history)) max(history$target_end_date) else as.Date(NA))
    }
  }
  peaks <- do.call(rbind, peaks)
  rownames(peaks) <- NULL
  list(settings = settings, peaks = peaks, cases = cases)
}

# Base forecasts per disease (prepare_settings$base); windowed_seasonal_forecasts.parquet
# is written by R/00_base_forecasts.R. Both tables have the same columns.
base_forecasts <- bind_rows(lapply(names(prepare_settings$base), function(disease) {
  name <- prepare_settings$base[[disease]]
  table <- if (name == "cmu_hub") forecasts else
    nanoparquet::read_parquet("windowed_seasonal_forecasts.parquet")
  table <- table[table$disease == disease, ]
  if (!nrow(table)) stop("No ", name, " forecasts for ", disease, ". See R/00_base_forecasts.R.")
  if (!"generation_date" %in% names(table)) table$generation_date <- as.Date(NA)
  table
}))
peak_aligned <- prepare_data(base_forecasts, vintages, prepare_settings)
dir.create("output/data", recursive = TRUE, showWarnings = FALSE)
saveRDS(peak_aligned, "output/data/peak_aligned.rds")


# 6. MISSING DATA: ALL DATES ----------------------------------------------------
# In every table: n_missing = n_absent + n_na; n_locations counts DISTINCT
# affected locations, and n_weeks counts DISTINCT affected weeks.
# Detail tables show only weeks with missing data; summaries include zeros.

# A. FORECASTS: count individual quantile values, separately for every horizon.
# Expect every observed location and horizon on every reference week between
# that disease's first and last forecasts, with all quantiles seen at that horizon.
# This describes coverage, including seasonal breaks or changes in horizons;
# an absent row does not establish that a forecast should have been submitted.
forecast_weeks <- forecasts %>%
  group_by(disease) %>%
  reframe(reference_date = seq(min(reference_date), max(reference_date), by = "week"))
forecast_expected <- forecasts %>%
  distinct(disease, location) %>%
  inner_join(
    forecasts %>% group_by(disease, horizon) %>%
      summarise(n_expected = n_distinct(quantile_level), .groups = "drop"),
    by = "disease", relationship = "many-to-many"
  )
forecast_missing_by_location <- forecast_weeks %>%
  inner_join(forecast_expected, by = "disease", relationship = "many-to-many") %>%
  left_join(
    forecasts %>% group_by(disease, horizon, reference_date, location) %>%
      summarise(
        n_present = n_distinct(quantile_level),
        n_na = sum(is.na(value)), .groups = "drop"
      ),
    by = c("disease", "horizon", "reference_date", "location")
  ) %>%
  mutate(
    n_absent = n_expected - coalesce(n_present, 0L),
    n_na = coalesce(n_na, 0L), n_missing = n_absent + n_na,
    is_missing = n_missing > 0, target_end_date = reference_date + 7 * horizon
  )
forecast_missing_summary <- forecast_missing_by_location %>%
  group_by(disease, horizon) %>%
  summarise(
    n_missing = sum(n_missing), n_absent = sum(n_absent), n_na = sum(n_na),
    n_locations = n_distinct(location[is_missing]),
    n_weeks = n_distinct(reference_date[is_missing]),
    .groups = "drop"
  )
forecast_missing_weeks <- forecast_missing_by_location %>%
  filter(n_missing > 0) %>%
  group_by(disease, horizon, reference_date, target_end_date) %>%
  summarise(
    n_missing = sum(n_missing), n_absent = sum(n_absent), n_na = sum(n_na),
    n_locations = n_distinct(location), .groups = "drop"
  )
cat("\nFORECASTS: missing quantile values by disease and horizon\n",
  "n_absent = absent quantile rows; n_na = stored NA values.\n",
  "For example, one completely absent 23-quantile forecast contributes 23 missing values.\n",
  sep = ""
)
print(forecast_missing_summary, n = Inf, width = 140)
cat("\nFORECASTS: affected reference weeks and corresponding target weeks\n")
print(forecast_missing_weeks, n = Inf, width = 140)

# B. REPORTS AND FINALIZED TRUTHS: weeks below are TARGET weeks, not issue dates.
# All reports: count NA rows across stored versions; one location/week may
# contribute several missing values. Do not invent expected revision records.
# Latest reports: one value per location/target week, including absent rows on
# the full weekly grid. Use geo_value so all 57 locations, including territories,
# are included. Finalized truths: count each forecast location/target week once,
# without multiplying by quantiles/horizons or inventing truths for absent forecasts.
latest_report_grid <- vintages %>%
  group_by(disease) %>%
  reframe(target_end_date = seq(min(target_end_date), max(target_end_date), by = "week")) %>%
  inner_join(vintages %>% distinct(disease, geo_value),
    by = "disease", relationship = "many-to-many"
  ) %>%
  left_join(
    as.data.frame(latest) %>%
      transmute(disease, geo_value,
        target_end_date = time_value,
        value, present = TRUE
      ),
    by = c("disease", "geo_value", "target_end_date")
  )
report_values <- bind_rows(
  vintages %>% transmute(
    source = "all_reports", disease, target_end_date,
    place = geo_value, n_absent = 0L, n_na = as.integer(is.na(value))
  ),
  latest_report_grid %>%
    transmute(
      source = "latest_reports", disease, target_end_date, place = geo_value,
      n_absent = as.integer(is.na(present)),
      n_na = as.integer(!is.na(present) & is.na(value))
    ),
  forecasts %>% distinct(disease, location, target_end_date, truth_finalized) %>%
    transmute(
      source = "truth_finalized", disease, target_end_date, place = location,
      n_absent = 0L, n_na = as.integer(is.na(truth_finalized))
    )
) %>% mutate(n_missing = n_absent + n_na, is_missing = n_missing > 0)
report_missing_summary <- report_values %>%
  group_by(source, disease) %>%
  summarise(
    n_missing = sum(n_missing), n_absent = sum(n_absent), n_na = sum(n_na),
    n_locations = n_distinct(place[is_missing]),
    n_weeks = n_distinct(target_end_date[is_missing]),
    .groups = "drop"
  )
report_missing_weeks <- report_values %>%
  filter(n_missing > 0) %>%
  group_by(source, disease, target_end_date) %>%
  summarise(
    n_missing = sum(n_missing), n_absent = sum(n_absent), n_na = sum(n_na),
    n_locations = n_distinct(place), .groups = "drop"
  )
cat("\nREPORTS AND FINALIZED TRUTHS: missing values by source and disease\n",
  "All reports count stored versions; latest reports and finalized truths count outcomes.\n",
  "Absence is checked only for latest reports; other sources count existing NA values.\n",
  sep = ""
)
print(report_missing_summary, n = Inf, width = 140)
cat("\nREPORTS AND FINALIZED TRUTHS: affected target weeks\n")
print(report_missing_weeks, n = Inf, width = 140)
