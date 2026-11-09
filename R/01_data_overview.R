# Purpose: understand the two data files and how reported outcomes are revised,
# and prepare the data for scripts 02 and 04.
# Run from the project root after R/00_base_forecasts.R: Rscript R/01_data_overview.R
# install.packages(c("nanoparquet", "pak", "epidatr", "epidatasets"))
# pak::pkg_install("cmu-delphi/epiprocess@main")
.libPaths(c("agent/tmp/multiqt-lib", .libPaths()))  # Same local package folder as script 02.
library(epiprocess)
library(dplyr)  # glimpse() used in this overview.
source("R/03_functions.R")  # feedback_schedule(): the delayed-feedback protocol (section 5).

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
# round of the design on. Seasons run from October to September and are set by the
# target week. Each row belongs to a period by its target week: "tuning" (R/05 chooses
# the methods' settings there, and every method warms up there), "evaluation" (scored
# by 02 and 04) or "other". Two designs:
#   "hub_first_months": the CMU hub forecasts of both diseases from 2024-11-23 (real
#     horizon -1 nowcasts, report versions throughout); tuning on the target weeks
#     through hub_tuning_through (the first six months), evaluation on the later
#     weeks of the evaluation seasons (June 2025 to September 2026).
#   "seasons": the base per disease (`base`) from 2023-10-14; tuning on the 2023/24
#     season, evaluation on the 2024/25 and 2025/26 seasons. Nothing is forecast from
#     May to November 2024, so the methods carry their 2023/24 state across that gap.
# Peaks only label weeks in hindsight; they do not affect online feedback. A peak is a
# local maximum of the latest-reported weekly counts (find_peaks() below): the largest
# value within peak_halfwidth weeks on either side, whose counts fall by at least
# peak_drop of its height on both sides before any higher week, by at least
# peak_min_prominence of the location's largest forecast-week count, on the centred
# peak_smooth-week mean (small locations are noisy). Peaks between the first and the
# last forecast target week are kept, in any period. Each target week belongs to its
# nearest peak.
#
# The output contains:
#   settings : the settings below, the feedback protocol (feedback_protocol() in
#              R/03_functions.R) and truth_version, the report version scored against;
#   peaks    : one row per disease, location and peak (find_peaks());
#   cases    : named disease_location_horizon lists, e.g. covid_18_0.
# Each case contains its raw forecast matrix, quantiles, the issue cutoff of each row
# (as_of, cutoff_source), every stored report of each row's target week (report_versions,
# for the feedback and the re-run design), the online feedback per delay with its audit
# (feedback_by_delay; all cases in one table: feedback_audit() in R/03_functions.R),
# latest-vintage truth for scoring, dates, seasons and periods, peak-relative week labels,
# the population (counts per 100k residents) and the historical counts for other
# scalings. Scripts 02/04 choose the scaling without reopening the data files.
# US remains a separate case and must not be counted as an additional state.
#
# base: the base forecasts per disease. "windowed_seasonal": R/00_base_forecasts.R,
#   Delphi's forecaster, made on Wednesdays (the generation date) from the reports as
#   of that day. "cmu_hub": the CMU-TimeSeries submissions (covid from 2024-11-23 only).
# Issue cutoff of a round (notes/notes.tex, Section 2): what is known just before its
#   forecasts are issued, the as-of date of its reports and of its feedback. The base's
#   actual generation date where known (windowed_seasonal: its generation_date; CMU hub:
#   hub_generation_dates), otherwise the hubs' Wednesday submission deadline, the
#   reference date - 3, also for a round without any submitted forecast (filled). The
#   archive dates reports but has no times: a report dated on the cutoff day counts as
#   known then (same-day convention). So the first report of a week, mostly published on
#   the Wednesday after it, is known at the issue of its horizon -1 nowcast.
# hub_generation_dates: the CMU hub rounds that Delphi's production pipeline generated
#   after the Wednesday deadline (flu's first round, the holidays), from the production
#   schedule in pipelines/<disease>_hosp_prod.R of the Delphi code pinned in
#   R/00_base_forecasts.R (the same dates as the windowed_seasonal generation dates).
# final_data_through: there are no report versions before 2024-11-19. Target weeks
#   ending by this date get a first report with their latest value on the Wednesday
#   after the week (target_end_date + 4), so in 2023/24 MultiQT learns from final
#   counts, as the windowed_seasonal forecasts were made from final counts. These
#   synthetic reports have source "synthetic_final"; archived versions "archive", and
#   the archive's first version (2024-11-19, a bulk load of all older weeks, without their
#   original publication dates) "archive_start".
# incomplete_weeks: mandatory hospital reporting paused from 2024-05-01 to 2024-10-31,
#   so the counts of these target weeks are incomplete (Texas covid: 643, 274, 8, 2 in
#   the weeks ending 2024-04-27 to 2024-05-18). Forecasts of them are dropped: never
#   learned from, never scored.
# invalid_submissions: CMU hub submissions broken by pipeline bugs or bad data, as
#   excluded in Delphi's calibration reports (location or horizon NA: all of them).
#   Covid 2024-11-23: made from a bad NHSN release corrected the next day (US median
#   21,614 at horizon 0, final count 8,104); covid 2026-04-25 horizon -1: a stale data
#   fetch (US 892, final 1,711); flu US 2024-12-14, 2024-12-21 and 2025-01-04: about half
#   the sum of the states. They are dropped and filled like missing rounds, as are
#   forecasts that are 0 at every quantile and location of a round and horizon.
# feedback_delays: the delays d prepared for the online feedback (feedback_schedule() in
#   R/03_functions.R): a forecast is learned once, at the first later round whose cutoff
#   is at least d calendar weeks after its first valid report and whose current report is
#   valid, from that report (missing rounds and invalid reports delay the update; later
#   revisions are ignored). feedback_delay in R/03_functions.R chooses one of them, and
#   feedback_revisions the re-run design, which re-learns the same rows in the same rounds
#   from their reports current at each later cutoff.
# fill_gaps_up_to: fill gaps of at most this many missing rounds (0: no filling),
#   see fill_missing_rounds() below; the gap from May to November 2024 stays.
design <- "hub_first_months"                  # "hub_first_months" or "seasons".
prepare_settings <- list(design = design,
                         base = if (design == "seasons")
                           c(covid = "windowed_seasonal", flu = "cmu_hub") else
                           c(covid = "cmu_hub", flu = "cmu_hub"),
                         start_date = as.Date(if (design == "seasons") "2023-10-14" else
                           "2024-11-23"), season_start_month = 10L,
                         tuning_seasons = "2023/24",                         # "seasons"
                         hub_tuning_through = as.Date("2025-05-31"),         # "hub_first_months"
                         evaluation_seasons = c("2024/25", "2025/26"),
                         peak_halfwidth = 8L, peak_drop = .3, peak_min_prominence = .1,
                         peak_smooth = 3L,
                         final_data_through = as.Date("2024-04-27"),
                         incomplete_weeks = as.Date(c("2024-05-04", "2024-11-02")),
                         invalid_submissions = data.frame(
                           disease = c("covid", "covid", "flu", "flu", "flu"),
                           reference_date = as.Date(c("2024-11-23", "2026-04-25",
                             "2024-12-14", "2024-12-21", "2025-01-04")),
                           location = c(NA, NA, "US", "US", "US"),
                           horizon = c(NA, -1L, NA, NA, NA)),
                         hub_generation_dates = data.frame(
                           disease = c("flu", rep(c("covid", "flu"), 3)),
                           reference_date = as.Date(c("2024-11-23", "2024-12-28",
                             "2024-12-28", "2025-01-04", "2025-01-04", "2025-12-27",
                             "2025-12-27")),
                           generation_date = as.Date(c("2024-11-21", "2024-12-26",
                             "2024-12-26", "2025-01-02", "2025-01-02", "2025-12-29",
                             "2025-12-29"))),
                         horizons = -1:3, feedback_delays = 0:4, fill_gaps_up_to = 6L)

# Season of a date (October to September, e.g. "2024/25") and its half season.
season_of <- function(date, start_month = 10L) {
  year <- as.integer(format(date, "%Y")) - (as.integer(format(date, "%m")) < start_month)
  paste0(year, "/", substr(year + 1L, 3L, 4L))
}
half_season_of <- function(date, start_month = 10L)
  ifelse((as.integer(format(date, "%m")) - start_month) %% 12L < 6L, "winter", "summer")
# Peaks of a weekly series x: local maxima of its centred `smooth`-week mean s that are
# the largest value of s within `halfwidth` weeks on either side (ties: the earliest) and
# whose prominence is at least `drop` times their height and `min_prominence` times
# `scale`. Prominence: the height above the higher of the two lowest points between the
# week and the nearest higher week on each side (or the end of the series), so a
# maximum at the end of the series has none. Each peak is placed at the largest x within
# one week. Returns the rows of x and the relative prominence.
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
  # Rounds with a submitted base forecast, broken or not (the others are filled).
  submitted <- unique(paste(selected$disease, selected$reference_date))
  # Broken CMU hub submissions (settings$invalid_submissions) are dropped and filled
  # like a missing round; a dropped first round stays missing.
  bad <- settings$invalid_submissions
  broken <- Reduce(`|`, lapply(seq_len(nrow(bad)), function(i)
    settings$base[selected$disease] == "cmu_hub" & selected$disease == bad$disease[i] &
      selected$reference_date == bad$reference_date[i] &
      (is.na(bad$location[i]) | selected$location == bad$location[i]) &
      (is.na(bad$horizon[i]) | selected$horizon == bad$horizon[i])), FALSE)
  if (any(broken)) message("Broken submissions dropped and filled: ",
    paste(unique(paste(selected$disease, selected$reference_date)[broken]), collapse = "; "))
  selected <- selected[!broken, , drop = FALSE]
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
  # Issue cutoff of each round (see the settings above), with its source: "generation_date"
  # (windowed_seasonal), "delphi_production_schedule" (hub_generation_dates),
  # "wednesday_deadline" (the hubs' deadline) or "wednesday_no_submission" (filled round).
  generation <- unique(as.data.frame(forecasts[!is.na(forecasts$generation_date),
    c("disease", "reference_date", "generation_date")]))
  generation$source <- rep("generation_date", nrow(generation))
  hub <- settings$hub_generation_dates
  hub <- hub[settings$base[hub$disease] == "cmu_hub", , drop = FALSE]
  generation <- rbind(generation, data.frame(hub,
    source = rep("delphi_production_schedule", nrow(hub))))
  stopifnot(!anyDuplicated(generation[c("disease", "reference_date")]),
    all(abs(as.numeric(generation$generation_date - (generation$reference_date - 3))) <= 7))
  round_cutoff <- function(disease, reference_date) {
    key <- paste(disease, reference_date)
    known <- match(key, paste(generation$disease, generation$reference_date))
    known[!key %in% submitted] <- NA
    date <- reference_date - 3
    date[!is.na(known)] <- generation$generation_date[known[!is.na(known)]]
    source <- ifelse(!is.na(known), generation$source[known], ifelse(key %in% submitted,
      "wednesday_deadline", "wednesday_no_submission"))
    data.frame(date, source)
  }
  reports <- as.data.frame(vintages[!is.na(vintages$location),
    c("disease", "location", "geo_value", "target_end_date", "version", "value")])
  reports$source <- ifelse(reports$version == min(reports$version), "archive_start", "archive")
  # Final counts before the report versions: a synthetic first report of the latest value
  # on the Wednesday after the week, before every stored version.
  reports <- reports[order(reports$disease, reports$location, reports$target_end_date,
    reports$version), , drop = FALSE]
  final <- reports[!duplicated(reports[c("disease", "location", "target_end_date")],
    fromLast = TRUE) & reports$target_end_date <= settings$final_data_through, , drop = FALSE]
  final$version <- final$target_end_date + 4
  final$source <- rep("synthetic_final", nrow(final))
  stopifnot(all(final$version < min(reports$version)))
  reports <- rbind(final, reports)
  # Scoring truth: the latest report of each week, the archive's last version.
  settings$truth_version <- max(reports$version)
  settings$feedback_protocol <- feedback_protocol()
  group_key <- function(x) paste(x$disease, x$location, sep = "_")
  forecast_groups <- split(seq_len(nrow(selected)), group_key(selected))
  report_groups <- split(seq_len(nrow(reports)), group_key(reports))
  # Period of a target week (see the designs above).
  period_of <- function(date) {
    season <- season_of(date, settings$season_start_month)
    tuning <- if (settings$design == "seasons") season %in% settings$tuning_seasons else
      date <= settings$hub_tuning_through
    ifelse(tuning, "tuning", ifelse(season %in% settings$evaluation_seasons, "evaluation", "other"))
  }
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

    # Peaks: local maxima of the latest reports of all weeks (find_peaks()), without the
    # incomplete weeks; kept between the first and the last forecast target week.
    series <- latest[is.finite(latest$value) &
      (latest$target_end_date < settings$incomplete_weeks[1] |
        latest$target_end_date > settings$incomplete_weeks[2]), , drop = FALSE]
    found <- find_peaks(series$value, settings$peak_halfwidth, settings$peak_drop,
      settings$peak_min_prominence, settings$peak_smooth,
      scale = max(series$value[series$target_end_date %in% f$target_end_date]))
    found <- found[series$target_end_date[found$row] >= min(f$target_end_date) &
      series$target_end_date[found$row] <= max(f$target_end_date), , drop = FALSE]
    peak_date <- series$target_end_date[found$row]
    location_peaks <- if (nrow(found)) data.frame(disease, location = r$location[1],
      geo_value = r$geo_value[1], location_name,
      season = season_of(peak_date, settings$season_start_month),
      half = half_season_of(peak_date, settings$season_start_month),
      period = period_of(peak_date), peak_date, peak_value = series$value[found$row],
      relative_prominence = found$relative_prominence,
      forecast_week = peak_date %in% f$target_end_date)
    if (is.null(location_peaks)) {
      warning("No forecast peak for ", key, ": no before/after windows.", call. = FALSE)
    } else {
      location_peaks <- location_peaks[order(location_peaks$peak_date), , drop = FALSE]
      rownames(location_peaks) <- NULL
      peaks[[key]] <- location_peaks
    }

    # The median identifies the issue/target rows; all quantiles enter base below.
    medians <- f[is.finite(f$truth_finalized) & f$quantile_level == .5, , drop = FALSE]
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

      # The issue cutoff of each row's round, and every stored report of its target week
      # (including invalid ones) for the feedback schedules and the re-run design.
      cutoff <- round_cutoff(disease, reference_date)
      versions <- r[r$target_end_date %in% target_date, , drop = FALSE]
      versions <- data.frame(row = match(versions$target_end_date, target_date),
        version = versions$version, value = versions$value, source = versions$source)
      versions <- versions[order(versions$row, versions$version), , drop = FALSE]
      rownames(versions) <- NULL
      # Each earlier forecast is learned once, at the first later round at least d weeks
      # after its first valid report whose current report is valid, from that report
      # (feedback_schedule() in R/03_functions.R, with its audit).
      feedback_by_delay <- lapply(setNames(settings$feedback_delays, settings$feedback_delays),
        function(delay) feedback_schedule(versions, cutoff$date, delay, reference_date))

      season <- season_of(target_date, settings$season_start_month)
      nearest <- if (is.null(location_peaks)) rep(NA_integer_, length(target_date)) else
        apply(abs(outer(as.numeric(target_date), as.numeric(location_peaks$peak_date), "-")),
          1, which.min)
      relative_week <- as.numeric(target_date - location_peaks$peak_date[nearest]) / 7

      # Save the exact pre-first-forecast history (reports by the first round's cutoff),
      # before filtering finalized truth. Scripts 02/03 can use its median, maximum, or no
      # scaling.
      scale_as_of <- round_cutoff(disease, min(fh$reference_date))$date
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
        # Online feedback per delay (get_case() in R/03_functions.R selects one) and
        # every stored report of each row's target week (row, version, value, source).
        feedback_by_delay = feedback_by_delay, report_versions = versions,
        # The base, and the issue cutoff of each round (the as-of date of its reports and
        # feedback) with its source.
        base_name = settings$base[[disease]], as_of = cutoff$date, cutoff_source = cutoff$source,
        # Season of the target week and its period: "tuning", "evaluation" or "other".
        season = season, period = period_of(target_date),
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
# Labels of the periods for the plots of 02, 04 and 05.
prepare_settings$tuning_label <- if (design == "seasons")
  paste("season", paste(prepare_settings$tuning_seasons, collapse = " and ")) else
  paste(format(prepare_settings$start_date, "%b %Y"), "to",
    format(prepare_settings$hub_tuning_through, "%b %Y"))
prepare_settings$evaluation_label <- if (design == "seasons")
  paste("seasons", paste(prepare_settings$evaluation_seasons, collapse = " and ")) else
  paste(format(prepare_settings$hub_tuning_through + 1, "%b %Y"), "to the end of season",
    tail(prepare_settings$evaluation_seasons, 1))
peak_aligned <- prepare_data(base_forecasts, vintages, prepare_settings)
dir.create("output/data", recursive = TRUE, showWarnings = FALSE)
saveRDS(peak_aligned, "output/data/peak_aligned.rds")
# The feedback schedules, all locations and horizons: rows with a valid report, consumed
# rows, rows learned from a synthetic final report, horizon -1 style rows whose first
# report is known at their own issue (learned at the next round at the earliest), consumed
# rows learned from a revision instead of their first report, consumed rows that waited
# at a round after the delay because their current report was invalid, rows consumed a
# week or more after the delay and their own issue (missing rounds, invalid reports, the
# off-season), and rows with a valid report never consumed (the last rounds; reports
# invalid ever since).
feedback_rows <- feedback_audit(peak_aligned)
cat("\nFEEDBACK SCHEDULES (", feedback_protocol()$id, ")\n", feedback_protocol()$cutoff, "\n",
  feedback_protocol()$vintage, "\n", sep = "")
print(do.call(rbind, lapply(split(feedback_rows, list(feedback_rows$delay,
  feedback_rows$disease), drop = TRUE), function(a) {
  consumed <- !is.na(a$consumed_row)
  data.frame(disease = a$disease[1], delay = a$delay[1], rows = nrow(a),
    with_valid_report = sum(!is.na(a$first_report)), consumed = sum(consumed),
    synthetic_final = sum(consumed & a$selected_source %in% "synthetic_final"),
    known_at_issue = sum(a$first_report <= a$issue_cutoff, na.rm = TRUE),
    learned_revision = sum(consumed & a$selected_version != a$first_report),
    waited_invalid = sum(consumed & a$invalid_cutoffs > 0),
    consumed_a_week_late = sum(consumed &
      a$consumed_cutoff - pmax(a$maturity, a$issue_cutoff + 1) >= 7),
    never_consumed = sum(!is.na(a$first_report) & !consumed))
})), row.names = FALSE)
cat("Issue cutoffs by source (rows of all cases):\n")
print(table(unlist(lapply(peak_aligned$cases, `[[`, "cutoff_source"))))


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
