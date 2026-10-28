# Purpose: understand the two data files and how reported outcomes are revised.
# install.packages(c("nanoparquet", "pak"))
# pak::pkg_install("cmu-delphi/epiprocess@main")
.libPaths(c("tmp/multiqt-lib", .libPaths()))  # Same local package folder as script 02.
library(epiprocess)

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
# truth_finalized   The later observed count supplied alongside forecasts for
#                   evaluation. It is the SAME outcome across the 23 quantile rows
#                   for the same disease, location and target week. Note that this #                   seems to be the final reported vintage for covid but not for flu.
#
forecasts <- nanoparquet::read_parquet("cmu_timeseries_submitted_forecasts.parquet")
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

# The report archive
# begins on Nov 19, 2024 with a bulk load of older target weeks; it does not
# show when those older counts were originally published.

vintages <- nanoparquet::read_parquet("nhsn_vintages.parquet")
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

