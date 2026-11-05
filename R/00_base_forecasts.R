# Purpose: download the hub forecasts and NHSN reports once.
# Run from the project root: Rscript R/00_base_forecasts.R
# Order of the workflow: 00 -> 01 -> 02 / 04.
# Writes to the project root:
#   cmu_timeseries_hub_forecasts.parquet  CMU-TimeSeries hub submissions;
#   nhsn_vintages_epidata.parquet         every stored NHSN report version.
library(dplyr)

# 0. SETTINGS -----------------------------------------------------------------
forecasts_path <- "cmu_timeseries_hub_forecasts.parquet"
vintages_path <- "nhsn_vintages_epidata.parquet"

# 1. DOWNLOAD THE HUB FORECASTS AND REPORTS ONCE ------------------------------
# Both tables are downloaded once and stored as Parquet in the project root. Delete
# either file to download both again, so the observed counts stored with the forecasts
# match the reports. R/01_data_overview.R describes their columns.
if (!all(file.exists(c(forecasts_path, vintages_path)))) {
  # Forecasts: every CMU-TimeSeries submission to the CDC forecast hubs. Keep the
  # hospital admission quantiles; the hubs also hold ED-visit targets.
  hubs <- c(covid = "CDCgov/covid19-forecast-hub", flu = "cdcepi/FluSight-forecast-hub")
  downloaded_forecasts <- bind_rows(lapply(names(hubs), function(disease) {
    files <- jsonlite::fromJSON(paste0("https://api.github.com/repos/", hubs[[disease]],
                                       "/contents/model-output/CMU-TimeSeries"))
    files <- files[grepl("\\.csv$", files$name), ]
    message("Fetching ", nrow(files), " ", disease, " forecast files")
    bind_rows(lapply(files$download_url, function(url) {
      read.csv(url, colClasses = c(location = "character", output_type_id = "character")) %>%
        filter(target == paste("wk inc", disease, "hosp"), output_type == "quantile")
    })) %>%
      transmute(disease, reference_date = as.Date(reference_date), target, horizon,
                target_end_date = as.Date(target_end_date), location,
                quantile_level = as.numeric(output_type_id), value)
  }))

  # Reports: every stored version of the weekly COVID and flu admissions for all
  # states, territories and the US total from the Delphi Epidata API (source
  # "nhsn") via epidatr. Needs DELPHI_EPIDATA_KEY for sensible rate limits.
  signals <- c(covid = "confirmed_admissions_covid_ew", flu = "confirmed_admissions_flu_ew")
  # One request per signal keeps each response small enough to retry cheaply.
  downloaded_reports <- bind_rows(lapply(names(signals), function(disease) {
    message("Fetching ", signals[[disease]])
    epidatr::epidata_archive(
      source = "nhsn",
      signals = signals[[disease]],
      geo_type = c("state", "nation"),
      geo_values = "*"
    ) %>% mutate(disease = disease)
  }))
  # reference_time is the Saturday ending the target week; report_time is the
  # version. Locations get the forecast-table codes (FIPS codes and "US");
  # territories without forecasts keep a missing location.
  census <- epidatasets::state_census
  downloaded_vintages <- downloaded_reports %>%
    transmute(disease, geo_value,
              location = if_else(geo_value == "us", "US",
                                 census$fips[match(geo_value, census$abbr)]),
              target_end_date = reference_time, version = as.Date(report_time), value) %>%
    mutate(location = if_else(location %in% downloaded_forecasts$location,
                              location, NA_character_)) %>%
    arrange(disease, location, geo_value, target_end_date, version)

  # The hubs supply no observed counts: truth_finalized is the latest stored
  # report of each target week, and missing while the week is unreported.
  latest_reports <- downloaded_vintages %>%
    filter(!is.na(location)) %>%
    group_by(disease, location, target_end_date) %>%
    slice_max(version, n = 1) %>%
    ungroup() %>%
    select(disease, location, target_end_date, truth_finalized = value)
  downloaded_forecasts %>%
    left_join(latest_reports, by = c("disease", "location", "target_end_date")) %>%
    arrange(disease, reference_date, horizon, location, quantile_level) %>%
    nanoparquet::write_parquet(forecasts_path)
  nanoparquet::write_parquet(downloaded_vintages, vintages_path)
}
