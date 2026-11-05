# Purpose: download the hub forecasts and NHSN reports, and make our own base
# forecasts with Delphi's windowed_seasonal forecaster (their code, unchanged).
# Run from the project root: Rscript R/00_base_forecasts.R
# Order of the workflow: 00 -> 01 -> 02 / 04. About 25 minutes per disease the
# first time; later runs read the cached forecasts (delete the parquet file below to
# rebuild it from the cache).
# Writes to the project root:
#   cmu_timeseries_hub_forecasts.parquet  CMU-TimeSeries hub submissions (section 1);
#   nhsn_vintages_epidata.parquet         every stored NHSN report version (section 1);
#   windowed_seasonal_forecasts.parquet   our windowed_seasonal forecasts (section 5).
# Delphi's code, its pinned packages and the caches go to external/delphi/ (not in Git).
library(dplyr)

# 0. SETTINGS -----------------------------------------------------------------
# cmu-delphi/delphi-hub-forecasting (formerly exploration-tooling): forecasters, as-of
# snapshots and the replay settings of its calibration harness
# (scripts/calibration/calibration_harness.R, CH_CONFIGS).
delphi_commit <- "9abebcfbcd4330aae5314fae51c5e8e5524675b5"  # 2026-11-04.
# The package versions in its renv.lock; other packages come from the usual library.
delphi_packages <- c(epiprocess = "889ab0304582a65e7ca1048babac1a5bfa5d9759",
                     epipredict = "43b352f18ccabbed68cf6783ee5669d875d64514",
                     epidatr = "b4360372e3319efb6d3e8b0d55674dd21e2a65b0")
delphi_dir <- "external/delphi"
diseases <- c("covid", "flu")
aheads <- -1:3                       # Horizons; Delphi's replays use 0:3.
# Forecast dates are Wednesdays (reference date = the Saturday after). 2023/24: Delphi's
# burn-in rounds. There are no NHSN report versions before 2024-11-19, so these rounds
# use the latest reports, cut at the forecast date (Delphi's "cheating" snapshot): no
# revisions. From 2024-11-20 on, every round uses the reports as of its generation date.
final_data_from <- as.Date("2023-10-11")
final_data_through <- as.Date("2024-04-24")
versioned_from <- as.Date("2024-11-20")
# Flu: Delphi's flu windowed_seasonal also trains on ILI+ and FluSurv rows, built with
# these files of their non-public aux_data/flusion_data folder. Put them in
# external/delphi/aux_data/flusion_data; without them, flu is skipped.
flu_aux_files <- c("flu_burden.csv", "us_pop.csv",
  paste0("WHO_NREVSS_Combined_prior_to_2015_16_", c("Nation", "HHS", "State"), ".csv"),
  paste0("WHO_NREVSS_Clinical_Labs_", c("Nation", "HHS", "State"), ".csv"))
forecasts_path <- "cmu_timeseries_hub_forecasts.parquet"
vintages_path <- "nhsn_vintages_epidata.parquet"
own_forecasts_path <- "windowed_seasonal_forecasts.parquet"

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

# 2. DELPHI'S CODE AND PACKAGES -----------------------------------------------
# The pinned packages are installed once into external/delphi/library, which comes
# first on the library path.
delphi_lib <- file.path(delphi_dir, "library")
dir.create(delphi_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(delphi_lib, .libPaths()))
for (package in names(delphi_packages)) {
  installed <- tryCatch(packageDescription(package, lib.loc = delphi_lib)$RemoteSha,
    error = function(e) NULL, warning = function(w) NULL)
  if (!identical(installed, delphi_packages[[package]]))
    remotes::install_github(paste0("cmu-delphi/", package, "@", delphi_packages[[package]]),
      lib = delphi_lib, dependencies = FALSE, upgrade = "never", force = TRUE)
}
# Checked without loading them: here must first load inside code_dir (below).
for (package in c("zeallot", "here"))
  if (system.file(package = package) == "")
    install.packages(package, lib = delphi_lib, repos = "https://cloud.r-project.org")

# Delphi's code at the pinned commit, downloaded once.
code_dir <- file.path(delphi_dir, paste0("delphi-hub-forecasting-", substr(delphi_commit, 1, 8)))
if (!dir.exists(code_dir)) {
  archive <- tempfile(fileext = ".tar.gz")
  download.file(paste0("https://codeload.github.com/cmu-delphi/delphi-hub-forecasting/tar.gz/",
    delphi_commit), archive, mode = "wb")
  top <- sub("/.*", "", untar(archive, list = TRUE)[1])
  untar(archive, exdir = delphi_dir)
  stopifnot(file.rename(file.path(delphi_dir, top), code_dir))
}
# Delphi's functions find aux_data/ through here::here(): make code_dir its root.
file.create(file.path(code_dir, ".here"))
owd <- setwd(code_dir)
here_root <- here::here()
setwd(owd)
if (normalizePath(here_root) != normalizePath(code_dir))
  stop("here::here() points to ", here_root, ". Run this script with Rscript from the project root.")

# Delphi's R/load_all.R attaches every package of R/imports.R and sources every file
# in R/. Attach the packages the forecasts use, in the same order, and source the files.
delphi_imports <- c("assertthat", "cli", "data.table", "dplyr", "epidatr", "epipredict",
  "epiprocess", "ggplot2", "glue", "httr2", "jsonlite", "lubridate", "magrittr", "MMWRweek",
  "parsnip", "purrr", "quantreg", "readr", "recipes", "rlang", "scales", "slider", "stringr",
  "tibble", "tidyr", "zeallot")
suppressPackageStartupMessages(for (package in delphi_imports)
  library(package, character.only = TRUE))
for (file in list.files(file.path(code_dir, "R"), "[.]R$", recursive = TRUE, full.names = TRUE))
  if (!basename(file) %in% c("imports.R", "load_all.R")) source(file)
# From the calibration harness: the replay settings and helpers only (the harness itself
# would source load_all.R and read Delphi's targets stores).
harness <- parse(file.path(code_dir, "scripts/calibration/calibration_harness.R"),
  keep.source = FALSE)
for (expression in harness)
  if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      as.character(expression[[2]]) %in% c("CH_SCHEDULE_DELAYS", "CH_CONFIGS", "ch_use",
        "ch_schedule", "ch_to_hub"))
    eval(expression, globalenv())

# 3. TRAINING DATA ------------------------------------------------------------
# Delphi's nhsn_prod_archive (pipelines/<disease>_hosp_prod.R): the NHSN reports with
# season columns, Wednesday time labels (week-ending Saturday - 3) and source "nhsn",
# without four territories. Delphi reads the same Epidata v5 archive with its own
# fetch; here it comes from nhsn_vintages_epidata.parquet, so the forecasts and
# R/01_data_overview.R use the same report versions.
g_insufficient_data_geos <- c("as", "mp", "vi", "gu")
g_time_value_adjust <- 3
vintages <- nanoparquet::read_parquet(vintages_path)
training_archive <- function(disease) {
  reports <- vintages[vintages$disease == disease, ]
  nhsn_archive_data <- tibble(geo_value = reports$geo_value, time_value = reports$target_end_date,
    version = reports$version, value = reports$value) %>%
    as_epi_archive(compactify = TRUE)
  nhsn_dt <- nhsn_archive_data$DT %>%
    add_season_info() %>%
    mutate(geo_value = ifelse(geo_value == "usa", "us", geo_value),
      time_value = floor_date(time_value, "week", week_start = 7) + 3, source = "nhsn") %>%
    filter(geo_value %nin% g_insufficient_data_geos)
  if (disease == "covid") return(as_epi_archive(nhsn_dt, other_keys = "source", compactify = TRUE))
  # Flu: plus Delphi's ILI+ and FluSurv rows, as its flusurv and ili_plus targets
  # (R/targets/flu_data_targets.R) and flu nhsn_prod_archive build them.
  flusurv <- generate_flusurv_adjusted()$DT %>%
    mutate(time_value = time_value + g_time_value_adjust, version = version + g_time_value_adjust) %>%
    mutate(source = "flusurv") %>%
    select(geo_value, time_value, version, hhs = adj_hosp_rate, source, agg_level, season, season_week)
  ili_plus <- gen_ili_data()$DT %>%
    tidyr::drop_na() %>%
    filter(hhs > 0.0001) %>%
    mutate(time_value = time_value + g_time_value_adjust, version = version + g_time_value_adjust) %>%
    select(geo_value, time_value, version, hhs, source, agg_level, season, season_week)
  to_keep <- ili_plus %>%
    group_by(geo_value, season) %>%
    summarize(total_count = n_distinct(time_value), .groups = "drop") %>%
    filter(total_count >= 20) %>%
    select(geo_value, season)
  ili_plus <- to_keep %>%
    left_join(ili_plus, by = join_by(geo_value, season), relationship = "many-to-many") %>%
    as_epi_archive(other_keys = "source", compactify = TRUE) %>%
    extract2("DT")
  extra_dt <- bind_rows(ili_plus, flusurv) %>%
    filter(agg_level %in% c("state", "nation"), geo_value %nin% g_insufficient_data_geos) %>%
    mutate(epiweek = epiweek(time_value), epiyear = epiyear(time_value)) %>%
    select(geo_value, source, time_value, version, value = hhs, season, season_week, epiweek, epiyear)
  bind_rows(nhsn_dt, extra_dt) %>%
    as_epi_archive(other_keys = "source", compactify = TRUE)
}

# 4. FORECASTS ----------------------------------------------------------------
# Delphi's replay of windowed_seasonal (CH_CONFIGS in the calibration harness): the
# grid row of pipelines/<disease>_hosp_prod.R, i.e. scaled_pop_seasonal with quantile
# regression (method "fn"), lags 0 and 7 days, whitening and sorted quantiles. Flu uses
# the seasonal training window; covid's windowed_seasonal uses none. At ahead -1 the
# outcome equals a predictor (the latest report), so those forecasts have zero width.
# The loop is that of ch_backfill(): one snapshot per forecast date, one
# run_forecaster() call per ahead, and a date that fails is reported and skipped.
# quantile_reg is deterministic, so ch_backfill()'s seeds are not needed. Each disease
# is cached in external/delphi/cache.
forecast_disease <- function(disease) {
  ch_use(paste0(disease, "_windowed_seasonal"))
  # The grid row's trainer is g_quantreg = quantile_reg(method = "fn"); CH_CONFIGS has
  # quantile_reg() with quantreg's default "br", which can stall when the outcome equals
  # a predictor, as at ahead -1 (Delphi's notes/quantreg-solver.md).
  CH_SPEC$params$trainer <<- epipredict::quantile_reg(method = "fn")
  archive <- training_archive(disease)
  schedule <- bind_rows(
    ch_schedule(final_data_through, final_data_from) %>% mutate(as_of_policy = "cheating"),
    ch_schedule(archive$versions_end, versioned_from) %>% mutate(as_of_policy = "asof"))
  cache_file <- file.path(delphi_dir, "cache", paste0(disease, "_windowed_seasonal_",
    rlang::hash(list(delphi_commit, CH_SPEC, aheads, schedule, archive$DT)), ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  failures <- list()
  forecasts <- purrr::map(seq_len(nrow(schedule)), .progress = paste(disease, "forecasts"),
    function(i) tryCatch({
      snapshot <- make_forecast_snapshot(archive, forecast_date = schedule$forecast_date[i],
        generation_date = schedule$generation_date[i], as_of_policy = schedule$as_of_policy[i])
      bind_rows(lapply(aheads, function(ahead) run_forecaster(snapshot = snapshot,
        forecaster = get(CH_SPEC$forecaster), aheads = ahead * CH_SPEC$ahead_multiplier,
        params = CH_SPEC$params, id = CH_FORECASTER_ID,
        target_date_shift = CH_SPEC$target_date_shift,
        join_extra_data = FALSE, extra_data = NULL, excluded_geos = NULL,
        sort_quantiles = CH_SPEC$sort_quantiles) %>% mutate(ahead = ahead)))
    }, error = function(e) {
      failures[[as.character(schedule$forecast_date[i])]] <<- conditionMessage(e)
      NULL
    })) %>%
    bind_rows() %>%
    left_join(schedule, by = "forecast_date")
  if (length(failures)) warning(disease, ": ", length(failures), " forecast date(s) failed and ",
    "were skipped: ", paste(names(failures), collapse = ", "), call. = FALSE, immediate. = TRUE)
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(forecasts, cache_file)
  forecasts
}

# 5. HUB FORMAT AND SAVE ------------------------------------------------------
# The columns of cmu_timeseries_hub_forecasts.parquet, plus forecast_date (the nominal
# Wednesday), generation_date (the as-of date of the reports) and data ("final" or
# "versions"). truth_finalized: the latest stored report, as for the hub forecasts.
if (!file.exists(own_forecasts_path)) {
  aux_dir <- file.path(delphi_dir, "aux_data", "flusion_data")
  if ("flu" %in% diseases && !all(file.exists(file.path(aux_dir, flu_aux_files)))) {
    message("Flu skipped: put Delphi's ", paste(flu_aux_files, collapse = ", "), " in ",
      aux_dir, ", delete ", own_forecasts_path, " and rerun.")
    diseases <- setdiff(diseases, "flu")
  }
  if ("flu" %in% diseases) {
    dir.create(file.path(code_dir, "aux_data", "flusion_data"), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(aux_dir, flu_aux_files), file.path(code_dir, "aux_data", "flusion_data"),
      overwrite = TRUE)
  }
  latest_reports <- vintages %>%
    filter(!is.na(location)) %>%
    group_by(disease, location, target_end_date) %>%
    slice_max(version, n = 1) %>%
    ungroup() %>%
    select(disease, location, target_end_date, truth_finalized = value)
  own_forecasts <- bind_rows(lapply(diseases, function(disease) {
    raw <- forecast_disease(disease)
    dates <- raw %>%
      distinct(forecast_date, generation_date, as_of_policy) %>%
      mutate(reference_date = get_forecast_reference_date(forecast_date))
    ch_to_hub(raw, disease) %>%
      left_join(dates, by = "reference_date") %>%
      transmute(disease, reference_date, target = paste("wk inc", disease, "hosp"), horizon,
        target_end_date, location, quantile_level = level, value, forecast_date, generation_date,
        data = if_else(as_of_policy == "cheating", "final", "versions"))
  })) %>%
    left_join(latest_reports, by = c("disease", "location", "target_end_date")) %>%
    arrange(disease, reference_date, horizon, location, quantile_level)
  nanoparquet::write_parquet(own_forecasts, own_forecasts_path)
}
own_forecasts <- nanoparquet::read_parquet(own_forecasts_path)
print(own_forecasts %>% count(disease, data, horizon) %>% as.data.frame())
