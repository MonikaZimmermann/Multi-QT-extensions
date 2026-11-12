# Purpose: download the hub forecasts and NHSN reports, rebuild the HHS reports before
# NHSN, and make our own base forecasts for covid and flu with Delphi's windowed_seasonal
# forecaster (their code, unchanged).
# Run from the project root: Rscript R/00_base_forecasts.R
# Order of the workflow: 00 -> 01 -> 05 -> 02 / 04. The first run downloads the daily HHS
# reports (section 4: minutes with DELPHI_EPIDATA_KEY, about two hours without, because
# Epidata allows anonymous clients about 60 requests per hour), builds their weekly
# archive (about an hour: Delphi's script recomputes all weeks at each of about 90 HHS
# issues) and makes the forecasts (about 25 minutes per disease); later runs read the
# caches. A change of any input of the forecasts (section 6: reports, exclusions, flu
# files, Delphi's code and packages, dates, horizons) rebuilds them; rounds that failed
# are tried again at the next run; the previous parquet is retained under agent/.
# Writes to the project root:
#   cmu_timeseries_hub_forecasts.parquet  CMU-TimeSeries hub submissions (section 1);
#   nhsn_vintages_epidata.parquet         every stored NHSN report version (section 1);
#   hhs_vintages_epidata.parquet          every weekly HHS report version (section 4);
#   windowed_seasonal_forecasts.parquet   our windowed_seasonal forecasts (section 7).
# Delphi's code, its pinned packages, its flu files and the caches go to external/delphi/
# (not in Git).
library(dplyr)
source("R/03_functions.R")

# 0. SETTINGS -----------------------------------------------------------------
# cmu-delphi/delphi-hub-forecasting (formerly exploration-tooling): forecasters, as-of
# snapshots, the replay settings of its calibration harness
# (scripts/calibration/calibration_harness.R, CH_CONFIGS) and its HHS archive
# (scripts/one_offs/hhs_2023_24_archive.R).
delphi_commit <- "9abebcfbcd4330aae5314fae51c5e8e5524675b5"  # 2026-11-04.
# The package versions in its renv.lock; other packages come from the usual library.
delphi_packages <- c(epiprocess = "889ab0304582a65e7ca1048babac1a5bfa5d9759",
                     epipredict = "43b352f18ccabbed68cf6783ee5669d875d64514",
                     epidatr = "b4360372e3319efb6d3e8b0d55674dd21e2a65b0")
delphi_dir <- "external/delphi"
diseases <- c("covid", "flu")
aheads <- 0:3                        # Horizons, as in Delphi's replays; no -1 (section 6).
# Forecast dates are Wednesdays (reference date = the Saturday after): Delphi's 2023/24
# burn-in rounds (BURN_IN_DATES in scripts/calibration/calibration_ws_replay.R), then
# every round from NHSN's first report versions on. Every round uses the reports as of
# its generation date: the HHS reports in 2023/24 (section 4), the NHSN reports from
# 2024-11-19 on. Nothing is forecast from May to November 2024.
burn_in_from <- as.Date("2023-10-11")
burn_in_through <- as.Date("2024-04-24")
versioned_from <- as.Date("2024-11-20")
# HHS reports: the daily admissions reported to HHS (Epidata source "hhs") from this
# admission day on, with every daily issue. Delphi runs its HHS archive from 2020-08 for
# the burn-in (hhs_weekly_archive_2020_2024.rds), so a 2023/24 snapshot also sees the
# earlier seasons.
hhs_start <- "20200801"
hhs_daily_dir <- file.path(delphi_dir, "cache", "hhs_daily")
hhs_archive_dir <- file.path(delphi_dir, "cache", "hhs_archive")
# Flu: Delphi's flu windowed_seasonal also trains on ILI+ and FluSurv rows, built with
# these files of their non-public aux_data/flusion_data folder (shared by Delphi as
# aux_data.zip). They are copied once from flu_aux_zip into
# external/delphi/aux_data/flusion_data; without them, flu is skipped.
flu_aux_zip <- "~/Downloads/aux_data.zip"
flu_aux_files <- c("flu_burden.csv", "us_pop.csv",
  paste0("WHO_NREVSS_Combined_prior_to_2015_16_", c("Nation", "HHS", "State"), ".csv"),
  paste0("WHO_NREVSS_Clinical_Labs_", c("Nation", "HHS", "State"), ".csv"))
forecasts_path <- "cmu_timeseries_hub_forecasts.parquet"
vintages_path <- "nhsn_vintages_epidata.parquet"
hhs_vintages_path <- "hhs_vintages_epidata.parquet"
own_forecasts_path <- "windowed_seasonal_forecasts.parquet"
# The same reviewed report exclusions are used here and in script 01. Record the
# filter on generated forecasts so script 01 can reject an older, unfiltered file.
reports_filter_id <- report_filter_id()

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
vintages <- filter_reports(nanoparquet::read_parquet(vintages_path))

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
# would source load_all.R and read Delphi's targets stores), among them
# ch_inputs_with_burn_in(), which puts the HHS reports before NHSN's (section 5).
harness <- parse(file.path(code_dir, "scripts/calibration/calibration_harness.R"),
  keep.source = FALSE)
for (expression in harness)
  if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      as.character(expression[[2]]) %in% c("CH_SCHEDULE_DELAYS", "CH_CONFIGS", "ch_use",
        "ch_schedule", "ch_to_hub", "ch_inputs_with_burn_in"))
    eval(expression, globalenv())

# 3. DELPHI'S FLU FILES -------------------------------------------------------
# Copied once from Delphi's aux_data.zip (flu_aux_zip) into
# external/delphi/aux_data/flusion_data, then into Delphi's code folder, where its
# functions read them. Without them, flu is skipped.
aux_dir <- file.path(delphi_dir, "aux_data", "flusion_data")
if (!all(file.exists(file.path(aux_dir, flu_aux_files))) && file.exists(path.expand(flu_aux_zip))) {
  members <- file.path("aux_data", "flusion_data", flu_aux_files)
  absent <- setdiff(members, unzip(path.expand(flu_aux_zip), list = TRUE)$Name)
  if (length(absent)) stop(flu_aux_zip, " lacks ", paste(absent, collapse = ", "))
  unzip(path.expand(flu_aux_zip), files = members, exdir = delphi_dir)
  message("Copied Delphi's flu files from ", flu_aux_zip, " to ", aux_dir)
}
if ("flu" %in% diseases && !all(file.exists(file.path(aux_dir, flu_aux_files)))) {
  message("Flu skipped: put Delphi's ", paste(flu_aux_files, collapse = ", "), " in ",
    aux_dir, " (or aux_data.zip at ", flu_aux_zip, ") and rerun.")
  diseases <- setdiff(diseases, "flu")
}
if ("flu" %in% diseases) {
  dir.create(file.path(code_dir, "aux_data", "flusion_data"), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file.path(aux_dir, flu_aux_files),
    file.path(code_dir, "aux_data", "flusion_data"), overwrite = TRUE))
}

# 4. HHS REPORTS BEFORE NHSN --------------------------------------------------
# NHSN's report versions begin on 2024-11-19 with a bulk load of all earlier weeks.
# Before May 2024, hospitals reported the same admissions daily to HHS, and Epidata keeps
# every daily issue. Delphi's scripts/one_offs/hhs_2023_24_archive.R turns them into
# weekly report versions, run here unchanged:
#   - week: the hhs time_value is the admission day; the week ending Saturday S is the
#     sum of the admission days S-7..S-1 (Delphi checked that other alignments match
#     NHSN far worse), labelled with the Wednesday S - 3 as in Delphi's NHSN archive;
#   - versions: every daily issue is a version (its actual date), and a week is reported
#     only once all its 7 days are; issues before 2023-07-01 collapse into one version
#     on that date;
#   - US: the hhs nation signal (the sum of all states and territories).
# The script asks for all locations at once, which needs DELPHI_EPIDATA_KEY. Here each
# location is asked for separately (the same rows), cached per location in
# external/delphi/cache/hhs_daily, so an interrupted download resumes; the script then
# reads these rows from its raw cache instead of downloading. Its validation against
# NHSN runs too, with the NHSN reports of nhsn_vintages_epidata.parquet (what its
# get_nhsn_data_archive() fetches from Epidata): it prints the agreement of the HHS and
# NHSN final counts in 2023/24 and the revisions, and writes them to hhs_archive_dir.
# hhs_vintages_epidata.parquet holds the weekly versions in the columns of
# nhsn_vintages_epidata.parquet (target_end_date: the Saturday), for script 01.
# Delete hhs_archive_dir to rebuild the archive, hhs_daily_dir to download again.
hhs_archive_path <- file.path(hhs_archive_dir, "hhs_2023_24_weekly_archive.rds")
if (!file.exists(hhs_archive_path)) {
  hhs_script <- parse(file.path(code_dir, "scripts/one_offs/hhs_2023_24_archive.R"),
    keep.source = FALSE)
  # Its first expression sources load_all.R (done above); its fetch is the `if` on raw_cache.
  stopifnot(grepl("load_all.R", paste(deparse(hhs_script[[1]]), collapse = ""), fixed = TRUE))
  fetch <- which(vapply(hhs_script, function(expression) is.call(expression) &&
    identical(expression[[1]], as.name("if")) &&
    identical(expression[[2]], quote(file.exists(raw_cache))), logical(1)))
  stopifnot(length(fetch) == 1L)
  hhs <- new.env()
  hhs$commandArgs <- function(trailingOnly = FALSE) c(hhs_archive_dir, hhs_start)
  hhs$get_nhsn_data_archive <- function(disease) {
    reports <- vintages[vintages$disease == disease, ]
    tibble(geo_value = reports$geo_value, time_value = reports$target_end_date,
      version = reports$version, value = reports$value) %>%
      as_epi_archive(compactify = TRUE)
  }
  for (expression in hhs_script[2:(fetch - 1L)]) eval(expression, hhs)
  # The download: hhs$signals and hhs$time_range are the script's own.
  hhs_geos <- sort(unique(vintages$geo_value))
  dir.create(hhs_daily_dir, recursive = TRUE, showWarnings = FALSE)
  if (!nzchar(Sys.getenv("DELPHI_EPIDATA_KEY")))
    message("No DELPHI_EPIDATA_KEY: the HHS download waits for Epidata's anonymous rate ",
      "limit (about 60 requests per hour, ", 2 * length(hhs_geos), " requests).")
  for (disease in names(hhs$signals)) for (geo in hhs_geos) {
    file <- file.path(hhs_daily_dir, paste0(disease, "_", geo, ".rds"))
    if (file.exists(file)) next
    for (attempt in 1:12) {
      rows <- tryCatch(suppressWarnings(epidatr::pub_covidcast("hhs", hhs$signals[[disease]],
        if (geo == "us") "nation" else "state", "day", geo, time_values = hhs$time_range,
        issues = "*")), error = function(e) e)
      if (!inherits(rows, "error")) break
      message("HHS ", disease, " ", geo, ": ", conditionMessage(rows), "; waiting 65 s")
      Sys.sleep(65)
    }
    if (inherits(rows, "error")) stop("HHS download failed for ", disease, " ", geo,
      ". Rerun later: the download resumes.")
    rows <- if (nrow(rows)) data.frame(geo_value = rows$geo_value, time_value = rows$time_value,
      version = rows$issue, value = rows$value) else data.frame(geo_value = character(),
      time_value = as.Date(character()), version = as.Date(character()), value = numeric())
    saveRDS(rows, file)
    message("HHS ", disease, " ", geo, ": ", nrow(rows), " daily reports")
  }
  # The script's raw cache: one table of both diseases, in the columns of its fetch.
  saveRDS(bind_rows(lapply(names(hhs$signals), function(disease)
    bind_rows(lapply(hhs_geos, function(geo)
      readRDS(file.path(hhs_daily_dir, paste0(disease, "_", geo, ".rds"))))) %>%
      transmute(disease = disease, geo_value, time_value, version, value))), hhs$raw_cache)
  for (expression in hhs_script[fetch:length(hhs_script)]) eval(expression, hhs)
  rm(hhs)
}
if (!file.exists(hhs_vintages_path) ||
    file.mtime(hhs_vintages_path) < file.mtime(hhs_archive_path)) {
  hhs_weekly <- readRDS(hhs_archive_path)
  bind_rows(lapply(names(hhs_weekly), function(disease) as_tibble(hhs_weekly[[disease]]) %>%
    transmute(disease, geo_value, target_end_date = time_value + 3, version, value))) %>%
    mutate(location = vintages$location[match(geo_value, vintages$geo_value)]) %>%
    select(disease, geo_value, location, target_end_date, version, value) %>%
    arrange(disease, location, geo_value, target_end_date, version) %>%
    nanoparquet::write_parquet(hhs_vintages_path)
}
# Recorded on the forecasts, so script 01 can check that its HHS reports are the ones
# the forecasts were trained on.
hhs_reports_id <- paste0("hhs-", unname(tools::md5sum(hhs_vintages_path)))

# 5. TRAINING DATA ------------------------------------------------------------
# Delphi's nhsn_prod_archive (pipelines/<disease>_hosp_prod.R): the NHSN reports with
# season columns, Wednesday time labels (week-ending Saturday - 3) and source "nhsn",
# without four territories. Delphi reads the same Epidata v5 archive with its own
# fetch; here it comes from nhsn_vintages_epidata.parquet, so the forecasts and
# R/01_data_overview.R use the same report versions. Then, as in Delphi's replay with the
# 2023/24 burn-in (ch_inputs_with_burn_in()), the HHS reports of section 4 come first:
# relabelled "nhsn", kept only at versions before NHSN's first, so a 2023/24 snapshot sees
# the HHS reports alone and a later one the NHSN reports alone; HHS levels are used as
# they are (2-4 % below NHSN in 2023/24, Delphi's ROADMAP 1b).
g_insufficient_data_geos <- c("as", "mp", "vi", "gu")
g_time_value_adjust <- 3
# Flu: Delphi's ILI+ and FluSurv rows, as its flusurv and ili_plus targets
# (R/targets/flu_data_targets.R) and flu nhsn_prod_archive build them. Their Epidata
# reports are downloaded once and cached with the flu files they used.
flu_extras <- function() {
  cache_file <- file.path(delphi_dir, "cache", paste0("flu_extras_",
    rlang::hash(unname(tools::md5sum(file.path(aux_dir, flu_aux_files)))), ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))
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
  saveRDS(extra_dt, cache_file)
  extra_dt
}
training_archive <- function(disease) {
  reports <- vintages[vintages$disease == disease, ]
  # Keep masked versions as NA in the archive. A snapshot must select the latest
  # version available then, even when excluded, rather than fall back to an older
  # value. A subsequent valid revision can make the same target week usable again.
  nhsn_archive_data <- tibble(geo_value = reports$geo_value, time_value = reports$target_end_date,
    version = reports$version, value = reports$value) %>%
    as_epi_archive(compactify = TRUE)
  nhsn_dt <- nhsn_archive_data$DT %>%
    add_season_info() %>%
    mutate(geo_value = ifelse(geo_value == "usa", "us", geo_value),
      time_value = floor_date(time_value, "week", week_start = 7) + 3, source = "nhsn") %>%
    filter(geo_value %nin% g_insufficient_data_geos)
  nhsn_train <- (if (disease == "covid") nhsn_dt else bind_rows(nhsn_dt, flu_extras())) %>%
    as_epi_archive(other_keys = "source", compactify = TRUE)
  ch_inputs_with_burn_in(hhs_path = hhs_archive_path,
    inputs = list(nhsn_train = nhsn_train))$nhsn_train
}

# 6. FORECASTS ----------------------------------------------------------------
# Delphi's replay of windowed_seasonal (CH_CONFIGS in the calibration harness): the
# grid row of pipelines/<disease>_hosp_prod.R, i.e. scaled_pop_seasonal with quantile
# regression (method "fn"), lags 0 and 7 days, whitening and sorted quantiles. Flu uses
# the seasonal training window; covid's windowed_seasonal uses none. Ahead -1 is not
# forecast: its outcome equals a predictor (the latest report), so its forecasts would
# have zero width whenever the week before is reported by the generation date, and
# Delphi's submissions take ahead -1 from other forecasters.
# The loop is that of ch_backfill(): one as-of snapshot per forecast date (the 2023/24
# rounds: ws_burn_in_forecasts() in calibration_ws_replay.R) and one run_forecaster() call
# per ahead. quantile_reg is deterministic, so ch_backfill()'s seeds are not needed.
# The identity of a disease's forecasts (inputs_id) covers everything they depend on:
# Delphi's code (the sourced files at the pinned commit) and the package versions, the
# replay settings and horizons, the reviewed report exclusions, the NHSN and HHS reports,
# the flu files and the training archive built from them (with flu's ILI+ and FluSurv
# reports). Each round is cached under it, in
# external/delphi/cache/<disease>_windowed_seasonal_<inputs_id>/, as soon as it succeeds.
# A round that fails is recorded (failures.csv in that folder) and tried again at the next
# run; a run with missing rounds is never taken as complete (section 7).
delphi_code_id <- rlang::hash(unname(tools::md5sum(c(sort(list.files(file.path(code_dir, "R"),
  "[.]R$", recursive = TRUE, full.names = TRUE)), file.path(code_dir, "scripts",
  c("calibration/calibration_harness.R", "one_offs/hhs_2023_24_archive.R"))))))
package_versions <- vapply(unique(c(names(delphi_packages), delphi_imports, "hardhat",
  "workflows")), function(package) {
    description <- packageDescription(package)
    paste(description$Version, if (is.null(description$RemoteSha)) "" else description$RemoteSha)
  }, character(1))
hhs_reports <- nanoparquet::read_parquet(hhs_vintages_path)
use_replay <- function(disease) {
  ch_use(paste0(disease, "_windowed_seasonal"))
  # The grid row's trainer is g_quantreg = quantile_reg(method = "fn"); CH_CONFIGS has
  # quantile_reg() with quantreg's default "br", which can stall when the outcome equals
  # a predictor, as at ahead -1 (not forecast here; Delphi's notes/quantreg-solver.md).
  CH_SPEC$params$trainer <<- epipredict::quantile_reg(method = "fn")
  # Flu: the grid row's keys_to_ignore is the cell g_very_latent_locations, which
  # make_forecaster_grid() unwraps to list("source", c("flusurv", "ILI+")) (the ILI+ and
  # FluSurv rows are not forecast). CH_CONFIGS has the cell still wrapped,
  # list(list(...)), which default_args_list() cannot read ("'names' attribute [2] must
  # be the same length as the vector [0]" in every round); use the grid row's value.
  keys <- CH_SPEC$params$keys_to_ignore
  if (is.list(keys) && length(keys) == 1L && is.list(keys[[1]]))
    CH_SPEC$params$keys_to_ignore <<- keys[[1]]
}
# The training archive, the forecast dates (both diseases through NHSN's last version;
# flu's ILI+ and FluSurv reports go on) and the identities: inputs_id, and forecast_id,
# which adds the forecast dates and is recorded in the parquet.
forecast_inputs <- function(disease) {
  use_replay(disease)
  archive <- training_archive(disease)
  schedule <- bind_rows(ch_schedule(burn_in_through, burn_in_from),
    ch_schedule(max(vintages$version), versioned_from))
  # Hashed in forms that are the same in every R session: the parsnip trainer and the
  # data.table hold environments and pointers, so the trainer enters by its class, mode,
  # engine and argument expressions, and the archive as a plain data frame.
  trainer <- CH_SPEC$params$trainer
  spec <- CH_SPEC
  spec$params$trainer <- list(class(trainer), trainer$mode, trainer$engine,
    lapply(trainer$args, rlang::quo_get_expr), lapply(trainer$eng_args, rlang::quo_get_expr))
  inputs_id <- rlang::hash(list(delphi_commit, delphi_code_id, package_versions,
    CH_FORECASTER_ID, spec, aheads, reviewed_report_exclusions(),
    vintages[vintages$disease == disease, ], hhs_reports[hhs_reports$disease == disease, ],
    if (disease == "flu") unname(tools::md5sum(file.path(aux_dir, flu_aux_files))),
    as.data.frame(archive$DT)))
  list(archive = archive, schedule = schedule, inputs_id = inputs_id,
    forecast_id = rlang::hash(list(inputs_id, schedule)),
    round_dir = file.path(delphi_dir, "cache", paste0(disease, "_windowed_seasonal_", inputs_id)))
}
# The scheduled rounds of one disease: cached rounds are read, the others forecast (and
# cached when they succeed). Returns the forecasts, the rounds still missing and the
# failures of this run.
forecast_disease <- function(disease, inputs) {
  use_replay(disease)
  schedule <- inputs$schedule
  dir.create(inputs$round_dir, recursive = TRUE, showWarnings = FALSE)
  files <- file.path(inputs$round_dir,
    paste0(schedule$forecast_date, "_", schedule$generation_date, ".rds"))
  failures <- data.frame(forecast_date = as.Date(character()),
    generation_date = as.Date(character()), error = character())
  todo <- which(!file.exists(files))
  message(disease, ": ", nrow(schedule) - length(todo), " of ", nrow(schedule),
    " rounds cached, ", length(todo), " to forecast")
  purrr::walk(todo, .progress = paste(disease, "forecasts"), function(i) {
    round <- tryCatch({
      snapshot <- make_forecast_snapshot(inputs$archive, forecast_date = schedule$forecast_date[i],
        generation_date = schedule$generation_date[i], as_of_policy = "asof")
      round <- bind_rows(lapply(aheads, function(ahead) run_forecaster(snapshot = snapshot,
        forecaster = get(CH_SPEC$forecaster), aheads = ahead * CH_SPEC$ahead_multiplier,
        params = CH_SPEC$params, id = CH_FORECASTER_ID,
        target_date_shift = CH_SPEC$target_date_shift,
        join_extra_data = FALSE, extra_data = NULL, excluded_geos = NULL,
        sort_quantiles = CH_SPEC$sort_quantiles) %>% mutate(ahead = ahead)))
      if (!nrow(round) || !setequal(unique(round$ahead), aheads))
        stop("no forecasts for some horizons")
      round
    }, error = function(e) e)
    if (inherits(round, "error")) {
      failures <<- rbind(failures, data.frame(forecast_date = schedule$forecast_date[i],
        generation_date = schedule$generation_date[i], error = conditionMessage(round)))
      return(invisible())
    }
    # Written under a temporary name first, so an interrupted run leaves no partial round.
    saveRDS(round, paste0(files[i], ".part"))
    stopifnot(file.rename(paste0(files[i], ".part"), files[i]))
  })
  write.csv(failures, file.path(inputs$round_dir, "failures.csv"), row.names = FALSE)
  done <- file.exists(files)
  raw <- if (any(done)) bind_rows(lapply(files[done], readRDS)) %>%
    left_join(schedule, by = "forecast_date") else NULL
  list(raw = raw, missing = schedule[!done, ], failures = failures)
}

# 7. HUB FORMAT AND SAVE ------------------------------------------------------
# The columns of cmu_timeseries_hub_forecasts.parquet, plus forecast_date (the nominal
# Wednesday), generation_date (the as-of date of the reports), data (the reports the round
# was trained on: "hhs" before NHSN's first version, "nhsn" after), forecast_id (section
# 6) and missing_rounds (scheduled rounds without forecasts; script 01 refuses forecasts
# with missing rounds). truth_finalized: the latest stored NHSN report after the reviewed
# exclusions. report_filter_id and hhs_reports_id record the report filter and the HHS
# reports used for training (script 01 checks both).
# The parquet is kept only if it has every disease to forecast with its current
# forecast_id, the current report filter and HHS reports and no missing rounds. Otherwise
# it is rebuilt from the cached rounds, and only the rounds not cached are forecast.
inputs <- lapply(setNames(diseases, diseases), forecast_inputs)
own_forecasts <- if (file.exists(own_forecasts_path))
  nanoparquet::read_parquet(own_forecasts_path) else NULL
is_current <- function(disease) {
  rows <- own_forecasts[own_forecasts$disease == disease, ]
  nrow(rows) > 0 && all(c("forecast_id", "missing_rounds", "report_filter_id",
    "hhs_reports_id") %in% names(rows)) &&
    all(rows$forecast_id %in% inputs[[disease]]$forecast_id) && all(rows$missing_rounds %in% 0) &&
    all(rows$report_filter_id %in% reports_filter_id) && all(rows$hhs_reports_id %in% hhs_reports_id)
}
if (is.null(own_forecasts) || !all(vapply(diseases, is_current, logical(1)))) {
  if (!is.null(own_forecasts)) {
    archive_dir <- "agent/outputs/archive/base_forecasts"
    dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
    archive_file <- tempfile(paste0("windowed_seasonal_",
      format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = archive_dir, fileext = ".parquet")
    if (!file.copy(own_forecasts_path, archive_file))
      stop("Could not retain the previous forecasts in ", archive_file)
    message("Own forecasts incomplete, without a disease or made from other inputs: ",
      "rebuilding; previous file retained in ", archive_file)
  }
  latest_reports <- vintages %>%
    filter(!is.na(location)) %>%
    group_by(disease, location, target_end_date) %>%
    slice_max(version, n = 1) %>%
    ungroup() %>%
    select(disease, location, target_end_date, truth_finalized = value)
  nhsn_first_version <- min(vintages$version)
  runs <- lapply(setNames(diseases, diseases), function(disease)
    forecast_disease(disease, inputs[[disease]]))
  own_forecasts <- bind_rows(lapply(diseases, function(disease) {
    run <- runs[[disease]]
    if (is.null(run$raw)) return(NULL)
    dates <- run$raw %>%
      distinct(forecast_date, generation_date) %>%
      mutate(reference_date = get_forecast_reference_date(forecast_date))
    forecast_id <- inputs[[disease]]$forecast_id
    missing_rounds <- nrow(run$missing)
    ch_to_hub(run$raw, disease) %>%
      left_join(dates, by = "reference_date") %>%
      transmute(disease, reference_date, target = paste("wk inc", disease, "hosp"), horizon,
        target_end_date, location, quantile_level = level, value, forecast_date, generation_date,
        data = if_else(generation_date < nhsn_first_version, "hhs", "nhsn"),
        forecast_id = .env$forecast_id, missing_rounds = .env$missing_rounds)
  }))
  if (!nrow(own_forecasts)) {
    print(bind_rows(lapply(runs, `[[`, "failures")))
    stop("No forecast round succeeded; see the failures above.")
  }
  own_forecasts <- own_forecasts %>%
    left_join(latest_reports, by = c("disease", "location", "target_end_date")) %>%
    mutate(report_filter_id = reports_filter_id, hhs_reports_id = hhs_reports_id) %>%
    arrange(disease, reference_date, horizon, location, quantile_level)
  nanoparquet::write_parquet(own_forecasts, own_forecasts_path)
  missing <- bind_rows(lapply(diseases, function(disease) runs[[disease]]$missing %>%
    mutate(disease = disease)))
  if (nrow(missing)) {
    failed <- bind_rows(lapply(diseases, function(disease) runs[[disease]]$failures %>%
      mutate(disease = disease)))
    print(as.data.frame(failed))
    warning(nrow(missing), " forecast round(s) missing (",
      paste(missing$disease, missing$forecast_date, collapse = ", "), "): the parquet records ",
      "them as missing_rounds and script 01 refuses it. Rerun Rscript R/00_base_forecasts.R ",
      "to try them again.", call. = FALSE, immediate. = TRUE)
  }
}
own_forecasts <- nanoparquet::read_parquet(own_forecasts_path)
print(own_forecasts %>% count(disease, data, missing_rounds, horizon) %>% as.data.frame())
