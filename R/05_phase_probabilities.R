# Phase probabilities for the methods phase_angle and phase_hidden (R/03; notes, "The angle from
# a phase probability", "Weather, holidays and the hidden-phase model"). For the base model of
# multiqt_settings() and its pooled locations, from what is known at each cutoff:
#   inputs   the reports (NHSN versions; before November 2024 the HHS versions for
#            windowed_seasonal), the last three corrected by their median revision at their age;
#            the base's projected growth (median at horizon 3 against 0); a seasonal timing
#            prior; NSSP emergency-department visits (Epidata, every version); the holidays;
#            specific humidity (Open-Meteo: ERA5 and the forecasts archived as issued); onset;
#   p_bias   the probability that the report falls below the base median: an online logistic
#            regression per disease, pooled over the pooled locations and the horizons, refitted
#            at every cutoff on the rows labelled by then (phase_angle);
#   hidden   the hidden-phase model (rise, overshoot weeks 1-4, tail): a hidden Markov chain
#            whose flu prior is a logistic hazard of the peak fitted on WHO/NREVSS laboratory
#            test positivity (2014/15-2019/20, 2022/23, 2023/24) with the season week, holidays
#            and humidity, and whose COVID prior is the same hazard fitted at each cutoff on the
#            COVID rises known then; data: the growth of the reports and of NSSP visits, Gaussian
#            per phase with a holiday shift. Predictive probabilities of each target week at its
#            cutoff, and the updated (smoothed) probabilities of earlier weeks at later cutoffs
#            (phase_hidden).
# Downloads once, into the project root: nssp_vintages_epidata.parquet (Delphi's Epidata, needs
# external/delphi/library from R/00) and weather_openmeteo.parquet (Open-Meteo, no key; raw
# responses cached in external/open_meteo/). NREVSS: Delphi's flu auxiliary data
# (external/delphi/aux_data/flusion_data, extracted by R/00).
# Run from the project root after R/01: [MULTIQT_BASE_MODEL=windowed_seasonal] Rscript
# R/05_phase_probabilities.R (about 15 minutes). Writes <output_dir>/data/phase_probabilities.rds.
source("R/03_functions.R")
lines <- readLines("R/01_data_overview.R")
start <- grep("^find_peaks <- function", lines)
eval(parse(text = lines[start:(start + which(lines[(start + 1):length(lines)] == "}")[1])]))
settings <- multiqt_settings()
prepared <- read_prepared(settings)
params <- method_parameters()
diseases <- c("covid", "flu")
states <- settings$pooled_locations
cores <- 10L
incomplete <- as.Date(c("2024-05-04", "2024-11-02"))   # R/01: incomplete reporting.
clamp <- function(p, lo = .02) pmin(pmax(p, lo), 1 - lo)
classes <- c("rise", "over", "tail")
phase_class <- function(rel) ifelse(rel <= 0, "rise", ifelse(rel <= 4, "over", "tail"))
season_week <- function(d) {                    # Weeks since the last 1 August, 0 to 52.
  y <- as.integer(format(d, "%Y")) - (as.integer(format(d, "%m")) < 8)
  pmin(as.integer(as.numeric(d - as.Date(paste0(y, "-08-01"))) %/% 7), 52L)
}

# 0. DOWNLOADS (ONCE) -------------------------------------------------------------------------
# NSSP: share of emergency-department visits with COVID / influenza, weekly, state level, every
# version (issue) published, for the pooled locations.
nssp_file <- "nssp_vintages_epidata.parquet"
if (!file.exists(nssp_file)) {
  .libPaths(c("external/delphi/library", .libPaths()))
  signals <- c(covid = "pct_ed_visits_covid", flu = "pct_ed_visits_influenza")
  abbr <- setNames(tolower(state.abb[match(c("California", "Texas", "New York", "Pennsylvania",
    "Illinois", "Ohio", "Georgia", "North Carolina", "Michigan", "New Jersey"), state.name)]),
    c("06", "48", "36", "42", "17", "39", "13", "37", "26", "34"))
  stopifnot(setequal(names(abbr), states))
  nssp <- do.call(rbind, lapply(names(signals), function(d) do.call(rbind, lapply(names(abbr),
    function(l) {
      x <- epidatr::pub_covidcast(source = "nssp", signals = signals[[d]], geo_type = "state",
        time_type = "week", geo_values = abbr[[l]], time_values = "*", issues = "*")
      data.frame(disease = d, location = l, week_end = as.Date(x$time_value) + 6,
        issue = as.Date(x$issue), value = x$value)
    }))))
  arrow::write_parquet(nssp, nssp_file)
}
# Weather: daily mean temperature and dew point of three large cities per state (ERA5,
# archive-api), and the hourly forecasts issued 1 to 7 days before each valid time
# (previous-runs-api, from January 2024), averaged per day; specific humidity (g/kg) from the
# dew point (Magnus formula, 1013.25 hPa); per state the population-weighted mean of its cities.
weather_file <- "weather_openmeteo.parquet"
if (!file.exists(weather_file)) {
  cities <- read.table(text = "
location city lat lon weight timezone
06 LosAngeles 34.05 -118.24 13.0 America/Los_Angeles
06 SanJose 37.34 -121.89 7.7 America/Los_Angeles
06 SanDiego 32.72 -117.16 3.3 America/Los_Angeles
48 Houston 29.76 -95.37 7.5 America/Chicago
48 Dallas 32.78 -96.80 8.1 America/Chicago
48 SanAntonio 29.42 -98.49 2.7 America/Chicago
36 NewYork 40.71 -74.01 12.0 America/New_York
36 Buffalo 42.89 -78.88 1.1 America/New_York
36 Rochester 43.16 -77.61 1.1 America/New_York
42 Philadelphia 39.95 -75.17 4.5 America/New_York
42 Pittsburgh 40.44 -80.00 2.4 America/New_York
42 Harrisburg 40.27 -76.88 0.6 America/New_York
17 Chicago 41.88 -87.63 8.5 America/Chicago
17 Peoria 40.69 -89.59 0.4 America/Chicago
17 Springfield 39.78 -89.65 0.2 America/Chicago
39 Columbus 39.96 -83.00 2.1 America/New_York
39 Cleveland 41.50 -81.69 2.1 America/New_York
39 Cincinnati 39.10 -84.51 1.5 America/New_York
13 Atlanta 33.75 -84.39 6.2 America/New_York
13 Augusta 33.47 -81.97 0.6 America/New_York
13 Savannah 32.08 -81.09 0.4 America/New_York
37 Charlotte 35.23 -80.84 2.0 America/New_York
37 Raleigh 35.78 -78.64 1.5 America/New_York
37 Greensboro 36.07 -79.79 0.8 America/New_York
26 Detroit 42.33 -83.05 4.3 America/Detroit
26 GrandRapids 42.96 -85.67 1.1 America/Detroit
26 Lansing 42.73 -84.56 0.5 America/Detroit
34 Newark 40.74 -74.17 6.0 America/New_York
34 Camden 39.93 -75.12 1.3 America/New_York
34 Trenton 40.22 -74.76 0.4 America/New_York
", header = TRUE, colClasses = c(location = "character"))
  raw <- "external/open_meteo"
  dir.create(raw, recursive = TRUE, showWarnings = FALSE)
  specific_humidity <- function(dew) {
    e <- 6.112 * exp(17.62 * dew / (243.12 + dew))
    622 * e / (1013.25 - .378 * e)
  }
  get_json <- function(url, file) {
    if (!file.exists(file)) for (attempt in 1:5) {
      if (tryCatch({ download.file(url, file, quiet = TRUE); TRUE }, error = function(e) FALSE)) break
      Sys.sleep(20 * attempt)
    }
    jsonlite::fromJSON(file)
  }
  end_date <- as.character(Sys.Date() - 7)
  tz <- function(c) utils::URLencode(c$timezone, reserved = TRUE)
  observed <- do.call(rbind, lapply(seq_len(nrow(cities)), function(i) {
    c <- cities[i, ]
    x <- get_json(sprintf(paste0("https://archive-api.open-meteo.com/v1/archive?latitude=%s",
      "&longitude=%s&start_date=2014-01-01&end_date=%s&daily=temperature_2m_mean,",
      "dew_point_2m_mean&timezone=%s"), c$lat, c$lon, end_date, tz(c)),
      file.path(raw, paste0("observed_", c$city, ".json")))$daily
    data.frame(location = c$location, weight = c$weight, kind = "observed", lead = NA_integer_,
      valid = as.Date(x$time), temperature = x$temperature_2m_mean,
      humidity = specific_humidity(x$dew_point_2m_mean))
  }))
  vars <- c(paste0("temperature_2m_previous_day", 1:7), paste0("dew_point_2m_previous_day", 1:7))
  chunks <- data.frame(start = c("2024-01-01", "2025-01-01", "2026-01-01"),
    end = c("2024-12-31", "2025-12-31", end_date))
  forecast <- do.call(rbind, lapply(seq_len(nrow(cities)), function(i) {
    c <- cities[i, ]
    do.call(rbind, lapply(seq_len(nrow(chunks)), function(k) {
      h <- get_json(sprintf(paste0("https://previous-runs-api.open-meteo.com/v1/forecast?",
        "latitude=%s&longitude=%s&start_date=%s&end_date=%s&hourly=%s&timezone=%s"), c$lat,
        c$lon, chunks$start[k], chunks$end[k], paste(vars, collapse = ","), tz(c)),
        file.path(raw, paste0("forecast_", c$city, "_", k, ".json")))$hourly
      day <- as.Date(substr(h$time, 1, 10))
      do.call(rbind, lapply(1:7, function(l) {
        a <- aggregate(cbind(temperature = h[[paste0("temperature_2m_previous_day", l)]],
          humidity = specific_humidity(h[[paste0("dew_point_2m_previous_day", l)]])) ~ day,
          FUN = mean, na.action = na.omit)
        data.frame(location = c$location, weight = c$weight, kind = "forecast", lead = l,
          valid = a$day, temperature = a$temperature, humidity = a$humidity)
      }))
    }))
  }))
  w <- rbind(observed, forecast)
  w$lead[is.na(w$lead)] <- 0L
  num <- aggregate(cbind(temperature = temperature * weight, humidity = humidity * weight,
    weight) ~ location + kind + lead + valid, w, sum)
  num$temperature <- num$temperature / num$weight; num$humidity <- num$humidity / num$weight
  num$issue <- num$valid - num$lead
  arrow::write_parquet(num[c("location", "kind", "lead", "valid", "issue", "temperature",
    "humidity")], weather_file)
}
weather_all <- as.data.frame(arrow::read_parquet(weather_file))
weather_all$location <- as.character(weather_all$location)
nssp_all <- as.data.frame(arrow::read_parquet(nssp_file))
nssp_all$location <- as.character(nssp_all$location)

# 1. REPORTS, NSSP AND WEATHER AS OF A CUTOFF ---------------------------------------------------------
vint <- as.data.frame(arrow::read_parquet("nhsn_vintages_epidata.parquet"))
# Before NHSN's first version (19 November 2024), the HHS report versions (R/00, as R/01
# does for this base): reports known at the cutoffs of 2023/24.
if (settings$base_model == "windowed_seasonal") vint <- rbind(as.data.frame(
  arrow::read_parquet("hhs_vintages_epidata.parquet"))[names(vint)], vint)
vint <- vint[vint$location %in% states & valid_report(vint$value), ]
vint <- vint[order(vint$disease, vint$location, vint$target_end_date, vint$version), ]
reports <- split(vint[c("target_end_date", "version", "value")],
  paste(vint$disease, vint$location))
# The reports of one disease and location known at `date`: the latest version per target.
known_at <- function(r, date) {
  r <- r[r$version <= date, ]
  r <- r[!duplicated(r$target_end_date, fromLast = TRUE), ]
  r[order(r$target_end_date), c("target_end_date", "value")]
}
# Value of targets `target` as of the dates `date` (one each).
value_at <- function(r, target, date) {
  idx <- split(seq_len(nrow(r)), as.numeric(r$target_end_date))
  vapply(seq_along(target), function(i) {
    j <- idx[[as.character(as.numeric(target[i]))]]
    if (is.null(j)) return(NA_real_)
    k <- findInterval(as.numeric(date[i]), as.numeric(r$version[j]))
    if (k == 0L) NA_real_ else r$value[j][k]
  }, numeric(1))
}
nssp <- nssp_all
weather <- weather_all[weather_all$kind == "observed", ]
weather$date <- weather$valid
weather_fc <- weather_all[weather_all$kind == "forecast", ]
weather <- split(weather, weather$location)
weather_fc <- split(weather_fc, weather_fc$location)
# Climatology of humidity per location and day of the year, 2014-2023, +-7 days.
climate <- lapply(weather, function(w) {
  w <- w[w$date < as.Date("2024-01-01"), ]
  doy <- as.integer(format(w$date, "%j"))
  vapply(1:366, function(d) mean(w$humidity[abs(((doy - d + 183) %% 366) - 183) <= 7]), numeric(1))
})
nssp <- nssp[is.finite(nssp$value) & nssp$value > 0, ]
nssp <- split(nssp[c("week_end", "issue", "value")], paste(nssp$disease, nssp$location))

# 2. INPUTS PER DISEASE, LOCATION AND CUTOFF ---------------------------------------------
# Median revision log1p(now) - log1p(at age a) of the reports of the 10 states whose first
# report is at least 9 weeks old at the cutoff (last 52 weeks), for ages 0, 1, 2 weeks after
# the newest report: the correction of the corrected labels, in its simplest form.
revision_shift <- function(disease, cutoff, lag_days) {
  vapply(0:2, function(a) {
    d <- unlist(lapply(states, function(l) {
      r <- reports[[paste(disease, l)]]
      targets <- unique(r$target_end_date[r$target_end_date <= cutoff - 63 &
        r$target_end_date > cutoff - 63 - 364 & r$target_end_date > incomplete[2]])
      if (!length(targets)) return(NULL)
      now <- value_at(r, targets, rep(cutoff, length(targets)))
      young <- value_at(r, targets, targets + lag_days + 7 * a)
      log1p(now) - log1p(young)
    }))
    d <- d[is.finite(d)]
    if (length(d) < 8L) 0 else median(d)
  }, numeric(1))
}
# Seasonal prior at a cutoff: per week of the year, the share of the 10 states' earlier weeks
# (at least 26 weeks before the cutoff, complete reporting) that lay after their nearest
# peak, +-1 week of the year.
season_table <- function(disease, cutoff) {
  after <- do.call(rbind, lapply(states, function(l) {
    s <- known_at(reports[[paste(disease, l)]], cutoff)
    s <- s[s$target_end_date <= cutoff - 182 & (s$target_end_date < incomplete[1] |
      s$target_end_date > incomplete[2]), ]
    if (nrow(s) < 52L) return(NULL)
    peaks <- find_peaks(s$value, 8L, .3, .1, 3L, scale = max(s$value))$row
    if (!length(peaks)) return(NULL)
    nearest <- peaks[apply(abs(outer(seq_len(nrow(s)), peaks, "-")), 1, which.min)]
    data.frame(week = as.integer(strftime(s$target_end_date, "%V")),
      after = as.numeric(seq_len(nrow(s)) > nearest))
  }))
  vapply(1:53, function(w) {
    near <- ((after$week - w + 26) %% 53) - 26
    x <- after$after[abs(near) <= 1]
    if (length(x) < 10L) .5 else mean(x)
  }, numeric(1))
}
cutoff_features <- function(disease, location, cutoff) {
  r <- known_at(reports[[paste(disease, location)]], cutoff)
  last <- nrow(r)
  lag_days <- as.numeric(cutoff - r$target_end_date[last])
  shift <- revision_shift_cache[[paste(disease, cutoff)]]
  y <- log1p(r$value[last - 0:3]) + c(shift, shift[3])
  n <- nssp[[paste(disease, location)]]
  n <- n[n$issue + 6 <= cutoff, ]
  n <- n[!duplicated(n$week_end, fromLast = TRUE), ]
  n <- n[order(n$week_end), ]
  nl <- if (nrow(n) >= 26L) log(n$value[nrow(n) - 0:3]) else rep(NA_real_, 4)
  w <- weather[[location]]
  past <- w$humidity[w$date >= cutoff - 14 & w$date < cutoff]
  days <- as.integer(format(cutoff - 14:1, "%j"))
  hum_level <- mean(past)
  f <- weather_fc[[location]]
  f <- f$humidity[f$issue == cutoff & f$lead %in% 1:7]
  recent <- tail(r, 20)
  trough <- which.min(recent$value)
  data.frame(disease, location, cutoff,
    hum_level, hum_anom = hum_level - mean(climate[[location]][days]),
    hum_fc = if (length(f) == 7L) mean(f) - hum_level else 0,
    trough_date = recent$target_end_date[trough],
    rising = as.numeric(tail(recent$value, 1) >= 1.5 * max(recent$value[trough], 1)),
    g1 = y[1] - y[2], g2 = (y[1] - y[3]) / 2, curv = (y[1] - y[2]) - (y[2] - y[3]),
    ng1 = nl[1] - nl[2], ng2 = (nl[1] - nl[3]) / 2, ncurv = (nl[1] - nl[2]) - (nl[2] - nl[3]),
    nratio = if (nrow(n) >= 26L) n$value[nrow(n)] / max(tail(n$value, 26)) else NA_real_,
    nssp_age = if (nrow(n)) as.numeric(cutoff - n$week_end[nrow(n)]) else NA_real_)
}
cutoffs <- lapply(setNames(nm = diseases), function(d)
  sort(unique(do.call(c, lapply(states, function(l)
    prepared$cases[[paste(d, l, 0, sep = "_")]]$as_of)))))
jobs <- do.call(rbind, lapply(diseases, function(d) data.frame(disease = d,
  cutoff = cutoffs[[d]])))
lag_of <- function(d, cutoff) {
  r <- known_at(reports[[paste(d, states[1])]], cutoff)
  as.numeric(cutoff - max(r$target_end_date))
}
shifts <- parallel::mclapply(seq_len(nrow(jobs)), function(i) with(jobs[i, ],
  revision_shift(disease, cutoff, lag_of(disease, cutoff))), mc.cores = cores)
revision_shift_cache <- setNames(shifts, paste(jobs$disease, jobs$cutoff))
seasons <- parallel::mclapply(seq_len(nrow(jobs)), function(i) with(jobs[i, ],
  season_table(disease, cutoff)), mc.cores = cores)
# Rise durations (weeks from the minimum of the 20 weeks before a peak to the peak) of the
# peaks of the reports known at the cutoff, at least 26 weeks old, 10 states.
rise_durations <- function(disease, cutoff) unlist(lapply(states, function(l) {
  s <- known_at(reports[[paste(disease, l)]], cutoff)
  s <- s[s$target_end_date <= cutoff - 182 & (s$target_end_date < incomplete[1] |
    s$target_end_date > incomplete[2]), ]
  if (nrow(s) < 52L) return(NULL)
  peaks <- find_peaks(s$value, 8L, .3, .1, 3L, scale = max(s$value))$row
  vapply(peaks[peaks > 20], function(i) {
    before <- (i - 20):i
    as.numeric(s$target_end_date[i] - s$target_end_date[before[which.min(s$value[before])]]) / 7
  }, numeric(1))
}))
durations <- parallel::mclapply(seq_len(nrow(jobs)), function(i) with(jobs[i, ],
  rise_durations(disease, cutoff)), mc.cores = cores)
duration_cache <- setNames(durations, paste(jobs$disease, jobs$cutoff))
season_cache <- setNames(seasons, paste(jobs$disease, jobs$cutoff))
features <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), function(i)
  do.call(rbind, lapply(states, function(l) cutoff_features(jobs$disease[i], l,
    jobs$cutoff[i]))), mc.cores = cores))

# 3. ROWS: INPUTS AND LABELS --------------------------------------------
med <- 12L
case_rows <- function(disease, location, horizon) {
  case <- get_case(prepared, disease, location, horizon, 0L, FALSE)
  stopifnot(case$alpha[med] == .5)
  h <- case$context_forecasts
  h <- h[!is.na(h$as_of) & is.finite(h$median) & h$median >= 0, ]
  m <- function(ref, hz) log1p(h$median[match(paste(ref, hz), paste(h$reference_date, h$horizon))])
  r <- reports[[paste(disease, location)]]
  label_date <- case$target_date + 11                     # One week after the first report.
  now <- value_at(r, case$target_date, label_date)
  before <- value_at(r, case$target_date - 7, label_date)
  week <- as.integer(strftime(case$target_date, "%V"))
  season <- vapply(seq_along(week), function(i) {
    table <- season_cache[[paste(disease, case$as_of[i])]]
    if (is.null(table)) .5 else table[week[i]]
  }, numeric(1))
  x <- data.frame(disease, location, horizon, row = seq_len(nrow(case$base)),
    reference_date = case$reference_date, target_date = case$target_date, cutoff = case$as_of,
    period = case$period, filled = case$filled, relative_week = case$relative_week,
    season = qlogis(clamp(season)), base_growth = m(case$reference_date, 3) -
      m(case$reference_date, 0), label_date,
    fall = as.numeric(now < before), bias = as.numeric(now < case$base[, med]),
    z = as.numeric(case$truth_latest < case$base[, med]),
    p_errors = decline_probability(case, params$phase_rows),
    p_season = clamp(season),
    p_oracle = ifelse(is.finite(case$relative_week), as.numeric(case$relative_week > 0), .5))
  x <- merge(x, features, by = c("disease", "location", "cutoff"), all.x = TRUE, sort = FALSE)
  # Holidays of the target week (its seven days up to target_date).
  t <- x$target_date
  in_break <- function(d) { md <- format(d, "%m-%d"); md >= "12-22" | md <= "01-03" }
  x$hol_break <- rowMeans(sapply(0:6, function(k) in_break(t - k)))
  jan3 <- as.Date(paste0(format(t, "%Y"), "-01-03"))
  x$hol_post <- as.numeric(as.numeric(t - jan3) %in% 1:21)
  thanksgiving <- vapply(as.integer(format(t, "%Y")), function(y) {
    nov <- as.Date(paste0(y, "-11-01")) + 0:29
    as.numeric(nov[format(nov, "%u") == "4"][4])
  }, numeric(1))
  x$hol_thanks <- as.numeric((as.numeric(t) - thanksgiving) %in% 0:6)
  # Onset: the share of earlier rises shorter than the rise to the target week so far.
  x$onset <- vapply(seq_len(nrow(x)), function(i) {
    D <- duration_cache[[paste(disease, x$cutoff[i])]]
    if (!isTRUE(x$rising[i] == 1) || length(D) < 5L) return(0)
    weeks <- as.numeric(x$target_date[i] - x$trough_date[i]) / 7
    qlogis(clamp((sum(D < weeks) + .5) / (length(D) + 1)))
  }, numeric(1))
  x
}
grid <- expand.grid(disease = diseases, location = states, horizon = settings$horizons,
  stringsAsFactors = FALSE)
rows <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), function(i)
  with(grid[i, ], case_rows(disease, location, horizon)), mc.cores = cores))
rows <- rows[order(rows$disease, rows$location, rows$horizon, rows$row), ]

# 4. P_BIAS: THE PROBABILITY THAT THE BASE IS TOO HIGH ----------------------------------------
# Online logistic regression per disease, pooled over the pooled locations and the horizons:
# refitted at each cutoff on the rows labelled by then (the report one week after its first,
# label_date <= cutoff; at least 40 rows with both labels); 1/2 before. Inputs: reports, base,
# NSSP, holidays, humidity and onset (the set "all" of the experiments of 16 November).
bias_terms <- c("season", "g1", "curv", "base_growth", "horizon", "horizon:g1", "horizon:curv",
  "ng1", "ncurv", "nratio", "horizon:ng1", "horizon:ncurv", "hol_break", "hol_post", "hol_thanks",
  "hum_level", "hum_anom", "hum_fc", "onset", "rising")
vars <- unique(unlist(strsplit(bias_terms, ":")))
rows$p_bias <- .5
for (d in diseases) {
  x <- rows$disease == d
  complete <- x & complete.cases(rows[vars])
  for (cut in sort(unique(rows$cutoff[x]))) {
    train <- rows[complete & rows$label_date <= cut & is.finite(rows$bias), ]
    now <- which(complete & rows$cutoff == cut)
    if (!length(now) || nrow(train) < 40L || length(unique(train$bias)) < 2L) next
    fit <- suppressWarnings(glm(reformulate(bias_terms, "bias"), binomial, train))
    rows$p_bias[now] <- clamp(predict(fit, rows[now, ], type = "response"))
  }
}

# 5. THE HIDDEN-PHASE MODEL ----------------------------------------------------------------
# Prior: WHO/NREVSS laboratory test positivity of the pooled states (New York with New York
# City; New Jersey has no series), peaks per season (the maximum of the 3-week mean, at least
# 5% positive), midpoints between consecutive seasons' peaks.
nrevss_dir <- "external/delphi/aux_data/flusion_data"
names_of <- c("06" = "California", "48" = "Texas", "36" = "New York", "42" = "Pennsylvania",
  "17" = "Illinois", "39" = "Ohio", "13" = "Georgia", "37" = "North Carolina",
  "26" = "Michigan", "34" = "New Jersey")
mmwr_end <- function(year, week) {
  jan4 <- as.Date(paste0(year, "-01-04"))
  first_sunday <- jan4 - as.integer(format(jan4, "%w"))
  first_sunday + 7 * (week - 1) + 6
}
num <- function(x) suppressWarnings(as.numeric(x))
old <- read.csv(file.path(nrevss_dir, "WHO_NREVSS_Combined_prior_to_2015_16_State.csv"), skip = 1,
  check.names = FALSE)
old$positives <- num(old[["TOTAL SPECIMENS"]]) * num(old[["PERCENT POSITIVE"]]) / 100
new <- read.csv(file.path(nrevss_dir, "WHO_NREVSS_Clinical_Labs_State.csv"), skip = 1,
  check.names = FALSE)
new$positives <- num(new[["TOTAL A"]]) + num(new[["TOTAL B"]])
lab <- rbind(old[c("REGION", "YEAR", "WEEK", "TOTAL SPECIMENS", "positives")],
  new[c("REGION", "YEAR", "WEEK", "TOTAL SPECIMENS", "positives")])
names(lab) <- c("region", "year", "week", "specimens", "positives")
lab$specimens <- num(lab$specimens)
lab$region[lab$region == "New York City"] <- "New York"
lab <- aggregate(cbind(specimens, positives) ~ region + year + week, lab, sum)
lab <- lab[lab$region %in% names_of, ]
lab$date <- mmwr_end(lab$year, lab$week)
lab <- lab[order(lab$region, lab$date), ]
lab$positivity <- ave(lab$positives / pmax(lab$specimens, 1), lab$region, FUN = function(x) {
  s <- as.numeric(stats::filter(x, rep(1 / 3, 3), sides = 2)); ifelse(is.na(s), x, s)
})
lab$season <- as.integer(format(lab$date, "%Y")) - (as.integer(format(lab$date, "%m")) < 8)
seasons_used <- c(2010:2019, 2022, 2023)
peaks <- do.call(rbind, lapply(split(lab[lab$season %in% seasons_used, ],
  list(lab$region[lab$season %in% seasons_used], lab$season[lab$season %in% seasons_used]),
  drop = TRUE), function(x) {
    if (max(x$positivity, na.rm = TRUE) < .05) return(NULL)
    data.frame(region = x$region[1], season = x$season[1],
      peak = x$date[which.max(x$positivity)])
  }))
peaks <- peaks[order(peaks$region, peaks$peak), ]
# Midpoints between consecutive peaks of a state in consecutive seasons.
peaks$next_peak <- ave(as.numeric(peaks$peak), peaks$region, FUN = function(p) c(p[-1], NA))
peaks$next_season <- ave(peaks$season, peaks$region, FUN = function(s) c(s[-1], NA))
peaks$midpoint <- as.Date(ifelse(!is.na(peaks$next_season) & peaks$next_season == peaks$season + 1,
  floor((as.numeric(peaks$peak) + peaks$next_peak) / 2), NA), origin = "1970-01-01")
peaks$prev_mid <- as.Date(ave(as.numeric(peaks$midpoint), peaks$region, FUN = function(m) c(NA, m[-length(m)])),
  origin = "1970-01-01")
# Hazard by season week: events / at risk, pooled over +-1 week.
hazard <- function(starts, events) {
  e <- tabulate(season_week(events) + 1L, 53)
  r <- numeric(53)
  for (i in seq_along(starts)) {
    a <- season_week(starts[i]); b <- season_week(events[i])
    if (b >= a) r[(a:b) + 1L] <- r[(a:b) + 1L] + 1 else r[c(a:52, 0:b) + 1L] <- r[c(a:52, 0:b) + 1L] + 1
  }
  h <- vapply(1:53, function(w) {
    near <- ((0:52) - (w - 1) + 26) %% 53 - 26
    sum(e[abs(near) <= 1]) / max(sum(r[abs(near) <= 1]), 1)
  }, numeric(1))
  h[r == 0] <- .5
  pmin(pmax(h, .005), .95)
}
ok <- !is.na(peaks$prev_mid)
flu_peak_hazard <- hazard(peaks$prev_mid[ok], peaks$peak[ok])
ok <- !is.na(peaks$midpoint)
flu_mid_hazard <- hazard(peaks$peak[ok] + 35, peaks$midpoint[ok])
# The reports known at `date`, without the weeks of incomplete reporting (the hidden-phase
# model's series).
known_complete <- function(r, date) {
  r <- known_at(r, date)
  r[r$target_end_date < incomplete[1] | r$target_end_date > incomplete[2], ]
}
# Weather and holidays.
obs <- weather                                  # Observed, per location (section 1).
fcs <- weather_fc                               # Forecasts as issued, per location.
holiday_covariates <- function(t) {
  in_break <- function(d) { md <- format(d, "%m-%d"); md >= "12-22" | md <= "01-03" }
  jan3 <- as.Date(paste0(format(t, "%Y"), "-01-03"))
  thanksgiving <- vapply(as.integer(format(t, "%Y")), function(y) {
    nov <- as.Date(paste0(y, "-11-01")) + 0:29
    as.numeric(nov[format(nov, "%u") == "4"][4])
  }, numeric(1))
  data.frame(hol_break = rowMeans(sapply(0:6, function(k) in_break(t - k))),
    hol_post = as.numeric(as.numeric(t - jan3) %in% 1:21),
    hol_thanks = as.numeric((as.numeric(t) - thanksgiving) %in% 0:6))
}
# Weekly humidity (14 days up to each week end) as known at `cutoff`: observed up to the
# cutoff, the forecast issued at the cutoff for 7 days, climatology after; cutoff NULL: observed.
weekly_humidity <- function(location, weeks, cutoff = NULL) {
  days <- seq(min(weeks) - 13, max(weeks), by = 1)
  o <- obs[[location]]
  v <- o$humidity[match(days, o$date)]
  clim <- climate[[location]][as.integer(format(days, "%j"))]
  if (!is.null(cutoff)) {
    v[days > cutoff] <- NA
    f <- fcs[[location]]
    f <- f[f$issue == cutoff, ]
    fc <- f$humidity[match(days, f$valid)]
    use <- days > cutoff & !is.na(fc)
    v[use] <- fc[use]
  }
  v[is.na(v)] <- clim[is.na(v)]
  cs <- cumsum(v); cc <- cumsum(v - clim)
  end <- match(weeks, days)
  roll <- function(x) (x[end] - c(0, x)[end - 13]) / 14
  data.frame(hum_level = roll(cs), hum_anom = roll(cc))
}
# Flu: discrete-time logistic hazard of the peak on the NREVSS seasons from 2014/15.
risk <- do.call(rbind, lapply(which(!is.na(peaks$prev_mid) & peaks$season >= 2014), function(i) {
  first <- peaks$prev_mid[i] + (6 - as.integer(format(peaks$prev_mid[i], "%w"))) %% 7
  wk <- seq(first, peaks$peak[i], by = 7)
  loc <- names(names_of)[names_of == peaks$region[i]]
  data.frame(location = loc, week = wk, event = as.numeric(wk == peaks$peak[i]))
}))
risk$sw <- season_week(risk$week)
risk <- risk[risk$sw >= 8, ]          # No peak before week 8 (late September) in the data.
risk <- cbind(risk, holiday_covariates(risk$week))
risk <- cbind(risk, do.call(rbind, lapply(seq_len(nrow(risk)), function(i)
  weekly_humidity(risk$location[i], risk$week[i]))))
hazard_formula <- event ~ splines::ns(sw, df = 4) + hol_break + hol_post + hol_thanks +
  hum_level + hum_anom
flu_hazard_fit <- glm(hazard_formula, binomial, risk)
flu_hazard_calendar <- glm(event ~ splines::ns(sw, df = 4), binomial, risk)
# Lag of the hospital-admission peaks behind the positivity peaks (2022/23, 2023/24).
latest <- function(disease, l) {
  r <- reports[[paste(disease, l)]]
  r <- r[!duplicated(r$target_end_date, fromLast = TRUE), ]
  r[order(r$target_end_date), ]
}
lags <- unlist(lapply(which(peaks$season %in% 2022:2023), function(i) {
  loc <- names(names_of)[names_of == peaks$region[i]]
  r <- latest("flu", loc)
  r <- r[r$target_end_date >= as.Date(paste0(peaks$season[i], "-08-01")) &
    r$target_end_date < as.Date(paste0(peaks$season[i] + 1, "-08-01")), ]
  if (!nrow(r)) return(NULL)
  as.numeric(r$target_end_date[which.max(r$value)] - peaks$peak[i]) / 7
}))
peak_lag <- if (length(lags)) round(median(lags)) else 0
hazard_summary <- capture.output({
cat("NREVSS logistic hazard of the flu peak, seasons from 2014/15 (rows", nrow(risk), ", peaks",
  sum(risk$event), ")\n")
print(round(summary(flu_hazard_fit)$coefficients, 3))
cat("AIC with holidays and humidity", round(AIC(flu_hazard_fit), 1), "| calendar only",
  round(AIC(flu_hazard_calendar), 1), "\n")
cat("Admission peak minus positivity peak (weeks), 2022/23-2023/24:", lags, "-> shift", peak_lag, "\n")
})
writeLines(hazard_summary)
flu_peak_hazard_at <- function(location, weeks, cutoff) {
  shifted <- weeks - 7 * peak_lag
  nd <- data.frame(sw = season_week(shifted), holiday_covariates(shifted),
    weekly_humidity(location, shifted, cutoff))
  h <- pmin(pmax(predict(flu_hazard_fit, nd, type = "response"), .002), .95)
  ifelse(nd$sw < 8, .002, h)
}
series_at <- function(disease, location, cutoff, shift) {
  s <- known_complete(reports[[paste(disease, location)]], cutoff)
  s$x <- log1p(s$value)
  n <- nrow(s)
  s$x[n - 0:2] <- s$x[n - 0:2] + shift
  s$growth <- c(NA, diff(s$x))
  s$growth[c(FALSE, diff(as.numeric(s$target_end_date)) != 7)] <- NA   # Across gaps.
  n_ <- nssp[[paste(disease, location)]]
  n_ <- n_[n_$issue + 6 <= cutoff, ]
  n_ <- n_[!duplicated(n_$week_end, fromLast = TRUE), ]
  ln <- log(n_$value[match(s$target_end_date, n_$week_end)])
  lp <- log(n_$value[match(s$target_end_date - 7, n_$week_end)])
  s$nssp_growth <- ln - lp
  cbind(s, holiday_covariates(s$target_end_date))
}
# Labels from the peaks confirmed by the cutoff (at least 8 weeks old), for weeks at least 12
# weeks old; peaks of every state.
labelled <- function(s, cutoff) {
  pk <- s$target_end_date[find_peaks(s$value, 8L, .3, .1, 3L, scale = max(s$value))$row]
  pk <- pk[pk <= cutoff - 56]
  if (!length(pk)) return(NULL)
  w <- s[s$target_end_date <= cutoff - 84 & is.finite(s$growth), ]
  rel <- vapply(as.numeric(w$target_end_date), function(t) {
    r <- (t - as.numeric(pk)) / 7; r[which.min(abs(r))]
  }, numeric(1))
  data.frame(growth = w$growth, nssp_growth = w$nssp_growth, hol_break = w$hol_break,
    hol_post = w$hol_post, phase = factor(phase_class(rel), classes), week = w$target_end_date,
    value = w$value)
}
# 6 hidden states: rise, over1..over4, tail. Transition matrix for the step into week w.
transition <- function(h_peak, h_mid) {
  T <- matrix(0, 6, 6)
  T[1, 1] <- 1 - h_peak; T[1, 2] <- h_peak
  T[2, 3] <- T[3, 4] <- T[4, 5] <- T[5, 6] <- 1
  T[6, 6] <- 1 - h_mid; T[6, 1] <- h_mid
  T
}
collapse <- function(a) cbind(rise = a[, 1], over = rowSums(a[, 2:5, drop = FALSE]), tail = a[, 6])
phase_probabilities <- function(disease, cutoff) {
  lag_days <- as.numeric(cutoff - max(known_complete(reports[[paste(disease, states[1])]], cutoff)$target_end_date))
  shift <- revision_shift(disease, cutoff, lag_days)
  ser <- lapply(setNames(nm = states), function(l) series_at(disease, l, cutoff, shift))
  # Emissions per phase from the labelled weeks of all 10 states: Gaussian, mean per phase
  # plus a holiday shift, sd per phase.
  lab_rows <- do.call(rbind, lapply(ser, labelled, cutoff = cutoff))
  if (is.null(lab_rows) || any(table(lab_rows$phase) < 20L)) return(NULL)
  emission <- function(column) {
    x <- lab_rows[is.finite(lab_rows[[column]]), ]
    if (any(table(x$phase) < 20L)) return(NULL)
    fit <- lm(reformulate(c("0", "phase", "hol_break", "hol_post"), column), x)
    sdv <- pmax(tapply(residuals(fit), fit$model$phase, sd), .05)[classes]
    list(fit = fit, sd = sdv)
  }
  em <- list(growth = emission("growth"), nssp_growth = emission("nssp_growth"))
  em <- Filter(Negate(is.null), em)
  # COVID: a logistic hazard with holidays and humidity on the rises known then, else constant.
  if (disease == "covid") {
    d <- do.call(rbind, lapply(states, function(l) {
      s <- ser[[l]]
      s <- s[s$target_end_date <= cutoff - 182, ]
      if (nrow(s) < 52L) return(NULL)
      pk <- s$target_end_date[find_peaks(s$value, 8L, .3, .1, 3L, scale = max(s$value))$row]
      if (length(pk) < 2L) return(NULL)
      mid <- as.Date(floor((as.numeric(pk[-1]) + as.numeric(pk[-length(pk)])) / 2), origin = "1970-01-01")
      do.call(rbind, lapply(seq_along(mid), function(j) {
        wk <- s$target_end_date[s$target_end_date > mid[j] & s$target_end_date <= pk[j + 1]]
        if (!length(wk)) return(NULL)
        data.frame(location = l, week = wk, event = as.numeric(wk == pk[j + 1]),
          tail_weeks = (as.numeric(mid[j]) - as.numeric(pk[j])) / 7 - 5)
      }))
    }))
    h_tail <- if (is.null(d) || !nrow(d)) 1 / 12 else
      1 / max(mean(unique(d[c("location", "tail_weeks")])$tail_weeks), 2)
    covid_fit <- NULL
    if (!is.null(d) && sum(d$event) >= 8L) {
      d <- cbind(d, holiday_covariates(d$week))
      d <- cbind(d, do.call(rbind, lapply(seq_len(nrow(d)), function(i)
        weekly_humidity(d$location[i], d$week[i]))))
      covid_fit <- tryCatch(suppressWarnings(glm(event ~ hol_break + hol_post + hol_thanks +
        hum_level + hum_anom, binomial, d)), error = function(e) NULL)
    }
    h_rise_const <- if (is.null(d) || !nrow(d)) 1 / 12 else sum(d$event) / nrow(d)
  }
  do.call(rbind, lapply(states, function(l) {
    s <- ser[[l]]
    w <- s$target_end_date
    n <- length(w)
    weeks_all <- c(w, w[n] + 7 * (1:6))
    h_peak <- if (disease == "flu") flu_peak_hazard_at(l, weeks_all, cutoff) else
      if (is.null(covid_fit)) rep(h_rise_const, length(weeks_all)) else {
        nd <- data.frame(holiday_covariates(weeks_all), weekly_humidity(l, weeks_all, cutoff))
        pmin(pmax(predict(covid_fit, nd, type = "response"), .002), .95)
      }
    h_mid <- if (disease == "flu") flu_mid_hazard[season_week(weeks_all) + 1L] else
      rep(h_tail, length(weeks_all))
    T_at <- function(date) { i <- match(date, weeks_all); transition(h_peak[i], h_mid[i]) }
    # Likelihood of week i: report growth and NSSP growth, each Gaussian per phase.
    P3 <- matrix(1, n, 3)
    for (column in names(em)) {
      cf <- coef(em[[column]]$fit)
      M <- outer(rep(1, n), cf[paste0("phase", classes)]) +
        (cf[["hol_break"]] * s$hol_break + cf[["hol_post"]] * s$hol_post)
      g <- s[[column]]
      D <- dnorm(g, M, matrix(em[[column]]$sd, n, 3, byrow = TRUE))
      D[!is.finite(g), ] <- 1
      P3 <- P3 * D
    }
    L <- cbind(P3[, 1], P3[, 2], P3[, 2], P3[, 2], P3[, 2], P3[, 3])
    lik <- function(i) L[i, ]
    alpha <- matrix(NA_real_, n, 6); Ts <- vector("list", n)
    a <- rep(1 / 6, 6) * lik(1); alpha[1, ] <- a / sum(a)
    for (i in 2:n) {
      Ts[[i]] <- T_at(w[i])
      a <- drop(alpha[i - 1, ] %*% Ts[[i]]) * lik(i)
      alpha[i, ] <- a / sum(a)
    }
    beta <- matrix(1, n, 6)
    for (i in (n - 1):1) {
      b <- drop(Ts[[i + 1]] %*% (lik(i + 1) * beta[i + 1, ]))
      beta[i, ] <- b / sum(b)
    }
    smooth <- alpha * beta; smooth <- smooth / rowSums(smooth)
    # Predictive for the 6 weeks after the latest report.
    ahead <- matrix(NA_real_, 6, 6); a <- alpha[n, ]
    for (m in 1:6) { a <- drop(a %*% T_at(w[n] + 7 * m)); ahead[m, ] <- a }
    keep <- w >= cutoff - 7 * 80
    rbind(data.frame(disease, location = l, cutoff, week = w[keep], kind = "smoothed",
        collapse(smooth[keep, , drop = FALSE])),
      data.frame(disease, location = l, cutoff, week = w[n] + 7 * (1:6), kind = "ahead",
        collapse(ahead)))
  }))
}
jobs <- unique(rows[c("disease", "cutoff")])
hidden <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), function(i)
  phase_probabilities(jobs$disease[i], jobs$cutoff[i]), mc.cores = cores))

# 6. SAVE -----------------------------------------------------------------------------------
# p_bias per case row; hidden: per disease, location and cutoff the phase probabilities of the
# weeks up to 80 weeks before the cutoff (smoothed) and of the 6 weeks after the latest report
# (predictive). get_case() (R/03) attaches them to the cases of the pooled locations.
out <- list(base_model = settings$base_model, created = Sys.time(), hazard = hazard_summary,
  bias = rows[c("disease", "location", "horizon", "row", "cutoff", "p_bias")],
  hidden = hidden[c("disease", "location", "cutoff", "week", "kind", classes)])
file <- file.path(settings$output_dir, "data", "phase_probabilities.rds")
saveRDS(out, file)
cat("Wrote", file, ":", nrow(out$bias), "rows,", nrow(out$hidden), "phase probabilities\n")
