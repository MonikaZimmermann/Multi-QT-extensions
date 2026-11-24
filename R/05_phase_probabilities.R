# Bias and phase probabilities for the methods phase_angle and phase_hidden (R/03; notes, "Bias
# and phase probabilities"). For the base of multiqt_settings() and its pooled locations, at every
# issuance cutoff, from what is known then (the functions are in R/03, section "PHASE AND BIAS
# PROBABILITIES"):
#   bias    p_t, the probability that the report of the target week falls below the base median:
#           one binomial GAM per disease (G0: standardized base growth, horizon and holidays with a
#           ridge penalty, smooths of the year fraction and of the newest hospital growth g1, a
#           state random intercept), refitted at every cutoff on the rows issued before it whose
#           label has arrived; label = the report as of target-week end + 11 days (COVID) or + 4
#           days (flu); else the smoothed label frequency;
#   hidden  P_t^k, the probabilities of rise, overshoot (weeks 1-4 after a peak) and tail: a
#           hidden semi-Markov chain over 84 states, reset at reporting gaps, smoothed within each
#           segment. COVID: exit hazards from the durations of the European waves (ERVISS) and
#           conditional-growth GAMs of the hospital and NSSP growth (O2) at the NSSP lead (0, 1 or
#           2 weeks) chosen per 13-week block by the predictive score of earlier blocks. Flu: exit-
#           hazard GAMs of duration and calendar on the completed episodes of the hospital reports
#           and NREVSS, and phase-mean Gaussians of the hospital and NSSP growth at lead 0 (O1).
# Inputs: R/01's prepared data (base, horizons, cutoffs); NHSN report versions
# (nhsn_vintages_epidata.parquet; for windowed_seasonal also the HHS versions,
# hhs_vintages_epidata.parquet, before NHSN's first one), reviewed exclusions masked (R/03); NSSP
# versions (nssp_vintages_epidata.parquet: Delphi's Epidata, downloaded once below, needs
# external/delphi/library from R/00); WHO/NREVSS laboratory files (Delphi's flu auxiliary data,
# external/delphi/aux_data/flusion_data, extracted by R/00; a snapshot); ECDC ERVISS hospital
# admissions (erviss_covid_admissions.parquet, downloaded once below; a snapshot).
# Run from the project root after R/01: [MULTIQT_BASE_MODEL=windowed_seasonal] [MULTIQT_CORES=n]
# Rscript R/05_phase_probabilities.R. Writes <output_dir>/data/phase_probabilities.rds (the
# tables bias and hidden read by phase_inputs() in R/03, and their details and records) and
# <output_dir>/data/phase_probability_models.rds (the fitted models at the season ends).
source("R/03_functions.R")
settings <- multiqt_settings()
prepared <- read_prepared(settings)
states <- settings$pooled_locations
cores <- as.integer(Sys.getenv("MULTIQT_CORES", 10L))
nhsn_file <- "nhsn_vintages_epidata.parquet"
hhs_file <- "hhs_vintages_epidata.parquet"
nssp_file <- "nssp_vintages_epidata.parquet"
erviss_file <- "erviss_covid_admissions.parquet"
nrevss_dir <- "external/delphi/aux_data/flusion_data"

# 0. DOWNLOADS (ONCE) -------------------------------------------------------------------------
# NSSP: share of emergency-department visits with COVID / influenza, weekly, state level, every
# version (issue) published, for the pooled locations.
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
# COVID waves of other countries: ECDC ERVISS (EU/EEA), weekly SARS-CoV-2 hospital admissions,
# all ages (a snapshot of the repository's current file, not versioned).
if (!file.exists(erviss_file)) {
  raw <- tempfile(fileext = ".csv")
  download.file(paste0("https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/",
    "main/data/nonSentinelSeverity.csv"), raw, quiet = TRUE)
  x <- read.csv(raw)
  x <- x[x$pathogen == "SARS-CoV-2" & x$indicator == "hospitaladmissions" & x$age == "total", ]
  arrow::write_parquet(data.frame(country = x$countryname, yearweek = x$yearweek,
    value = x$value), erviss_file)
}

# 1. DATA AS USED -------------------------------------------------------------------------------
# Hospital report versions of the pooled states; before NHSN's first version (19 November 2024)
# the HHS report versions for windowed_seasonal (R/00, as R/01 does for this base). Reviewed
# exclusions masked (value NA, the version kept: never replaced by an older one).
vint <- as.data.frame(arrow::read_parquet(nhsn_file))
with_hhs <- settings$base_model == "windowed_seasonal"
if (with_hhs) vint <- rbind(as.data.frame(arrow::read_parquet(hhs_file))[names(vint)], vint)
vint <- filter_reports(vint[vint$location %in% states, ])
vint <- vint[valid_report(vint$value) | vint$report_excluded, ]
vint <- vint[order(vint$disease, vint$location, vint$target_end_date, vint$version), ]
nssp <- as.data.frame(arrow::read_parquet(nssp_file))
nssp$location <- as.character(nssp$location)
nssp <- nssp[nssp$location %in% states & is.finite(nssp$value) & nssp$value > 0, ]
pp <- list(states = states, incomplete = as.Date(c("2024-05-04", "2024-11-02")),
  reports = split(vint[c("target_end_date", "version", "value")], paste(vint$disease,
    vint$location)),
  nssp = split(nssp[c("week_end", "issue", "value")], paste(nssp$disease, nssp$location)),
  nrevss = pp_nrevss_series(nrevss_dir, states), erviss = pp_erviss_series(erviss_file))
input_files <- c(prepared = file.path(settings$output_dir, "data", "peak_aligned.rds"),
  nhsn = nhsn_file, hhs = if (with_hhs) hhs_file, nssp = nssp_file, erviss = erviss_file,
  nrevss_old = file.path(nrevss_dir, "WHO_NREVSS_Combined_prior_to_2015_16_State.csv"),
  nrevss_new = file.path(nrevss_dir, "WHO_NREVSS_Clinical_Labs_State.csv"))
sources <- data.frame(input = names(input_files), file = unname(input_files),
  md5 = unname(tools::md5sum(input_files)), bytes = file.size(input_files),
  modified = file.mtime(input_files))
coverage <- data.frame(
  input = c("nhsn", "nssp", "nrevss", "erviss"),
  first_week = c(min(vint$target_end_date), min(nssp$week_end),
    min(do.call(c, lapply(pp$nrevss, function(x) min(x$target_end_date)))),
    min(do.call(c, lapply(pp$erviss, function(x) min(x$target_end_date))))),
  last_week = c(max(vint$target_end_date), max(nssp$week_end),
    max(do.call(c, lapply(pp$nrevss, function(x) max(x$target_end_date)))),
    max(do.call(c, lapply(pp$erviss, function(x) max(x$target_end_date))))),
  last_version = c(max(vint$version), max(nssp$issue), NA, NA),
  series = c(length(pp$reports), length(pp$nssp), length(pp$nrevss), length(pp$erviss)))
season_end <- function(d) as.Date(paste0(as.integer(format(d, "%Y")) -
  (format(d, "%m") < "08") + 1, "-08-01"))

# 2. BIAS PROBABILITY -----------------------------------------------------------------------------
rows <- pp_bias_rows(pp, prepared, settings$horizons, cores = cores)$rows
jobs <- unique(rows[c("disease", "cutoff")])
jobs <- jobs[order(jobs$disease, jobs$cutoff), ]
jobs$keep <- FALSE
for (d in unique(jobs$disease)) {
  k <- jobs$disease == d
  jobs$keep[k] <- as.numeric(jobs$cutoff[k]) %in% tapply(as.numeric(jobs$cutoff[k]),
    season_end(jobs$cutoff[k]), max)
}
fits <- parallel::mclapply(seq_len(nrow(jobs)), function(i) pp_bias_fit_at(rows,
  jobs$disease[i], jobs$cutoff[i], states, keep = jobs$keep[i]), mc.cores = cores,
  mc.preschedule = FALSE)
failed <- vapply(fits, inherits, logical(1), "try-error")
if (any(failed)) stop(fits[[which(failed)[1]]])
rows$p_raw <- NA_real_; rows$p_bias <- NA_real_
rows$source <- NA_character_; rows$fallback_reason <- NA_character_
for (f in fits) {
  rows$p_raw[f$rows] <- f$p_raw; rows$p_bias[f$rows] <- f$p
  rows$source[f$rows] <- f$source; rows$fallback_reason[f$rows] <- f$reason
}
stopifnot(all(is.finite(rows$p_bias)))
bias_records <- do.call(rbind, lapply(fits, `[[`, "record"))
cat("Bias:", nrow(rows), "rows,", nrow(jobs), "fits;", sum(bias_records$fitted), "GAMs, sources",
  paste(names(table(rows$source)), table(rows$source), collapse = ", "), "\n")

# 3. PHASE PROBABILITIES ----------------------------------------------------------------------------
# COVID: O2 at leads 0, 1, 2 at every cutoff; flu: O1 at lead 0.
variants_of <- list(covid = pp_settings$covid_variants, flu = pp_settings$flu_variant)
hid <- parallel::mclapply(seq_len(nrow(jobs)), function(i) pp_hidden_at(pp, jobs$disease[i],
  jobs$cutoff[i], variants_of[[jobs$disease[i]]], keep = jobs$keep[i]), mc.cores = cores,
  mc.preschedule = FALSE)
failed <- vapply(hid, inherits, logical(1), "try-error")
if (any(failed)) stop(hid[[which(failed)[1]]])
# The COVID lead of each 13-week block (blocks from the first cutoff of the base), from the
# predictive scores of the weeks reported before the block; held for the whole block.
first <- min(rows$cutoff)
scores <- pp_predictive_scores(do.call(rbind, lapply(hid, `[[`, "scores")))
leads <- pp_select_leads(scores[scores$disease == "covid", ], "O2", first,
  sort(unique(jobs$cutoff)))
hidden <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(i) {
  if (jobs$disease[i] == "flu") return(cbind(hid[[i]]$probabilities[[pp_settings$flu_variant]],
    lead = 0L))
  L <- leads$lead[leads$block == pp_block_of(jobs$cutoff[i], first)]
  cbind(hid[[i]]$probabilities[[paste0("O2_L", L)]], lead = L)
}))
stopifnot(max(abs(rowSums(hidden[pp_classes]) - 1)) < 1e-9)
hidden_records <- do.call(rbind, lapply(hid, `[[`, "records"))
transitions <- do.call(rbind, lapply(hid, `[[`, "transition"))
cat("Phases:", nrow(hidden), "probabilities; COVID leads by block:",
  paste(leads$lead, collapse = ""), "\n")

# 4. SAVE ---------------------------------------------------------------------------------------------
# bias: p_bias (operational, clipped) per case row; hidden: per disease, location and cutoff the
# probabilities of the weeks up to 80 weeks back (smoothed) and of the 6 weeks after the newest
# hospital report (nssp-informed or predicted). get_case() (R/03) attaches them to the cases of the
# pooled locations; phase_inputs() floors each component at .005 without renormalizing.
out <- list(base_model = settings$base_model, created = Sys.time(), sources = sources,
  coverage = coverage, settings = pp_settings,
  bias = rows[c("disease", "location", "horizon", "row", "cutoff", "p_bias")],
  bias_details = rows[c("disease", "location", "horizon", "row", "reference_date", "target_date",
    "cutoff", "period", "label_date", "label_report", "base_median", "bias", "bias_mature",
    "base_growth", "g1", "year_phase", "hol_break", "hol_post", "hol_thanks", "p_raw", "p_bias",
    "source", "fallback_reason")],
  bias_records = bias_records,
  hidden = hidden[c("disease", "location", "cutoff", "week", "kind", pp_classes, "lead")],
  hidden_records = hidden_records, transitions = transitions, lead_selection = leads,
  predictive_scores = scores)
file <- file.path(settings$output_dir, "data", "phase_probabilities.rds")
saveRDS(out, file)
models <- c(lapply(which(jobs$keep), function(i) list(kind = "bias", disease = jobs$disease[i],
  cutoff = jobs$cutoff[i], object = fits[[i]]$object)),
  lapply(which(jobs$keep), function(i) list(kind = "observation", disease = jobs$disease[i],
    cutoff = jobs$cutoff[i], object = hid[[i]]$models)))
saveRDS(models, file.path(settings$output_dir, "data", "phase_probability_models.rds"))
cat("Wrote", file, ":", nrow(out$bias), "bias rows,", nrow(out$hidden), "phase probabilities\n")
