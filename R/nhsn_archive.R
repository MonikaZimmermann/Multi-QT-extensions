# Pull the full NHSN state-level history (all report versions) via epidatr
# and save it as Parquet. Needs DELPHI_EPIDATA_KEY for sensible rate limits.
#
# All these packages are on CRAN.
suppressPackageStartupMessages({
  library(epidatr)
  library(dplyr)
  library(nanoparquet)
})

out_path <- "nhsn_archive.parquet"
signals <- c(
  "confirmed_admissions_covid_ew",
  "confirmed_admissions_flu_ew",
  "confirmed_admissions_rsv_ew"
)

# One request per signal keeps each response small enough to retry cheaply.
archive <- bind_rows(lapply(signals, function(signal) {
  message("Fetching ", signal)
  epidata_archive(
    source = "nhsn",
    signals = signal,
    geo_type = "state",
    geo_values = "*"
  ) |> mutate(signal = signal)
}))

write_parquet(archive, out_path)
message("Wrote ", nrow(archive), " rows to ", out_path)
