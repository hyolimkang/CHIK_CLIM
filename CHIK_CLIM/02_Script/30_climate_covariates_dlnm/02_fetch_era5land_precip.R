# ===========================================================================
# 02_fetch_era5land_precip.R
#
# Purpose
# -------
# Fetch Brazil-wide precipitation (total_precipitation) from Copernicus
# CDS's `reanalysis-era5-land` (raw hourly) dataset.
#
# Why raw hourly
# -----------------
# `derived-era5-land-daily-statistics` (the "convenience" dataset CDS
# computes on-demand per request) excludes accumulated variables
# (precipitation included) entirely, and timing showed it's slow anyway
# since it's an on-demand computation (script 01 also switched to raw
# hourly for this reason). Precipitation only exists in the raw hourly
# dataset in the first place, so this is the only option.
#
# ERA5-Land's accumulation convention and the "00:00 snapshot" trick
# --------------------------------------------------------------------
# total_precipitation resets its accumulation at UTC 00:00 every day and
# builds up from there. So the value at "day D, 00:00" is the accumulation
# from "day D-1 00:00" to "day D 00:00" — i.e. it's **day D-1's full daily
# total** (confirmed by testing: the 2020-01-01 00:00 value was the right
# order of magnitude for a full day's total on 2019-12-31). This means:
#   - only the "00:00" time needs to be requested (not all 24 hours),
#     cutting data volume to 1/24
#   - to get day D's real precipitation, you need to request "day D+1, 00:00"
#   - this script loops years from YEAR_START through (YEAR_END+1), pulling
#     the 00:00 snapshot for every (month, day) of each year. That way a
#     year boundary like Dec 31 -> Jan 1 is naturally covered by "next
#     year's January request" — the extraction step
#     (03_build_national_climate_weekly.R) just shifts the label back by
#     one day (date - 1) and drops any duplicate dates.
#
# Resumability
# ------------
# Caches to 01_Data/era5_cache/precip_<year>.nc per year; skips if it already exists.
#
# Output
# ------
#   01_Data/era5_cache/precip_<year>.nc  (year = YEAR_START .. YEAR_END+1)
# ===========================================================================

for (p in c("here", "ecmwfr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({ library(here); library(ecmwfr) })

BRAZIL_AREA <- c(6, -75, -35, -31)  # N, W, S, E — same as 01_fetch_era5land_temperature.R

fetch_era5land_precip_year <- function(year, cache_dir) {
  out_file <- file.path(cache_dir, sprintf("precip_%d.nc", year))
  if (file.exists(out_file)) {
    message("[cache] ", basename(out_file), " already exists, skipping")
    return(out_file)
  }

  request <- list(
    dataset_short_name = "reanalysis-era5-land",
    variable            = "total_precipitation",
    year                = as.character(year),
    month               = sprintf("%02d", 1:12),
    day                 = sprintf("%02d", 1:31),
    time                = "00:00",
    area                = BRAZIL_AREA,
    data_format         = "netcdf",
    download_format     = "unarchived"
  )

  message(sprintf("[cds] requesting precip %d (00:00 snapshots) ...", year))
  f <- ecmwfr::wf_request(
    request  = request,
    transfer = TRUE,
    path     = cache_dir,
    time_out = 3600 * 3
  )
  file.rename(f, out_file)
  message("[cds] saved ", out_file)
  out_file
}

if (sys.nframe() == 0) {
  # ============================================================
  CDS_KEY    <- "bd2e4ae0-4af4-4559-8def-c501fb4d673f"
  YEAR_START <- 2014
  YEAR_END   <- 2024
  # ============================================================

  ecmwfr::wf_set_key(key = CDS_KEY)

  cache_dir <- here::here("01_Data/era5_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # +1 year so December 31 of YEAR_END is captured via Jan-1-next-year's 00:00 snapshot.
  for (yr in YEAR_START:(YEAR_END + 1)) {
    fetch_era5land_precip_year(yr, cache_dir)
  }

  message("[done] all precip year files fetched into ", cache_dir)
}
