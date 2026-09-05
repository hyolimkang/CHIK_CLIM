# ===========================================================================
# 01_fetch_era5land_temperature.R
#
# Purpose
# -------
# Fetch Brazil-wide (bbox) temperature (2m_temperature) from Copernicus
# CDS's `reanalysis-era5-land` (raw hourly) dataset, at just 4 times of day
# (00/06/12/18 UTC). The muni-week aggregation step
# (03_build_national_climate_weekly.R) approximates each day's Tmin/Tmax/Tmean
# from these 4 values (not a true 24-hour min/max, but a standard-enough
# approximation for a national weekly DLNM covariate).
#
# Why raw hourly instead of "derived-era5-land-daily-statistics"
# ------------------------------------------------------------------
# The daily-statistics dataset (which CDS computes daily stats for you) was
# tried first, but timing showed it recomputes on-demand from raw hourly
# data on every request, making it slow (a small area + 1 year took nearly
# 20 minutes). Switching to pulling just a few raw hourly timesteps cut the
# same-sized request (small area + 1 year, 4 timesteps/day) down to a
# measured 4-5 minutes of actual processing.
#
# Why a national grid download instead of per-muni API calls
# --------------------------------------------------------
# An earlier attempt (archive/open_meteo_national_climate_weekly_ratelimited.R)
# called the API once per municipality and hit rate limiting (HTTP 429).
# CDS instead serves a whole grid for a geographic area (bbox) in one go, so
# the number of requests is independent of how many municipalities there
# are — fetching all of Brazil in one request works fine.
#
# CDS allows only one concurrent request per account (confirmed by testing:
# a new request sits in "accepted" indefinitely if another job is still
# active). So multiple years can't be fetched in parallel and must be
# requested sequentially — and if the script is killed mid-run, the CDS-side
# job must be cancelled separately too (killing only the local process
# leaves the server-side job running, which then blocks the next request) —
# DELETE https://cds.climate.copernicus.eu/api/retrieve/v1/jobs/<jobID>
#
# One file per year (unified with 02_fetch_era5land_precip.R's cache unit)
# --------------------------------------------------------------------
# Requests were originally batched 2 years at a time (temp_2014_2015.nc
# style), but at full-Brazil scale a single chunk took over an hour, so this
# was split down to one file per year instead. That gives finer-grained
# resumability and matches the precip cache's naming convention, which is
# easier to manage.
#
# Resumability
# ------------
# Caches to 01_Data/era5_cache/temp_<year>.nc per year; skips if it already exists.
#
# Output
# ------
#   01_Data/era5_cache/temp_<year>.nc  (year = YEAR_START .. YEAR_END)
#     variable t2m, layer order = ascending (day, time) within that year
# ===========================================================================

for (p in c("here", "ecmwfr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({ library(here); library(ecmwfr) })

# Brazil bounding box, [North, West, South, East] — from
# 01_Data/ibge_muni_centroids.rds lon/lat range, padded by ~1 degree.
BRAZIL_AREA <- c(6, -75, -35, -31)
DAILY_TIMES <- c("00:00", "06:00", "12:00", "18:00")

fetch_era5land_temp_year <- function(year, cache_dir) {
  out_file <- file.path(cache_dir, sprintf("temp_%d.nc", year))
  if (file.exists(out_file)) {
    message("[cache] ", basename(out_file), " already exists, skipping")
    return(out_file)
  }

  request <- list(
    dataset_short_name = "reanalysis-era5-land",
    variable            = "2m_temperature",
    year                = as.character(year),
    month               = sprintf("%02d", 1:12),
    day                 = sprintf("%02d", 1:31),
    time                = DAILY_TIMES,
    area                = BRAZIL_AREA,
    data_format         = "netcdf",
    download_format     = "unarchived"
  )

  message(sprintf("[cds] requesting temperature %d ...", year))
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
  # Case panel's first week (2014-12-28) needs late-Dec-2014 daily data,
  # so YEAR_START stays 2014 even though most analysis years are 2015+.
  YEAR_START <- 2014
  YEAR_END   <- 2024
  # ============================================================

  ecmwfr::wf_set_key(key = CDS_KEY)

  cache_dir <- here::here("01_Data/era5_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  for (yr in YEAR_START:YEAR_END) {
    fetch_era5land_temp_year(yr, cache_dir)
  }

  message("[done] all temperature year-files fetched into ", cache_dir)
}
