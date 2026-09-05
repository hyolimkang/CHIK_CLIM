# ===========================================================================
# 03_build_national_climate_weekly.R
#
# Purpose
# -------
# From the Brazil-wide ERA5-Land grid (raw hourly) that scripts 01/02
# fetched into 01_Data/era5_cache/*.nc, extract an area-weighted mean per
# municipality **polygon** and aggregate to muni6 x week_start (Sunday
# start). The final output shape matches the earlier Open-Meteo version
# (archive/open_meteo_national_climate_weekly_ratelimited.R), so
# 04_build_dlnm_panel.R can read this file as-is.
#
# Why polygon area-weighted mean instead of centroid
# ----------------------------------------------------
# Sampling a single centroid point per municipality can fail to represent
# the actual local climate for large, irregularly-shaped municipalities
# (especially in the Amazon). So `exactextractr::exact_extract(raster,
# polygons, "mean")` is used to weight grid cells by how much of their area
# overlaps the municipality boundary (coverage fraction):
#   T_{m,t} = sum(A_mg * T_{g,t}) / sum(A_mg)
# exact_extract's "mean" already uses this coverage-fraction weighting by
# default (no separate weight raster like population needed), so it matches
# this formula directly. Polygons come from
# 01_Data/ibge_muni_polygons.rds, saved by
# 00_data_prep/06_fetch_muni_centroids.R (EPSG:4326, same CRS as the
# ERA5-Land grid).
#
# Temperature: approximating min/max/mean from 4 timesteps (00/06/12/18 UTC)
# ----------------------------------------------------------
# Since 01_fetch_era5land_temperature.R only pulls 4 timesteps/day from raw
# hourly data, each NetCDF file has the same date appearing 4 times (at
# different times). So after extracting per muni6, the timestamp (date +
# time) is floored to date and grouped by (muni6, date) to compute that
# day's min/max/mean here directly — CDS isn't pre-computing this for us
# (which is exactly why the raw hourly dataset was used instead of a daily
# statistics one).
#
# Handling NetCDF's valid_time
# ------------------------
# The NetCDF files CDS returns encode time using CF-standard units like
# "seconds since 1970-01-01". terra attaches this raw number directly to
# the layer name, so instead of relying on that, ncdf4 is used to read the
# valid_time variable's units attribute directly and convert it to a real
# timestamp (parse_cf_time()).
#
# Precipitation's -1-day shift
# ---------------
# As explained in 02_fetch_era5land_precip.R, the "day D 00:00" snapshot
# value is actually day (D-1)'s full daily total precipitation. The label
# is shifted back here (date - 1). The same date can end up extracted from
# two different yearly files (at year boundaries), so distinct() drops the
# duplicate (same reanalysis value either way, so it doesn't matter which
# copy is kept).
#
# Output
# ------
#   01_Data/chik_muni_week_climate_2015_2024.rds
#     columns: muni6, week_start, PRCP, Tmin, Tmax, Tmean, n_days
# ===========================================================================

for (p in c("here", "terra", "ncdf4", "exactextractr", "sf", "dplyr", "lubridate", "tibble", "readr", "purrr", "tidyr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(terra); library(ncdf4); library(exactextractr); library(sf)
  library(dplyr); library(lubridate); library(tibble); library(readr); library(purrr); library(tidyr)
})

# CF time units strings vary in whether a time-of-day is present at all
# (observed: "seconds since 1970-01-01" for reanalysis-era5-land, "days
# since 2020-01-01 00:00:00" for the daily-statistics dataset) — so the
# unit and the optional time-of-day are parsed as two independent regexes
# rather than one combined pattern.
parse_cf_time <- function(nc, varname = "valid_time") {
  vt <- ncdf4::ncvar_get(nc, varname)
  units_str <- ncdf4::ncatt_get(nc, varname, "units")$value

  m_unit <- regmatches(units_str, regexec("^(days|hours|seconds) since (\\d{4}-\\d{2}-\\d{2})", units_str))[[1]]
  if (length(m_unit) == 0) stop("Could not parse CF time units: ", units_str)
  unit_secs <- c(days = 86400, hours = 3600, seconds = 1)[[m_unit[2]]]

  m_time <- regmatches(units_str, regexec("(\\d{2}:\\d{2}:\\d{2})", units_str))[[1]]
  time_part <- if (length(m_time) > 0) m_time[1] else "00:00:00"

  origin <- as.POSIXct(paste(m_unit[3], time_part), tz = "UTC")
  origin + vt * unit_secs
}

# Returns one row per (muni6, timestamp) — NOT collapsed to date yet, since
# temperature files have 4 timestamps/day (00/06/12/18 UTC) that must stay
# distinct until the daily min/max/mean aggregation step.
#
# exact_extract() is called ONCE on the whole multi-layer raster (not once
# per layer/timestep) — the expensive part (polygon/grid coverage-fraction
# overlap) is computed once per municipality and reused across all layers,
# which is what makes this fast enough for ~1,460 layers x 5,570 polygons.
extract_nc_zonal <- function(nc_path, varname, polygons) {
  nc <- ncdf4::nc_open(nc_path)
  timestamps <- parse_cf_time(nc)
  ncdf4::nc_close(nc)

  ras <- terra::rast(nc_path)
  vals <- exactextractr::exact_extract(ras, polygons, "mean", progress = FALSE)

  colnames(vals) <- as.character(as.numeric(timestamps))  # unique even within the same day
  vals$muni6 <- polygons$muni6

  vals |>
    tidyr::pivot_longer(-muni6, names_to = "timestamp", values_to = varname) |>
    dplyr::mutate(timestamp = as.POSIXct(as.numeric(timestamp), origin = "1970-01-01", tz = "UTC"))
}

build_temperature_daily <- function(cache_dir, polygons) {
  files <- list.files(cache_dir, pattern = "^temp_\\d{4}\\.nc$", full.names = TRUE)
  if (length(files) == 0) stop("No cached temperature files found.")
  message("[extract] temperature: ", length(files), " file(s)")

  # A handful of tiny/island municipalities (e.g. Fernando de Noronha) have
  # no overlapping ERA5-Land land pixels at all, so exact_extract() returns
  # NA for every timestamp there. max()/min() on an all-NA vector emit a
  # warning and silently return -Inf (not NA) — must be converted back to
  # NA explicitly, or these rows would pass downstream NA-checks while
  # holding garbage -Inf values. Affected munis are dropped in
  # 04_build_dlnm_panel.R the same way the unmatched SINAN codes are.
  purrr::map_dfr(files, extract_nc_zonal, varname = "t2m", polygons = polygons) |>
    dplyr::mutate(date = as.Date(timestamp)) |>
    dplyr::group_by(muni6, date) |>
    dplyr::summarise(
      Tmean = mean(t2m, na.rm = TRUE) - 273.15,
      Tmax  = suppressWarnings(max(t2m, na.rm = TRUE)) - 273.15,
      Tmin  = suppressWarnings(min(t2m, na.rm = TRUE)) - 273.15,
      n_obs = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(c(Tmean, Tmax, Tmin), ~ dplyr::if_else(is.finite(.x), .x, NA_real_)))
}

build_precip_daily <- function(cache_dir, polygons) {
  files <- list.files(cache_dir, pattern = "^precip_\\d{4}\\.nc$", full.names = TRUE)
  if (length(files) == 0) stop("No cached precip files found.")
  message("[extract] PRCP: ", length(files), " file(s)")

  purrr::map_dfr(files, extract_nc_zonal, varname = "tp", polygons = polygons) |>
    dplyr::mutate(
      date = as.Date(timestamp) - 1,   # 00:00 snapshot on `date` is the total for `date - 1`
      PRCP = tp * 1000                 # m -> mm
    ) |>
    dplyr::distinct(muni6, date, .keep_all = TRUE) |>
    dplyr::select(muni6, date, PRCP)
}

build_national_climate_weekly <- function(cache_dir, polygons, week_start_vec) {
  temp_daily   <- build_temperature_daily(cache_dir, polygons)
  precip_daily <- build_precip_daily(cache_dir, polygons)

  daily <- dplyr::inner_join(temp_daily, precip_daily, by = c("muni6", "date")) |>
    dplyr::mutate(week_start = lubridate::floor_date(date, unit = "week", week_start = 7))

  week_start_vec <- as.Date(unique(week_start_vec))

  # Same -Inf-on-all-NA-input issue as build_temperature_daily() (see its
  # comment) — the 3 no-land-coverage munis have all-NA daily Tmin/Tmax for
  # every day, so min()/max() here need the same is.finite() cleanup.
  #
  # PRCP has a quieter version of the same bug: sum(x, na.rm=TRUE) on an
  # all-NA week doesn't warn or return Inf, it silently returns 0 — which
  # would misrepresent "no data" as "confirmed zero rainfall". Guard with
  # an explicit all(is.na(...)) check instead of relying on is.finite()
  # after the fact.
  daily |>
    dplyr::filter(week_start %in% week_start_vec) |>
    dplyr::group_by(muni6, week_start) |>
    dplyr::summarise(
      PRCP   = if (all(is.na(PRCP))) NA_real_ else sum(PRCP, na.rm = TRUE),
      Tmin   = suppressWarnings(min(Tmin, na.rm = TRUE)),
      Tmax   = suppressWarnings(max(Tmax, na.rm = TRUE)),
      Tmean  = mean(Tmean, na.rm = TRUE),
      n_days = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(c(Tmin, Tmax), ~ dplyr::if_else(is.finite(.x), .x, NA_real_))) |>
    dplyr::arrange(muni6, week_start)
}

if (sys.nframe() == 0) {
  cache_dir <- here::here("01_Data/era5_cache")
  polygons  <- readRDS(here::here("01_Data/ibge_muni_polygons.rds"))
  panel     <- readRDS(here::here("01_Data/chik_brazil_muni_week_2015_2024.rds"))

  weekly_climate <- build_national_climate_weekly(cache_dir, polygons, panel$week_start)

  out_path <- here::here("01_Data/chik_muni_week_climate_2015_2024.rds")
  saveRDS(weekly_climate, out_path)
  readr::write_csv(weekly_climate, sub("\\.rds$", ".csv", out_path))

  message(sprintf(
    "[save] %s (%d rows, %d municipalities, %d weeks)",
    out_path, nrow(weekly_climate),
    dplyr::n_distinct(weekly_climate$muni6),
    dplyr::n_distinct(weekly_climate$week_start)
  ))
}
