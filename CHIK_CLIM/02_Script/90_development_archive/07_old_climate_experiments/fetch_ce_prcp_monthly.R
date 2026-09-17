# fetch_ce_prcp_monthly.R
# Ceará monthly precipitation time series (2015-2024).
# Source: Open-Meteo Historical API (ERA5-Land), daily sum -> monthly total.
# Location: CE state centroid (Fortaleza region).

fetch_ce_prcp_monthly <- function(
    year_start = 2015,
    year_end   = 2024,
    lat        = -5.2,
    lon        = -39.5,
    out_csv    = "01_Data/ce_prcp_monthly.csv"
) {
  if (!requireNamespace("httr2", quietly = TRUE)) install.packages("httr2")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate")
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
  if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")

  start_date <- sprintf("%04d-01-01", year_start)
  end_date   <- sprintf("%04d-12-31", year_end)

  url <- paste0(
    "https://archive-api.open-meteo.com/v1/archive?",
    "latitude=", lat, "&longitude=", lon,
    "&start_date=", start_date, "&end_date=", end_date,
    "&daily=precipitation_sum&timezone=auto"
  )

  message("Fetching Open-Meteo daily precipitation ...")
  resp <- httr2::request(url) |>
    httr2::req_timeout(180) |>
    httr2::req_perform()
  dat <- httr2::resp_body_json(resp)

  daily <- tibble::tibble(
    date = as.Date(unlist(dat$daily$time)),
    prcp_day = as.numeric(unlist(dat$daily$precipitation_sum))
  )

  prcp_monthly <- daily |>
    dplyr::mutate(
      year  = lubridate::year(date),
      month = lubridate::month(date),
      year_month = lubridate::floor_date(date, "month")
    ) |>
    dplyr::group_by(year, month, year_month) |>
    dplyr::summarise(
      PRCP = sum(prcp_day, na.rm = TRUE),
      n_days = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(year_month)

  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(prcp_monthly, out_csv)
  saveRDS(prcp_monthly, sub("\\.csv$", ".rds", out_csv))

  message("Saved: ", out_csv, " (", nrow(prcp_monthly), " rows)")
  prcp_monthly
}

if (sys.nframe() == 0) {
  print(fetch_ce_prcp_monthly())
}
