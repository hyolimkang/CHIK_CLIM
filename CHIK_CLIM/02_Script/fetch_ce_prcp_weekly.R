# fetch_ce_prcp_weekly.R
# Weekly precipitation aligned to ce_weekly$week_start (SINAN epi weeks).
# Source: Open-Meteo daily precipitation_sum -> sum within each week_start.

fetch_ce_prcp_weekly <- function(
    week_start_vec,
    lat = -5.2,
    lon = -39.5,
    out_csv = "01_Data/ce_prcp_weekly.csv"
) {
  if (!requireNamespace("httr2", quietly = TRUE)) install.packages("httr2")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
  if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")

  week_start_vec <- sort(unique(as.Date(week_start_vec)))
  start_date <- min(week_start_vec)
  end_date   <- max(week_start_vec) + 6L   # cover last epi week

  url <- paste0(
    "https://archive-api.open-meteo.com/v1/archive?",
    "latitude=", lat, "&longitude=", lon,
    "&start_date=", start_date, "&end_date=", end_date,
    "&daily=precipitation_sum&timezone=auto"
  )

  message("Fetching Open-Meteo daily precipitation ...")
  dat <- httr2::request(url) |>
    httr2::req_timeout(180) |>
    httr2::req_perform() |>
    httr2::resp_body_json()

  daily <- tibble::tibble(
    date = as.Date(unlist(dat$daily$time)),
    prcp_day = as.numeric(unlist(dat$daily$precipitation_sum))
  )

  weeks <- tibble::tibble(week_start = week_start_vec) |>
    dplyr::mutate(week_end = week_start + 6L)

  prcp_weekly <- weeks |>
    dplyr::rowwise() |>
    dplyr::mutate(
      PRCP = sum(
        daily$prcp_day[daily$date >= week_start & daily$date <= week_end],
        na.rm = TRUE
      ),
      n_days = sum(daily$date >= week_start & daily$date <= week_end)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(week_start, PRCP, n_days)

  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(prcp_weekly, out_csv)
  saveRDS(prcp_weekly, sub("\\.csv$", ".rds", out_csv))

  message("Saved: ", out_csv, " (", nrow(prcp_weekly), " rows)")
  prcp_weekly
}

if (sys.nframe() == 0) {
  ce_path <- "01_Data/ce_weekly_2015_2024.rds"
  if (!file.exists(ce_path)) {
    stop("Need ce_weekly RDS at ", ce_path)
  }
  ce <- readRDS(ce_path)
  print(fetch_ce_prcp_weekly(ce$week_start))
}
