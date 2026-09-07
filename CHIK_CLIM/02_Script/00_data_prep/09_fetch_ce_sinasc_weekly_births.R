# ===========================================================================
# 09_fetch_ce_sinasc_weekly_births.R
#
# Purpose
# -------
# Real, date-level live births for Ceara from SINASC (Sistema de Informacao
# sobre Nascidos Vivos) individual records, aggregated directly to
# epidemiological weeks -- the "best" tier in the roadmap's births/deaths
# hierarchy (docs/model_development_roadmap.md section 4.2), ahead of
# distributing an annual total across days.
#
# Source
# ------
# National, year-level SINASC microdata CSVs (all of Brazil, one row per
# live birth), confirmed reachable at:
#   https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SINASC/csv/SINASC_<year>_csv.zip
# for every year 2014-2025. Each file is ~100-900MB uncompressed, so only
# DTNASC (birth date) and CODMUNRES (6-digit residence municipality code,
# same convention as this repo's muni6) are read, and rows are filtered to
# Ceara (CODMUNRES starting "23") immediately after reading.
#
# The legacy ftp.datasus.gov.br DBC-per-UF files (which would avoid the
# national-file overhead) were unreachable from this environment (connection
# timeout on FTP/HTTP/HTTPS) as of 2026-09-07 -- apparently decommissioned in
# favor of this S3-hosted CSV distribution, the same pattern this repo
# already uses for SINAN (00_data_prep/01_fetch_chik_sinan_brazil.R).
#
# Output
# ------
#   01_Data/ce_sinasc_births_daily.rds   (date, births)
#   01_Data/ce_sinasc_births_weekly.rds  (week_start, births, n_days)
# ===========================================================================

for (p in c("here", "data.table", "dplyr", "lubridate", "readr", "purrr", "tibble")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(data.table); library(dplyr); library(lubridate); library(readr)
})

SOURCE_URL_TEMPLATE <- paste0(
  "https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/",
  "SINASC/csv/SINASC_%d_csv.zip"
)

fetch_ce_births_one_year <- function(year, cache_dir) {
  zip_path <- file.path(cache_dir, sprintf("SINASC_%d_csv.zip", year))
  if (!file.exists(zip_path)) {
    url <- sprintf(SOURCE_URL_TEMPLATE, year)
    message("[download] ", url)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    # These are ~100-900MB national files; R's default 60s download.file()
    # timeout is too short even on a healthy connection.
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = 1800)
    utils::download.file(url, zip_path, mode = "wb", quiet = FALSE)
  } else {
    message("[cache] using existing download: ", zip_path)
  }

  csv_name <- sprintf("SINASC_%d.csv", year)
  extract_cmd <- sprintf('unzip -p "%s" "%s"', zip_path, csv_name)

  # Column naming is not consistent across release years (case, quoting and
  # an optional leading "contador" row-index column all vary), so the
  # header is read first and DTNASC/CODMUNRES are matched case-insensitively
  # rather than assuming a fixed casing. nrows=0 lets fread's own pipe
  # handling read just the header -- unlike a hand-rolled
  # `system("... | head -1", intern = TRUE)`, which hangs on Windows because
  # R's system() shells out via cmd.exe, and unzip never receives a SIGPIPE
  # (a POSIX-only signal) when head closes the pipe early, so it blocks
  # forever trying to write the rest of a ~600-900MB stream into a closed pipe.
  header_names <- names(data.table::fread(cmd = extract_cmd, nrows = 0))
  find_col <- function(pattern) {
    hit <- header_names[toupper(header_names) == pattern]
    if (length(hit) == 0) stop("Column ", pattern, " not found in SINASC ", year, " header")
    hit[1]
  }
  col_dtnasc <- find_col("DTNASC")
  col_codmunres <- find_col("CODMUNRES")

  message("[read] ", csv_name, " (national file; filtering to Ceara only)")
  dt <- data.table::fread(
    cmd = extract_cmd,
    select = c(col_dtnasc, col_codmunres),
    colClasses = "character",
    sep = ";", quote = "\"", encoding = "Latin-1"
  )
  data.table::setnames(dt, c(col_dtnasc, col_codmunres), c("DTNASC", "CODMUNRES"))

  dt <- dt[substr(CODMUNRES, 1, 2) == "23"]
  dt <- dt[!is.na(DTNASC) & nchar(DTNASC) == 8]

  parsed_date <- as.Date(dt$DTNASC, format = "%d%m%Y")
  n_unparsed <- sum(is.na(parsed_date))
  if (n_unparsed > 0) {
    message(sprintf(
      "[warn] %d/%d Ceara SINASC %d records had an unparseable DTNASC and were dropped",
      n_unparsed, nrow(dt), year
    ))
  }

  tibble::tibble(date = parsed_date[!is.na(parsed_date)]) |>
    dplyr::filter(lubridate::year(date) == year) |>
    dplyr::count(date, name = "births")
}

if (sys.nframe() == 0) {
  # Scoped to cover the v2.2 short-period fit window (2015-2019). YEAR_END
  # is one year past the last fit year because a Sunday-start week
  # containing December 2019 dates (e.g. week_start 2019-12-29) needs
  # January-2020 days to be complete -- rerun with a wider range later to
  # extend to the full 2015-2025 period; already-downloaded years are
  # cached and skipped.
  YEAR_START <- 2014L
  YEAR_END <- 2020L

  cache_dir <- here::here("01_Data/sinasc_cache")
  daily <- purrr::map_dfr(YEAR_START:YEAR_END, fetch_ce_births_one_year, cache_dir = cache_dir)

  # Fill in zero-birth calendar days explicitly rather than leaving gaps --
  # SINASC only has rows for days with at least one recorded birth.
  full_days <- tibble::tibble(date = seq(as.Date(sprintf("%d-01-01", YEAR_START)),
                                          as.Date(sprintf("%d-12-31", YEAR_END)), by = "day"))
  daily <- full_days |>
    dplyr::left_join(daily, by = "date") |>
    dplyr::mutate(births = dplyr::if_else(is.na(births), 0L, as.integer(births))) |>
    dplyr::arrange(date)

  weekly <- daily |>
    dplyr::mutate(week_start = lubridate::floor_date(date, unit = "week", week_start = 7)) |>
    dplyr::group_by(week_start) |>
    dplyr::summarise(births = sum(births), n_days = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(week_start)

  daily_path <- here::here("01_Data/ce_sinasc_births_daily.rds")
  weekly_path <- here::here("01_Data/ce_sinasc_births_weekly.rds")
  saveRDS(daily, daily_path)
  saveRDS(weekly, weekly_path)
  readr::write_csv(weekly, sub("\\.rds$", ".csv", weekly_path))

  message(sprintf(
    "[save] %s (%s rows) and %s (%s rows, only weeks with n_days=7 are complete)",
    daily_path, format(nrow(daily), big.mark = ","),
    weekly_path, format(nrow(weekly), big.mark = ",")
  ))
}
