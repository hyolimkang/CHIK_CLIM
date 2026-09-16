# ===========================================================================
# 09c_fetch_pernambuco_sinasc_weekly_births.R
#
# Pernambuco equivalent of 09b_fetch_bahia_sinasc_weekly_births.R, for the
# Pernambuco external-replication of the frozen v4.9 renewal model.
# SINASC_<year>_csv.zip files are NATIONAL, already cached under
# 01_Data/sinasc_cache/ -- this script reuses that same cache and only
# changes the CODMUNRES filter to "26" (Pernambuco).
#
# Output
# ------
#   01_Data/pernambuco_sinasc_births_daily.rds
#   01_Data/pernambuco_sinasc_births_weekly.rds
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
UF_PREFIX <- "26" # Pernambuco

fetch_uf_births_one_year <- function(year, cache_dir, uf_prefix) {
  zip_path <- file.path(cache_dir, sprintf("SINASC_%d_csv.zip", year))
  if (!file.exists(zip_path)) {
    url <- sprintf(SOURCE_URL_TEMPLATE, year)
    message("[download] ", url)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = 1800)
    utils::download.file(url, zip_path, mode = "wb", quiet = FALSE)
  } else {
    message("[cache] reusing existing national download: ", zip_path)
  }

  csv_name <- sprintf("SINASC_%d.csv", year)
  extract_cmd <- sprintf('unzip -p "%s" "%s"', zip_path, csv_name)

  header_names <- names(data.table::fread(cmd = extract_cmd, nrows = 0))
  find_col <- function(pattern) {
    hit <- header_names[toupper(header_names) == pattern]
    if (length(hit) == 0) stop("Column ", pattern, " not found in SINASC ", year, " header")
    hit[1]
  }
  col_dtnasc <- find_col("DTNASC")
  col_codmunres <- find_col("CODMUNRES")

  message("[read] ", csv_name, " (national file; filtering to UF prefix ", uf_prefix, ")")
  dt <- data.table::fread(
    cmd = extract_cmd,
    select = c(col_dtnasc, col_codmunres),
    colClasses = "character",
    sep = ";", quote = "\"", encoding = "Latin-1"
  )
  data.table::setnames(dt, c(col_dtnasc, col_codmunres), c("DTNASC", "CODMUNRES"))

  dt <- dt[substr(CODMUNRES, 1, 2) == uf_prefix]
  dt <- dt[!is.na(DTNASC) & nchar(DTNASC) == 8]

  parsed_date <- as.Date(dt$DTNASC, format = "%d%m%Y")
  n_unparsed <- sum(is.na(parsed_date))
  if (n_unparsed > 0) {
    message(sprintf(
      "[warn] %d/%d Pernambuco SINASC %d records had an unparseable DTNASC and were dropped",
      n_unparsed, nrow(dt), year
    ))
  }

  tibble::tibble(date = parsed_date[!is.na(parsed_date)]) |>
    dplyr::filter(lubridate::year(date) == year) |>
    dplyr::count(date, name = "births")
}

if (sys.nframe() == 0) {
  YEAR_START <- 2014L
  YEAR_END <- 2025L

  cache_dir <- here::here("01_Data/sinasc_cache")
  daily <- purrr::map_dfr(YEAR_START:YEAR_END, fetch_uf_births_one_year, cache_dir = cache_dir, uf_prefix = UF_PREFIX)

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

  daily_path <- here::here("01_Data/pernambuco_sinasc_births_daily.rds")
  weekly_path <- here::here("01_Data/pernambuco_sinasc_births_weekly.rds")
  saveRDS(daily, daily_path)
  saveRDS(weekly, weekly_path)
  readr::write_csv(weekly, sub("\\.rds$", ".csv", weekly_path))

  message(sprintf(
    "[save] %s (%s rows) and %s (%s rows, only weeks with n_days=7 are complete)",
    daily_path, format(nrow(daily), big.mark = ","),
    weekly_path, format(nrow(weekly), big.mark = ",")
  ))
}
