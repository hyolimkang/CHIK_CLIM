# ===========================================================================
# 10b_build_bahia_weekly_demography.R
#
# Bahia equivalent of 10_build_ceara_weekly_demography.R, for the Bahia
# external-replication of the frozen Ceara v4.9 renewal model. Reuses the
# exact same generic interpolation/reconciliation logic (copied verbatim,
# not reimplemented) with UF_CODE=29 (Bahia) and the Bahia births file from
# 09b. Does NOT modify the Ceara script or its outputs.
#
# Output
# ------
#   01_Data/bahia_weekly_demography_2015_2025.rds (+ .csv)
# ===========================================================================

for (p in c("here", "dplyr", "lubridate", "readr", "purrr", "tibble")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(dplyr); library(lubridate); library(readr)
})

interpolate_population <- function(annual_pop, dates) {
  reference_dates <- as.Date(sprintf("%d-07-01", annual_pop$year))
  stats::approx(x = reference_dates, y = annual_pop$pop_total, xout = dates, rule = 2)$y
}

distribute_annual_flow_to_weekly <- function(annual_flow, week_start_vec) {
  daily_rate <- setNames(
    annual_flow$total / ifelse(lubridate::leap_year(annual_flow$year), 366, 365),
    annual_flow$year
  )
  vapply(week_start_vec, function(ws) {
    days <- seq(as.Date(ws), as.Date(ws) + 6, by = "day")
    sum(daily_rate[as.character(lubridate::year(days))])
  }, numeric(1))
}

if (sys.nframe() == 0) {
  UF_CODE <- 29L # Bahia

  demography <- readRDS(here::here("01_Data/ibge_population_projection_uf_2024revision.rds")) |>
    dplyr::filter(uf_code == UF_CODE)

  births_weekly <- readRDS(here::here("01_Data/bahia_sinasc_births_weekly.rds"))
  incomplete <- births_weekly |> dplyr::filter(n_days != 7)
  if (nrow(incomplete) > 0) {
    message(sprintf(
      "[drop] %d incomplete boundary week(s) from bahia_sinasc_births_weekly.rds: %s",
      nrow(incomplete), paste(format(incomplete$week_start), collapse = ", ")
    ))
    births_weekly <- births_weekly |> dplyr::filter(n_days == 7)
  }

  week_start_vec <- births_weekly$week_start
  week_end_vec <- week_start_vec + 6

  boundary_dates <- sort(unique(c(week_start_vec, week_end_vec + 1)))
  n_pop <- interpolate_population(
    demography |> dplyr::select(year, pop_total), boundary_dates
  )
  pop_lookup <- setNames(n_pop, as.character(boundary_dates))

  deaths_annual <- demography |> dplyr::transmute(year, total = deaths_total)
  deaths_weekly <- distribute_annual_flow_to_weekly(deaths_annual, week_start_vec)

  weekly_demography <- tibble::tibble(
    week_start = week_start_vec,
    week_end = week_end_vec,
    N_start = pop_lookup[as.character(week_start_vec)],
    N_end = pop_lookup[as.character(week_end_vec + 1)],
    births = births_weekly$births,
    all_cause_deaths = deaths_weekly
  ) |>
    dplyr::mutate(
      net_population_reconciliation = N_end - N_start - births + all_cause_deaths,
      population_revision = "IBGE Projecoes da Populacao, Revisao 2024 (3a edicao)",
      population_reference_date_rule = "1 July each year; linear interpolation between years",
      birth_source = "SINASC individual records (DTNASC, CODMUNRES), aggregated to epidemiological week",
      death_source = "IBGE Revisao 2024 annual OBT_T (all-cause deaths), constant-daily-rate distribution within year",
      death_method = "annual_total_day_weighted",
      provisional_data_flag = week_start >= as.Date("2025-01-01")
    )

  if (anyNA(weekly_demography)) {
    stop("Weekly demographic contract has missing values -- check input coverage")
  }
  if (any(weekly_demography$N_start <= 0) || any(weekly_demography$N_end <= 0)) {
    stop("Interpolated population must be strictly positive")
  }

  out_path <- here::here("01_Data/bahia_weekly_demography_2015_2025.rds")
  saveRDS(weekly_demography, out_path)
  readr::write_csv(weekly_demography, sub("\\.rds$", ".csv", out_path))

  message(sprintf(
    "[save] %s\n  %s weeks | %s to %s | N range %s to %s",
    out_path, format(nrow(weekly_demography), big.mark = ","),
    format(min(weekly_demography$week_start)), format(max(weekly_demography$week_start)),
    format(round(min(weekly_demography$N_start)), big.mark = ","),
    format(round(max(weekly_demography$N_end)), big.mark = ",")
  ))
  message(sprintf(
    "[check] mean weekly births=%.1f | mean weekly deaths=%.1f | mean |reconciliation|=%.1f",
    mean(weekly_demography$births), mean(weekly_demography$all_cause_deaths),
    mean(abs(weekly_demography$net_population_reconciliation))
  ))
}
