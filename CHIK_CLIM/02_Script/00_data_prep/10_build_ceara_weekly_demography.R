# ===========================================================================
# 10_build_ceara_weekly_demography.R
#
# Purpose
# -------
# Build the weekly demographic input contract described in
# docs/model_development_roadmap.md section 4.5 for Ceara:
#   week_start, week_end, N_start, N_end, births, all_cause_deaths,
#   net_population_reconciliation, plus source/vintage metadata.
#
# Population stock (N)
# ---------------------
# IBGE's 2024-revision UF population series (already fetched by
# 07_fetch_ibge_population_projection_uf.R) has an official reference date
# of 1 July each year (confirmed against IBGE's published "Estimativas da
# Populacao" methodological notes, not assumed). N is linearly interpolated
# between consecutive 1-July values onto each week's start/end boundary
# date, per roadmap section 4.1 -- not step-assigned by calendar year the
# way the muni-level DLNM panel does it.
#
# Births
# ------
# Real, date-level SINASC records aggregated directly to epidemiological
# weeks (01_Data/ce_sinasc_births_weekly.rds, from
# 09_fetch_ce_sinasc_weekly_births.R) -- the roadmap's preferred tier over
# distributing an annual total.
#
# Deaths
# ------
# SIM (Sistema de Informacao sobre Mortalidade) individual-record access was
# attempted and could not be resolved to a working download path from this
# environment as of 2026-09-07 (the legacy ftp.datasus.gov.br host timed out
# entirely, and no working S3/CKAN path analogous to SINASC's was found
# after a reasonably thorough search). Deaths therefore fall back to the
# roadmap's second-tier method: IBGE 2024-revision annual all-cause deaths
# (OBT_T, from the same fetch as the population stock -- internally
# coherent with N, per roadmap section 4.2's stated preference for a single
# documented source over silently mixing systems) distributed across the
# year's days at a constant daily rate, then summed to weeks. This
# preserves annual totals exactly, including leap days and cross-year
# weeks (roadmap section 4.2). Revisit if a working SIM source is found.
#
# Reconciliation
# --------------
# net_population_reconciliation[t] = N_end[t] - N_start[t] - births[t] +
# deaths[t], the exact identity in roadmap section 4.3 -- not an estimated
# migration flow, just the residual that makes the stock and flows add up
# given the externally fixed N series.
#
# Output
# ------
#   01_Data/ceara_weekly_demography.rds (+ .csv)
#     week_start, week_end, N_start, N_end, births, all_cause_deaths,
#     net_population_reconciliation, population_revision,
#     population_reference_date_rule, birth_source, death_source,
#     death_method, provisional_data_flag
# ===========================================================================

for (p in c("here", "dplyr", "lubridate", "readr", "purrr", "tibble")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(dplyr); library(lubridate); library(readr)
})

interpolate_population <- function(annual_pop, dates) {
  # annual_pop: tibble(year, pop_total), reference date = 1 July that year.
  reference_dates <- as.Date(sprintf("%d-07-01", annual_pop$year))
  stats::approx(x = reference_dates, y = annual_pop$pop_total, xout = dates, rule = 2)$y
}

distribute_annual_flow_to_weekly <- function(annual_flow, week_start_vec) {
  # annual_flow: tibble(year, total). Constant daily rate within each
  # calendar year (total / 365 or 366), summed over each week's 7 days --
  # exact for leap years and for weeks that straddle a year boundary.
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
  UF_CODE <- 23L

  demography <- readRDS(here::here("01_Data/ibge_population_projection_uf_2024revision.rds")) |>
    dplyr::filter(uf_code == UF_CODE)

  births_weekly <- readRDS(here::here("01_Data/ce_sinasc_births_weekly.rds"))
  incomplete <- births_weekly |> dplyr::filter(n_days != 7)
  if (nrow(incomplete) > 0) {
    # Only the very first/last calendar week of the fetched SINASC range can
    # be incomplete (a Sunday-start week straddling the fetch boundary) --
    # drop them rather than fail, since the fit script filters to its own
    # DATE_START/DATE_END window anyway. Extend 09's year range if a week
    # actually needed by a fit ends up among those dropped.
    message(sprintf(
      "[drop] %d incomplete boundary week(s) from ce_sinasc_births_weekly.rds: %s",
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

  out_path <- here::here("01_Data/ceara_weekly_demography.rds")
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
