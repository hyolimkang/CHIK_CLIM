# ===========================================================================
# 10e_build_riodejaneiro_city_weekly_demography.R
#
# Rio de Janeiro CITY (municipality, muni6=330455) demographic input,
# analogous to 10d (state). Reuses the exact same interpolation logic.
#
# IMPORTANT DISCLOSED APPROXIMATION: no municipality-level annual
# all-cause-death count is cached in this project (only UF-level, via
# ibge_population_projection_uf_2024revision.rds). City deaths are
# therefore obtained by allocating the RJ STATE's annual all-cause deaths
# proportional to the city's share of state population each year
# (deaths_city[y] = deaths_state[y] * pop_city[y] / pop_state[y]) -- the
# SAME population-scaled-allocation convention already used and disclosed
# elsewhere in this project (e.g. Pernambuco spatial-model importation
# allocation). City BIRTHS are genuine (SINASC filtered to CODMUNRES=
# 330455, script 09e) and city POPULATION is genuine (IBGE municipality
# population estimates, ibge_pop_muni_year_2015_2025.rds) -- only deaths
# are a proxy. This is disclosed in RJ_CITY_input_audit.md.
#
# Output
# ------
#   01_Data/riodejaneiro_city_weekly_demography_2015_2025.rds (+ .csv)
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
  MUNI6 <- "330455"
  UF_CODE <- 33L

  city_pop_annual <- readRDS(here::here("01_Data/ibge_pop_muni_year_2015_2025.rds")) |>
    dplyr::mutate(muni6 = as.character(muni6)) |>
    dplyr::filter(muni6 == MUNI6) |>
    dplyr::transmute(year, pop_total = population)
  # City births data starts late Dec 2014 (a boundary week under
  # week_start=Sunday floor_date), one year before the cached municipality
  # population series (2015-2025) begins. Back-fill 2014 with the 2015
  # population figure -- the same constant-extrapolation convention
  # `interpolate_population()`'s rule=2 already applies elsewhere in this
  # pipeline, just made explicit here so the death-share join also covers it.
  if (!2014 %in% city_pop_annual$year) {
    city_pop_annual <- dplyr::bind_rows(
      tibble::tibble(year = 2014, pop_total = city_pop_annual$pop_total[city_pop_annual$year == min(city_pop_annual$year)]),
      city_pop_annual
    )
  }
  state_demography <- readRDS(here::here("01_Data/ibge_population_projection_uf_2024revision.rds")) |>
    dplyr::filter(uf_code == UF_CODE) |>
    dplyr::transmute(year, pop_state = pop_total, deaths_state = deaths_total)

  births_weekly <- readRDS(here::here("01_Data/riodejaneiro_city_sinasc_births_weekly.rds"))
  incomplete <- births_weekly |> dplyr::filter(n_days != 7)
  if (nrow(incomplete) > 0) {
    message(sprintf(
      "[drop] %d incomplete boundary week(s) from riodejaneiro_city_sinasc_births_weekly.rds: %s",
      nrow(incomplete), paste(format(incomplete$week_start), collapse = ", ")
    ))
    births_weekly <- births_weekly |> dplyr::filter(n_days == 7)
  }

  week_start_vec <- births_weekly$week_start
  week_end_vec <- week_start_vec + 6

  boundary_dates <- sort(unique(c(week_start_vec, week_end_vec + 1)))
  n_pop_city <- interpolate_population(city_pop_annual, boundary_dates)
  pop_lookup <- setNames(n_pop_city, as.character(boundary_dates))

  # City deaths: allocate state deaths proportional to city/state population share (documented proxy)
  city_share_annual <- city_pop_annual |>
    dplyr::inner_join(state_demography, by = "year") |>
    dplyr::mutate(deaths_city_annual = deaths_state * pop_total / pop_state) |>
    dplyr::transmute(year, total = deaths_city_annual)
  deaths_weekly <- distribute_annual_flow_to_weekly(city_share_annual, week_start_vec)

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
      population_revision = "IBGE municipality population estimates (ibge_pop_muni_year_2015_2025)",
      population_reference_date_rule = "1 July each year; linear interpolation between years",
      birth_source = "SINASC individual records (DTNASC, CODMUNRES=330455), aggregated to epidemiological week",
      death_source = "PROXY: RJ-state IBGE Revisao 2024 annual all-cause deaths, allocated to the city proportional to city/state population share (no cached municipality-level death count) -- disclosed approximation, see RJ_CITY_input_audit.md",
      death_method = "state_annual_total_population_share_allocated",
      provisional_data_flag = week_start >= as.Date("2025-01-01")
    )

  if (anyNA(weekly_demography)) {
    stop("Weekly demographic contract has missing values -- check input coverage")
  }
  if (any(weekly_demography$N_start <= 0) || any(weekly_demography$N_end <= 0)) {
    stop("Interpolated population must be strictly positive")
  }

  out_path <- here::here("01_Data/riodejaneiro_city_weekly_demography_2015_2025.rds")
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
    "[check] mean weekly births=%.1f | mean weekly deaths (proxy)=%.1f | mean |reconciliation|=%.1f",
    mean(weekly_demography$births), mean(weekly_demography$all_cause_deaths),
    mean(abs(weekly_demography$net_population_reconciliation))
  ))
}
