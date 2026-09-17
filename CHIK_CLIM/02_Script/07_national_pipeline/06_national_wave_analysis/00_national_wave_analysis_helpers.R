# Shared helpers for the nationwide early-Re and retrospective susceptibility
# branches.  This file deliberately does not alter the frozen wave census or
# either Stan model; it only standardises their inputs and output locations.

required_national_packages <- c("dplyr", "readr", "tidyr", "tibble", "lubridate")
missing_national_packages <- required_national_packages[
  !vapply(required_national_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_national_packages)) {
  stop("Missing package(s): ", paste(missing_national_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(lubridate)
})

national_project_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM R-project root.")
}

national_path <- function(...) file.path(national_project_root(), ...)

national_output_paths <- function() {
  list(
    table = national_path("03_Output", "tables", "national_wave_analysis"),
    figure = national_path("03_Output", "figures", "national_wave_analysis"),
    fit = national_path("02_Script", "stan", "national_wave_analysis")
  )
}

ensure_national_output_dirs <- function() {
  paths <- national_output_paths()
  invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
  paths
}

read_national_wave_census <- function() {
  census_path <- national_path(
    "03_Output", "tables", "national_wave_census_v2", "brazil_chik_wave_census_v2.csv"
  )
  if (!file.exists(census_path)) stop("Frozen V2 wave census is missing: ", census_path)
  parse_census_date <- function(x) {
    # The frozen V2 CSV was written under a locale that rendered Date values
    # as m/d/Y rather than ISO-8601.  Parse explicitly so the census itself is
    # not rewritten or otherwise altered.
    as.Date(x, format = "%m/%d/%Y")
  }
  read_csv(census_path, show_col_types = FALSE) |>
    mutate(across(c(onset_week, peak_week, end_week, preceding_trough_week), parse_census_date)) |>
    arrange(state, onset_week)
}

read_national_state_week <- function() {
  panel_path <- national_path("01_Data", "chik_brazil_muni_week_2015_2025.rds")
  population_path <- national_path("01_Data", "ibge_population_projection_uf_2024revision.rds")
  if (!file.exists(panel_path) || !file.exists(population_path)) {
    stop("The national confirmed-case panel or IBGE UF population projection is missing.")
  }

  panel <- readRDS(panel_path) |>
    transmute(
      muni6 = sprintf("%06d", as.integer(muni6)),
      week_start = as.Date(week_start),
      reported_cases = as.integer(cases_confirmed),
      state = substr(muni6, 1L, 2L)
    ) |>
    dplyr::filter(week_start >= as.Date("2015-01-04"), week_start <= as.Date("2025-12-28"))

  state_lookup <- readRDS(population_path) |>
    transmute(state = sprintf("%02d", as.integer(uf_code)), state_abbr = uf_sigla) |>
    distinct()
  dates <- seq(as.Date("2015-01-04"), as.Date("2025-12-28"), by = "week")

  expand_grid(state = state_lookup$state, week_start = dates) |>
    left_join(state_lookup, by = "state") |>
    left_join(
      panel |>
        group_by(state, week_start) |>
        summarise(reported_cases = sum(reported_cases), .groups = "drop"),
      by = c("state", "week_start")
    ) |>
    mutate(reported_cases = coalesce(reported_cases, 0L)) |>
    arrange(state_abbr, week_start) |>
    select(state = state_abbr, week_start, reported_cases)
}

read_national_demography <- function() {
  population_path <- national_path("01_Data", "ibge_population_projection_uf_2024revision.rds")
  population <- readRDS(population_path) |>
    dplyr::filter(year %in% 2014:2026) |>
    transmute(
      state = uf_sigla,
      year = as.integer(year),
      pop_total = as.numeric(pop_total),
      births = as.numeric(births_total),
      deaths = as.numeric(deaths_total)
    ) |>
    arrange(state, year)
  if (n_distinct(population$state) != 27L || anyNA(population) || any(population$pop_total <= 0)) {
    stop("IBGE UF demographic input is incomplete or invalid for 2014-2026.")
  }
  population
}

interpolate_uf_population <- function(demography, dates) {
  reference_dates <- as.Date(sprintf("%d-07-01", demography$year))
  stats::approx(reference_dates, demography$pop_total, xout = dates, rule = 2)$y
}

distribute_annual_flow_to_weeks <- function(annual_flow, week_start) {
  daily_rate <- setNames(
    annual_flow$total / ifelse(leap_year(annual_flow$year), 366, 365),
    annual_flow$year
  )
  vapply(week_start, function(start) {
    days <- seq(as.Date(start), as.Date(start) + 6, by = "day")
    sum(daily_rate[as.character(year(days))])
  }, numeric(1))
}

centered_mean_3 <- function(x) {
  smoothed <- as.numeric(stats::filter(x, rep(1 / 3, 3), sides = 2))
  smoothed[is.na(smoothed)] <- x[is.na(smoothed)]
  smoothed
}

summarise_draws <- function(x, prefix) {
  tibble(
    !!paste0(prefix, "_q025") := as.numeric(quantile(x, .025, na.rm = TRUE)),
    !!paste0(prefix, "_median") := median(x, na.rm = TRUE),
    !!paste0(prefix, "_q975") := as.numeric(quantile(x, .975, na.rm = TRUE))
  )
}

hmc_gate <- function(summary_table, sampler_params, max_treedepth, fit = NULL, pars = NULL) {
  relevant <- summary_table |>
    dplyr::filter(!grepl("__$", parameter))
  divergences <- sum(vapply(sampler_params, function(x) sum(x[, "divergent__"]), numeric(1)))
  treedepth_hits <- sum(vapply(sampler_params, function(x) sum(x[, "treedepth__"] >= max_treedepth), numeric(1)))
  energy_bfmi <- vapply(sampler_params, function(x) {
    energy <- x[, "energy__"]
    mean(diff(energy)^2) / stats::var(energy)
  }, numeric(1))
  tail_ess <- NA_real_
  if (!is.null(fit)) {
    if (!requireNamespace("posterior", quietly = TRUE)) stop("The posterior package is required for tail ESS.")
    # RStan retains post-warmup iterations in iteration x chain x parameter
    # order. posterior::ess_tail then computes the requested rank-normalised
    # tail ESS without pooling chains first.
    chain_draws <- rstan::extract(fit, pars = pars, permuted = FALSE, inc_warmup = FALSE)
    tail_table <- posterior::summarise_draws(
      posterior::as_draws_array(chain_draws), posterior::ess_tail
    )
    tail_ess <- min(tail_table[[2]], na.rm = TRUE)
  }
  tibble(
    divergences = divergences,
    treedepth_hits = treedepth_hits,
    max_rhat = max(relevant$Rhat, na.rm = TRUE),
    min_bulk_ess = min(relevant$n_eff, na.rm = TRUE),
    min_tail_ess = tail_ess,
    min_bfmi = min(energy_bfmi),
    bfmi_by_chain = paste(round(energy_bfmi, 3), collapse = ";"),
    hmc_pass = divergences == 0L && treedepth_hits == 0L &&
      max(relevant$Rhat, na.rm = TRUE) <= 1.01 &&
      min(relevant$n_eff, na.rm = TRUE) >= 400 &&
      (is.na(tail_ess) || tail_ess >= 400) &&
      min(energy_bfmi) >= .30
  )
}
