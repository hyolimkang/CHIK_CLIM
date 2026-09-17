# =============================================================================
# Prepare the Ceará 2015-2025 annual FOI / susceptibility pilot input.
#
# Outcome: SINAN CLASSI_FIN == 13 (confirmed chikungunya), assigned to the
# calendar year of the cleaned event date (onset preferred, notification
# fallback). This is not the notified/suspected case count.
# =============================================================================

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr", "lubridate")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tibble); library(lubridate) })

project_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM R-project root.")
}
project_path <- function(...) file.path(project_root(), ...)

make_q_scenarios <- function(n_draws = 500000L, seed = 20260910L) {
  # Existing renewal priors are used only to construct an externally informed
  # *overall confirmed-case detection* anchor: q = p_symp * rho_sym_Brazil.
  # q itself remains one constant logit-scale Stan parameter in this model.
  set.seed(seed)
  q_reference <- stats::rbeta(n_draws, 30, 28) * stats::rbeta(n_draws, 20, 60)
  logit_reference <- stats::qlogis(q_reference)
  primary_mean <- mean(logit_reference)
  primary_sd <- stats::sd(logit_reference)
  primary_median <- stats::median(q_reference)
  tibble(
    scenario = c("Q1_primary_product_anchor", "Q2_broader_product_anchor", "Q3_lower_detection"),
    q_logit_prior_mean = c(primary_mean, primary_mean, stats::qlogis(primary_median / 2)),
    q_logit_prior_sd = c(primary_sd, 2 * primary_sd, primary_sd),
    prior_description = c(
      "Moment-matched logit-normal approximation to p_symp * rho_sym_Brazil",
      "Same centre as Q1 with twice the logit-scale SD",
      "Half the Q1 median detection probability, with Q1 logit-scale SD"
    )
  )
}

interpolate_population <- function(population_ce, boundary_dates) {
  ref_dates <- as.Date(sprintf("%d-07-01", population_ce$year))
  stats::approx(ref_dates, population_ce$pop_total, xout = boundary_dates, rule = 2)$y
}

build_dynamic_annual_foi_data <- function(write_outputs = TRUE) {
  years <- 2015:2025
  uf_code <- "23"
  uf_numeric <- 23L
  date_start <- as.Date("2015-01-01")
  date_end <- as.Date("2025-12-31")

  foi_ensemble <- readr::read_csv(project_path("01_Data", "brazil_uf_foi_ensemble.csv"), show_col_types = FALSE) |>
    filter(uf == "CE")
  foi_summary <- readr::read_csv(project_path("01_Data", "brazil_uf_foi_longterm.csv"), show_col_types = FALSE) |>
    filter(uf == "CE")
  if (nrow(foi_ensemble) != 100L || nrow(foi_summary) != 1L || any(foi_ensemble$foi_equiv <= 0)) {
    stop("Expected one Ceará summary and 100 positive Ceará foi_equiv ensemble members.")
  }
  log_foi <- log(foi_ensemble$foi_equiv)

  # The individual file is the authoritative calendar-date case source.
  individual <- readRDS(project_path("01_Data", "chik_sinan_individual_2015_2025.rds"))
  annual_cases <- individual |>
    transmute(
      muni6 = coalesce(muni_residence6, muni_notif6),
      event_date = as.Date(event_date),
      is_confirmed_chik = is_confirmed_chik
    ) |>
    filter(substr(muni6, 1, 2) == uf_code, is_confirmed_chik,
           event_date >= date_start, event_date <= date_end) |>
    mutate(year = lubridate::year(event_date)) |>
    count(year, name = "observed_cases") |>
    tidyr::complete(year = years, fill = list(observed_cases = 0L)) |>
    arrange(year)

  # Rebuild the audited state weekly series from the zero-filled source panel
  # for the diagnostic figure and verify all weekly counts against the same
  # cleaned individual records before calendar-year aggregation.
  weekly_source <- readRDS(project_path("01_Data", "chik_brazil_muni_week_2015_2025.rds")) |>
    mutate(week_start = as.Date(week_start), muni6 = sprintf("%06d", as.integer(muni6))) |>
    filter(substr(muni6, 1, 2) == uf_code) |>
    group_by(week_start) |>
    summarise(cases_confirmed = sum(cases_confirmed), .groups = "drop") |>
    arrange(week_start)
  individual_weekly <- individual |>
    transmute(
      muni6 = coalesce(muni_residence6, muni_notif6), event_date = as.Date(event_date),
      is_confirmed_chik = is_confirmed_chik
    ) |>
    filter(substr(muni6, 1, 2) == uf_code, is_confirmed_chik,
           event_date >= date_start, event_date <= date_end) |>
    mutate(week_start = lubridate::floor_date(event_date, "week", week_start = 7)) |>
    count(week_start, name = "individual_confirmed")
  weekly_check <- full_join(weekly_source, individual_weekly, by = "week_start") |>
    mutate(across(c(cases_confirmed, individual_confirmed), ~ coalesce(.x, 0L)))
  if (sum(weekly_check$cases_confirmed) != sum(weekly_check$individual_confirmed)) {
    stop("Confirmed-case total differs between weekly source and individual source.")
  }

  population <- readRDS(project_path("01_Data", "ibge_population_projection_uf_2024revision.rds")) |>
    filter(.data$uf_code == uf_numeric) |>
    arrange(year)
  required_population_columns <- c("year", "pop_total", "deaths_total")
  if (!all(required_population_columns %in% names(population))) stop("Population projection columns are incomplete.")
  boundary_dates <- as.Date(sprintf("%d-01-01", 2015:2026))
  boundary_population <- interpolate_population(population, boundary_dates)

  births_daily <- readRDS(project_path("01_Data", "ce_sinasc_births_daily.rds")) |>
    mutate(date = as.Date(date), year = lubridate::year(date)) |>
    filter(year %in% years)
  expected_birth_dates <- seq(date_start, date_end, by = "day")
  if (!identical(births_daily$date, expected_birth_dates)) stop("SINASC daily Ceará births are not complete for 2015-2025.")
  annual_births <- births_daily |>
    group_by(year) |>
    summarise(births = sum(births), .groups = "drop")
  annual_deaths <- population |>
    filter(year %in% years) |>
    transmute(year, deaths = deaths_total)

  annual <- tibble(
    year = years,
    N_start = boundary_population[-length(boundary_population)],
    N_end = boundary_population[-1]
  ) |>
    left_join(annual_births, by = "year") |>
    left_join(annual_deaths, by = "year") |>
    left_join(annual_cases, by = "year") |>
    mutate(reconciliation = N_end - N_start - births + deaths)
  if (anyNA(annual) || any(annual$N_start <= 0 | annual$N_end <= 0 | annual$births < 0 | annual$deaths < 0)) {
    stop("Annual demographic input contains missing or invalid values.")
  }
  if (any(abs(annual$N_end - (annual$N_start + annual$births - annual$deaths + annual$reconciliation)) > 1e-5)) {
    stop("Annual demographic reconciliation identity failed.")
  }

  q_scenarios <- make_q_scenarios()
  stan_base_data <- list(
    Y = nrow(annual), C = as.integer(annual$observed_cases),
    N_start = annual$N_start, N_end = annual$N_end, births = annual$births,
    deaths = annual$deaths, reconciliation = annual$reconciliation,
    log_foi_prior_mean = mean(log_foi), log_foi_prior_sd = stats::sd(log_foi)
  )
  result <- list(
    annual = annual, weekly_cases = weekly_source, weekly_case_audit = weekly_check,
    stan_base_data = stan_base_data, q_scenarios = q_scenarios,
    foi_anchor = list(
      definition = "foi_equiv", ensemble_members = 100L,
      log_mean = mean(log_foi), log_sd = stats::sd(log_foi),
      median = stats::median(foi_ensemble$foi_equiv),
      q025 = stats::quantile(foi_ensemble$foi_equiv, .025),
      q975 = stats::quantile(foi_ensemble$foi_equiv, .975),
      population_coverage = foi_summary$population_coverage,
      population_coverage_min_member = foi_summary$population_coverage_min_member,
      population_coverage_max_member = foi_summary$population_coverage_max_member
    ),
    metadata = list(
      state = "Ceara", uf = "CE", years = years,
      case_definition = "SINAN CLASSI_FIN == 13: confirmed chikungunya",
      initial_susceptibility = "S_start[2015] = N_start[2015]; U_start[2015] = 0"
    )
  )

  if (write_outputs) {
    out_dir <- project_path("03_Output", "tables", "dynamic_annual_foi")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, project_path("01_Data", "dynamic_annual_foi_ceara_input.rds"))
    write_csv(annual, file.path(out_dir, "dynamic_annual_foi_ceara_input_audit.csv"))
    write_csv(q_scenarios, file.path(out_dir, "dynamic_annual_foi_q_prior_scenarios.csv"))
    write_csv(tibble(
      state = "Ceara", foi_definition = "foi_equiv", population_coverage = foi_summary$population_coverage,
      population_coverage_min_member = foi_summary$population_coverage_min_member,
      population_coverage_max_member = foi_summary$population_coverage_max_member
    ), file.path(out_dir, "dynamic_annual_foi_foi_coverage_audit.csv"))
    write_csv(weekly_check, file.path(out_dir, "dynamic_annual_foi_weekly_case_audit.csv"))
  }
  result
}

if (sys.nframe() == 0L) {
  prepared <- build_dynamic_annual_foi_data(write_outputs = TRUE)
  message("[prepare] Ceará annual dynamic-FOI input: ", min(prepared$annual$year), "-", max(prepared$annual$year))
}
