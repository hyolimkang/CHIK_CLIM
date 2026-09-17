# v2.3_state full-period fit (2015-2025).
#
# This script is a pre-climate state-level benchmark. It preserves the v2.3
# model specification and saves one self-contained fit bundle; it never
# overwrites v2.2 sources or results.
## 1. Packages and project helper ------------------------------------------------

required_packages <- c("here", "rstan", "dplyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(dplyr)
})

rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

source(here::here(
  "02_Script", "00_shared", "functions", "01_build_ceara_state_weekly.R"
))

## 2. Utilities and settings -----------------------------------------------------

gi_weights <- function() {
  weights <- diff(pgamma(0:8, shape = 4, rate = 2))
  weights / sum(weights)
}

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }

  value <- as.integer(value)
  if (is.na(value) || value < 1L) stop(name, " must be a positive integer.")
  value
}
run_v2_3_state <- function() {
  ## Fit configuration ----------------------------------------------------------
  start_date <- as.Date("2015-01-04")
  end_date <- as.Date("2025-12-21")
  generation_interval_weeks <- 8L
  seed_weeks <- 8L

  iterations <- env_int("RENEWAL_V2_3_ITER", 2000L)
  warmup <- env_int("RENEWAL_V2_3_WARMUP", 1000L)
  chains <- env_int("RENEWAL_V2_3_CHAINS", 4L)
  if (warmup >= iterations) stop("Warmup must be smaller than iterations.")

  # The primary serology likelihood is deliberately tempered in v2.3.
  sero_power <- as.numeric(Sys.getenv("RENEWAL_V2_3_SERO_POWER", "0.5"))
  if (!is.finite(sero_power) || sero_power <= 0 || sero_power > 1) {
    stop("RENEWAL_V2_3_SERO_POWER must be in (0, 1].")
  }

  # This Beta prior is on q = overall infection-to-report probability.
  q_prior_alpha <- 16.17486
  q_prior_beta <- 108.896

  ## Read and validate the weekly inputs ---------------------------------------
  case_panel <- readRDS(here::here(
    "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"
  ))
  weekly_cases <- build_state_weekly(case_panel, "23") |>
    filter(week_start >= start_date, week_start <= end_date)

  demography_path <- here::here(
    "01_Data", "ceara_weekly_demography_2015_2025.rds"
  )
  weekly_demography <- readRDS(demography_path) |>
    filter(week_start >= start_date, week_start <= end_date)

  # Use only the corrected demographic series for population accounting.
  weekly <- weekly_cases |>
    select(-population) |>
    inner_join(weekly_demography, by = "week_start")

  expected_dates <- seq(start_date, end_date, by = "week")
  if (
    !identical(as.Date(weekly$week_start), expected_dates) ||
      nrow(weekly) != nrow(weekly_cases) ||
      nrow(weekly) != nrow(weekly_demography)
  ) {
    stop("Incomplete weekly input sequence.")
  }

  demographic_accounting_error <- max(abs(
    weekly$N_end - weekly$N_start - weekly$births +
      weekly$all_cause_deaths - weekly$net_population_reconciliation
  ))
  if (demographic_accounting_error > 1e-6) {
    stop("Demographic accounting error.")
  }

  ## Serology windows -----------------------------------------------------------
  serology_start <- match(
    as.Date(c("2018-06-03", "2018-06-03")),
    weekly$week_start
  )
  serology_end <- match(
    as.Date(c("2018-12-30", "2019-12-29")), weekly$week_start
  )
  if (anyNA(c(serology_start, serology_end))) {
    stop("Serology windows unavailable in the fitting data.")
  }

  ## Stan data: this block mirrors renewal_ceara_v2_3_state.stan --------------
  n_weeks <- nrow(weekly)
  stan_data <- list(
    N = n_weeks,
    G = generation_interval_weeks,
    seed_weeks = seed_weeks,
    C = as.integer(weekly$cases),
    w = gi_weights(),
    N_start = weekly$N_start,
    N_end = weekly$N_end,
    births = weekly$births,
    deaths = weekly$all_cause_deaths,
    seasonal_sin = sin(2 * pi * seq_len(n_weeks) / 52.1775),
    seasonal_cos = cos(2 * pi * seq_len(n_weeks) / 52.1775),
    log_seed_prior_mean = rep(log(10), seed_weeks),
    seed_prior_sd = 1.5,
    q_prior_a = q_prior_alpha,
    q_prior_b = q_prior_beta,
    J = 2L,
    sero_window_start = as.integer(serology_start),
    sero_window_end = as.integer(serology_end),
    sero_positive = c(103L, 289L),
    sero_n = c(404L, 409L),
    sero_power = sero_power
  )

  initial_values <- function() {
    list(
      alpha_R = log(1.2),
      beta_sin = 0,
      beta_cos = 0,
      phi_R = 0.5,
      sigma_R = 0.05,
      z_R = rep(0, stan_data$N),
      log_seed_hazard = rep(log(10 / stan_data$N_start[1]), seed_weeks),
      q = q_prior_alpha / (q_prior_alpha + q_prior_beta),
      sigma_geo = 0.5,
      z_site = rep(0, stan_data$J),
      phi_obs = 20
    )
  }

  ## Sample and save ------------------------------------------------------------
  stan_source <- here::here("02_Script", "stan", "renewal_ceara_v2_3_state.stan")
  output_path <- here::here("02_Script", "stan", "renewal_ceara_v2_3_state_fit.rds")

  message(sprintf(
    "[data] %s to %s | N = %d | cases = %s | serology power = %.2f",
    start_date, end_date, stan_data$N,
    format(sum(stan_data$C), big.mark = ","), sero_power
  ))

  elapsed_seconds <- system.time({
    fit <- rstan::stan(
      file = stan_source,
      data = stan_data,
      chains = chains,
      iter = iterations,
      warmup = warmup,
      seed = 24052018,
      init = initial_values,
      control = list(adapt_delta = 0.97, max_treedepth = 13),
      refresh = max(1L, floor(iterations / 10L))
    )
  })[["elapsed"]]

  result_bundle <- list(
    fit = fit,
    stan_data = stan_data,
    weekly_data = weekly,
    config = list(
      model_version = "v2.3_state",
      date_start = start_date,
      date_end = end_date,
      demography_source = normalizePath(demography_path, winslash = "/", mustWork = TRUE),
      demographic_accounting_error = demographic_accounting_error,
      sero_power = sero_power,
      q_prior = c(alpha = q_prior_alpha, beta = q_prior_beta),
      recurrent_introduction = FALSE,
      iterations = iterations,
      warmup = warmup,
      chains = chains,
      elapsed_seconds = elapsed_seconds,
      stan_source = normalizePath(stan_source, winslash = "/", mustWork = TRUE)
    )
  )
  saveRDS(result_bundle, output_path)
  message("[save] ", output_path)
}

if (sys.nframe() == 0L) run_v2_3_state()
