# Unchanged-model v2.2 full-period stress test (2015--2025).
# This script extends only the validated data window. It never alters
# renewal_ceara_v2_2.stan, its priors, likelihoods, or demographic recursion.
required_packages <- c("here", "rstan", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(dplyr)
})
rstan_options(auto_write = TRUE)
options(mc.cores = min(4, parallel::detectCores()))
source(here::here("02_Script/40_renewal_model/00_shared/01_build_ceara_state_weekly.R"))

discretize_gamma_generation_interval <- function(mean_weeks, sd_weeks, G) {
  shape <- (mean_weeks / sd_weeks)^2
  rate <- mean_weeks / sd_weeks^2
  x <- diff(stats::pgamma(0:G, shape = shape, rate = rate))
  x / sum(x)
}
env_integer <- function(name, default) {
  x <- Sys.getenv(name, "")
  if (!nzchar(x)) {
    return(default)
  }
  x <- suppressWarnings(as.integer(x))
  if (is.na(x) || x < 1L) stop(name, " must be a positive integer")
  x
}

run_v2_2_full_stress <- function() {
  date_start <- as.Date("2015-01-04")
  date_end <- as.Date("2025-12-21")
  G <- 8L
  seed_weeks <- 8L
  n_iter <- env_integer("RENEWAL_V2_2_FULL_STRESS_ITER", 2000L)
  n_warmup <- env_integer("RENEWAL_V2_2_FULL_STRESS_WARMUP", 1000L)
  n_chains <- env_integer("RENEWAL_V2_2_FULL_STRESS_CHAINS", 4L)
  if (n_warmup >= n_iter) stop("Warmup must be smaller than iterations")
  panel <- readRDS(here::here("01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
  cases <- build_state_weekly(panel, "23") |> filter(week_start >= date_start, week_start <= date_end)
  demography_path <- here::here("01_Data/ceara_weekly_demography_2015_2025.rds")
  demography <- readRDS(demography_path) |> filter(week_start >= date_start, week_start <= date_end)
  weekly <- cases |>
    select(-population) |>
    inner_join(demography, by = "week_start")
  expected_dates <- seq(date_start, date_end, by = "week")
  if (!identical(as.Date(weekly$week_start), expected_dates) || nrow(weekly) != nrow(cases) || nrow(weekly) != nrow(demography)) stop("Cases and validated demography do not form the same complete 2015--2025 weekly sequence")
  if (anyNA(weekly$cases) || any(weekly$cases < 0) || any(abs(weekly$cases - round(weekly$cases)) > 1e-8)) stop("Cases must be non-negative integers")
  accounting_error <- max(abs(weekly$N_end - weekly$N_start - weekly$births + weekly$all_cause_deaths - weekly$net_population_reconciliation))
  if (accounting_error > 1e-6 || any(weekly$N_start <= 0) || any(weekly$N_end <= 0) || any(weekly$births < 0) || any(weekly$all_cause_deaths < 0)) stop("Validated demographic accounting contract failed")
  sero_start <- match(as.Date(c("2018-06-03", "2018-06-03")), weekly$week_start)
  sero_end <- match(as.Date(c("2018-12-30", "2019-12-29")), weekly$week_start)
  if (anyNA(sero_start) || anyNA(sero_end)) stop("Serology windows are outside the full-period data")
  w <- discretize_gamma_generation_interval(2, 1, G)
  stan_data <- list(
    N = nrow(weekly), G = G, seed_weeks = seed_weeks, C = as.integer(weekly$cases), w = as.vector(w),
    N_start = as.vector(weekly$N_start), N_end = as.vector(weekly$N_end), births = as.vector(weekly$births), deaths = as.vector(weekly$all_cause_deaths),
    time_scaled = seq(-.5, .5, length.out = nrow(weekly)), log_seed_prior_mean = rep(log(10), seed_weeks), seed_prior_sd = 1.5, reporting_region_sd = .75,
    J = 2L, sero_window_start = as.integer(sero_start), sero_window_end = as.integer(sero_end), sero_positive = c(103L, 289L), sero_n = c(404L, 409L)
  )
  # Identical initialization construction to v2.2, extended only in time.
  smooth_cases <- vapply(seq_len(stan_data$N), function(t) mean(stan_data$C[max(1L, t - 2L):min(stan_data$N, t + 2L)] + .5), numeric(1))
  rough_X <- smooth_cases * (.25 * stan_data$N_start[1] / sum(smooth_cases))
  rough_S <- stan_data$N_start[1]
  rough_U <- 0
  rough_hazard <- numeric(stan_data$N)
  net_migration <- stan_data$N_end - stan_data$N_start - stan_data$births + stan_data$deaths
  for (t in seq_len(stan_data$N)) {
    rough_X[t] <- min(rough_X[t], .2 * rough_S)
    rough_hazard[t] <- log(-log1p(-rough_X[t] / rough_S))
    s_after <- rough_S - rough_X[t]
    u_after <- rough_U + rough_X[t]
    frac_s <- s_after / (s_after + u_after)
    rough_S <- s_after + stan_data$births[t] - stan_data$deaths[t] * frac_s + net_migration[t] * frac_s
    rough_U <- u_after - stan_data$deaths[t] * (1 - frac_s) + net_migration[t] * (1 - frac_s)
  }
  init_fn <- function() {
    list(
      sigma_R = .05, log_infection_scale = rough_hazard[1], log_hazard_relative = rough_hazard[-1] - rough_hazard[1],
      p_symp = .52, rho_sym_brazil = .25, logit_rho_sym_ceara_mid = qlogis(.116), reporting_trend = .3, sigma_geo = .5, z_site = rep(0, 2), phi_obs = 20
    )
  }
  stan_file <- here::here("02_Script/stan/renewal_ceara_v2_2.stan")
  message(sprintf("[data] %s to %s | %d weeks | %s cases | demographic accounting max error %.3g", date_start, date_end, stan_data$N, format(sum(stan_data$C), big.mark = ","), accounting_error))
  elapsed <- system.time(fit <- rstan::stan(file = stan_file, data = stan_data, iter = n_iter, warmup = n_warmup, chains = n_chains, seed = 24052018L, init = init_fn, control = list(adapt_delta = .97, max_treedepth = 13), refresh = max(1L, floor(n_iter / 10L))))[["elapsed"]]
  output <- here::here("02_Script/stan/renewal_ceara_v2_2_full_stress_fit.rds")
  saveRDS(list(fit = fit, stan_data = stan_data, weekly_data = weekly, config = list(stress_test = TRUE, model_version = "v2.2 unchanged", date_start = date_start, date_end = date_end, demography_source = normalizePath(demography_path, winslash = "/", mustWork = TRUE), demographic_accounting_error = accounting_error, population_revision = unique(weekly$population_revision), serology_sites = c("Juazeiro do Norte", "Quixada"), serology_positive = c(103L, 289L), serology_n = c(404L, 409L), iterations = n_iter, warmup = n_warmup, chains = n_chains, elapsed_seconds = elapsed, stan_source = normalizePath(stan_file, winslash = "/", mustWork = TRUE))), output)
  message("[save] ", output)
}
if (sys.nframe() == 0) run_v2_2_full_stress()
