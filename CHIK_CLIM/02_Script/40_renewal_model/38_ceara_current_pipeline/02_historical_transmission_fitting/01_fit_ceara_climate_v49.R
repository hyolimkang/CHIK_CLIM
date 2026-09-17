# Climate-forced v4.9 canary -- CE @ q=0.05, the cleanest of the 4
# candidate reconstructions (Phase 0 audit). Reuses CE's OWN existing
# make_v4_9_data()/make_v4_9_init_fn() (from
# 28_v4_9_hierarchical_seasonality/scripts/01_fit_v4_9.R) UNCHANGED --
# the only new step is appending the Phase-1 z_T_anom/z_P_anom vectors
# (primary lag spec) to the exact same stan_data list, matched by week.
#
# NA handling (disclosed convention, smallest possible choice): the first
# few weeks of the series lack a complete lag window (climate coverage
# starts 2014-12-28, before some early lag windows can be filled) and
# have NA z_T_anom/z_P_anom. These are set to 0 (i.e. "assume no known
# anomaly / average climate") for those boundary weeks only -- NOT
# imputed via any new method. This affects a handful of weeks in
# 2015-01/02, far from every epidemic wave used elsewhere in this project.
#
# This script performs ONLY a smoke test (1 chain, short iterations) to
# verify the new Stan file compiles and runs without dimension/parsing
# errors. It does NOT constitute the Phase 2 full canary fit.

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script/40_renewal_model/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script/40_renewal_model/22_v4_3_short_2015_2019_q_calibration/scripts/02_fit_v4_3.R"))
source(file.path(root, "02_Script/40_renewal_model/27_v4_8_seeded_recurrence/scripts/01_fit_v4_8.R"))
source(file.path(root, "02_Script/40_renewal_model/28_v4_9_hierarchical_seasonality/scripts/01_fit_v4_9.R"))

base_dir <- file.path(root, "02_Script/40_renewal_model/36_climate_forced_v4_9")
table_dir <- file.path(root, "03_Output/tables/climate_forced_v4_9")
out_dir <- file.path(base_dir, "outputs", "ce_canary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

QVAL <- 0.05 # CE's established reference scenario -- unchanged, fixed q, not re-privileged
SEED_START_DATE <- as.Date("2022-01-02"); SEED_LENGTH <- 8L
FIXED <- list(date_start = as.Date("2015-01-04"), date_end = as.Date("2025-12-21"),
              imports_per_week = 1, year_effect_prior_sd = 0.40,
              sero_pos = 103L, sero_n = 404L, kappa_sero = 50)

build_climate_forced_ce_data <- function() {
  panel <- readRDS(file.path(root, "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |> filter(week_start >= FIXED$date_start, week_start <= FIXED$date_end)
  weekly_demography <- readRDS(file.path(root, "01_Data/ceara_weekly_demography_2015_2025.rds")) |>
    filter(week_start >= FIXED$date_start, week_start <= FIXED$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch.")

  sero_week <- find_sero_week_index(as.Date(weekly$week_start))
  prepared <- make_v4_9_data(weekly, FIXED$imports_per_week, FIXED$year_effect_prior_sd,
                              sero_week$index, FIXED$sero_pos, FIXED$sero_n, FIXED$kappa_sero, QVAL,
                              SEED_START_DATE, SEED_LENGTH)

  anom <- read_csv(file.path(table_dir, "climate_anomaly_covariates.csv"), show_col_types = FALSE) |>
    filter(state == "CE", lag_spec == "primary") |> mutate(week_start = as.Date(week_start))
  joined <- tibble::tibble(week_start = as.Date(weekly$week_start)) |> left_join(anom |> select(week_start, z_T_anom, z_P_anom), by = "week_start")
  if (nrow(joined) != nrow(weekly)) stop("Climate anomaly join produced wrong row count.")
  n_na <- sum(is.na(joined$z_T_anom) | is.na(joined$z_P_anom))
  message(sprintf("Climate-forced CE data: %d weeks, %d with NA anomaly (boundary weeks, set to 0)", nrow(joined), n_na))
  joined$z_T_anom[is.na(joined$z_T_anom)] <- 0
  joined$z_P_anom[is.na(joined$z_P_anom)] <- 0

  prepared$stan_data$z_T_anom <- joined$z_T_anom
  prepared$stan_data$z_P_anom <- joined$z_P_anom
  list(stan_data = prepared$stan_data, weekly = weekly, sero_week = sero_week, seed_idx = prepared$seed_idx, seed_dates = prepared$seed_dates)
}

make_climate_init_fn <- function(Y, td14_scalar_inits) {
  base_init_fn <- make_v4_9_init_fn(Y, td14_scalar_inits)
  function(chain_id = 1L) {
    v <- base_init_fn(chain_id)
    v$beta_T <- 0 # neutral start: no assumed climate effect
    v$beta_P <- 0
    v
  }
}

run_smoke_test <- function() {
  prepared <- build_climate_forced_ce_data()
  message(sprintf("Stan data ready: N=%d, Y=%d, z_T_anom range=[%.2f,%.2f], z_P_anom range=[%.2f,%.2f]",
                   prepared$stan_data$N, prepared$stan_data$Y, min(prepared$stan_data$z_T_anom), max(prepared$stan_data$z_T_anom),
                   min(prepared$stan_data$z_P_anom), max(prepared$stan_data$z_P_anom)))

  stan_path <- file.path(root, "02_Script/stan/renewal_ceara_v4_9_climate_forced_canary.stan")
  message("Compiling ", stan_path, " ...")
  mod <- rstan::stan_model(file = stan_path)
  message("Compiled OK.")

  td14_scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_climate_init_fn(prepared$stan_data$Y, td14_scalar_inits)

  message("Running smoke test: 1 chain, 20 warmup + 20 sampling ...")
  fit <- rstan::sampling(mod, data = prepared$stan_data, chains = 1L, iter = 40L, warmup = 20L,
                          init = init_fn, seed = 1L,
                          control = list(adapt_delta = 0.95, max_treedepth = 10, metric = "dense_e"), refresh = 5)
  message("SMOKE TEST COMPLETE -- no dimension/parsing errors.")
  print(fit, pars = c("alpha_R", "beta_T", "beta_P", "phi_obs"))
  saveRDS(list(fit = fit, stan_data = prepared$stan_data), file.path(out_dir, "ce_climate_canary_smoketest.rds"))
  invisible(fit)
}

run_full_canary_fit <- function() {
  settings <- list(warmup = 800L, sampling = 400L, chains = 4L, adapt_delta = 0.95, max_treedepth = 14L, metric = "dense_e")
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  prepared <- build_climate_forced_ce_data()
  message(sprintf("[CE climate canary] Stan data: N=%d, Y=%d, q=%.2f (fixed, unchanged from frozen v4.9)",
                   prepared$stan_data$N, prepared$stan_data$Y, QVAL))

  stan_path <- file.path(root, "02_Script/stan/renewal_ceara_v4_9_climate_forced_canary.stan")
  mod <- rstan::stan_model(file = stan_path)
  td14_scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_climate_init_fn(prepared$stan_data$Y, td14_scalar_inits)

  sample_file <- file.path(progress_dir, "chain")
  started <- Sys.time()
  fit <- rstan::sampling(
    mod, data = prepared$stan_data, chains = settings$chains,
    iter = settings$warmup + settings$sampling, warmup = settings$warmup,
    seed = 20260916L, init = init_fn,
    control = list(adapt_delta = settings$adapt_delta, max_treedepth = settings$max_treedepth, metric = settings$metric),
    refresh = 25, sample_file = sample_file
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  gate_pars <- c("alpha_R", "z_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                 "sigma_season_year", "z_season_year", "phi_obs", "beta_T", "beta_P")
  hmc <- compute_hmc_gate(fit, gate_pars, settings$max_treedepth)

  bundle <- list(fit = fit, stan_data = prepared$stan_data, weekly_data = prepared$weekly,
                  sero_week = prepared$sero_week, seed_idx = prepared$seed_idx, seed_dates = prepared$seed_dates,
                  config = c(settings, list(model_version = "ceara_v4_9_climate_forced_canary", q = QVAL,
                                             elapsed_seconds = elapsed_seconds, seed = 20260916L)),
                  hmc = hmc)
  fit_path <- file.path(out_dir, "ce_climate_forced_canary_q0.05.rds")
  saveRDS(bundle, fit_path)
  message(sprintf("[CE climate canary] HMC gate: %s | divergences=%d | max Rhat=%.4f | %.0fs | %s",
                   if (hmc$hmc_pass) "PASS" else "FAIL", hmc$divergences, hmc$maximum_rhat, elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_full_canary_fit()
