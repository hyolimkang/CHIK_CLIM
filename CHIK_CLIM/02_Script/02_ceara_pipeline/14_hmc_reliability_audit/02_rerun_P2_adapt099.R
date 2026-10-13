# STEP 4 of the HMC reliability audit: rerun ONLY the primary P2 model
# (HalfNormal(0,0.40) prior on tau_year, full 2015-2025 period) with
# adapt_delta=0.99 (was 0.95) and increased warmup (1500, was 800) +
# sampling (500, was 400), to test whether P2's 2/1600 divergences and
# borderline Rhat=1.011 resolve with more conservative HMC settings.
# Does NOT touch/overwrite the original ce_tau_year_P2.rds -- saved
# separately as ce_tau_year_P2_adapt099.rds under this module's own
# DIR_FITS_HA. P1/P3 are explicitly NOT rerun per the task's instruction.

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

source(file.path(ROOT, "02_Script/02_ceara_pipeline/14_hmc_reliability_audit/00_config_hmc_audit.R"))
source(file.path(ROOT, "02_Script/02_ceara_pipeline/12_future_year_effect_validation/00_config_year_effect_validation.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/22_v4_3_short_2015_2019_q_calibration/scripts/02_fit_v4_3.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/27_v4_8_seeded_recurrence/scripts/01_fit_v4_8.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/28_v4_9_hierarchical_seasonality/scripts/01_fit_v4_9.R"))

# Duplicated from 12_future_year_effect_validation/02_phase1_fit_tau_year.R
# (byte-identical data construction) -- NOT re-sourcing that script
# directly, since its own top level re-fits P1+P2+P3 sequentially.
build_full_period_tau_data <- function(tau_prior_scale) {
  panel <- readRDS(file.path(ROOT, "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |> dplyr::filter(week_start >= FIXED_FULL$date_start, week_start <= FIXED_FULL$date_end)
  weekly_demography <- readRDS(file.path(ROOT, "01_Data/ceara_weekly_demography_2015_2025.rds")) |>
    dplyr::filter(week_start >= FIXED_FULL$date_start, week_start <= FIXED_FULL$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch.")
  sero_week <- find_sero_week_index(as.Date(weekly$week_start))
  prepared <- make_v4_9_data(weekly, FIXED_FULL$imports_per_week, tau_prior_scale,
                              sero_week$index, FIXED_FULL$sero_pos, FIXED_FULL$sero_n, FIXED_FULL$kappa_sero, QVAL,
                              SEED_START_DATE_FULL, SEED_LENGTH_FULL)
  prepared$stan_data$tau_prior_scale <- prepared$stan_data$year_effect_prior_sd
  prepared$stan_data$year_effect_prior_sd <- NULL
  anom <- read_csv(ACCEPTED_ANOMALY_CSV, show_col_types = FALSE) |>
    dplyr::filter(state == "CE", lag_spec == "primary") |> mutate(week_start = as.Date(week_start))
  joined <- tibble::tibble(week_start = as.Date(weekly$week_start)) |> left_join(anom |> select(week_start, z_T_anom, z_P_anom), by = "week_start")
  if (nrow(joined) != nrow(weekly)) stop("Climate anomaly join produced wrong row count.")
  joined$z_T_anom[is.na(joined$z_T_anom)] <- 0; joined$z_P_anom[is.na(joined$z_P_anom)] <- 0
  prepared$stan_data$z_T_anom <- joined$z_T_anom; prepared$stan_data$z_P_anom <- joined$z_P_anom
  list(stan_data = prepared$stan_data, weekly = weekly, sero_week = sero_week, seed_idx = prepared$seed_idx, seed_dates = prepared$seed_dates)
}
make_tau_init_fn <- function(Y, td14_scalar_inits, tau_init = 0.40) {
  base_init_fn <- make_v4_9_init_fn(Y, td14_scalar_inits)
  function(chain_id = 1L) { v <- base_init_fn(chain_id); v$beta_T <- 0; v$beta_P <- 0; v$tau_year <- tau_init; v }
}

HIGH_ADAPT_SETTINGS <- list(warmup = 1500L, sampling = 500L, chains = 4L, adapt_delta = 0.99, max_treedepth = 14L, metric = "dense_e")

message("[HMC audit Step4] Compiling tau_year Stan model (reused source, same model as P2)...")
mod <- rstan::stan_model(file = STAN_SOURCE_TAU_YEAR)

message("[HMC audit Step4] Rebuilding P2's own data (tau_prior_scale=0.40, full 2015-2025, 2022 seed retained)...")
prepared <- build_full_period_tau_data(TAU_PRIORS$P2)
td14_scalar_inits <- load_td14_scalar_inits(ROOT)
init_fn <- make_tau_init_fn(prepared$stan_data$Y, td14_scalar_inits, tau_init = TAU_PRIORS$P2)

progress_dir <- file.path(DIR_FITS_HA, "chains_P2_adapt099"); dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
fit <- rstan::sampling(
  mod, data = prepared$stan_data, chains = HIGH_ADAPT_SETTINGS$chains,
  iter = HIGH_ADAPT_SETTINGS$warmup + HIGH_ADAPT_SETTINGS$sampling, warmup = HIGH_ADAPT_SETTINGS$warmup,
  seed = 20261013L, init = init_fn,
  control = list(adapt_delta = HIGH_ADAPT_SETTINGS$adapt_delta, max_treedepth = HIGH_ADAPT_SETTINGS$max_treedepth, metric = HIGH_ADAPT_SETTINGS$metric),
  refresh = 50, sample_file = file.path(progress_dir, "chain")
)
elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

gate_pars <- c("alpha_R", "z_year", "tau_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2", "sigma_season_year", "z_season_year", "phi_obs", "beta_T", "beta_P")
hmc <- compute_hmc_gate(fit, gate_pars, HIGH_ADAPT_SETTINGS$max_treedepth)

bundle <- list(fit = fit, stan_data = prepared$stan_data, weekly_data = prepared$weekly, sero_week = prepared$sero_week,
                seed_idx = prepared$seed_idx, seed_dates = prepared$seed_dates,
                config = c(HIGH_ADAPT_SETTINGS, list(prior_id = "P2_adapt099", model_version = "ceara_v4_9_climate_forced_tau_year",
                                                      tau_prior_scale = TAU_PRIORS$P2, q = QVAL, elapsed_seconds = elapsed_seconds, seed = 20261013L)),
                hmc = hmc)
saveRDS(bundle, file.path(DIR_FITS_HA, "ce_tau_year_P2_adapt099.rds"))
message(sprintf("[HMC audit Step4] P2 rerun (adapt_delta=0.99, warmup=1500): HMC gate: %s | divergences=%d | max Rhat=%.4f | min ESS=%.0f | %.0fs",
                 if (hmc$hmc_pass) "PASS" else "FAIL", hmc$divergences, hmc$maximum_rhat, hmc$minimum_bulk_ess, elapsed_seconds))
message("[HMC audit Step4] Saved (NEW file, original ce_tau_year_P2.rds untouched): ", file.path(DIR_FITS_HA, "ce_tau_year_P2_adapt099.rds"))
