# PHASE 1 -- Estimate year-effect variance (tau_year), full 2015-2025
# period, RETAINING the 2022 conditioned seed exactly as in the accepted
# fit (explicitly disclosed -- this is a conditional historical fit, NOT
# the seed-free holdout validation of Phase 3).
#
# Fits P1 (HalfNormal(0,0.20)), P2 (HalfNormal(0,0.40)), P3 (HalfNormal(0,0.80))
# using the SAME compiled tau_year Stan model and the SAME HMC settings
# as the accepted M0 fit (4 chains, 800 warmup + 400 sampling,
# adapt_delta=0.95, max_treedepth=14, dense_e).

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

source(file.path(ROOT, "02_Script/02_ceara_pipeline/12_future_year_effect_validation/00_config_year_effect_validation.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/22_v4_3_short_2015_2019_q_calibration/scripts/02_fit_v4_3.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/27_v4_8_seeded_recurrence/scripts/01_fit_v4_8.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/28_v4_9_hierarchical_seasonality/scripts/01_fit_v4_9.R"))

build_full_period_tau_data <- function(tau_prior_scale) {
  panel <- readRDS(file.path(ROOT, "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |> dplyr::filter(week_start >= FIXED_FULL$date_start, week_start <= FIXED_FULL$date_end)
  weekly_demography <- readRDS(file.path(ROOT, "01_Data/ceara_weekly_demography_2015_2025.rds")) |>
    dplyr::filter(week_start >= FIXED_FULL$date_start, week_start <= FIXED_FULL$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch.")

  sero_week <- find_sero_week_index(as.Date(weekly$week_start))
  # make_v4_9_data() expects year_effect_prior_sd as its 2nd arg -- reused unmodified for everything EXCEPT
  # that we immediately rename the field below (tau_prior_scale replaces year_effect_prior_sd in the tau_year Stan model's data block).
  prepared <- make_v4_9_data(weekly, FIXED_FULL$imports_per_week, tau_prior_scale,
                              sero_week$index, FIXED_FULL$sero_pos, FIXED_FULL$sero_n, FIXED_FULL$kappa_sero, QVAL,
                              SEED_START_DATE_FULL, SEED_LENGTH_FULL)
  prepared$stan_data$tau_prior_scale <- prepared$stan_data$year_effect_prior_sd
  prepared$stan_data$year_effect_prior_sd <- NULL

  anom <- read_csv(ACCEPTED_ANOMALY_CSV, show_col_types = FALSE) |>
    dplyr::filter(state == "CE", lag_spec == "primary") |> mutate(week_start = as.Date(week_start))
  joined <- tibble::tibble(week_start = as.Date(weekly$week_start)) |> left_join(anom |> select(week_start, z_T_anom, z_P_anom), by = "week_start")
  if (nrow(joined) != nrow(weekly)) stop("Climate anomaly join produced wrong row count.")
  joined$z_T_anom[is.na(joined$z_T_anom)] <- 0
  joined$z_P_anom[is.na(joined$z_P_anom)] <- 0
  prepared$stan_data$z_T_anom <- joined$z_T_anom
  prepared$stan_data$z_P_anom <- joined$z_P_anom
  list(stan_data = prepared$stan_data, weekly = weekly, sero_week = sero_week, seed_idx = prepared$seed_idx, seed_dates = prepared$seed_dates)
}

make_tau_init_fn <- function(Y, td14_scalar_inits, tau_init = 0.40) {
  base_init_fn <- make_v4_9_init_fn(Y, td14_scalar_inits)
  function(chain_id = 1L) {
    v <- base_init_fn(chain_id)
    v$beta_T <- 0; v$beta_P <- 0; v$tau_year <- tau_init
    v
  }
}

message("[Phase1] Compiling the tau_year Stan model ONCE, reused for P1-P3...")
mod <- rstan::stan_model(file = STAN_SOURCE_TAU_YEAR)
message("[Phase1] Compiled OK.")

run_full_fit <- function(prior_id) {
  tau_prior_scale <- TAU_PRIORS[[prior_id]]
  out_path <- file.path(DIR_FITS_YV, sprintf("ce_tau_year_%s.rds", prior_id))
  message(sprintf("\n[Phase1 %s] FULL FIT: HalfNormal(0, %.2f) prior on tau_year, full 2015-2025 period, 2022 seed RETAINED", prior_id, tau_prior_scale))
  prepared <- build_full_period_tau_data(tau_prior_scale)
  td14_scalar_inits <- load_td14_scalar_inits(ROOT)
  init_fn <- make_tau_init_fn(prepared$stan_data$Y, td14_scalar_inits, tau_init = tau_prior_scale)

  progress_dir <- file.path(DIR_FITS_YV, paste0("chains_", prior_id)); dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
  started <- Sys.time()
  fit <- rstan::sampling(
    mod, data = prepared$stan_data, chains = HMC_SETTINGS$chains,
    iter = HMC_SETTINGS$warmup + HMC_SETTINGS$sampling, warmup = HMC_SETTINGS$warmup,
    seed = 20261013L, init = init_fn,
    control = list(adapt_delta = HMC_SETTINGS$adapt_delta, max_treedepth = HMC_SETTINGS$max_treedepth, metric = HMC_SETTINGS$metric),
    refresh = 50, sample_file = file.path(progress_dir, "chain")
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  gate_pars <- c("alpha_R", "z_year", "tau_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                 "sigma_season_year", "z_season_year", "phi_obs", "beta_T", "beta_P")
  hmc <- compute_hmc_gate(fit, gate_pars, HMC_SETTINGS$max_treedepth)

  bundle <- list(fit = fit, stan_data = prepared$stan_data, weekly_data = prepared$weekly,
                  sero_week = prepared$sero_week, seed_idx = prepared$seed_idx, seed_dates = prepared$seed_dates,
                  config = c(HMC_SETTINGS, list(prior_id = prior_id, model_version = "ceara_v4_9_climate_forced_tau_year",
                                                 tau_prior_scale = tau_prior_scale, q = QVAL, elapsed_seconds = elapsed_seconds, seed = 20261013L)),
                  hmc = hmc)
  saveRDS(bundle, out_path)
  message(sprintf("[Phase1 %s] HMC gate: %s | divergences=%d | max Rhat=%.4f | min ESS=%.0f | %.0fs | saved: %s",
                   prior_id, if (hmc$hmc_pass) "PASS" else "FAIL", hmc$divergences, hmc$maximum_rhat, hmc$minimum_bulk_ess, elapsed_seconds, out_path))
  if (!hmc$hmc_pass) warning(sprintf("[Phase1 %s] HMC gate FAILED -- reported, not hidden.", prior_id))
  invisible(bundle)
}

fit_results <- list()
for (p in names(TAU_PRIORS)) fit_results[[p]] <- run_full_fit(p)

hmc_overview <- bind_rows(lapply(names(TAU_PRIORS), function(p) {
  h <- fit_results[[p]]$hmc
  tibble::tibble(prior_id = p, tau_prior_scale = TAU_PRIORS[[p]], divergences = h$divergences, max_treedepth_hits = h$max_treedepth_hits,
                 maximum_rhat = h$maximum_rhat, minimum_bulk_ess = h$minimum_bulk_ess, hmc_pass = h$hmc_pass,
                 elapsed_seconds = fit_results[[p]]$config$elapsed_seconds)
}))
write_csv(hmc_overview, file.path(DIR_TABLES_YV, "PHASE1_hmc_diagnostics_P1_P3.csv"))
print(as.data.frame(hmc_overview), digits = 4)
message("\n[Phase1] All 3 tau_year fits complete.")
