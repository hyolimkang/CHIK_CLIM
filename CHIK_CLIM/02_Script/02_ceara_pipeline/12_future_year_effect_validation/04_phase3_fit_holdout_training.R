# PHASE 3 -- Fit the tau_year model on TRAINING-ONLY data for H1 and H2.
#
# Uses the dedicated seed-free builder (03_build_training_only_data.R) --
# NEVER the accepted seed-requiring builder. The 2022 seed is not present
# in either training window by construction (H1 ends 2019, H2 ends 2021).
#
# Efficiency/scope decision (disclosed in the final report): rather than
# fitting FOUR separate training models (fixed-SD x {H1,H2} and
# tau_year x {H1,H2}), each holdout trains ONE tau_year model; the
# "fixed-SD generator" comparison arm in Phase 4 is produced by sampling
# FUTURE year effects with sigma=0.40 directly on top of this SAME
# training posterior's other parameters (alpha_R, seasonal, climate,
# sigma_season_year, phi_obs, terminal state) -- i.e. the transmission-
# parameter posterior is shared and matched between the two generator
# arms, isolating the year-effect-generation-mechanism contrast exactly
# as Phase 2 does for the full-period comparison. This halves the
# training-fit cost (2 fits instead of 4) without weakening the
# comparison requested.

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

source(file.path(ROOT, "02_Script/02_ceara_pipeline/12_future_year_effect_validation/00_config_year_effect_validation.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/22_v4_3_short_2015_2019_q_calibration/scripts/02_fit_v4_3.R"))
source(file.path(ROOT, "02_Script/02_ceara_pipeline/12_future_year_effect_validation/03_build_training_only_data.R"))

TRAINING_TAU_PRIOR_SCALE <- TAU_PRIORS$P2 # 0.40 -- matches the original fixed-SD scale, the natural default; revisit after Phase 1 results if warranted

message("[Phase3 fit] Compiling the tau_year Stan model (reused from Phase 1 cache if available)...")
mod <- rstan::stan_model(file = STAN_SOURCE_TAU_YEAR)
message("[Phase3 fit] Compiled OK.")

make_training_init_fn <- function(Y, td14_scalar_inits, tau_init = 0.40) {
  offsets <- c(-0.06, -0.02, 0.02, 0.06)
  function(chain_id = 1L) {
    cid <- as.integer(chain_id)
    base <- td14_scalar_inits[[(cid - 1L) %% length(td14_scalar_inits) + 1L]]
    offset <- offsets[(cid - 1L) %% 4L + 1L]
    list(alpha_R = base$alpha_R, phi_obs = base$phi_obs, beta_sin1 = base$beta_sin, beta_cos1 = base$beta_cos,
         beta_sin2 = 0, beta_cos2 = 0, sigma_season_year = 0.05, z_season_year = rep(0, Y),
         z_year = rep(offset / tau_init, Y), tau_year = tau_init, beta_T = 0, beta_P = 0)
  }
}

run_training_fit <- function(holdout_id) {
  h <- HOLDOUTS[[holdout_id]]
  out_path <- file.path(DIR_FITS_YV, sprintf("ce_tau_year_training_%s.rds", holdout_id))
  message(sprintf("\n[Phase3 fit %s] TRAINING-ONLY FIT: %s | tau prior HalfNormal(0,%.2f) | NO 2022 seed (seed-free builder)",
                   holdout_id, h$label, TRAINING_TAU_PRIOR_SCALE))
  prepared <- build_training_only_tau_data(h$train_start, h$train_end, TRAINING_TAU_PRIOR_SCALE)
  message(sprintf("  N=%d weeks, Y=%d years, serology included=%s, is_seed sum=%d (MUST be 0)",
                   prepared$n_weeks, prepared$stan_data$Y, prepared$serology_included, sum(prepared$stan_data$is_seed)))
  stopifnot(sum(prepared$stan_data$is_seed) == 0L, prepared$stan_data$N_fit == prepared$n_weeks)

  td14_scalar_inits <- load_td14_scalar_inits(ROOT)
  init_fn <- make_training_init_fn(prepared$stan_data$Y, td14_scalar_inits, tau_init = TRAINING_TAU_PRIOR_SCALE)

  progress_dir <- file.path(DIR_FITS_YV, paste0("chains_training_", holdout_id)); dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
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

  bundle <- list(fit = fit, stan_data = prepared$stan_data, weekly_data = prepared$weekly, sero_week = prepared$sero_week,
                  config = c(HMC_SETTINGS, list(holdout_id = holdout_id, model_version = "ceara_v4_9_climate_forced_tau_year_TRAINING_ONLY",
                                                 tau_prior_scale = TRAINING_TAU_PRIOR_SCALE, q = QVAL, train_start = as.character(h$train_start),
                                                 train_end = as.character(h$train_end), n_weeks = prepared$n_weeks, elapsed_seconds = elapsed_seconds, seed = 20261013L)),
                  hmc = hmc)
  saveRDS(bundle, out_path)
  message(sprintf("[Phase3 fit %s] HMC gate: %s | divergences=%d | max Rhat=%.4f | min ESS=%.0f | %.0fs | saved: %s",
                   holdout_id, if (hmc$hmc_pass) "PASS" else "FAIL", hmc$divergences, hmc$maximum_rhat, hmc$minimum_bulk_ess, elapsed_seconds, out_path))
  if (!hmc$hmc_pass) warning(sprintf("[Phase3 fit %s] HMC gate FAILED -- reported, not hidden.", holdout_id))
  invisible(bundle)
}

fit_results <- list()
for (h in names(HOLDOUTS)) fit_results[[h]] <- run_training_fit(h)

hmc_overview <- bind_rows(lapply(names(HOLDOUTS), function(h) {
  hh <- fit_results[[h]]$hmc
  tibble::tibble(holdout_id = h, divergences = hh$divergences, max_treedepth_hits = hh$max_treedepth_hits,
                 maximum_rhat = hh$maximum_rhat, minimum_bulk_ess = hh$minimum_bulk_ess, hmc_pass = hh$hmc_pass,
                 elapsed_seconds = fit_results[[h]]$config$elapsed_seconds)
}))
write_csv(hmc_overview, file.path(DIR_TABLES_YV, "PHASE3_hmc_diagnostics_training_fits.csv"))
print(as.data.frame(hmc_overview), digits = 4)
message("\n[Phase3 fit] Both training-only fits complete.")
