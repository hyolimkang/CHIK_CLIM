# PHASE 2 -- Climate robustness fits M1-M4.
#
# M0 (accepted, primary lag, prior SD=0.40) is READ-ONLY and never
# refit. M1-M4 reuse the SAME accepted Stan source
# (renewal_ceara_v4_9_climate_forced_canary.stan, byte-identical, no
# changes) -- only the DATA differs (z_T_anom/z_P_anom from a different
# lag spec for M1/M2, or year_effect_prior_sd for M3/M4). Reuses the
# exact accepted data-prep chain (make_v4_9_data() etc.) unmodified.
# Same HMC settings as M0: 4 chains, 800 warmup + 400 sampling,
# adapt_delta=0.95, max_treedepth=14, dense_e metric.
#
# Run Phase 1 first and confirm PASS before running this script.

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

source(file.path(ROOT, "02_Script/02_ceara_pipeline/11_climate_robustness/00_config_climate_robustness.R"))

# ---- Reuse the EXACT accepted data-prep chain (unmodified source files) ----
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/22_v4_3_short_2015_2019_q_calibration/scripts/02_fit_v4_3.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/27_v4_8_seeded_recurrence/scripts/01_fit_v4_8.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/28_v4_9_hierarchical_seasonality/scripts/01_fit_v4_9.R"))

build_climate_data_for_spec <- function(lag_spec, year_effect_prior_sd) {
  panel <- readRDS(file.path(ROOT, "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |> dplyr::filter(week_start >= FIXED$date_start, week_start <= FIXED$date_end)
  weekly_demography <- readRDS(file.path(ROOT, "01_Data/ceara_weekly_demography_2015_2025.rds")) |>
    dplyr::filter(week_start >= FIXED$date_start, week_start <= FIXED$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch.")

  sero_week <- find_sero_week_index(as.Date(weekly$week_start))
  prepared <- make_v4_9_data(weekly, FIXED$imports_per_week, year_effect_prior_sd,
                              sero_week$index, FIXED$sero_pos, FIXED$sero_n, FIXED$kappa_sero, QVAL,
                              SEED_START_DATE, SEED_LENGTH)

  anom <- read_csv(ACCEPTED_ANOMALY_CSV, show_col_types = FALSE) |>
    dplyr::filter(state == "CE", lag_spec == !!lag_spec) |> mutate(week_start = as.Date(week_start))
  joined <- tibble::tibble(week_start = as.Date(weekly$week_start)) |> left_join(anom |> select(week_start, z_T_anom, z_P_anom), by = "week_start")
  if (nrow(joined) != nrow(weekly)) stop("Climate anomaly join produced wrong row count.")
  n_na <- sum(is.na(joined$z_T_anom) | is.na(joined$z_P_anom))
  message(sprintf("  [%s] %d weeks, %d with NA anomaly (boundary weeks, set to 0)", lag_spec, nrow(joined), n_na))
  joined$z_T_anom[is.na(joined$z_T_anom)] <- 0
  joined$z_P_anom[is.na(joined$z_P_anom)] <- 0

  prepared$stan_data$z_T_anom <- joined$z_T_anom
  prepared$stan_data$z_P_anom <- joined$z_P_anom
  list(stan_data = prepared$stan_data, weekly = weekly, sero_week = sero_week, seed_idx = prepared$seed_idx, seed_dates = prepared$seed_dates)
}

make_climate_init_fn <- function(Y, td14_scalar_inits) {
  base_init_fn <- make_v4_9_init_fn(Y, td14_scalar_inits)
  function(chain_id = 1L) { v <- base_init_fn(chain_id); v$beta_T <- 0; v$beta_P <- 0; v }
}

message("[Phase2] Compiling the accepted (unmodified) climate-forced Stan model ONCE, reused for M1-M4...")
mod <- rstan::stan_model(file = STAN_SOURCE_CLIMATE_M0)
message("[Phase2] Compiled OK.")

run_smoke_test <- function(model_id) {
  spec <- MODEL_SPECS[[model_id]]
  message(sprintf("\n[Phase2 %s] Smoke test: lag_spec=%s, year_effect_prior_sd=%.2f", model_id, spec$lag_spec, spec$year_effect_prior_sd))
  prepared <- build_climate_data_for_spec(spec$lag_spec, spec$year_effect_prior_sd)
  td14_scalar_inits <- load_td14_scalar_inits(ROOT)
  init_fn <- make_climate_init_fn(prepared$stan_data$Y, td14_scalar_inits)
  fit <- rstan::sampling(mod, data = prepared$stan_data, chains = 1L, iter = 40L, warmup = 20L, init = init_fn, seed = 1L,
                          control = list(adapt_delta = 0.95, max_treedepth = 10, metric = "dense_e"), refresh = 0)
  message(sprintf("[Phase2 %s] Smoke test OK -- no dimension/parsing errors.", model_id))
  invisible(TRUE)
}

run_full_fit <- function(model_id) {
  spec <- MODEL_SPECS[[model_id]]
  out_path <- file.path(DIR_FITS_CR, sprintf("ce_climate_robustness_%s.rds", model_id))
  message(sprintf("\n[Phase2 %s] FULL FIT: %s", model_id, spec$label))
  prepared <- build_climate_data_for_spec(spec$lag_spec, spec$year_effect_prior_sd)
  td14_scalar_inits <- load_td14_scalar_inits(ROOT)
  init_fn <- make_climate_init_fn(prepared$stan_data$Y, td14_scalar_inits)

  progress_dir <- file.path(DIR_FITS_CR, paste0("chains_", model_id)); dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
  started <- Sys.time()
  fit <- rstan::sampling(
    mod, data = prepared$stan_data, chains = HMC_SETTINGS$chains,
    iter = HMC_SETTINGS$warmup + HMC_SETTINGS$sampling, warmup = HMC_SETTINGS$warmup,
    seed = 20261013L, init = init_fn,
    control = list(adapt_delta = HMC_SETTINGS$adapt_delta, max_treedepth = HMC_SETTINGS$max_treedepth, metric = HMC_SETTINGS$metric),
    refresh = 50, sample_file = file.path(progress_dir, "chain")
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  gate_pars <- c("alpha_R", "z_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                 "sigma_season_year", "z_season_year", "phi_obs", "beta_T", "beta_P")
  hmc <- compute_hmc_gate(fit, gate_pars, HMC_SETTINGS$max_treedepth)

  bundle <- list(fit = fit, stan_data = prepared$stan_data, weekly_data = prepared$weekly,
                  sero_week = prepared$sero_week, seed_idx = prepared$seed_idx, seed_dates = prepared$seed_dates,
                  config = c(HMC_SETTINGS, list(model_id = model_id, model_version = "ceara_v4_9_climate_forced_canary",
                                                 lag_spec = spec$lag_spec, year_effect_prior_sd = spec$year_effect_prior_sd,
                                                 q = QVAL, elapsed_seconds = elapsed_seconds, seed = 20261013L)),
                  hmc = hmc)
  saveRDS(bundle, out_path)
  message(sprintf("[Phase2 %s] HMC gate: %s | divergences=%d | max Rhat=%.4f | min ESS=%.0f | %.0fs | saved: %s",
                   model_id, if (hmc$hmc_pass) "PASS" else "FAIL", hmc$divergences, hmc$maximum_rhat, hmc$minimum_bulk_ess, elapsed_seconds, out_path))
  if (!hmc$hmc_pass) warning(sprintf("[Phase2 %s] HMC gate FAILED -- do not treat this fit as an accepted result. Reported, not hidden.", model_id))
  invisible(bundle)
}

NEW_MODEL_IDS <- setdiff(names(MODEL_SPECS), "M0")

message("\n[Phase2] === SMOKE TESTS (all 4 new specs) ===")
for (m in NEW_MODEL_IDS) run_smoke_test(m)
message("\n[Phase2] All smoke tests passed -- no dimension/parsing errors across M1-M4. Proceeding to full fits.")

message("\n[Phase2] === FULL FITS (sequential, ~14 min each based on M0's precedent) ===")
fit_results <- list()
for (m in NEW_MODEL_IDS) fit_results[[m]] <- run_full_fit(m)

message("\n[Phase2] All full fits complete.")
hmc_overview <- bind_rows(lapply(NEW_MODEL_IDS, function(m) {
  h <- fit_results[[m]]$hmc
  tibble::tibble(model_id = m, divergences = h$divergences, max_treedepth_hits = h$max_treedepth_hits,
                 maximum_rhat = h$maximum_rhat, minimum_bulk_ess = h$minimum_bulk_ess, hmc_pass = h$hmc_pass,
                 elapsed_seconds = fit_results[[m]]$config$elapsed_seconds)
}))
write_csv(hmc_overview, file.path(DIR_TABLES_CR, "PHASE2_hmc_diagnostics_M1_M4.csv"))
print(as.data.frame(hmc_overview), digits = 4)
message("\n[Phase2] HMC diagnostics saved -> ", file.path(DIR_TABLES_CR, "PHASE2_hmc_diagnostics_M1_M4.csv"))
