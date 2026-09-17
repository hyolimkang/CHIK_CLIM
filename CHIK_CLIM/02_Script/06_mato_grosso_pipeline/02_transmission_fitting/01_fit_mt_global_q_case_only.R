# Mato Grosso Phase 1 -- global-q CASE-ONLY frozen v4.9 fit.
# Reuses the EXACT SAME generic Stan file already validated for Bahia/PE/RJ
# (renewal_bahia_v4_9_global_q_multisite_serology.stan is state-agnostic:
# all state-specific content enters only via `data`) -- NOT copied or
# modified. Same q prior, same HMC settings, same transmission priors as
# CE/BA/PE/RJ. J_sero=0 (case-only, Phase 1 -- no serology; MT_SEROLOGY_AUDIT.md
# found no independent MT serological anchor anyway). No seed (Section 4-5
# input audit found no qualifying MT inter-wave gap: max 184 weeks vs the
# 399-week Ceara-2022 standard). NO Pernambuco ridge/curved
# reparameterisation imported.

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script/00_shared/legacy_model_functions/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script/00_shared/legacy_model_functions/22_v4_3_short_2015_2019_q_calibration/scripts/02_fit_v4_3.R"))

base_dir <- file.path(root, "03_Output/model_fits/mato_grosso/v4_9_replication")
out_dir <- file.path(base_dir, "outputs", "caseonly")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

LOGIT_Q_PRIOR_MEAN <- qlogis(0.10) # SAME prior as CE/BA/PE/RJ -- not tuned for MT
LOGIT_Q_PRIOR_SD <- 1.0

make_mt_global_q_data <- function(weekly, imports_per_week = 1, year_effect_prior_sd = 0.40) {
  prepared <- make_v4_3_data(weekly, imports_per_week, year_effect_prior_sd,
                              t_sero = 1L, sero_pos = 0L, sero_n = 1L, kappa_sero = 50)
  n <- nrow(weekly)
  prepared$stan_data$seasonal_sin2 <- sin(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$seasonal_cos2 <- cos(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$t_sero <- NULL; prepared$stan_data$sero_pos <- NULL; prepared$stan_data$sero_n <- NULL
  prepared$stan_data$is_seed <- integer(n) # no seed -- Section 4-5 found no qualifying MT gap
  prepared$stan_data$X_seed <- numeric(n)
  prepared$stan_data$N_fit <- n
  prepared$stan_data$fit_index <- seq_len(n)
  prepared$stan_data$logit_q_prior_mean <- LOGIT_Q_PRIOR_MEAN
  prepared$stan_data$logit_q_prior_sd <- LOGIT_Q_PRIOR_SD
  # Case-only: J_sero=0 (Phase 1, no serology -- and none identified for MT anyway)
  prepared$stan_data$J_sero <- 0L
  prepared$stan_data$sero_n_positive <- integer(0)
  prepared$stan_data$sero_n_tested <- integer(0)
  prepared$stan_data$sero_window_start_idx <- integer(0)
  prepared$stan_data$sero_window_n_weeks <- integer(0)
  prepared$stan_data$sero_geographic_sd <- 1.0
  prepared
}

make_mt_init_fn <- function(Y, td14_scalar_inits) {
  # Dispersed but prior-valid: q centred on the prior mean with modest
  # offsets (same convention validated for CE/BA/PE/RJ's global-q models).
  offsets <- c(-0.06, -0.02, 0.02, 0.06)
  function(chain_id = 1L) {
    cid <- as.integer(chain_id)
    base <- td14_scalar_inits[[(cid - 1L) %% length(td14_scalar_inits) + 1L]]
    offset <- offsets[(cid - 1L) %% 4L + 1L]
    list(
      alpha_R = base$alpha_R, phi_obs = base$phi_obs,
      beta_sin1 = base$beta_sin, beta_cos1 = base$beta_cos,
      beta_sin2 = 0, beta_cos2 = 0,
      sigma_season_year = 0.05, z_season_year = rep(0, Y),
      z_year = rep(offset / 0.40, Y),
      z_geo = numeric(0),
      logit_q = LOGIT_Q_PRIOR_MEAN + offset
    )
  }
}

run_fit_mt_global_q_case_only <- function() {
  settings <- list(warmup = 1200L, sampling = 750L, chains = 4L, adapt_delta = 0.95, max_treedepth = 14L, metric = "dense_e")
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  weekly <- readRDS(file.path(root, "03_Output/tables/mato_grosso_v4_9_replication/mato_grosso_weekly_input.rds"))
  prepared <- make_mt_global_q_data(weekly)
  Y <- prepared$stan_data$Y
  message("MT case-only Stan data: N=", prepared$stan_data$N, " Y=", Y, " J_sero=", prepared$stan_data$J_sero)

  scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_mt_init_fn(Y, scalar_inits)
  message("Initial values (per chain):")
  for (i in 1:4) {
    v <- init_fn(i)
    message(sprintf("  chain %d: logit_q=%.4f (q=%.4f), alpha_R=%.4f, phi_obs=%.4f",
                     i, v$logit_q, plogis(v$logit_q), v$alpha_R, v$phi_obs))
  }

  stan_path <- file.path(root, "02_Script/stan/renewal_bahia_v4_9_global_q_multisite_serology.stan")
  sample_file <- file.path(progress_dir, "chain")

  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data, chains = settings$chains,
    iter = settings$warmup + settings$sampling, warmup = settings$warmup,
    seed = 20260915L, init = init_fn,
    control = list(adapt_delta = settings$adapt_delta, max_treedepth = settings$max_treedepth, metric = settings$metric),
    refresh = 25, sample_file = sample_file
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  gate_pars <- c("alpha_R", "z_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                 "sigma_season_year", "z_season_year", "phi_obs", "logit_q")
  hmc <- compute_hmc_gate(fit, gate_pars, settings$max_treedepth)

  bundle <- list(fit = fit, stan_data = prepared$stan_data, weekly_data = weekly,
                  config = c(settings, list(model_version = "mato_grosso_v4_9_global_q_case_only",
                                             elapsed_seconds = elapsed_seconds, seed = 20260915L,
                                             logit_q_prior_mean = LOGIT_Q_PRIOR_MEAN, logit_q_prior_sd = LOGIT_Q_PRIOR_SD)),
                  hmc = hmc)
  fit_path <- file.path(out_dir, "mt_global_q_case_only.rds")
  saveRDS(bundle, fit_path)
  message(sprintf("[mt case-only] HMC gate: %s | %.0fs | %s", if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_mt_global_q_case_only()
