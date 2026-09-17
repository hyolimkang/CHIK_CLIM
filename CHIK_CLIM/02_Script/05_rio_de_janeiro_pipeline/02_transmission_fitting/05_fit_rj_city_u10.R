# Rio de Janeiro CITY auxiliary -- Model B: global-q CASE + U10 serology,
# frozen v4.9 transmission structure, DESIGN-AWARE logit-normal serology
# likelihood on U10's WEIGHTED prevalence (not naive beta-binomial on raw
# counts), geography matched (NO eta_geo). Uses the new
# renewal_riodejaneiro_city_v4_9_global_q_logitnormal_serology.stan file.
#
# Dispersed prior-valid inits (0.05/0.10/0.15/0.20), explicitly NOT
# centred at q=0.010, and NOT using the state posterior as an informative
# prior (same broad logit_q ~ Normal(logit(0.10), 1) prior family).

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

base_dir <- file.path(root, "03_Output/05_rio_de_janeiro_pipeline/model_fits/v4_9_replication")
out_dir <- file.path(base_dir, "outputs", "city_u10")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(root, "03_Output/05_rio_de_janeiro_pipeline/tables/rio_de_janeiro_v4_9_replication/city")

LOGIT_Q_PRIOR_MEAN <- qlogis(0.10)
LOGIT_Q_PRIOR_SD <- 1.0

U10_WEIGHTED_PREV <- 0.180
U10_CI_LO <- 0.148
U10_CI_HI <- 0.212
U10_SURVEY_START <- as.Date("2018-07-01")
U10_SURVEY_END <- as.Date("2018-10-31")

# SE on the logit scale, derived from the reported 95% CI via the standard
# delta-method approximation: SE_logit = (logit(hi) - logit(lo)) / (2*1.96).
se_logit_u10 <- (qlogis(U10_CI_HI) - qlogis(U10_CI_LO)) / (2 * 1.96)
message(sprintf("Derived SE_logit_U10 = (logit(%.3f) - logit(%.3f)) / 3.92 = (%.4f - %.4f) / 3.92 = %.4f",
                 U10_CI_HI, U10_CI_LO, qlogis(U10_CI_HI), qlogis(U10_CI_LO), se_logit_u10))

make_rj_city_u10_data <- function(weekly, imports_per_week = 1, year_effect_prior_sd = 0.40) {
  prepared <- make_v4_3_data(weekly, imports_per_week, year_effect_prior_sd,
                              t_sero = 1L, sero_pos = 0L, sero_n = 1L, kappa_sero = 50)
  n <- nrow(weekly)
  prepared$stan_data$seasonal_sin2 <- sin(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$seasonal_cos2 <- cos(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$t_sero <- NULL; prepared$stan_data$sero_pos <- NULL; prepared$stan_data$sero_n <- NULL
  prepared$stan_data$is_seed <- integer(n)
  prepared$stan_data$X_seed <- numeric(n)
  prepared$stan_data$N_fit <- n
  prepared$stan_data$fit_index <- seq_len(n)
  prepared$stan_data$logit_q_prior_mean <- LOGIT_Q_PRIOR_MEAN
  prepared$stan_data$logit_q_prior_sd <- LOGIT_Q_PRIOR_SD

  dates <- as.Date(weekly$week_start)
  window_idx <- which(dates >= U10_SURVEY_START & dates <= U10_SURVEY_END)
  stopifnot(length(window_idx) > 0, all(diff(window_idx) == 1L))
  prepared$stan_data$J_sero <- 1L
  prepared$stan_data$p_hat_sero <- as.array(U10_WEIGHTED_PREV)
  prepared$stan_data$se_logit_sero <- as.array(se_logit_u10)
  prepared$stan_data$sero_window_start_idx <- as.array(min(window_idx))
  prepared$stan_data$sero_window_n_weeks <- as.array(length(window_idx))
  prepared$window_idx <- window_idx
  prepared
}

make_dispersed_init_fn <- function(Y, td14_scalar_inits) {
  q_starts <- c(0.05, 0.10, 0.15, 0.20) # SAME dispersion as Model A -- NOT centred at 0.010
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
      logit_q = qlogis(q_starts[cid])
    )
  }
}

run_fit_rj_city_u10 <- function() {
  settings <- list(warmup = 1200L, sampling = 1000L, chains = 4L, adapt_delta = 0.98, max_treedepth = 14L, metric = "dense_e")
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  weekly <- readRDS(file.path(table_dir, "riodejaneiro_city_weekly_input.rds"))
  prepared <- make_rj_city_u10_data(weekly)
  Y <- prepared$stan_data$Y
  message("RJ CITY + U10 Stan data: N=", prepared$stan_data$N, " Y=", Y, " J_sero=", prepared$stan_data$J_sero,
          " sero_window=[", prepared$stan_data$sero_window_start_idx, ",",
          prepared$stan_data$sero_window_start_idx + prepared$stan_data$sero_window_n_weeks - 1, "]")

  scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_dispersed_init_fn(Y, scalar_inits)
  message("Initial values (per chain, dispersed 0.05/0.10/0.15/0.20):")
  for (i in 1:4) {
    v <- init_fn(i)
    message(sprintf("  chain %d: logit_q=%.4f (q=%.4f), alpha_R=%.4f", i, v$logit_q, plogis(v$logit_q), v$alpha_R))
  }

  stan_path <- file.path(root, "02_Script/00_shared/stan/current/rio_de_janeiro/renewal_riodejaneiro_city_v4_9_global_q_logitnormal_serology.stan")
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

  bundle <- list(fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, window_idx = prepared$window_idx,
                  u10 = list(weighted_prev = U10_WEIGHTED_PREV, ci_lo = U10_CI_LO, ci_hi = U10_CI_HI, se_logit = se_logit_u10),
                  config = c(settings, list(model_version = "rio_de_janeiro_CITY_v4_9_global_q_u10_logitnormal",
                                             elapsed_seconds = elapsed_seconds, seed = 20260915L)),
                  hmc = hmc)
  fit_path <- file.path(out_dir, "rj_city_u10.rds")
  saveRDS(bundle, fit_path)
  message(sprintf("[rj city + U10] HMC gate: %s | %.0fs | %s", if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_rj_city_u10()
