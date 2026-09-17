# Mato Grosso -- targeted fixed-q identification profile. NOT a new model:
# reuses the EXACT SAME unmodified global-q Stan file
# (renewal_bahia_v4_9_global_q_multisite_serology.stan, J_sero=0, same as
# the Phase-1 fit). "Fixed q" is achieved with ZERO Stan-file changes by
# setting an extremely tight logit-normal prior (sd=0.02) centred at each
# target q -- a continuous relaxation of a hard constraint, not a new
# model (same convention as
# 34_rio_de_janeiro_v4_9_replication/scripts/06_fit_rj_targeted_fixed_q.R).
# Grid is FIXED per instruction: 0.020, 0.030, 0.040, 0.050, 0.070 -- not
# to be changed after inspecting results.

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
source(file.path(root, "02_Script/06_mato_grosso_pipeline/02_transmission_fitting/01_fit_mt_global_q_case_only.R"))

base_dir <- file.path(root, "03_Output/model_fits/mato_grosso/v4_9_replication")
out_dir <- file.path(base_dir, "outputs", "fixedq")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(root, "03_Output/tables/mato_grosso_v4_9_replication")

Q_GRID <- c(0.020, 0.030, 0.040, 0.050, 0.070) # FIXED per instruction -- do not change
FIXED_Q_PRIOR_SD <- 0.02 # tight enough to pin q to within ~+/-2% multiplicatively

weekly <- readRDS(file.path(table_dir, "mato_grosso_weekly_input.rds"))
stan_path <- file.path(root, "02_Script/stan/renewal_bahia_v4_9_global_q_multisite_serology.stan")

fit_one_q <- function(q_target) {
  prepared <- make_mt_global_q_data(weekly)
  prepared$stan_data$logit_q_prior_mean <- qlogis(q_target)
  prepared$stan_data$logit_q_prior_sd <- FIXED_Q_PRIOR_SD
  Y <- prepared$stan_data$Y

  scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- function(chain_id = 1L) {
    cid <- as.integer(chain_id)
    base <- scalar_inits[[(cid - 1L) %% length(scalar_inits) + 1L]]
    offsets <- c(-0.06, -0.02, 0.02, 0.06)
    offset <- offsets[(cid - 1L) %% 4L + 1L]
    list(alpha_R = base$alpha_R, phi_obs = base$phi_obs, beta_sin1 = base$beta_sin, beta_cos1 = base$beta_cos,
         beta_sin2 = 0, beta_cos2 = 0, sigma_season_year = 0.05, z_season_year = rep(0, Y),
         z_year = rep(offset / 0.40, Y), z_geo = numeric(0), logit_q = qlogis(q_target))
  }

  tag <- sprintf("q%.3f", q_target)
  progress_dir <- file.path(out_dir, "chains", tag)
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("\n===== Fitting MT fixed-q (targeted) q=%.3f =====", q_target))
  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data, chains = 4, iter = 1200 + 750, warmup = 1200,
    seed = 20260915L, init = init_fn,
    control = list(adapt_delta = 0.98, max_treedepth = 14, metric = "dense_e"),
    refresh = 0, sample_file = file.path(progress_dir, "chain")
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                                  "sigma_season_year", "z_season_year", "phi_obs", "logit_q"), 14L)
  bundle <- list(fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, q_target = q_target,
                  config = list(elapsed_seconds = elapsed_seconds, adapt_delta = 0.98, warmup = 1200L, sampling = 750L, chains = 4L, max_treedepth = 14L),
                  hmc = hmc)
  saveRDS(bundle, file.path(out_dir, sprintf("mt_fixedq_%s.rds", tag)))
  message(sprintf("[q=%.3f] HMC gate: %s | divergences=%d | max Rhat=%.4f | %.0fs",
                   q_target, if (hmc$hmc_pass) "PASS" else "FAIL", hmc$divergences, hmc$maximum_rhat, elapsed_seconds))
  invisible(bundle)
}

if (sys.nframe() == 0L) {
  for (q in Q_GRID) fit_one_q(q)
  message("\n===== MT fixed-q profile: all fits complete =====")
}
