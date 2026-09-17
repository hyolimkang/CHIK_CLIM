# Rio de Janeiro Phase 1 -- rescue attempt for the 4-divergence marginal
# HMC fail: SAME warmup/sampling (1200/750), SAME model/priors/inits, only
# adapt_delta raised 0.95 -> 0.98. Justified because the original fit's
# divergences were few (4/3000, 0.13%), with zero treedepth-14 hits,
# max Rhat=1.003, min ESS~1000, healthy BFMI, no chain-specific collapse --
# the established "benign noise, adaptation-resolvable" signature in this
# project, not a structural geometry failure.

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
source(file.path(root, "02_Script/05_rio_de_janeiro_pipeline/02_transmission_fitting/01_fit_rj_global_q_case_only.R"))

base_dir <- file.path(root, "03_Output/model_fits/rio_de_janeiro/v4_9_replication")
out_dir <- file.path(base_dir, "outputs", "caseonly_ad098")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_rescue <- function() {
  settings <- list(warmup = 1200L, sampling = 750L, chains = 4L, adapt_delta = 0.98, max_treedepth = 14L, metric = "dense_e")
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  weekly <- readRDS(file.path(root, "03_Output/tables/rio_de_janeiro_v4_9_replication/rio_de_janeiro_weekly_input.rds"))
  prepared <- make_rj_global_q_data(weekly)
  Y <- prepared$stan_data$Y
  scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_rj_init_fn(Y, scalar_inits)

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
                  config = c(settings, list(model_version = "rio_de_janeiro_v4_9_global_q_case_only_ad098",
                                             elapsed_seconds = elapsed_seconds, seed = 20260915L)),
                  hmc = hmc)
  fit_path <- file.path(out_dir, "rj_global_q_case_only_ad098.rds")
  saveRDS(bundle, fit_path)
  message(sprintf("[rj case-only, adapt_delta=0.98] HMC gate: %s | %.0fs | %s", if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_rescue()
