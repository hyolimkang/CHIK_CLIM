# Pernambuco Spatial-0 -- branch-connectivity HMC test (revised spec).
#
# 6 chains, 1000 warmup, 500 sampling, adapt_delta=0.95, max_treedepth=14,
# dense_e. Spatial-1 revealed TWO distinct branches (Chain 1: low-q/
# high-depletion/high-R0; Chains 2-4: high-q/near-fully-susceptible/
# low-R0). Rather than assume either is right, Spatial-0 deliberately
# starts 2 chains in each of three classes -- LOW (near Spatial-1 Chain 1's
# mode), INTERMEDIATE (central prior-supported region), HIGH (near
# Spatial-1 Chains 2-4's mode) -- to test whether Spatial-0's simpler
# (no sigma_region) geometry lets all chains converge to ONE posterior
# region, or whether the branch separation persists.

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(6L, parallel::detectCores()))

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "40_renewal_model", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script", "40_renewal_model", "22_v4_3_short_2015_2019_q_calibration", "scripts", "02_fit_v4_3.R"))
source(file.path(root, "02_Script", "40_renewal_model", "33_pernambuco_v4_9_replication", "scripts", "36_build_pe_spatial0_stan_data.R"))

base_dir <- file.path(root, "02_Script", "40_renewal_model", "33_pernambuco_v4_9_replication")
out_dir <- file.path(base_dir, "outputs", "spatial0")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

LOGIT_Q_PRIOR_MEAN <- qlogis(0.10)

built <- build_pe_spatial0_stan_data()
stan_data <- built$stan_data
Y <- stan_data$Y

# Branch-representative starting points, taken DIRECTLY from the Spatial-1
# geometry diagnostic (PE_spatial1_geometry_summary.csv), not guessed:
#   LOW branch  <- Chain 1's mode (median logit_q=-3.199, alpha_global=0.185)
#   HIGH branch <- Chains 2-4's mode (median logit_q=-0.462, alpha_global=-0.009)
#   INTERMEDIATE <- central prior-supported region (prior means)
BRANCH_STARTS <- list(
  low  = list(logit_q = -3.199, alpha_R = 0.185),
  mid  = list(logit_q = LOGIT_Q_PRIOR_MEAN, alpha_R = log(1.2)),
  high = list(logit_q = -0.462, alpha_R = -0.009)
)
CHAIN_BRANCH <- c("low", "low", "mid", "mid", "high", "high")

make_init_fn <- function() {
  within_branch_offset <- c(-0.02, 0.02) # small jitter so the 2 chains per branch are not identical
  scalar_inits <- load_td14_scalar_inits(root, n_chains = 6L)
  function(chain_id = 1L) {
    cid <- as.integer(chain_id)
    branch <- CHAIN_BRANCH[cid]
    start <- BRANCH_STARTS[[branch]]
    jitter <- within_branch_offset[(cid - 1L) %% 2L + 1L]
    base <- scalar_inits[[cid]]
    list(
      alpha_R = start$alpha_R + jitter,
      z_year_contrast = rep(jitter / 0.40, Y - 1L),
      beta_sin1 = base$beta_sin, beta_cos1 = base$beta_cos,
      beta_sin2 = 0, beta_cos2 = 0,
      sigma_season_year = 0.05, z_season_year = rep(0, Y),
      phi_obs = base$phi_obs,
      logit_q = start$logit_q + jitter
    )
  }
}

run_spatial0_pilot <- function() {
  settings <- list(warmup = 1000L, sampling = 500L, chains = 6L, adapt_delta = 0.95, max_treedepth = 14L, metric = "dense_e")
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
  stan_path <- file.path(root, "02_Script", "stan", "renewal_pernambuco_v4_9_spatial0.stan")
  sample_file <- file.path(progress_dir, "chain")

  init_fn <- make_init_fn()
  message("Initial values (per chain):")
  for (i in 1:6) {
    v <- init_fn(i)
    message(sprintf("  chain %d [%s branch]: logit_q=%.4f (q=%.4f), alpha_R=%.4f, phi_obs=%.4f",
                     i, CHAIN_BRANCH[i], v$logit_q, plogis(v$logit_q), v$alpha_R, v$phi_obs))
  }

  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = stan_data, chains = settings$chains,
    iter = settings$warmup + settings$sampling, warmup = settings$warmup,
    seed = 20260914L, init = init_fn,
    control = list(adapt_delta = settings$adapt_delta, max_treedepth = settings$max_treedepth, metric = settings$metric),
    refresh = 25, sample_file = sample_file
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  gate_pars <- c("alpha_R", "z_year_contrast", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                 "sigma_season_year", "z_season_year", "phi_obs", "logit_q")
  hmc <- compute_hmc_gate(fit, gate_pars, settings$max_treedepth)

  bundle <- list(fit = fit, stan_data = stan_data, stratum_levels = built$stratum_levels,
                  recife_index = built$recife_index, sero_window = built$sero_window,
                  config = c(settings, list(model_version = "pernambuco_v4_9_spatial0", elapsed_seconds = elapsed_seconds)),
                  hmc = hmc)
  fit_path <- file.path(out_dir, "pe_spatial0_pilot.rds")
  saveRDS(bundle, fit_path)
  message(sprintf("[spatial0 pilot] HMC gate: %s | %.0fs | %s", if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_spatial0_pilot()
