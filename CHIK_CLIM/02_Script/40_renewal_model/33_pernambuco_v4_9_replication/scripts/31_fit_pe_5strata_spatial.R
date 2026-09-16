# Pernambuco Spatial Model 1 -- Section 12: initial HMC pilot.
#
# 5-stratum extension of the frozen v4.9 global-q architecture. Fresh
# parameterisation (Section 11): does NOT import the homogeneous model's
# affine/quadratic ridge coordinates -- those characterised a different
# (unstratified) posterior geometry.
#
# 4 chains, 1000 warmup, 500 sampling, adapt_delta=0.95, max_treedepth=14,
# dense_e, dispersed valid inits. Per instruction: NOT a 300-warmup pilot.

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "40_renewal_model", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script", "40_renewal_model", "22_v4_3_short_2015_2019_q_calibration", "scripts", "02_fit_v4_3.R"))
source(file.path(root, "02_Script", "40_renewal_model", "33_pernambuco_v4_9_replication", "scripts", "30_build_pe_5strata_stan_data.R"))

base_dir <- file.path(root, "02_Script", "40_renewal_model", "33_pernambuco_v4_9_replication")
out_dir <- file.path(base_dir, "outputs", "spatial5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

built <- build_pe_5strata_stan_data()
stan_data <- built$stan_data
R <- stan_data$R; Y <- stan_data$Y

make_init_fn <- function() {
  q_starts <- c(0.02, 0.05, 0.10, 0.20)
  offsets <- c(-0.06, -0.02, 0.02, 0.06)
  scalar_inits <- load_td14_scalar_inits(root)
  function(chain_id = 1L) {
    cid <- (as.integer(chain_id) - 1L) %% 4L + 1L
    base <- scalar_inits[[cid]]
    list(
      alpha_global = log(1.2) + offsets[cid],
      z_region_contrast = rep(0, R - 1L),
      sigma_region = 0.05,
      z_year_contrast = rep(0, Y - 1L),
      beta_sin1 = base$beta_sin, beta_cos1 = base$beta_cos,
      beta_sin2 = 0, beta_cos2 = 0,
      sigma_season_year = 0.05, z_season_year = rep(0, Y),
      phi_obs = base$phi_obs,
      logit_q = qlogis(q_starts[cid])
    )
  }
}

run_pilot <- function() {
  settings <- list(warmup = 1000L, sampling = 500L, chains = 4L, adapt_delta = 0.95, max_treedepth = 14L, metric = "dense_e")
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
  stan_path <- file.path(root, "02_Script", "stan", "renewal_pernambuco_v4_9_5strata_spatial.stan")
  sample_file <- file.path(progress_dir, "chain")

  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = stan_data, chains = settings$chains,
    iter = settings$warmup + settings$sampling, warmup = settings$warmup,
    seed = 20260914L, init = make_init_fn(),
    control = list(adapt_delta = settings$adapt_delta, max_treedepth = settings$max_treedepth, metric = settings$metric),
    refresh = 25, sample_file = sample_file
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  gate_pars <- c("alpha_global", "delta_region", "sigma_region", "z_year_contrast", "beta_sin1", "beta_cos1",
                 "beta_sin2", "beta_cos2", "sigma_season_year", "z_season_year", "phi_obs", "logit_q")
  hmc <- compute_hmc_gate(fit, gate_pars, settings$max_treedepth)

  bundle <- list(fit = fit, stan_data = stan_data, stratum_levels = built$stratum_levels,
                  recife_index = built$recife_index, sero_window = built$sero_window,
                  config = c(settings, list(model_version = "pernambuco_v4_9_5strata_spatial",
                                             elapsed_seconds = elapsed_seconds)),
                  hmc = hmc)
  fit_path <- file.path(out_dir, "pe_5strata_pilot.rds")
  saveRDS(bundle, fit_path)
  message(sprintf("[5strata pilot] HMC gate: %s | %.0fs | %s", if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_pilot()
