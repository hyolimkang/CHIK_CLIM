# Rio de Janeiro CITY auxiliary -- Model B rescue: same data, same priors,
# same inits, SAME iteration count (1200 warmup / 1000 sampling) as
# 12_fit_rj_city_u10.R. The only change is adapt_delta 0.98 -> 0.99, to
# resolve the single stray divergent transition observed in the original
# run (otherwise clean: max Rhat=1.0011, min ESS=1661, BFMI>=0.84) -- the
# same benign near-pass pattern and rescue approach already used and
# approved for the RJ STATE global-q case-only model (04b_..._ad098.R).

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script/40_renewal_model/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script/40_renewal_model/22_v4_3_short_2015_2019_q_calibration/scripts/02_fit_v4_3.R"))
source(file.path(root, "02_Script/40_renewal_model/34_rio_de_janeiro_v4_9_replication/scripts/12_fit_rj_city_u10.R"))

base_dir <- file.path(root, "02_Script/40_renewal_model/34_rio_de_janeiro_v4_9_replication")
out_dir <- file.path(base_dir, "outputs", "city_u10_ad099")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(root, "03_Output/tables/rio_de_janeiro_v4_9_replication/city")

run_fit_rj_city_u10_ad099 <- function() {
  settings <- list(warmup = 1200L, sampling = 1000L, chains = 4L, adapt_delta = 0.99, max_treedepth = 14L, metric = "dense_e")
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  weekly <- readRDS(file.path(table_dir, "riodejaneiro_city_weekly_input.rds"))
  prepared <- make_rj_city_u10_data(weekly)
  Y <- prepared$stan_data$Y
  message("RJ CITY + U10 (ad099 rescue) Stan data: N=", prepared$stan_data$N, " Y=", Y, " J_sero=", prepared$stan_data$J_sero,
          " sero_window=[", prepared$stan_data$sero_window_start_idx, ",",
          prepared$stan_data$sero_window_start_idx + prepared$stan_data$sero_window_n_weeks - 1, "]")

  scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_dispersed_init_fn(Y, scalar_inits)

  stan_path <- file.path(root, "02_Script/stan/renewal_riodejaneiro_city_v4_9_global_q_logitnormal_serology.stan")
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
                  config = c(settings, list(model_version = "rio_de_janeiro_CITY_v4_9_global_q_u10_logitnormal_ad099",
                                             elapsed_seconds = elapsed_seconds, seed = 20260915L)),
                  hmc = hmc)
  fit_path <- file.path(out_dir, "rj_city_u10_ad099.rds")
  saveRDS(bundle, fit_path)
  message(sprintf("[rj city + U10 ad099] HMC gate: %s | %.0fs | %s", if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_rj_city_u10_ad099()
