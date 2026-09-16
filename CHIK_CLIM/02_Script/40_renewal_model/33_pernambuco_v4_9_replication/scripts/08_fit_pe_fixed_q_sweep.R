# Pernambuco external replication -- Section 13 fixed-q sensitivity sweep.
# Same Bahia-style q grid (0.05-0.30), U14 Recife serology ON, SAME
# geo-adjusted observation model (sero_geographic_sd=1.0 fixed, same
# kappa_sero). Reuses the EXACT same generic fixed-q Stan file already
# validated for Bahia (renewal_ceara_v4_9_bahia_multisite_serology.stan)
# -- state-agnostic, not modified. This is a mandatory SENSITIVITY
# analysis, not the primary PE fit (Models A/B remain primary).
#
# Usage: QVAL=0.05 Rscript 08_fit_pe_fixed_q_sweep.R

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
setwd(root)
source(file.path(root, "02_Script/40_renewal_model/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script/40_renewal_model/22_v4_3_short_2015_2019_q_calibration/scripts/02_fit_v4_3.R"))
source(file.path(root, "02_Script/40_renewal_model/33_pernambuco_v4_9_replication/scripts/02_build_pe_serology_window.R"))

base_dir <- file.path(root, "02_Script/40_renewal_model/33_pernambuco_v4_9_replication")
SERO_GEOGRAPHIC_SD <- 1.0 # SAME fixed value as Bahia -- not tuned for PE

make_pe_fixedq_data <- function(weekly, imports_per_week, year_effect_prior_sd, q) {
  prepared <- make_v4_3_data(weekly, imports_per_week, year_effect_prior_sd,
                              t_sero = 1L, sero_pos = 0L, sero_n = 1L, kappa_sero = 50)
  n <- nrow(weekly)
  prepared$stan_data$seasonal_sin2 <- sin(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$seasonal_cos2 <- cos(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$q <- q
  prepared$stan_data$t_sero <- NULL; prepared$stan_data$sero_pos <- NULL; prepared$stan_data$sero_n <- NULL
  prepared$stan_data$is_seed <- integer(n)
  prepared$stan_data$X_seed <- numeric(n)
  prepared$stan_data$N_fit <- n
  prepared$stan_data$fit_index <- seq_len(n)

  sero_built <- build_pe_sero_stan_fields(as.Date(weekly$week_start))
  prepared$stan_data <- c(prepared$stan_data, sero_built$stan_fields)
  prepared$stan_data$sero_geographic_sd <- SERO_GEOGRAPHIC_SD
  prepared$sero_audit <- sero_built$audit
  prepared
}

make_pe_fixedq_init_fn <- function(Y, J_sero, td14_scalar_inits) {
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
      z_geo = as.array(rep(0, J_sero))
    )
  }
}

run_fit_pe_fixedq <- function(qval = as.numeric(Sys.getenv("QVAL", "0.05")), stage = Sys.getenv("STAGE", "full")) {
  tag <- paste0("q", sprintf("%.2f", qval))
  message("[pe-fixedq-", tag, "] frozen v4.9 structure + fixed q=", qval, " + U14 Recife serology | stage=", stage)

  out_dir <- file.path(base_dir, "outputs", "fixed_q_sweep", tag)
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  fixed <- list(
    date_start = as.Date("2015-01-04"), date_end = as.Date("2025-12-21"),
    imports_per_week = 1, year_effect_prior_sd = 0.40,
    chains = 4L, max_treedepth = 14L, adapt_delta = .95, metric = "dense_e"
  )
  default_warmup <- if (stage == "pilot") "150" else "800"
  default_sampling <- if (stage == "pilot") "100" else "400"
  settings <- list(warmup = as.integer(Sys.getenv("WARMUP", default_warmup)),
                    iter_sampling = as.integer(Sys.getenv("ITER_SAMPLING", default_sampling)))

  weekly <- readRDS(file.path(root, "03_Output/tables/pernambuco_v4_9_replication/pernambuco_weekly_input.rds"))
  prepared <- make_pe_fixedq_data(weekly, fixed$imports_per_week, fixed$year_effect_prior_sd, qval)

  iter_total <- settings$warmup + settings$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  td14_scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_pe_fixedq_init_fn(prepared$stan_data$Y, prepared$stan_data$J_sero, td14_scalar_inits)

  stan_path <- file.path(root, "02_Script/stan/renewal_ceara_v4_9_bahia_multisite_serology.stan")
  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data,
    chains = fixed$chains, iter = iter_total, warmup = settings$warmup,
    seed = 20260913L, init = init_fn,
    control = list(adapt_delta = fixed$adapt_delta, max_treedepth = fixed$max_treedepth, metric = fixed$metric),
    refresh = 25, sample_file = sample_file_base
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                                  "sigma_season_year", "z_season_year", "phi_obs", "z_geo"), fixed$max_treedepth)

  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, sero_audit = prepared$sero_audit,
    config = c(fixed, settings, list(model_version = "pernambuco_fixedq_U14_serology", q = qval, tag = tag, stage = stage,
                                      sero_geographic_sd = SERO_GEOGRAPHIC_SD,
                                      years = prepared$years, elapsed_seconds = elapsed_seconds,
                                      stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE))),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_pe_fixedq_fit_", tag, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[pe-fixedq-%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", tag, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_pe_fixedq()
