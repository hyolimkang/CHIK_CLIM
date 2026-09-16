# Bahia identification experiment: frozen v4.9 transmission structure + the
# SAME six-site geographic-offset serology model, but with q replaced by ONE
# global, time-constant, logit-scale ESTIMATED parameter
# (renewal_bahia_v4_9_global_q_multisite_serology.stan). Prior:
# logit_q ~ normal(logit(0.10), 1.0), deliberately broad, NOT tightened
# after seeing the posterior. Also supports the J_sero=0 identification
# ablation (SERO_ON=FALSE) using the exact same Stan file/prior.
#
# Usage:
#   Rscript 11_fit_bahia_global_q_multisite.R              # pilot, serology ON
#   SERO_ON=FALSE Rscript 11_fit_bahia_global_q_multisite.R # pilot, serology OFF (ablation)
#   STAGE=full Rscript 11_fit_bahia_global_q_multisite.R

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

bahia_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

root <- bahia_root()
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "40_renewal_model", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script", "40_renewal_model", "22_v4_3_short_2015_2019_q_calibration", "scripts", "02_fit_v4_3.R"))
source(file.path(root, "02_Script", "40_renewal_model", "30_bahia_v4_9_replication", "scripts", "06_build_bahia_serology_windows.R"))

base_dir <- file.path(root, "02_Script", "40_renewal_model", "30_bahia_v4_9_replication")
SERO_GEOGRAPHIC_SD <- 1.0 # FIXED, not estimated
LOGIT_Q_PRIOR_MEAN <- qlogis(0.10) # logit(0.10)
LOGIT_Q_PRIOR_SD <- 1.0 # deliberately broad regularising prior; NOT tightened post-hoc

make_bahia_global_q_data <- function(weekly, imports_per_week, year_effect_prior_sd, sero_on = TRUE) {
  prepared <- make_v4_3_data(weekly, imports_per_week, year_effect_prior_sd,
                              t_sero = 1L, sero_pos = 0L, sero_n = 1L, kappa_sero = 50)
  n <- nrow(weekly)
  prepared$stan_data$seasonal_sin2 <- sin(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$seasonal_cos2 <- cos(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$t_sero <- NULL; prepared$stan_data$sero_pos <- NULL; prepared$stan_data$sero_n <- NULL
  prepared$stan_data$q_fixed <- NULL # this model estimates q; no fixed-q field in the data block
  prepared$stan_data$is_seed <- integer(n)
  prepared$stan_data$X_seed <- numeric(n)
  prepared$stan_data$N_fit <- n
  prepared$stan_data$fit_index <- seq_len(n)
  prepared$stan_data$logit_q_prior_mean <- LOGIT_Q_PRIOR_MEAN
  prepared$stan_data$logit_q_prior_sd <- LOGIT_Q_PRIOR_SD

  if (sero_on) {
    sero_built <- build_bahia_sero_stan_fields(as.Date(weekly$week_start))
    prepared$stan_data <- c(prepared$stan_data, sero_built$stan_fields)
    prepared$sero_audit <- sero_built$audit
  } else {
    # J_sero=0 identification ablation: same file, empty serology arrays.
    prepared$stan_data$J_sero <- 0L
    prepared$stan_data$sero_n_positive <- integer(0)
    prepared$stan_data$sero_n_tested <- integer(0)
    prepared$stan_data$sero_window_start_idx <- integer(0)
    prepared$stan_data$sero_window_n_weeks <- integer(0)
    prepared$sero_audit <- NULL
  }
  prepared$stan_data$sero_geographic_sd <- SERO_GEOGRAPHIC_SD
  prepared
}

make_bahia_global_q_init_fn <- function(Y, J_sero, td14_scalar_inits) {
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
      z_geo = if (J_sero > 0) rep(0, J_sero) else numeric(0),
      logit_q = LOGIT_Q_PRIOR_MEAN + offset # start near prior centre, small between-chain spread
    )
  }
}

run_fit_bahia_global_q <- function(qval_tag = Sys.getenv("QTAG", "on"),
                                    stage = Sys.getenv("STAGE", "pilot"),
                                    sero_on = as.logical(Sys.getenv("SERO_ON", "TRUE"))) {
  # Short tags only ("on"/"off") -- this OneDrive path is already close to
  # Windows' 260-char MAX_PATH; a longer tag (e.g. "global_q_sero_ablation")
  # pushed saveRDS()'s target path to 274 chars and failed with
  # "cannot open the connection" only at the very last (save) step, AFTER
  # all 4 chains had already sampled successfully.
  tag <- if (sero_on) qval_tag else "off"
  message("[bahia-", tag, "] frozen v4.9 structure + global-q", if (sero_on) " + 6-site serology" else " (SEROLOGY ABLATION, J_sero=0)",
          " | stage=", stage)

  out_dir <- file.path(base_dir, "outputs_global_q_multisite_serology", tag)
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

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "29") |> filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  demography_path <- file.path(root, "01_Data", "bahia_weekly_demography_2015_2025.rds")
  weekly_demography <- readRDS(demography_path) |> filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch.")

  prepared <- make_bahia_global_q_data(weekly, fixed$imports_per_week, fixed$year_effect_prior_sd, sero_on = sero_on)

  iter_total <- settings$warmup + settings$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  td14_scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_bahia_global_q_init_fn(prepared$stan_data$Y, prepared$stan_data$J_sero, td14_scalar_inits)

  stan_path <- file.path(root, "02_Script", "stan", "renewal_bahia_v4_9_global_q_multisite_serology.stan")
  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data,
    chains = fixed$chains, iter = iter_total, warmup = settings$warmup,
    seed = 20260913L, init = init_fn,
    control = list(adapt_delta = fixed$adapt_delta, max_treedepth = fixed$max_treedepth, metric = fixed$metric),
    refresh = 25, sample_file = sample_file_base
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  gate_pars <- c("alpha_R", "z_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                 "sigma_season_year", "z_season_year", "phi_obs", "logit_q")
  if (sero_on) gate_pars <- c(gate_pars, "z_geo")
  hmc <- compute_hmc_gate(fit, gate_pars, fixed$max_treedepth)

  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, sero_audit = prepared$sero_audit,
    config = c(fixed, settings, list(model_version = "bahia_v4_9_global_q_multisite_serology", tag = tag, stage = stage,
                                      sero_on = sero_on, sero_geographic_sd = SERO_GEOGRAPHIC_SD,
                                      logit_q_prior_mean = LOGIT_Q_PRIOR_MEAN, logit_q_prior_sd = LOGIT_Q_PRIOR_SD,
                                      years = prepared$years, elapsed_seconds = elapsed_seconds,
                                      stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE))),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_bahia_global_q_fit_", tag, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[bahia-%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", tag, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_bahia_global_q()
