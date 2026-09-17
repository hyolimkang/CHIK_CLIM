# v5.0b = v4.9 + a low-rank (K=8 cyclic basis coefficients/year) dynamic-R
# deviation. See the .stan file's header for full rationale; supersedes
# v5.0a (weekly AR(1), archived as computationally infeasible).
#
# STAGE "pilot": short run (150 warmup / 100 sampling) whose ONLY purpose
# is to diagnose sampling geometry (seconds/iteration, divergences,
# treedepth hits, BFMI) BEFORE committing to a full production run, per
# explicit instruction ("Return the computational pilot before launching
# the full production fit").
# STAGE "full": production run (800 warmup / 400 sampling, matching v4.9's
# longer-run settings), only run after the pilot shows clean/practical geometry.
#
# Usage: STAGE=pilot Rscript 05_fit_v5_0b.R
#        STAGE=full  Rscript 05_fit_v5_0b.R

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

v5_0_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

root <- v5_0_root()
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "22_v4_3_short_2015_2019_q_calibration", "scripts", "02_fit_v4_3.R"))
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "27_v4_8_seeded_recurrence", "scripts", "01_fit_v4_8.R"))
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "28_v4_9_hierarchical_seasonality", "scripts", "01_fit_v4_9.R"))
source(file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "29_v5_0_dynamic_R", "scripts", "04_build_lowrank_basis_and_prior_predictive.R"))

base_dir <- file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "29_v5_0_dynamic_R")
K_KNOTS <- 8L

stage_settings <- function(stage) {
  switch(stage,
    pilot = list(warmup = 150L, iter_sampling = 100L, label = "PILOT (geometry check only)"),
    full  = list(warmup = 800L, iter_sampling = 400L, label = "FULL (matches v4.9 longer-run settings)"),
    stop("Unknown STAGE '", stage, "' -- must be 'pilot' or 'full'.")
  )
}

make_v5_0b_init_fn <- function(Y, K, td14_scalar_inits) {
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
      tau_dynamic = 0.05, z_dynamic = matrix(0, Y, K)
    )
  }
}

run_fit_v5_0b <- function(qval = 0.05, stage = Sys.getenv("STAGE", "pilot")) {
  settings <- stage_settings(stage)
  tag <- paste0("q", sprintf("%.2f", qval), "_", stage)
  message("[", tag, "] v5.0b low-rank dynamic-R | q=", qval, " | ", settings$label)

  out_dir <- file.path(base_dir, "outputs", paste0("v5_0b_", tag))
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  fixed <- list(
    date_start = as.Date("2015-01-04"), date_end = as.Date("2025-12-21"),
    imports_per_week = 1, year_effect_prior_sd = 0.40,
    chains = 4L, max_treedepth = 14L, adapt_delta = .95, metric = "dense_e",
    sero_pos = 103L, sero_n = 404L, kappa_sero = 50
  )
  SEED_START_DATE <- as.Date("2022-01-02")
  SEED_LENGTH <- 8L

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |> filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  demography_path <- file.path(root, "01_Data", "ceara_weekly_demography_2015_2025.rds")
  weekly_demography <- readRDS(demography_path) |> filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch.")

  sero_week <- find_sero_week_index(as.Date(weekly$week_start))
  prepared <- make_v4_9_data(weekly, fixed$imports_per_week, fixed$year_effect_prior_sd,
                              sero_week$index, fixed$sero_pos, fixed$sero_n, fixed$kappa_sero, qval,
                              SEED_START_DATE, SEED_LENGTH)

  basis <- build_dynamic_basis(as.Date(weekly$week_start), prepared$stan_data$year_id, K = K_KNOTS)
  prepared$stan_data$K <- K_KNOTS
  prepared$stan_data$B_dynamic <- basis$B_centred

  iter_total <- settings$warmup + settings$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  td14_scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_v5_0b_init_fn(prepared$stan_data$Y, K_KNOTS, td14_scalar_inits)

  stan_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v5_0b_lowrank_dynamic_R.stan")
  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data,
    chains = fixed$chains, iter = iter_total, warmup = settings$warmup,
    seed = 20260912L, init = init_fn,
    control = list(adapt_delta = fixed$adapt_delta, max_treedepth = fixed$max_treedepth, metric = fixed$metric),
    refresh = 10, sample_file = sample_file_base
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                                  "sigma_season_year", "z_season_year", "tau_dynamic", "phi_obs"), fixed$max_treedepth)

  seconds_per_iteration <- elapsed_seconds / (iter_total * fixed$chains)
  message(sprintf("[%s] seconds/iteration (wall, across %d chains): %.3f", tag, fixed$chains, seconds_per_iteration))

  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, sero_week = sero_week,
    seed_idx = prepared$seed_idx, seed_dates = prepared$seed_dates, basis = basis,
    config = c(fixed, settings, list(model_version = "v5_0b_lowrank_dynamic_R", q = qval, tag = tag, stage = stage,
                                      K = K_KNOTS, seed_start_date = as.character(SEED_START_DATE), seed_length = SEED_LENGTH,
                                      years = prepared$years, elapsed_seconds = elapsed_seconds, seconds_per_iteration = seconds_per_iteration,
                                      stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE))),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v5_0b_fit_", tag, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", tag, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_v5_0b()
