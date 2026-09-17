# v4.9 = v4.8 + a parsimonious hierarchical seasonal-shape model, replacing
# ONLY the temporal R0(t) parameterisation (see the .stan file's header
# comment for the full rationale). Everything else -- fixed q sweep, q
# grid, demographic S/U recursion, lifelong immunity, weekly NB
# observation likelihood, phi_obs prior, Juazeiro serology treatment, the
# 2022 conditional seed mechanism and its exclusion from the likelihood,
# imports_per_week, generation interval, HMC configuration -- is retained
# unchanged from v4.8 by construction (make_v4_9_data() calls v4.8's
# make_v4_8_data() unmodified and only appends the second-harmonic
# covariates).
#
# Usage: QVAL=0.10 Rscript 01_fit_v4_9.R

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

v4_9_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

root <- v4_9_root()
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "22_v4_3_short_2015_2019_q_calibration", "scripts", "02_fit_v4_3.R"))
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "27_v4_8_seeded_recurrence", "scripts", "01_fit_v4_8.R"))

base_dir <- file.path(root, "03_Output", "model_fits", "ceara", "v4_9_hierarchical_seasonality")

Q_GRID <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30) # identical grid to v4.7/v4.8

# v4.9's data prep is v4.8's make_v4_8_data() (is_seed/X_seed/fit_index
# construction is byte-for-byte unchanged) plus the second annual harmonic.
make_v4_9_data <- function(weekly, imports_per_week, year_effect_prior_sd, t_sero, sero_pos, sero_n, kappa_sero, q,
                            seed_start_date, seed_length) {
  prepared <- make_v4_8_data(weekly, imports_per_week, year_effect_prior_sd, t_sero, sero_pos, sero_n, kappa_sero, q,
                              seed_start_date, seed_length)
  n <- nrow(weekly)
  prepared$stan_data$seasonal_sin2 <- sin(4 * pi * seq_len(n) / 52.1775) # second harmonic: double the first-harmonic frequency
  prepared$stan_data$seasonal_cos2 <- cos(4 * pi * seq_len(n) / 52.1775)
  prepared
}

make_v4_9_init_fn <- function(Y, td14_scalar_inits) {
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
      z_year = rep(offset / 0.40, Y)
    )
  }
}

run_fit_v4_9 <- function(qval = as.numeric(Sys.getenv("QVAL", "0.10"))) {
  tag <- paste0("q", sprintf("%.2f", qval))
  message("[", tag, "] v4.9 hierarchical-seasonality FULL PERIOD (2015-2025) | q=", qval)

  out_dir <- file.path(base_dir, "outputs", tag)
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  fixed <- list(
    date_start = as.Date("2015-01-04"), date_end = as.Date("2025-12-21"), # identical to v4.7/v4.8
    imports_per_week = 1, year_effect_prior_sd = 0.40,
    chains = 4L, max_treedepth = 14L, adapt_delta = .95, metric = "dense_e", # identical HMC config
    sero_pos = 103L, sero_n = 404L, kappa_sero = 50
  )
  settings <- list(warmup = as.integer(Sys.getenv("WARMUP", "400")), iter_sampling = as.integer(Sys.getenv("ITER_SAMPLING", "200")))
  SEED_START_DATE <- as.Date("2022-01-02") # identical seed window to v4.8
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
  message(sprintf("[%s] seed window: %s to %s (%d weeks) | N_fit=%d of N=%d",
                   tag, min(prepared$seed_dates), max(prepared$seed_dates), SEED_LENGTH,
                   prepared$stan_data$N_fit, prepared$stan_data$N))

  iter_total <- settings$warmup + settings$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  td14_scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_v4_9_init_fn(prepared$stan_data$Y, td14_scalar_inits)

  stan_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_9_hierarchical_seasonality.stan")
  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data,
    chains = fixed$chains, iter = iter_total, warmup = settings$warmup,
    seed = 20260912L, init = init_fn,
    control = list(adapt_delta = fixed$adapt_delta, max_treedepth = fixed$max_treedepth, metric = fixed$metric),
    refresh = 25, sample_file = sample_file_base
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                                  "sigma_season_year", "z_season_year", "phi_obs"), fixed$max_treedepth)

  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, sero_week = sero_week,
    seed_idx = prepared$seed_idx, seed_dates = prepared$seed_dates,
    config = c(fixed, settings, list(model_version = "v4_9_hierarchical_seasonality", q = qval, tag = tag,
                                      seed_start_date = as.character(SEED_START_DATE), seed_length = SEED_LENGTH,
                                      years = prepared$years, elapsed_seconds = elapsed_seconds,
                                      stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE))),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_9_fit_", tag, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", tag, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_v4_9()
