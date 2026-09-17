# Fits v4.3 (2015-2019 q-identifiability calibration): Fit A (cases only),
# Fit B (cases + Juazeiro), or the legacy-prior sensitivity variant. Uses
# EXACTLY the v4.2/Stage-B4 sampler configuration (dense_e, adapt_delta=0.95,
# max_treedepth=14) -- Section 6 explicitly forbids returning to diag_e.
#
# Usage: FIT=A STAGE=canary Rscript 02_fit_v4_3.R
#        FIT=B STAGE=full   Rscript 02_fit_v4_3.R
#        FIT=legacy STAGE=canary Rscript 02_fit_v4_3.R

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

v4_3_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

root <- v4_3_root()
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))

base_dir <- file.path(root, "03_Output", "02_ceara_pipeline", "model_fits", "baseline", "v4_3_short_2015_2019_q_calibration")

stage_settings <- function(stage) {
  switch(stage,
    canary = list(warmup = 400L, iter_sampling = 200L, label = "canary (geometry check)"),
    full   = list(warmup = 1000L, iter_sampling = 1000L, label = "full (matches v4.0/td14 iteration count)"),
    stop("Unknown STAGE '", stage, "' -- must be 'canary' or 'full'.")
  )
}

make_v4_3_data <- function(weekly, imports_per_week, year_effect_prior_sd, t_sero, sero_pos, sero_n, kappa_sero, use_serology = NULL) {
  dates <- as.Date(weekly$week_start)
  years <- sort(unique(as.integer(format(dates, "%Y"))))
  expected_dates <- seq(min(dates), max(dates), by = "week")
  if (!identical(dates, expected_dates)) stop("Weekly data must be complete and consecutive.")
  accounting_error <- max(abs(
    weekly$N_end - weekly$N_start - weekly$births + weekly$all_cause_deaths - weekly$net_population_reconciliation
  ))
  if (accounting_error > 1e-6) stop("Demographic accounting failed before Stan fitting.")
  stan_data <- list(
    N = nrow(weekly), G = 8L, C = as.integer(weekly$cases), w = generation_weights(),
    N_start = as.numeric(weekly$N_start), N_end = as.numeric(weekly$N_end),
    births = as.numeric(weekly$births), deaths = as.numeric(weekly$all_cause_deaths),
    Y = length(years), year_id = match(as.integer(format(dates, "%Y")), years),
    seasonal_sin = sin(2 * pi * seq_len(nrow(weekly)) / 52.1775),
    seasonal_cos = cos(2 * pi * seq_len(nrow(weekly)) / 52.1775),
    imports_per_week = imports_per_week, year_effect_prior_sd = year_effect_prior_sd,
    t_sero = t_sero, sero_pos = sero_pos, sero_n = sero_n, kappa_sero = kappa_sero
  )
  if (!is.null(use_serology)) stan_data$use_serology <- use_serology
  list(stan_data = stan_data, years = years, demographic_accounting_error = accounting_error)
}

find_sero_week_index <- function(dates) {
  survey_start <- as.Date("2018-06-01"); survey_end <- as.Date("2018-12-31")
  midpoint <- survey_start + floor(as.numeric(survey_end - survey_start) / 2)
  idx <- which.min(abs(dates - midpoint))
  list(index = idx, midpoint_date = midpoint, matched_week_start = dates[idx])
}

# Section 2/3 dispersed inits: eta_q at 4 dispersed values spanning the weak
# prior's support (q ~ 0.02/0.05/0.10/0.20), transmission scalars reused
# from 4 distinct v4.0/td14 posterior draws (z_year truncated/re-drawn to
# the shorter Y=5 span with small deterministic offsets, since td14's
# Y=11-length z_year cannot be reused directly for a 5-year window).
eta_q_init_values <- qlogis(c(0.02, 0.05, 0.10, 0.20))
q_init_values_legacy <- qbeta(c(.20, .40, .60, .80), 16.2, 108.7)

load_td14_scalar_inits <- function(root, n_chains = 4L, seed = 20260912L) {
  td14_path <- file.path(root, "03_Output", "02_ceara_pipeline", "model_fits", "initialization", "renewal_ceara_v4_0_minimal_no_vaccine_fit_td14.rds")
  td14 <- readRDS(td14_path)
  draws <- rstan::extract(td14$fit, pars = c("alpha_R", "beta_sin", "beta_cos", "phi_obs"), permuted = TRUE)
  set.seed(seed)
  idx <- sample.int(length(draws$alpha_R), n_chains)
  lapply(idx, function(i) list(alpha_R = draws$alpha_R[i], beta_sin = draws$beta_sin[i], beta_cos = draws$beta_cos[i], phi_obs = draws$phi_obs[i]))
}

make_init_function <- function(Y, td14_scalar_inits, variant) {
  offsets <- c(-0.06, -0.02, 0.02, 0.06)
  function(chain_id = 1L) {
    cid <- as.integer(chain_id)
    base <- td14_scalar_inits[[(cid - 1L) %% length(td14_scalar_inits) + 1L]]
    offset <- offsets[(cid - 1L) %% 4L + 1L]
    z_year_init <- rep(offset / 0.40, Y)
    q_part <- if (identical(variant, "legacy")) list(q = q_init_values_legacy[(cid - 1L) %% 4L + 1L]) else list(eta_q = eta_q_init_values[(cid - 1L) %% 4L + 1L])
    c(base, list(z_year = z_year_init), q_part)
  }
}

run_fit_v4_3 <- function(fit_id = Sys.getenv("FIT", "A"), stage = Sys.getenv("STAGE", "canary")) {
  is_legacy <- identical(fit_id, "legacy")
  use_serology <- if (is_legacy) NULL else if (identical(fit_id, "A")) 0L else 1L
  settings <- stage_settings(stage)
  tag <- paste0("fit", fit_id, "_", stage)
  message("[", tag, "] ", if (is_legacy) "LEGACY PRIOR SENSITIVITY (cases+Juazeiro, Beta(16.2,108.7) prior)" else sprintf("Fit %s (use_serology=%d) -- %s", fit_id, use_serology, settings$label))

  out_dir <- file.path(base_dir, "outputs", tag)
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  fixed <- list(
    date_start = as.Date("2015-01-04"), date_end = as.Date("2019-12-29"),
    imports_per_week = 1, year_effect_prior_sd = 0.40,
    chains = 4L, max_treedepth = 14L, adapt_delta = .95, metric = "dense_e",
    sero_pos = 103L, sero_n = 404L, kappa_sero = 50
  )

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |> filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  demography_path <- file.path(root, "01_Data", "ceara_weekly_demography_2015_2025.rds")
  weekly_demography <- readRDS(demography_path) |> filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch.")
  message(sprintf("[%s] %s to %s | %d weeks", tag, fixed$date_start, fixed$date_end, nrow(weekly)))

  sero_week <- find_sero_week_index(as.Date(weekly$week_start))
  prepared <- make_v4_3_data(weekly, fixed$imports_per_week, fixed$year_effect_prior_sd,
                              sero_week$index, fixed$sero_pos, fixed$sero_n, fixed$kappa_sero,
                              use_serology = use_serology)
  iter_total <- settings$warmup + settings$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  td14_scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_init_function(prepared$stan_data$Y, td14_scalar_inits, if (is_legacy) "legacy" else "weak")

  stan_path <- file.path(root, "02_Script", "00_shared", "stan", "current", "ceara", if (is_legacy) "renewal_ceara_v4_3_legacy_prior.stan" else "renewal_ceara_v4_3_weak_q.stan")
  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data,
    chains = fixed$chains, iter = iter_total, warmup = settings$warmup,
    seed = 20260912L, init = init_fn,
    control = list(adapt_delta = fixed$adapt_delta, max_treedepth = fixed$max_treedepth, metric = fixed$metric),
    refresh = 25, sample_file = sample_file_base
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  q_par_name <- if (is_legacy) "q" else "eta_q"
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs", q_par_name), fixed$max_treedepth)

  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, sero_week = sero_week,
    config = c(fixed, settings, list(
      model_version = "v4_3_short_2015_2019", fit_id = fit_id, stage = stage, is_legacy = is_legacy,
      use_serology = use_serology, years = prepared$years, elapsed_seconds = elapsed_seconds,
      stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE)
    )),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_3_fit_", tag, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", tag, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_v4_3()
