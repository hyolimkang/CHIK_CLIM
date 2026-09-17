# Bahia external replication of the FROZEN Ceara v4.9 hierarchical-
# seasonality renewal model. See BAHIA_V4_9_FROZEN_MODEL_ASSESSMENT.md.
#
# Uses renewal_ceara_v4_9_hierarchical_seasonality_optional_serology.stan
# (the ONLY change from frozen v4.9: a use_serology gate) with
# use_serology=0 for Bahia (no serology anchor exists; not invented/
# borrowed from Ceara). Reuses make_v4_3_data()/generation_weights()/
# compute_hmc_gate()/load_td14_scalar_inits() UNCHANGED from the Ceara
# pipeline; only the input data (Bahia case/demography series, UF="29")
# and the serology/seed fields differ.
#
# NO SEED MECHANISM: per BAHIA_V4_9_FROZEN_MODEL_ASSESSMENT.md Section 6,
# applying the SAME objective long-gap standard implied by Ceara's own
# 2022 event (weeks_since_previous_wave=399, ~7.7 years -- the only Ceara
# transition ever requiring re-seeding; all other Ceara transitions, up to
# 56 weeks, were handled by the continuous renewal process alone) to
# Bahia's 8 inter-wave gaps (14-49 weeks, from
# 03_Output/07_national_pipeline/tables/national_wave_analysis/brazil_chik_wave_analysis_master.csv)
# finds NONE qualify for external re-seeding. is_seed is therefore all-zero
# and every week is included in the likelihood (fit_index = 1:N).
#
# Usage: QVAL=0.05 Rscript 02_fit_bahia_v4_9.R

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

root <- ROOT # from 00_project_setup.R -- bahia_root() left defined above but unused, no scientific change
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "22_v4_3_short_2015_2019_q_calibration", "scripts", "02_fit_v4_3.R"))

base_dir <- file.path(root, "03_Output", "03_bahia_pipeline", "model_fits", "v4_9_replication")

Q_GRID <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30) # identical grid to the final Ceara sensitivity analysis; q=0.05 not privileged

make_bahia_v4_9_data <- function(weekly, imports_per_week, year_effect_prior_sd, q) {
  # Dummy serology fields (t_sero=1, sero_pos=0, sero_n=1): never enter the
  # likelihood because use_serology=0, but make_v4_3_data()/the Stan data
  # block still require a valid-shaped placeholder.
  prepared <- make_v4_3_data(weekly, imports_per_week, year_effect_prior_sd,
                              t_sero = 1L, sero_pos = 0L, sero_n = 1L, kappa_sero = 50)
  n <- nrow(weekly)
  prepared$stan_data$seasonal_sin2 <- sin(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$seasonal_cos2 <- cos(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$q <- q
  prepared$stan_data$use_serology <- 0L
  # No seed mechanism (see header comment): every week is fit.
  prepared$stan_data$is_seed <- integer(n)
  prepared$stan_data$X_seed <- numeric(n)
  prepared$stan_data$N_fit <- n
  prepared$stan_data$fit_index <- seq_len(n)
  prepared
}

make_bahia_init_fn <- function(Y, td14_scalar_inits) {
  offsets <- c(-0.06, -0.02, 0.02, 0.06)
  function(chain_id = 1L) {
    cid <- as.integer(chain_id)
    base <- td14_scalar_inits[[(cid - 1L) %% length(td14_scalar_inits) + 1L]] # Ceara td14 posterior draws used ONLY as a computational starting point, not scientific borrowing
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

run_fit_bahia_v4_9 <- function(qval = as.numeric(Sys.getenv("QVAL", "0.05")), stage = Sys.getenv("STAGE", "pilot")) {
  tag <- paste0("q", sprintf("%.2f", qval))
  message("[bahia-", tag, "] frozen v4.9 structure, use_serology=0, no seed | q=", qval, " | stage=", stage)

  out_dir <- file.path(base_dir, "outputs", tag)
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  fixed <- list(
    date_start = as.Date("2015-01-04"), date_end = as.Date("2025-12-21"), # identical window to Ceara v4.9
    imports_per_week = 1, year_effect_prior_sd = 0.40,
    chains = 4L, max_treedepth = 14L, adapt_delta = .95, metric = "dense_e" # identical HMC config to Ceara v4.9
  )
  default_warmup <- if (stage == "pilot") "150" else "800"
  default_sampling <- if (stage == "pilot") "100" else "400"
  settings <- list(warmup = as.integer(Sys.getenv("WARMUP", default_warmup)),
                    iter_sampling = as.integer(Sys.getenv("ITER_SAMPLING", default_sampling)))

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "29") |> dplyr::filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  demography_path <- file.path(root, "01_Data", "bahia_weekly_demography_2015_2025.rds")
  weekly_demography <- readRDS(demography_path) |> dplyr::filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch.")

  prepared <- make_bahia_v4_9_data(weekly, fixed$imports_per_week, fixed$year_effect_prior_sd, qval)

  iter_total <- settings$warmup + settings$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  td14_scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_bahia_init_fn(prepared$stan_data$Y, td14_scalar_inits)

  stan_path <- file.path(root, "02_Script", "00_shared", "stan", "current", "multistate", "renewal_ceara_v4_9_hierarchical_seasonality_optional_serology.stan")
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
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly,
    config = c(fixed, settings, list(model_version = "bahia_v4_9_frozen_replication", q = qval, tag = tag, stage = stage,
                                      years = prepared$years, elapsed_seconds = elapsed_seconds,
                                      stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE))),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_bahia_v4_9_fit_", tag, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[bahia-%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", tag, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_bahia_v4_9()
