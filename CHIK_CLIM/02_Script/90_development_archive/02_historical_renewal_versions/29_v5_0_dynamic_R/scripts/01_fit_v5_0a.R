# v5.0a = v4.9 + a strongly regularized weekly AR(1) residual (delta_R) on
# log R0(t), keeping the direct Juazeiro serology anchor unchanged. See the
# .stan file's header comment for the full rationale (v4.9's smooth
# two-harmonic seasonal curve, even with year-specific amplitude, cannot
# sharpen the 2017 peak; the hypothesis is that a regularized short-
# timescale departure -- centred to zero within each year so it cannot
# absorb year_effect's or A_year's role -- can, without free weekly hazards
# or a return to v2's flexible R process).
#
# Reference q=0.05 only (per the current development-scenario decision --
# q is not treated as an estimated ascertainment fraction here). Everything
# except the log_R0(t) formula (mu_R + delta_R replacing v4.9's log_R0) is
# retained unchanged: fixed q, demographic S/U recursion, lifelong
# immunity, GI kernel, 2022 seed mechanism + likelihood exclusion, NB2
# observation model, phi_obs prior, harmonic/A_year structure, imports,
# Juazeiro serology treatment, HMC configuration.
#
# Usage: Rscript 01_fit_v5_0a.R

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

base_dir <- file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "29_v5_0_dynamic_R")

make_v5_0_init_fn <- function(Y, N, td14_scalar_inits) {
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
      rho_R = 0.8, sigma_R_dynamic = 0.02, z_delta = rep(0, N) # near-zero dynamic residual at init
    )
  }
}

run_fit_v5_0a <- function(qval = 0.05) {
  tag <- paste0("q", sprintf("%.2f", qval))
  message("[", tag, "] v5.0a dynamic-R (direct serology) | q=", qval)

  out_dir <- file.path(base_dir, "outputs", paste0("v5_0a_", tag))
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  fixed <- list(
    date_start = as.Date("2015-01-04"), date_end = as.Date("2025-12-21"),
    imports_per_week = 1, year_effect_prior_sd = 0.40,
    chains = 4L, max_treedepth = 14L, adapt_delta = .95, metric = "dense_e",
    sero_pos = 103L, sero_n = 404L, kappa_sero = 50
  )
  settings <- list(warmup = as.integer(Sys.getenv("WARMUP", "800")), iter_sampling = as.integer(Sys.getenv("ITER_SAMPLING", "400")))
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

  iter_total <- settings$warmup + settings$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  td14_scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_v5_0_init_fn(prepared$stan_data$Y, prepared$stan_data$N, td14_scalar_inits)

  stan_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v5_0a_dynamic_R_direct_serology.stan")
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
                                  "sigma_season_year", "z_season_year", "rho_R", "sigma_R_dynamic", "phi_obs"), fixed$max_treedepth)

  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, sero_week = sero_week,
    seed_idx = prepared$seed_idx, seed_dates = prepared$seed_dates,
    config = c(fixed, settings, list(model_version = "v5_0a_dynamic_R_direct_serology", q = qval, tag = tag,
                                      seed_start_date = as.character(SEED_START_DATE), seed_length = SEED_LENGTH,
                                      years = prepared$years, elapsed_seconds = elapsed_seconds,
                                      stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE))),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v5_0a_fit_", tag, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", tag, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_v5_0a()
