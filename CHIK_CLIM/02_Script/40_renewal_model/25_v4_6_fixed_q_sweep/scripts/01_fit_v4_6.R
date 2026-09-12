# v4.6: q reverts to FIXED DATA (as in v4.0/td14), swept across a
# pre-registered grid, instead of jointly estimated (which v4.1-v4.5 showed
# is not resolvable by observation-model engineering alone). Same weekly
# renewal process and Juazeiro serology treatment as v4.3-v4.5, restricted
# to the same 2015-2019 window. Each q value is an independent fit; q's
# effect on case fit and serology compatibility is read off by comparing
# fits across the grid, not by estimating q inside one model.
#
# Q_GRID pre-registered BEFORE fitting: 0.05, 0.10, 0.15, 0.20, 0.25, 0.30
# (spans the low-q mode found throughout the audit, the legacy Beta(16.2,108.7)
# prior mean ~0.13, and a plausible upper range). q=0.10 is run first as a
# canary (matches td14's original fixed value) to confirm the reverted
# model reproduces clean HMC behavior before sweeping the rest.
#
# Usage: QVAL=0.10 Rscript 01_fit_v4_6.R

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

v4_6_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

root <- v4_6_root()
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "40_renewal_model", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script", "40_renewal_model", "22_v4_3_short_2015_2019_q_calibration", "scripts", "02_fit_v4_3.R"))

base_dir <- file.path(root, "02_Script", "40_renewal_model", "25_v4_6_fixed_q_sweep")

Q_GRID <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

make_v4_6_data <- function(weekly, imports_per_week, year_effect_prior_sd, t_sero, sero_pos, sero_n, kappa_sero, q) {
  prepared <- make_v4_3_data(weekly, imports_per_week, year_effect_prior_sd, t_sero, sero_pos, sero_n, kappa_sero)
  prepared$stan_data$q <- q
  prepared
}

# Fixed transmission-only inits (no eta_q/q dimension to disperse across
# chains -- q is data now, so a single dispersed-scalar init set suffices).
make_v4_6_init_fn <- function(Y, td14_scalar_inits) {
  offsets <- c(-0.06, -0.02, 0.02, 0.06)
  function(chain_id = 1L) {
    cid <- as.integer(chain_id)
    base <- td14_scalar_inits[[(cid - 1L) %% length(td14_scalar_inits) + 1L]]
    offset <- offsets[(cid - 1L) %% 4L + 1L]
    c(base, list(z_year = rep(offset / 0.40, Y)))
  }
}

run_fit_v4_6 <- function(qval = as.numeric(Sys.getenv("QVAL", "0.10"))) {
  tag <- paste0("q", sprintf("%.2f", qval))
  message("[", tag, "] v4.6 fixed-q-sweep | q=", qval)

  out_dir <- file.path(base_dir, "outputs", tag)
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  fixed <- list(
    date_start = as.Date("2015-01-04"), date_end = as.Date("2019-12-29"),
    imports_per_week = 1, year_effect_prior_sd = 0.40,
    chains = 4L, max_treedepth = 14L, adapt_delta = .95, metric = "dense_e",
    sero_pos = 103L, sero_n = 404L, kappa_sero = 50
  )
  settings <- list(warmup = 400L, iter_sampling = 200L)

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |> filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  demography_path <- file.path(root, "01_Data", "ceara_weekly_demography_2015_2025.rds")
  weekly_demography <- readRDS(demography_path) |> filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch.")

  sero_week <- find_sero_week_index(as.Date(weekly$week_start))
  prepared <- make_v4_6_data(weekly, fixed$imports_per_week, fixed$year_effect_prior_sd,
                              sero_week$index, fixed$sero_pos, fixed$sero_n, fixed$kappa_sero, qval)

  iter_total <- settings$warmup + settings$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  td14_scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_v4_6_init_fn(prepared$stan_data$Y, td14_scalar_inits)

  stan_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_6_fixed_q_sweep.stan")
  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data,
    chains = fixed$chains, iter = iter_total, warmup = settings$warmup,
    seed = 20260912L, init = init_fn,
    control = list(adapt_delta = fixed$adapt_delta, max_treedepth = fixed$max_treedepth, metric = fixed$metric),
    refresh = 25, sample_file = sample_file_base
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs"), fixed$max_treedepth)

  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, sero_week = sero_week,
    config = c(fixed, settings, list(model_version = "v4_6_fixed_q_sweep", q = qval, tag = tag,
                                      years = prepared$years, elapsed_seconds = elapsed_seconds,
                                      stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE))),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_6_fit_", tag, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", tag, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_v4_6()
