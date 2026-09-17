# v4.8 = v4.7 + externally conditioned 2022 epidemic re-seeding.
#
# v4.7 (full 2015-2025 fixed-q sweep) systematically underpredicts the 2022
# recurrence (observed 50,986 vs. predicted medians ~23,900-26,400 across
# q=0.05-0.30) because the renewal kernel's convolution term has decayed to
# near zero after ~4 years of very low incidence and cannot reignite growth
# without an external seed. There is independent genomic evidence of a
# newly introduced ECSA lineage in 2022, so v4.8 CONDITIONS ON a 2022
# reintroduction (does not try to infer whether one occurred) and asks:
# given that seed, can the EXISTING (never-reset) susceptibility and the
# EXISTING renewal/transmission model reproduce the subsequent wave?
#
# Retained unchanged from v4.7: fixed q, weekly NB observation likelihood
# (outside the seed window), GI kernel, R0(t) parameterisation, all
# priors, S/U demographic bookkeeping, lifelong immunity, Juazeiro
# serology (beta-binomial, kappa_sero=50), HMC configuration.
#
# 2022 SEED WINDOW (pre-registered here in R, not hard-coded in Stan):
#   SEED_START_DATE = 2022-01-02, SEED_LENGTH = 8 weeks (= G, the
#   generation-interval kernel length -- the minimum window needed to
#   prime the renewal convolution before handing control back to it).
#   Chosen from the raw weekly series: sustained exponential growth is
#   already visible within this window (44 -> 75 -> 112 -> 153 -> 309 ->
#   388 -> 561 -> 647 cases/week, 2022-01-02 to 2022-02-20) and matches
#   the calendar-year boundary used everywhere else in this analysis for
#   "the 2022 recurrence". The remaining 44 weeks of 2022 (2022-02-27
#   onward, including the entire peak of ~3,500/week in May and the full
#   decline) are left ENTIRELY to the renewal process and the fitted
#   likelihood -- this is the actual test.
#
# Seed weeks are conditioned via X_seed[t] = C[t]/q (the SAME q-mapping
# v4.7 already assumes) and excluded from the case likelihood via
# fit_index, so the same 8 observations are never used both to define the
# latent state and to score it (see the .stan file's header comment).
#
# Usage: QVAL=0.10 Rscript 01_fit_v4_8.R

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

v4_8_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

root <- v4_8_root()
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "22_v4_3_short_2015_2019_q_calibration", "scripts", "02_fit_v4_3.R"))

base_dir <- file.path(root, "03_Output", "02_ceara_pipeline", "model_fits", "baseline", "v4_8_seeded_recurrence")

Q_GRID <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30) # identical grid to v4.7
SEED_START_DATE <- as.Date("2022-01-02")
SEED_LENGTH <- 8L # = G (generation-interval kernel length)

make_v4_8_data <- function(weekly, imports_per_week, year_effect_prior_sd, t_sero, sero_pos, sero_n, kappa_sero, q,
                            seed_start_date, seed_length) {
  prepared <- make_v4_3_data(weekly, imports_per_week, year_effect_prior_sd, t_sero, sero_pos, sero_n, kappa_sero)
  prepared$stan_data$q <- q

  dates <- as.Date(weekly$week_start)
  seed_start_idx <- which(dates == seed_start_date)
  if (length(seed_start_idx) != 1L) stop("SEED_START_DATE not found (or not unique) in the weekly series.")
  seed_idx <- seq.int(seed_start_idx, seed_start_idx + seed_length - 1L)
  if (max(seed_idx) > nrow(weekly)) stop("Seed window extends beyond the available weekly series.")

  is_seed <- integer(nrow(weekly))
  is_seed[seed_idx] <- 1L
  X_seed <- numeric(nrow(weekly))
  X_seed[seed_idx] <- weekly$cases[seed_idx] / q # same q*X mapping as the rest of the observation model
  fit_index <- which(is_seed == 0L)

  prepared$stan_data$is_seed <- is_seed
  prepared$stan_data$X_seed <- X_seed
  prepared$stan_data$N_fit <- length(fit_index)
  prepared$stan_data$fit_index <- fit_index
  prepared$seed_idx <- seed_idx
  prepared$seed_dates <- dates[seed_idx]
  prepared
}

make_v4_8_init_fn <- function(Y, td14_scalar_inits) {
  offsets <- c(-0.06, -0.02, 0.02, 0.06)
  function(chain_id = 1L) {
    cid <- as.integer(chain_id)
    base <- td14_scalar_inits[[(cid - 1L) %% length(td14_scalar_inits) + 1L]]
    offset <- offsets[(cid - 1L) %% 4L + 1L]
    c(base, list(z_year = rep(offset / 0.40, Y)))
  }
}

run_fit_v4_8 <- function(qval = as.numeric(Sys.getenv("QVAL", "0.10"))) {
  tag <- paste0("q", sprintf("%.2f", qval))
  message("[", tag, "] v4.8 seeded-recurrence FULL PERIOD (2015-2025) | q=", qval)

  out_dir <- file.path(base_dir, "outputs", tag)
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  fixed <- list(
    date_start = as.Date("2015-01-04"), date_end = as.Date("2025-12-21"), # identical to v4.7
    imports_per_week = 1, year_effect_prior_sd = 0.40,
    chains = 4L, max_treedepth = 14L, adapt_delta = .95, metric = "dense_e", # identical HMC config to v4.7
    sero_pos = 103L, sero_n = 404L, kappa_sero = 50
  )
  settings <- list(warmup = as.integer(Sys.getenv("WARMUP", "400")), iter_sampling = as.integer(Sys.getenv("ITER_SAMPLING", "200")))

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |> filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  demography_path <- file.path(root, "01_Data", "ceara_weekly_demography_2015_2025.rds")
  weekly_demography <- readRDS(demography_path) |> filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch.")

  sero_week <- find_sero_week_index(as.Date(weekly$week_start))
  prepared <- make_v4_8_data(weekly, fixed$imports_per_week, fixed$year_effect_prior_sd,
                              sero_week$index, fixed$sero_pos, fixed$sero_n, fixed$kappa_sero, qval,
                              SEED_START_DATE, SEED_LENGTH)
  message(sprintf("[%s] seed window: %s to %s (%d weeks) | N_fit=%d of N=%d",
                   tag, min(prepared$seed_dates), max(prepared$seed_dates), SEED_LENGTH,
                   prepared$stan_data$N_fit, prepared$stan_data$N))

  iter_total <- settings$warmup + settings$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  td14_scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_v4_8_init_fn(prepared$stan_data$Y, td14_scalar_inits)

  stan_path <- file.path(root, "02_Script", "00_shared", "stan", "current", "ceara", "renewal_ceara_v4_8_seeded_recurrence.stan")
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
    seed_idx = prepared$seed_idx, seed_dates = prepared$seed_dates,
    config = c(fixed, settings, list(model_version = "v4_8_seeded_recurrence", q = qval, tag = tag,
                                      seed_start_date = as.character(SEED_START_DATE), seed_length = SEED_LENGTH,
                                      years = prepared$years, elapsed_seconds = elapsed_seconds,
                                      stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE))),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_8_fit_", tag, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", tag, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_v4_8()
