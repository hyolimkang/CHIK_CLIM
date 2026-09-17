# Staged fitting for the v4.1 q-only ablation. Reuses v4.0's exact helper
# functions (generation_weights, compute_hmc_gate, build_state_weekly) by
# sourcing the archived v4.0 fit script's logic path (the SAME source file
# referenced in archive/CHECKSUM_MANIFEST.md), not a reimplementation.
#
# STAGE controls warmup/sampling length only -- everything else (data,
# priors, initial-value strategy, adapt_delta, max_treedepth) is identical
# across stages, matching the archived v4.0/td14 numerical settings except
# for the deliberately shortened iteration counts in Stages A-C:
#   A: 5 warmup + 5 sampling   (implementation validation only, not for gate)
#   B: 300 warmup + 150 sampling (geometry canary)
#   C: 500 warmup + 300 sampling (intermediate)
#   D: 1000 warmup + 1000 sampling (production, matches v4.0/td14 exactly)
#
# Per-chain CmdStan-style CSVs are written to a PERMANENT path from the
# start of sampling (outputs/stage<X>/chains/chain_<id>.csv) so
# 03_monitor_v4_1_q_only.R can inspect live progress in another terminal.
#
# Usage: STAGE=B Rscript 02_fit_v4_1_q_only.R

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

q_only_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

root <- q_only_root()
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
# Reuse v4.0's exact generation_weights()/compute_hmc_gate()/
# build_state_weekly() -- sourcing the archived-checksum-verified fit
# script's own definitions, not reimplementing them.
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))

q_only_dir <- function(root = q_only_root()) file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "20_v4_1_q_only")

stage_settings <- function(stage) {
  switch(stage,
    A = list(warmup = 5L,    iter_sampling = 5L,    label = "Stage A: implementation validation only -- NOT for gate assessment"),
    B = list(warmup = 300L,  iter_sampling = 150L,  label = "Stage B: geometry canary"),
    C = list(warmup = 500L,  iter_sampling = 300L,  label = "Stage C: intermediate run"),
    D = list(warmup = 1000L, iter_sampling = 1000L, label = "Stage D: production run (matches v4.0/td14 exactly)"),
    stop("Unknown STAGE '", stage, "' -- must be one of A, B, C, D.")
  )
}

# Mirrors v4.0's make_v4_data() exactly, except q_fixed is NOT included in
# stan_data (q is now a sampled parameter, not data).
make_v4_1_q_only_data <- function(weekly, imports_per_week, year_effect_prior_sd) {
  dates <- as.Date(weekly$week_start)
  years <- sort(unique(as.integer(format(dates, "%Y"))))
  expected_dates <- seq(min(dates), max(dates), by = "week")
  if (!identical(dates, expected_dates)) stop("Weekly data must be complete and consecutive.")
  accounting_error <- max(abs(
    weekly$N_end - weekly$N_start - weekly$births + weekly$all_cause_deaths - weekly$net_population_reconciliation
  ))
  if (accounting_error > 1e-6) stop("Demographic accounting failed before Stan fitting.")
  list(
    stan_data = list(
      N = nrow(weekly), G = 8L, C = as.integer(weekly$cases), w = generation_weights(),
      N_start = as.numeric(weekly$N_start), N_end = as.numeric(weekly$N_end),
      births = as.numeric(weekly$births), deaths = as.numeric(weekly$all_cause_deaths),
      Y = length(years), year_id = match(as.integer(format(dates, "%Y")), years),
      seasonal_sin = sin(2 * pi * seq_len(nrow(weekly)) / 52.1775),
      seasonal_cos = cos(2 * pi * seq_len(nrow(weekly)) / 52.1775),
      imports_per_week = imports_per_week, year_effect_prior_sd = year_effect_prior_sd
    ),
    years = years, demographic_accounting_error = accounting_error
  )
}

# Dispersed-but-plausible q initial values across chains: the 10th, 35th,
# 65th, and 90th percentiles of the SAME Beta(16.2, 108.7) prior used in the
# model -- not identical across chains, but not extreme tail values either.
q_init_values <- qbeta(c(.10, .35, .65, .90), 16.2, 108.7)

v4_1_q_only_initial_values <- function(stan_data) {
  function(chain_id = 1L) {
    offsets <- c(-0.06, -0.02, 0.02, 0.06)
    offset <- offsets[(as.integer(chain_id) - 1L) %% length(offsets) + 1L]
    list(
      alpha_R = log(1.2) + offset,
      z_year = rep(offset / stan_data$year_effect_prior_sd, stan_data$Y),
      beta_sin = offset / 2, beta_cos = -offset / 2, phi_obs = 20,
      q = q_init_values[(as.integer(chain_id) - 1L) %% 4L + 1L]
    )
  }
}

run_fit_v4_1_q_only <- function(stage = Sys.getenv("STAGE", "B")) {
  base_dir <- q_only_dir(root)
  settings <- stage_settings(stage)
  message("[", stage, "] ", settings$label)

  out_dir <- file.path(base_dir, "outputs", paste0("stage", stage))
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  fixed <- list(
    date_start = as.Date("2015-01-04"), date_end = as.Date("2025-12-21"),
    imports_per_week = 1, year_effect_prior_sd = 0.40,
    chains = 4L, adapt_delta = .95, max_treedepth = 14L
  )

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |>
    filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  demography_path <- file.path(root, "01_Data", "ceara_weekly_demography_2015_2025.rds")
  weekly_demography <- readRDS(demography_path) |>
    filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) {
    stop("Case and demographic panels do not have the same full weekly sequence.")
  }

  prepared <- make_v4_1_q_only_data(weekly, fixed$imports_per_week, fixed$year_effect_prior_sd)
  iter_total <- settings$warmup + settings$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  message(sprintf(
    "[%s] %d weeks | warmup=%d sampling=%d | adapt_delta=%.2f max_treedepth=%d | q inits=%s | sample_file base=%s",
    stage, prepared$stan_data$N, settings$warmup, settings$iter_sampling, fixed$adapt_delta, fixed$max_treedepth,
    paste(round(q_init_values, 4), collapse = ","), sample_file_base
  ))

  stan_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_1_q_only.stan")
  refresh <- if (stage == "A") 1L else 25L
  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data,
    chains = fixed$chains, iter = iter_total, warmup = settings$warmup,
    seed = 20260911L,
    init = v4_1_q_only_initial_values(prepared$stan_data),
    control = list(adapt_delta = fixed$adapt_delta, max_treedepth = fixed$max_treedepth),
    refresh = refresh, sample_file = sample_file_base
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs", "q"), fixed$max_treedepth)

  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly,
    config = c(fixed, settings, list(
      model_version = "v4_1_q_only", stage = stage, years = prepared$years,
      demographic_accounting_error = prepared$demographic_accounting_error,
      demography_source = normalizePath(demography_path, winslash = "/", mustWork = TRUE),
      elapsed_seconds = elapsed_seconds, stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE),
      q_role = "estimated single constant reporting fraction, Beta(16.2,108.7) prior; q-ONLY ablation of the HMC-passing v4.0/td14 backbone"
    )),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_1_q_only_fit_stage", stage, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", stage, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  if (stage != "A" && !hmc$hmc_pass) {
    message("[", stage, "] Gate not passed at this stage. Per Section 7, do NOT proceed automatically to the next stage -- report first.")
  }
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_v4_1_q_only()
