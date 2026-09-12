# Staged fitting for v4.2 (v4.0/td14 backbone + estimated q + Juazeiro do
# Norte 2018 serology anchor). Reuses v4.0's exact helper functions
# (generation_weights, compute_hmc_gate, build_state_weekly) by sourcing the
# archived-checksum-verified fit script, not reimplementing them.
#
# STAGE controls warmup/sampling length only (Section 12):
#   A: 5+5     (implementation validation only, not for gate assessment)
#   B: 300+150 (geometry canary)
#   C: 500+300 (intermediate)
#   D: 1000+1000 (production, matches v4.0/td14 exactly)
#
# Usage: STAGE=B Rscript 02_fit_v4_2.R

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

v4_2_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

root <- v4_2_root()
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "40_renewal_model", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))

v4_2_dir <- function(root = v4_2_root()) file.path(root, "02_Script", "40_renewal_model", "21_v4_2_q_juazeiro_serology")

stage_settings <- function(stage) {
  switch(stage,
    A = list(warmup = 5L,    iter_sampling = 5L,    adapt_delta = .95, metric = "diag_e", label = "Stage A: implementation validation only -- NOT for gate assessment"),
    B = list(warmup = 300L,  iter_sampling = 150L,  adapt_delta = .95, metric = "diag_e", label = "Stage B: geometry canary"),
    B2 = list(warmup = 400L, iter_sampling = 200L,  adapt_delta = .99, metric = "diag_e", label = "Stage B2: adapt_delta=0.99 only -- FAILED (worse treedepth; see FAILED_STAGE_B2_ADAPT099_TREEDEPTH)"),
    B4 = list(warmup = 400L, iter_sampling = 200L,  adapt_delta = .95, metric = "dense_e", label = "Stage B4: adapt_delta=0.95 (back to Stage B's value) + dense_e metric -- isolates whether strong posterior correlations (not step size) are the difficulty"),
    C = list(warmup = 500L,  iter_sampling = 300L,  adapt_delta = .95, metric = "diag_e", label = "Stage C: intermediate run"),
    D = list(warmup = 1000L, iter_sampling = 1000L, adapt_delta = .95, metric = "diag_e", label = "Stage D: production run (matches v4.0/td14 exactly)"),
    stop("Unknown STAGE '", stage, "' -- must be one of A, B, B2, B4, C, D.")
  )
}

make_v4_2_data <- function(weekly, imports_per_week, year_effect_prior_sd, t_sero, sero_pos, sero_n, kappa_sero) {
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
      imports_per_week = imports_per_week, year_effect_prior_sd = year_effect_prior_sd,
      t_sero = t_sero, sero_pos = sero_pos, sero_n = sero_n, kappa_sero = kappa_sero
    ),
    years = years, demographic_accounting_error = accounting_error
  )
}

# Section 5: survey window June-December 2018; midpoint used as the primary
# anchor date (no individual sample dates available in this repository).
find_sero_week_index <- function(dates) {
  survey_start <- as.Date("2018-06-01"); survey_end <- as.Date("2018-12-31")
  midpoint <- survey_start + floor(as.numeric(survey_end - survey_start) / 2)
  idx <- which.min(abs(dates - midpoint))
  list(index = idx, midpoint_date = midpoint, matched_week_start = dates[idx])
}

# Stage B used q inits at the prior's 10th/35th/65th/90th percentiles.
# Stage B2 onward (Section 4) uses 20th/40th/60th/80th -- dispersed but
# further from the tails, since B already showed the extreme-offset
# transmission-parameter inits were not themselves the source of trouble.
q_percentiles_for_stage <- function(stage) if (stage %in% c("A", "B")) c(.10, .35, .65, .90) else c(.20, .40, .60, .80)

# Section 4: prefer initialising transmission parameters from DIFFERENT
# posterior draws of the successful v4.0/td14 fit (one distinct draw per
# chain) rather than deterministic offsets, from Stage B2 onward. Falls
# back to the offset scheme (used for Stage A/B, and if td14 draws cannot
# be loaded) so this remains robust.
load_td14_chain_inits <- function(root, n_chains = 4L, seed = 20260911L) {
  td14_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_0_minimal_no_vaccine_fit_td14.rds")
  if (!file.exists(td14_path)) return(NULL)
  td14 <- readRDS(td14_path)
  draws <- rstan::extract(td14$fit, pars = c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs"), permuted = TRUE)
  set.seed(seed)
  idx <- sample.int(nrow(draws$alpha_R), n_chains)
  lapply(idx, function(i) list(
    alpha_R = draws$alpha_R[i], z_year = as.numeric(draws$z_year[i, ]),
    beta_sin = draws$beta_sin[i], beta_cos = draws$beta_cos[i], phi_obs = draws$phi_obs[i]
  ))
}

v4_2_initial_values <- function(stan_data, stage, td14_chain_inits = NULL) {
  q_init_values <- qbeta(q_percentiles_for_stage(stage), 16.2, 108.7)
  use_td14_inits <- stage != "A" && stage != "B" && !is.null(td14_chain_inits)
  function(chain_id = 1L) {
    cid <- as.integer(chain_id)
    q_val <- q_init_values[(cid - 1L) %% 4L + 1L]
    if (use_td14_inits) {
      base <- td14_chain_inits[[(cid - 1L) %% length(td14_chain_inits) + 1L]]
      return(c(base, list(q = q_val)))
    }
    offsets <- c(-0.06, -0.02, 0.02, 0.06)
    offset <- offsets[(cid - 1L) %% length(offsets) + 1L]
    list(
      alpha_R = log(1.2) + offset,
      z_year = rep(offset / stan_data$year_effect_prior_sd, stan_data$Y),
      beta_sin = offset / 2, beta_cos = -offset / 2, phi_obs = 20,
      q = q_val
    )
  }
}

run_fit_v4_2 <- function(stage = Sys.getenv("STAGE", "B")) {
  base_dir <- v4_2_dir(root)
  settings <- stage_settings(stage)
  message("[", stage, "] ", settings$label)

  out_dir <- file.path(base_dir, "outputs", paste0("stage", stage))
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  fixed <- list(
    date_start = as.Date("2015-01-04"), date_end = as.Date("2025-12-21"),
    imports_per_week = 1, year_effect_prior_sd = 0.40,
    chains = 4L, max_treedepth = 14L,
    sero_pos = 103L, sero_n = 404L, kappa_sero = 50
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

  sero_week <- find_sero_week_index(as.Date(weekly$week_start))
  message(sprintf("[%s] serology anchor: survey midpoint=%s -> matched week_start=%s (index %d of %d)",
                   stage, sero_week$midpoint_date, sero_week$matched_week_start, sero_week$index, nrow(weekly)))

  prepared <- make_v4_2_data(weekly, fixed$imports_per_week, fixed$year_effect_prior_sd,
                              sero_week$index, fixed$sero_pos, fixed$sero_n, fixed$kappa_sero)
  iter_total <- settings$warmup + settings$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  q_init_values <- qbeta(q_percentiles_for_stage(stage), 16.2, 108.7)
  td14_chain_inits <- if (stage %in% c("A", "B")) NULL else load_td14_chain_inits(root)
  message(sprintf(
    "[%s] %d weeks | warmup=%d sampling=%d | adapt_delta=%.2f metric=%s | q inits=%s | init source=%s | sero=%d/%d kappa=%g | sample_file base=%s",
    stage, prepared$stan_data$N, settings$warmup, settings$iter_sampling, settings$adapt_delta, settings$metric,
    paste(round(q_init_values, 4), collapse = ","), if (!is.null(td14_chain_inits)) "v4.0/td14 posterior draws" else "deterministic offsets",
    fixed$sero_pos, fixed$sero_n, fixed$kappa_sero, sample_file_base
  ))

  stan_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_2_q_juazeiro_serology.stan")
  refresh <- if (stage == "A") 1L else 25L
  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data,
    chains = fixed$chains, iter = iter_total, warmup = settings$warmup,
    seed = 20260911L,
    init = v4_2_initial_values(prepared$stan_data, stage, td14_chain_inits),
    control = list(adapt_delta = settings$adapt_delta, max_treedepth = fixed$max_treedepth, metric = settings$metric),
    refresh = refresh, sample_file = sample_file_base
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs", "q"), fixed$max_treedepth)

  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, sero_week = sero_week,
    config = c(fixed, settings, list(
      model_version = "v4_2_q_juazeiro_serology", stage = stage, years = prepared$years,
      demographic_accounting_error = prepared$demographic_accounting_error,
      demography_source = normalizePath(demography_path, winslash = "/", mustWork = TRUE),
      elapsed_seconds = elapsed_seconds, stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE),
      q_role = "estimated single constant reporting fraction, Beta(16.2,108.7) prior, stabilised (if it works) by ONE Juazeiro do Norte 2018 serology anchor via beta-binomial(kappa_sero=50)"
    )),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_2_fit_stage", stage, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", stage, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  if (stage != "A" && !hmc$hmc_pass) {
    message("[", stage, "] Gate not passed at this stage. Per Section 13, do NOT proceed automatically to the next stage -- report first.")
  }
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_v4_2()
