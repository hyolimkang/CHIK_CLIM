# Serology leverage/conflict analysis (Section 2): fits the SAME
# successful Stage B4 configuration (v4.0/td14 backbone + estimated q +
# Juazeiro do Norte 2018 serology; adapt_delta=0.95, max_treedepth=14,
# dense_e metric, td14-posterior-draw transmission inits, q inits at prior
# 20th/40th/60th/80th percentiles; warmup=400, sampling=200, 4 chains) at
# each kappa_sero on a fixed ladder, changing NOTHING else. kappa=50 is
# Stage B4 itself and is not refit -- reused directly. A final "binomial"
# run uses the literal binomial(sero_n, p_state_sero) limiting model as a
# stress test only (NOT a proposed production model -- Section 2).
#
# Usage: KAPPA=100 Rscript 05_fit_kappa_ladder.R
#        VARIANT=binomial Rscript 05_fit_kappa_ladder.R

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

kappa_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

root <- kappa_root()
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))

base_dir <- file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "21_v4_2_q_juazeiro_serology")

make_v4_2_data <- function(weekly, imports_per_week, year_effect_prior_sd, t_sero, sero_pos, sero_n, kappa_sero = NULL) {
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
    t_sero = t_sero, sero_pos = sero_pos, sero_n = sero_n
  )
  if (!is.null(kappa_sero)) stan_data$kappa_sero <- kappa_sero
  list(stan_data = stan_data, years = years, demographic_accounting_error = accounting_error)
}

find_sero_week_index <- function(dates) {
  survey_start <- as.Date("2018-06-01"); survey_end <- as.Date("2018-12-31")
  midpoint <- survey_start + floor(as.numeric(survey_end - survey_start) / 2)
  idx <- which.min(abs(dates - midpoint))
  list(index = idx, midpoint_date = midpoint, matched_week_start = dates[idx])
}

q_init_values <- qbeta(c(.20, .40, .60, .80), 16.2, 108.7)

load_td14_chain_inits <- function(root, n_chains = 4L, seed = 20260911L) {
  td14_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_0_minimal_no_vaccine_fit_td14.rds")
  td14 <- readRDS(td14_path)
  draws <- rstan::extract(td14$fit, pars = c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs"), permuted = TRUE)
  set.seed(seed)
  idx <- sample.int(nrow(draws$alpha_R), n_chains)
  lapply(idx, function(i) list(
    alpha_R = draws$alpha_R[i], z_year = as.numeric(draws$z_year[i, ]),
    beta_sin = draws$beta_sin[i], beta_cos = draws$beta_cos[i], phi_obs = draws$phi_obs[i]
  ))
}

v4_2_initial_values <- function(td14_chain_inits) {
  function(chain_id = 1L) {
    cid <- as.integer(chain_id)
    base <- td14_chain_inits[[(cid - 1L) %% length(td14_chain_inits) + 1L]]
    c(base, list(q = q_init_values[(cid - 1L) %% 4L + 1L]))
  }
}

run_kappa_ladder <- function(kappa = as.numeric(Sys.getenv("KAPPA", "100")), variant = Sys.getenv("VARIANT", "betabinomial")) {
  is_binomial <- identical(variant, "binomial")
  tag <- if (is_binomial) "binomial" else paste0("kappa", kappa)
  out_dir <- file.path(base_dir, "outputs", paste0("stageB4_", tag))
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  fixed <- list(
    date_start = as.Date("2015-01-04"), date_end = as.Date("2025-12-21"),
    imports_per_week = 1, year_effect_prior_sd = 0.40,
    chains = 4L, max_treedepth = 14L, adapt_delta = .95, metric = "dense_e",
    warmup = 400L, iter_sampling = 200L,
    sero_pos = 103L, sero_n = 404L
  )

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |> filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  demography_path <- file.path(root, "01_Data", "ceara_weekly_demography_2015_2025.rds")
  weekly_demography <- readRDS(demography_path) |> filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch.")

  sero_week <- find_sero_week_index(as.Date(weekly$week_start))
  prepared <- make_v4_2_data(weekly, fixed$imports_per_week, fixed$year_effect_prior_sd,
                              sero_week$index, fixed$sero_pos, fixed$sero_n,
                              kappa_sero = if (is_binomial) NULL else kappa)
  iter_total <- fixed$warmup + fixed$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  td14_chain_inits <- load_td14_chain_inits(root)

  message(sprintf("[%s] variant=%s kappa=%s | %d weeks | warmup=%d sampling=%d adapt_delta=%.2f metric=%s",
                   tag, variant, if (is_binomial) "Inf(binomial)" else kappa, prepared$stan_data$N,
                   fixed$warmup, fixed$iter_sampling, fixed$adapt_delta, fixed$metric))

  stan_path <- file.path(root, "02_Script", "stan",
                          if (is_binomial) "renewal_ceara_v4_2_q_juazeiro_serology_binomial.stan" else "renewal_ceara_v4_2_q_juazeiro_serology.stan")
  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data,
    chains = fixed$chains, iter = iter_total, warmup = fixed$warmup,
    seed = 20260911L,
    init = v4_2_initial_values(td14_chain_inits),
    control = list(adapt_delta = fixed$adapt_delta, max_treedepth = fixed$max_treedepth, metric = fixed$metric),
    refresh = 25, sample_file = sample_file_base
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs", "q"), fixed$max_treedepth)

  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, sero_week = sero_week,
    config = c(fixed, list(
      model_version = "v4_2_kappa_ladder", tag = tag, variant = variant, kappa_sero = if (is_binomial) Inf else kappa,
      years = prepared$years, elapsed_seconds = elapsed_seconds, stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE)
    )),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_2_fit_", tag, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", tag, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_kappa_ladder()
