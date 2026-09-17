# ===========================================================================
# 03_simulation_recovery_v4_0_minimal_no_vaccine.R
#
# Purpose: prespecified parameter-recovery check for the v4.0 model.  The
# synthetic data use the same weekly Ceará demography and the exact model
# recurrence, but not the observed case data.  This is a computational
# validation check, not a calibration exercise.
# ===========================================================================

recovery_project_root <- function() {
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(current, "CHIK_CLIM.Rproj"))) return(current)
  candidate <- file.path(current, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Run from the outer repository or inner CHIK_CLIM R-project root.")
}

source(file.path(
  recovery_project_root(), "02_Script", "00_shared", "legacy_model_functions", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"
))

simulate_v4_trajectory <- function(stan_data, truth) {
  N <- stan_data$N
  G <- stan_data$G
  S <- numeric(N)
  U <- numeric(N)
  X <- numeric(N)
  R0_t <- numeric(N)
  expected_cases <- numeric(N)
  reconciliation <- stan_data$N_end - stan_data$N_start - stan_data$births + stan_data$deaths

  year_effect <- stan_data$year_effect_prior_sd * (truth$z_year - mean(truth$z_year))
  S[1] <- stan_data$N_start[1]

  for (t in seq_len(N)) {
    total <- S[t] + U[t]
    infectiousness <- sum(vapply(seq_len(G), function(g) {
      if (t > g) stan_data$w[g] * X[t - g] else 0
    }, numeric(1)))
    R0_t[t] <- exp(
      truth$alpha_R + year_effect[stan_data$year_id[t]] +
        truth$beta_sin * stan_data$seasonal_sin[t] +
        truth$beta_cos * stan_data$seasonal_cos[t]
    )
    force <- R0_t[t] * infectiousness / total + stan_data$imports_per_week / total
    X[t] <- S[t] * (-expm1(-force))
    expected_cases[t] <- stan_data$q_fixed * X[t] + 1e-9

    if (t < N) {
      susceptible_after <- S[t] - X[t]
      immune_after <- U[t] + X[t]
      susceptible_fraction <- susceptible_after / (susceptible_after + immune_after)
      S[t + 1] <- susceptible_after + stan_data$births[t] -
        stan_data$deaths[t] * susceptible_fraction + reconciliation[t] * susceptible_fraction
      U[t + 1] <- immune_after - stan_data$deaths[t] * (1 - susceptible_fraction) +
        reconciliation[t] * (1 - susceptible_fraction)
    }
  }

  list(X = X, R0_t = R0_t, expected_cases = expected_cases)
}

run_v4_0_simulation_recovery <- function() {
  root <- v4_root()
  paths <- v4_paths(root)
  dir.create(paths$tables, recursive = TRUE, showWarnings = FALSE)

  settings <- list(
    q_fixed = env_number("RENEWAL_V4_Q_FIXED", 0.10, lower = 1e-6, upper = 1 - 1e-6),
    imports_per_week = env_number("RENEWAL_V4_IMPORTS_PER_WEEK", 1, lower = 0),
    year_effect_prior_sd = env_number("RENEWAL_V4_YEAR_EFFECT_PRIOR_SD", 0.40, lower = 1e-6),
    iter = env_integer("RENEWAL_V4_RECOVERY_ITER", 1600L),
    warmup = env_integer("RENEWAL_V4_RECOVERY_WARMUP", 800L),
    chains = env_integer("RENEWAL_V4_RECOVERY_CHAINS", 4L),
    adapt_delta = 0.95,
    max_treedepth = 12L
  )
  if (settings$warmup >= settings$iter) stop("Recovery warmup must be smaller than iterations.")

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- dplyr::filter(build_state_weekly(panel, "23"),
    week_start >= as.Date("2015-01-04"), week_start <= as.Date("2025-12-21")
  )
  weekly_demography <- dplyr::filter(
    readRDS(file.path(root, "01_Data", "ceara_weekly_demography_2015_2025.rds")),
    week_start >= as.Date("2015-01-04"), week_start <= as.Date("2025-12-21")
  )
  weekly <- dplyr::arrange(dplyr::inner_join(
    dplyr::select(weekly_cases, -population), weekly_demography, by = "week_start"
  ), week_start)
  prepared <- make_v4_data(weekly, settings$q_fixed, settings$imports_per_week, settings$year_effect_prior_sd)

  # Fixed before fitting: moderate annual variation and seasonal effects on
  # the scale allowed by the fitted priors.
  truth <- list(
    alpha_R = log(1.45),
    z_year = seq(-0.8, 0.8, length.out = prepared$stan_data$Y),
    beta_sin = 0.12,
    beta_cos = -0.08,
    phi_obs = 25
  )
  simulated <- simulate_v4_trajectory(prepared$stan_data, truth)
  synthetic_data <- prepared$stan_data
  set.seed(20260912L)
  synthetic_data$C <- as.integer(stats::rnbinom(
    synthetic_data$N, mu = simulated$expected_cases, size = truth$phi_obs
  ))

  message("[v4.0 recovery] fitting one pre-specified synthetic Ceará trajectory")
  started <- Sys.time()
  recovery_fit <- rstan::stan(
    file = paths$stan, data = synthetic_data,
    chains = settings$chains, iter = settings$iter, warmup = settings$warmup,
    seed = 20260912L,
    init = v4_initial_values(synthetic_data),
    control = list(adapt_delta = settings$adapt_delta, max_treedepth = settings$max_treedepth),
    refresh = max(1L, floor(settings$iter / 10L))
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  pars <- c("alpha_R", "beta_sin", "beta_cos", "phi_obs", "year_effect")
  fit_summary <- as.data.frame(rstan::summary(recovery_fit, pars = pars)$summary)
  fit_summary$parameter <- rownames(fit_summary)
  rownames(fit_summary) <- NULL

  truth_table <- tibble::tibble(
    parameter = c("alpha_R", "beta_sin", "beta_cos", "phi_obs", paste0("year_effect[", seq_len(prepared$stan_data$Y), "]")),
    truth = c(truth$alpha_R, truth$beta_sin, truth$beta_cos, truth$phi_obs,
      settings$year_effect_prior_sd * (truth$z_year - mean(truth$z_year)))
  )
  recovery <- dplyr::left_join(truth_table, fit_summary, by = "parameter") |>
    dplyr::transmute(
      parameter, truth, posterior_mean = mean, posterior_median = `50%`,
      lower_95 = `2.5%`, upper_95 = `97.5%`, rhat = Rhat, bulk_ESS = n_eff,
      covered_by_95pct_interval = truth >= `2.5%` & truth <= `97.5%`
    )
  hmc <- compute_hmc_gate(recovery_fit, pars, settings$max_treedepth)
  recovery_summary <- tibble::tibble(
    metric = c("all_parameters_95pct_coverage", "core_parameters_95pct_coverage", "HMC_gate_pass", "elapsed_seconds"),
    value = c(
      mean(recovery$covered_by_95pct_interval),
      mean(recovery$covered_by_95pct_interval[recovery$parameter %in% c("alpha_R", "beta_sin", "beta_cos", "phi_obs")]),
      as.numeric(hmc$hmc_pass), elapsed_seconds
    )
  )
  readr::write_csv(recovery, file.path(paths$tables, "simulation_recovery.csv"))
  readr::write_csv(recovery_summary, file.path(paths$tables, "simulation_recovery_summary.csv"))
  saveRDS(list(
    fit = recovery_fit, synthetic_data = synthetic_data, truth = truth,
    simulated = simulated, hmc = hmc, config = settings, elapsed_seconds = elapsed_seconds
  ), file.path(root, "03_Output", "02_ceara_pipeline", "model_fits", "initialization", "renewal_ceara_v4_0_minimal_no_vaccine_recovery_fit.rds"))
  message("[v4.0 recovery] HMC gate: ", if (hmc$hmc_pass) "PASS" else "FAIL")
  invisible(list(recovery = recovery, hmc = hmc))
}

if (sys.nframe() == 0L) run_v4_0_simulation_recovery()
