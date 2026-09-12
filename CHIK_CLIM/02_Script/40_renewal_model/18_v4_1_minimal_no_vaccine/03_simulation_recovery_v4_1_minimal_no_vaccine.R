# ===========================================================================
# 03_simulation_recovery_v4_1_minimal_no_vaccine.R
#
# Purpose: prespecified parameter-recovery check for v4.1, mandatory before
# calling v4.1 a stabilized backbone. Simulates from known q, alpha_R,
# sigma_year, beta_sin, beta_cos, phi_obs using the exact v4.1 recurrence and
# the real weekly Ceará demography, then refits and checks recovery of those
# six parameters plus S_2022, S_2025, and cumulative infections (2015-2025).
# Particular attention: q and S are the two quantities newly identified (or
# newly uncertain) relative to v4.0's fixed-q design.
# ===========================================================================

recovery_project_root <- function() {
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(current, "CHIK_CLIM.Rproj"))) return(current)
  candidate <- file.path(current, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Run from the outer repository or inner CHIK_CLIM R-project root.")
}

source(file.path(
  recovery_project_root(), "02_Script", "40_renewal_model",
  "18_v4_1_minimal_no_vaccine", "01_fit_v4_1_minimal_no_vaccine.R"
))

simulate_v4_1_trajectory <- function(stan_data, truth) {
  N <- stan_data$N; G <- stan_data$G
  S <- numeric(N); U <- numeric(N); X <- numeric(N); R0_t <- numeric(N); expected_cases <- numeric(N)
  reconciliation <- stan_data$N_end - stan_data$N_start - stan_data$births + stan_data$deaths
  year_effect <- truth$sigma_year * as.numeric(stan_data$B_year %*% truth$z_year_free)
  S[1] <- stan_data$N_start[1]

  for (t in seq_len(N)) {
    total <- S[t] + U[t]
    infectiousness <- sum(vapply(seq_len(G), function(g) if (t > g) stan_data$w[g] * X[t - g] else 0, numeric(1)))
    R0_t[t] <- exp(truth$alpha_R + year_effect[stan_data$year_id[t]] + truth$beta_sin * stan_data$seasonal_sin[t] + truth$beta_cos * stan_data$seasonal_cos[t])
    force <- R0_t[t] * infectiousness / total + stan_data$imports_per_week / total
    X[t] <- S[t] * (-expm1(-force))
    expected_cases[t] <- truth$q * X[t] + 1e-9
    if (t < N) {
      susceptible_after <- S[t] - X[t]; immune_after <- U[t] + X[t]
      susceptible_fraction <- susceptible_after / (susceptible_after + immune_after)
      S[t + 1] <- susceptible_after + stan_data$births[t] - stan_data$deaths[t] * susceptible_fraction + reconciliation[t] * susceptible_fraction
      U[t + 1] <- immune_after - stan_data$deaths[t] * (1 - susceptible_fraction) + reconciliation[t] * (1 - susceptible_fraction)
    }
  }
  list(X = X, S = S, R0_t = R0_t, expected_cases = expected_cases, year_effect = year_effect)
}

run_v4_1_simulation_recovery <- function() {
  root <- v4_1_root()
  paths <- v4_1_paths(root)
  dir.create(paths$tables, recursive = TRUE, showWarnings = FALSE)

  settings <- list(
    q_prior_a = env_number("RENEWAL_V4_1_Q_PRIOR_A", 16.2, lower = 1e-6),
    q_prior_b = env_number("RENEWAL_V4_1_Q_PRIOR_B", 108.7, lower = 1e-6),
    imports_per_week = env_number("RENEWAL_V4_1_IMPORTS_PER_WEEK", 1, lower = 0),
    iter = env_integer("RENEWAL_V4_1_RECOVERY_ITER", 1600L),
    warmup = env_integer("RENEWAL_V4_1_RECOVERY_WARMUP", 800L),
    chains = env_integer("RENEWAL_V4_1_RECOVERY_CHAINS", 4L),
    adapt_delta = 0.95, max_treedepth = 14L
  )
  if (settings$warmup >= settings$iter) stop("Recovery warmup must be smaller than iterations.")

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- dplyr::filter(build_state_weekly(panel, "23"), week_start >= as.Date("2015-01-04"), week_start <= as.Date("2025-12-21"))
  weekly_demography <- dplyr::filter(readRDS(file.path(root, "01_Data", "ceara_weekly_demography_2015_2025.rds")), week_start >= as.Date("2015-01-04"), week_start <= as.Date("2025-12-21"))
  weekly <- dplyr::arrange(dplyr::inner_join(dplyr::select(weekly_cases, -population), weekly_demography, by = "week_start"), week_start)
  prepared <- make_v4_1_data(weekly, settings$q_prior_a, settings$q_prior_b, settings$imports_per_week)
  years <- as.integer(format(as.Date(weekly$week_start), "%Y"))

  # Fixed before fitting: q near the prior mean, moderate annual/seasonal
  # variation on the scale allowed by the priors.
  q_prior_mean <- settings$q_prior_a / (settings$q_prior_a + settings$q_prior_b)
  truth <- list(
    q = q_prior_mean, alpha_R = log(1.45), sigma_year = 0.30,
    z_year_free = seq(-1.1, 1.1, length.out = prepared$stan_data$Y - 1L),
    beta_sin = 0.12, beta_cos = -0.08, phi_obs = 25
  )
  simulated <- simulate_v4_1_trajectory(prepared$stan_data, truth)
  synthetic_data <- prepared$stan_data
  set.seed(20260912L)
  synthetic_data$C <- as.integer(stats::rnbinom(synthetic_data$N, mu = simulated$expected_cases, size = truth$phi_obs))

  message("[v4.1 recovery] fitting one pre-specified synthetic Ceará trajectory (q=", round(truth$q, 4), ", sigma_year=", truth$sigma_year, ")")
  started <- Sys.time()
  recovery_fit <- rstan::stan(
    file = paths$stan, data = synthetic_data,
    chains = settings$chains, iter = settings$iter, warmup = settings$warmup,
    seed = 20260912L, init = v4_1_initial_values(synthetic_data),
    control = list(adapt_delta = settings$adapt_delta, max_treedepth = settings$max_treedepth),
    refresh = max(1L, floor(settings$iter / 10L))
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  core_pars <- c("alpha_R", "sigma_year", "beta_sin", "beta_cos", "phi_obs", "q")
  fit_summary <- as.data.frame(rstan::summary(recovery_fit, pars = c(core_pars, "year_effect"))$summary)
  fit_summary$parameter <- rownames(fit_summary); rownames(fit_summary) <- NULL
  truth_table <- tibble::tibble(
    parameter = c(core_pars, paste0("year_effect[", seq_len(prepared$stan_data$Y), "]")),
    truth = c(truth$q, truth$alpha_R, truth$sigma_year, truth$beta_sin, truth$beta_cos, truth$phi_obs, simulated$year_effect)
  )
  recovery <- dplyr::left_join(truth_table, fit_summary, by = "parameter") |>
    dplyr::transmute(parameter, truth, posterior_mean = mean, posterior_median = `50%`, lower_95 = `2.5%`, upper_95 = `97.5%`,
                      rhat = Rhat, bulk_ESS = n_eff, covered_by_95pct_interval = truth >= `2.5%` & truth <= `97.5%`)

  # Derived-quantity recovery: S_2022, S_2025, cumulative infections.
  X_draws <- rstan::extract(recovery_fit, pars = "X")$X
  S_draws <- rstan::extract(recovery_fit, pars = "S")$S
  idx_2022 <- which(years == 2022L)[1]; idx_2025 <- which(years == 2025L)[1]
  cumulative_truth <- sum(simulated$X)
  derived_recovery <- tibble::tibble(
    parameter = c("S_2022", "S_2025", "cumulative_infections_2015_2025"),
    truth = c(simulated$S[idx_2022], simulated$S[idx_2025], cumulative_truth),
    posterior_mean = c(mean(S_draws[, idx_2022]), mean(S_draws[, idx_2025]), mean(rowSums(X_draws))),
    posterior_median = c(median(S_draws[, idx_2022]), median(S_draws[, idx_2025]), median(rowSums(X_draws))),
    lower_95 = c(quantile(S_draws[, idx_2022], .025), quantile(S_draws[, idx_2025], .025), quantile(rowSums(X_draws), .025)),
    upper_95 = c(quantile(S_draws[, idx_2022], .975), quantile(S_draws[, idx_2025], .975), quantile(rowSums(X_draws), .975)),
    rhat = NA_real_, bulk_ESS = NA_real_
  ) |> dplyr::mutate(covered_by_95pct_interval = truth >= lower_95 & truth <= upper_95)
  recovery <- dplyr::bind_rows(recovery, derived_recovery)

  hmc <- compute_hmc_gate(recovery_fit, c(core_pars, "z_year_free"), settings$max_treedepth)
  recovery_summary <- tibble::tibble(
    metric = c("all_parameters_95pct_coverage", "core_parameters_95pct_coverage", "q_and_S_95pct_coverage", "HMC_gate_pass", "elapsed_seconds"),
    value = c(
      mean(recovery$covered_by_95pct_interval),
      mean(recovery$covered_by_95pct_interval[recovery$parameter %in% core_pars]),
      mean(recovery$covered_by_95pct_interval[recovery$parameter %in% c("q", "S_2022", "S_2025")]),
      as.numeric(hmc$hmc_pass), elapsed_seconds
    )
  )
  readr::write_csv(recovery, file.path(paths$tables, "simulation_recovery.csv"))
  readr::write_csv(recovery_summary, file.path(paths$tables, "simulation_recovery_summary.csv"))
  saveRDS(list(fit = recovery_fit, synthetic_data = synthetic_data, truth = truth, simulated = simulated, hmc = hmc, config = settings, elapsed_seconds = elapsed_seconds),
          file.path(root, "02_Script", "stan", "renewal_ceara_v4_1_minimal_no_vaccine_recovery_fit.rds"))
  message("[v4.1 recovery] HMC gate: ", if (hmc$hmc_pass) "PASS" else "FAIL", "; q coverage: ", recovery$covered_by_95pct_interval[recovery$parameter == "q"])
  invisible(list(recovery = recovery, hmc = hmc))
}

if (sys.nframe() == 0L) run_v4_1_simulation_recovery()
