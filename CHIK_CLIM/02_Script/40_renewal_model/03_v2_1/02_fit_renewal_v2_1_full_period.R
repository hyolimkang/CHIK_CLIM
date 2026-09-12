# =============================================================================
# 02_fit_renewal_v2_1_full_period.R
#
# Diagnostic stress-test extension of renewal v2.1 through December 2025.
# This script uses the unchanged v2.1 Stan source, priors, and likelihood. It
# is not the final full-history model and should not overwrite the v2.1 fit.
#
# Environment overrides:
#   RENEWAL_V2_1_FULL_ITER, RENEWAL_V2_1_FULL_WARMUP,
#   RENEWAL_V2_1_FULL_CHAINS, RENEWAL_V2_1_FULL_TAG
# =============================================================================

# Source v2.1 for its package checks and shared non-model helpers. Its main
# block is guarded, so sourcing does not run a fit.
source(here::here("02_Script/40_renewal_model/03_v2_1/02_fit_renewal_v2_1.R"))

if (sys.nframe() == 0) {
  UF_CODE <- "23"
  DATE_START <- as.Date("2015-01-04")
  DATE_END <- as.Date("2025-12-28")
  G <- 8L
  SEED_WEEKS <- 8L
  GI_MEAN <- 2
  GI_SD <- 1

  # These are exactly the v2.1 data assumptions, retained for the stress test.
  REPORTING_REGION_SD <- 0.75
  SERO_SITE <- c("Juazeiro do Norte", "Quixada")
  SERO_WINDOW_START <- as.Date(c("2018-06-03", "2018-06-03"))
  SERO_WINDOW_END <- as.Date(c("2018-12-30", "2019-12-29"))
  SERO_POSITIVE <- c(103L, 289L)
  SERO_N <- c(404L, 409L)

  N_ITER <- env_integer("RENEWAL_V2_1_FULL_ITER", 2000L)
  N_WARMUP <- env_integer("RENEWAL_V2_1_FULL_WARMUP", 1000L)
  N_CHAINS <- env_integer("RENEWAL_V2_1_FULL_CHAINS", 4L)
  FIT_SEED <- 20251228L
  OUTPUT_TAG <- Sys.getenv("RENEWAL_V2_1_FULL_TAG", unset = "")
  if (N_WARMUP >= N_ITER) {
    stop("RENEWAL_V2_1_FULL_WARMUP must be smaller than RENEWAL_V2_1_FULL_ITER")
  }

  panel <- readRDS(here::here("01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
  ceara_weekly <- build_state_weekly(panel, UF_CODE) |>
    dplyr::filter(week_start >= DATE_START, week_start <= DATE_END)
  expected_dates <- seq(DATE_START, DATE_END, by = "week")
  if (!identical(ceara_weekly$week_start, expected_dates)) {
    stop("Ceara data do not form the expected complete, ordered weekly sequence")
  }
  if (anyNA(ceara_weekly$cases) || any(ceara_weekly$cases < 0) ||
    any(abs(ceara_weekly$cases - round(ceara_weekly$cases)) > 1e-8)) {
    stop("Weekly case counts must be complete, non-negative integers")
  }
  if (SEED_WEEKS < G || SEED_WEEKS >= nrow(ceara_weekly)) {
    stop("Require G <= SEED_WEEKS < number of weeks")
  }

  sero_window_start <- match(SERO_WINDOW_START, ceara_weekly$week_start)
  sero_window_end <- match(SERO_WINDOW_END, ceara_weekly$week_start)
  if (anyNA(sero_window_start) || anyNA(sero_window_end) ||
    any(sero_window_start > sero_window_end)) {
    stop("Serology collection windows must be complete, ordered weekly intervals")
  }

  N_pop_t <- as.numeric(ceara_weekly$population)
  if (length(N_pop_t) != nrow(ceara_weekly) ||
    any(!is.finite(N_pop_t)) || any(N_pop_t < 1)) {
    stop("Weekly state population must be finite and at least one in every week")
  }
  w <- discretize_gamma_generation_interval(GI_MEAN, GI_SD, G)
  time_scaled <- seq(-0.5, 0.5, length.out = nrow(ceara_weekly))
  log_seed_prior_mean <- rep(log(10), SEED_WEEKS)

  stan_data_v2_1_full_period <- list(
    N = nrow(ceara_weekly),
    G = G,
    seed_weeks = SEED_WEEKS,
    C = as.integer(ceara_weekly$cases),
    w = as.vector(w),
    N_pop_t = as.vector(N_pop_t),
    time_scaled = as.vector(time_scaled),
    log_seed_prior_mean = log_seed_prior_mean,
    seed_prior_sd = 1.5,
    reporting_region_sd = REPORTING_REGION_SD,
    J = length(SERO_SITE),
    sero_window_start = as.integer(sero_window_start),
    sero_window_end = as.integer(sero_window_end),
    sero_positive = SERO_POSITIVE,
    sero_n = SERO_N
  )

  # The initialization follows v2.1 exactly, including its demographic update.
  smooth_cases <- vapply(seq_len(stan_data_v2_1_full_period$N), function(t) {
    index <- max(1L, t - 2L):min(stan_data_v2_1_full_period$N, t + 2L)
    mean(stan_data_v2_1_full_period$C[index] + 0.5)
  }, numeric(1))
  rough_X <- smooth_cases * (0.25 * N_pop_t[1] / sum(smooth_cases))
  rough_log_hazard <- numeric(stan_data_v2_1_full_period$N)
  rough_susceptible <- N_pop_t[1]
  for (t in seq_len(stan_data_v2_1_full_period$N)) {
    rough_X[t] <- min(rough_X[t], 0.2 * rough_susceptible)
    rough_log_hazard[t] <- log(-log1p(-rough_X[t] / rough_susceptible))
    rough_susceptible <- rough_susceptible - rough_X[t]
    if (t < stan_data_v2_1_full_period$N) {
      population_change <- N_pop_t[t + 1L] - N_pop_t[t]
      if (population_change >= 0) {
        rough_susceptible <- rough_susceptible + population_change
      } else {
        rough_susceptible <- rough_susceptible * N_pop_t[t + 1L] / N_pop_t[t]
      }
    }
  }

  init_fn <- function() {
    list(
      sigma_R = 0.05,
      log_infection_scale = rough_log_hazard[1],
      log_hazard_relative = rough_log_hazard[-1] - rough_log_hazard[1],
      p_symp = 0.52,
      rho_sym_brazil = 0.25,
      logit_rho_sym_ceara_mid = stats::qlogis(0.116),
      reporting_trend = 0.3,
      sigma_geo = 0.5,
      z_site = rep(0, length(SERO_SITE)),
      phi_obs = 20
    )
  }

  message(sprintf(
    "[stress test data] %s to %s | N=%d | cases=%s | population range=%s to %s",
    DATE_START, DATE_END, nrow(ceara_weekly),
    format(sum(ceara_weekly$cases), big.mark = ","),
    format(min(N_pop_t), big.mark = ",", scientific = FALSE),
    format(max(N_pop_t), big.mark = ",", scientific = FALSE)
  ))

  stan_file <- here::here("02_Script/stan/renewal_ceara_v2_1.stan")
  fit_renewal_v2_1_full_period <- rstan::stan(
    file = stan_file,
    data = stan_data_v2_1_full_period,
    iter = N_ITER,
    warmup = N_WARMUP,
    chains = N_CHAINS,
    seed = FIT_SEED,
    init = init_fn,
    control = list(adapt_delta = 0.97, max_treedepth = 13),
    refresh = max(1L, floor(N_ITER / 10L))
  )

  posterior_draws <- rstan::extract(
    fit_renewal_v2_1_full_period,
    pars = c("p_symp", "rho_sym_t", "overall_detection")
  )
  p_symp_posterior <- summarise_scalar_draws(posterior_draws$p_symp, "p_symp")
  symptomatic_reporting_posterior <- summarise_weekly_draws(
    posterior_draws$rho_sym_t, ceara_weekly$week_start
  )
  overall_detection_posterior <- summarise_weekly_draws(
    posterior_draws$overall_detection, ceara_weekly$week_start
  )
  posterior_correlations <- data.frame(
    week_start = ceara_weekly$week_start,
    p_symp_vs_symptomatic_reporting = vapply(
      seq_len(stan_data_v2_1_full_period$N),
      function(t) stats::cor(posterior_draws$p_symp, posterior_draws$rho_sym_t[, t]), numeric(1)
    ),
    p_symp_vs_overall_detection = vapply(
      seq_len(stan_data_v2_1_full_period$N),
      function(t) stats::cor(posterior_draws$p_symp, posterior_draws$overall_detection[, t]), numeric(1)
    ),
    symptomatic_reporting_vs_overall_detection = vapply(
      seq_len(stan_data_v2_1_full_period$N),
      function(t) stats::cor(posterior_draws$rho_sym_t[, t], posterior_draws$overall_detection[, t]), numeric(1)
    )
  )

  fit_summary <- summary(fit_renewal_v2_1_full_period)$summary
  key_parameters <- c(
    "sigma_R", "p_symp", "rho_sym_brazil", "rho_sym_ceara_mid",
    "reporting_ceara_offset", "reporting_trend", "sigma_geo", "phi_obs",
    "state_attack_window[1]", "state_attack_window[2]", "p_site[1]", "p_site[2]"
  )
  suffix <- if (nzchar(OUTPUT_TAG)) paste0("_", OUTPUT_TAG) else ""
  out_path <- here::here(
    "02_Script/stan", paste0("renewal_ceara_v2_1_full_period_fit", suffix, ".rds")
  )
  fit_bundle <- list(
    fit = fit_renewal_v2_1_full_period,
    stan_data = stan_data_v2_1_full_period,
    weekly_data = ceara_weekly,
    posterior_summaries = list(
      p_symp = p_symp_posterior,
      symptomatic_reporting = symptomatic_reporting_posterior,
      overall_detection = overall_detection_posterior,
      correlations = posterior_correlations
    ),
    config = list(
      stress_test = TRUE,
      purpose = "Diagnostic full-period extension; not the final full-history model",
      uf_code = UF_CODE,
      date_start = DATE_START,
      date_end = DATE_END,
      generation_interval_mean_weeks = GI_MEAN,
      generation_interval_sd_weeks = GI_SD,
      p_symp_prior = c(alpha = 30, beta = 28),
      reporting_region_sd = REPORTING_REGION_SD,
      serology_sites = SERO_SITE,
      serology_window_start = SERO_WINDOW_START,
      serology_window_end = SERO_WINDOW_END,
      serology_positive = SERO_POSITIVE,
      serology_n = SERO_N,
      closed_population = FALSE,
      population_process = "weekly observed population; positive changes enter susceptible population and negative changes remove people proportionally",
      iterations = N_ITER,
      warmup = N_WARMUP,
      chains = N_CHAINS,
      seed = FIT_SEED,
      stan_source = normalizePath(stan_file, winslash = "/", mustWork = TRUE)
    )
  )
  saveRDS(fit_bundle, out_path)
  message("[save] ", out_path)
  print(round(fit_summary[key_parameters, c("mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat")], 3))
  print(p_symp_posterior)
  message(sprintf(
    "[diagnostics] divergences=%d | max_treedepth hits=%d | min E-BFMI=%.3f",
    rstan::get_num_divergent(fit_renewal_v2_1_full_period),
    rstan::get_num_max_treedepth(fit_renewal_v2_1_full_period),
    min(rstan::get_bfmi(fit_renewal_v2_1_full_period))
  ))
}
