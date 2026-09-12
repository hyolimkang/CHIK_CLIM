# =============================================================================
# 02_fit_renewal_v2_1.R
#
# Fit the dynamic susceptible-depletion renewal model for Ceara, 2015-2019.
# Version 2.1 uses the state-weekly population series to update susceptibility
# between weeks; v2 sources and results are preserved.
#
# Main assumptions
# ----------------
# * X[t] is true weekly incident infection count.
# * State population can change weekly: positive net changes enter susceptible
#   population and negative changes are removed proportionally.
# * Infection-derived immunity does not wane.
# * R0[t] follows a smooth first-order Gaussian random walk.
# * p_symp is estimated with a Beta(30, 28) prior.
# * Symptomatic reporting changes monotonically on the logit scale.
# * Beta(20, 60) describes Brazil-wide long-run symptomatic reporting;
#   Ceara can differ through a region-level log-odds offset.
# * The Juazeiro do Norte and Quixada surveys are linked to the mean Ceara
#   cumulative proportion over their respective collection windows.
# * Geographic heterogeneity is a non-centered site-level log-odds offset with
#   a half-normal(0, 1) SD; assay performance is not modelled separately.
#
# Environment overrides for a computational smoke test:
#   RENEWAL_V2_1_ITER, RENEWAL_V2_1_WARMUP, RENEWAL_V2_1_CHAINS,
#   RENEWAL_V2_1_TAG
# =============================================================================

required_packages <- c("here", "rstan", "posterior", "dplyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(dplyr)
})

rstan_options(auto_write = TRUE)
options(mc.cores = min(4, parallel::detectCores()))

source(here::here("02_Script/40_renewal_model/00_shared/01_build_ceara_state_weekly.R"))

discretize_gamma_generation_interval <- function(mean_weeks, sd_weeks, G) {
  shape <- (mean_weeks / sd_weeks)^2
  rate <- mean_weeks / sd_weeks^2
  interval_mass <- diff(stats::pgamma(0:G, shape = shape, rate = rate))
  interval_mass / sum(interval_mass)
}

env_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1) stop(name, " must be a positive integer")
  parsed
}

summarise_scalar_draws <- function(draws, parameter) {
  quantiles <- stats::quantile(draws, probs = c(0.025, 0.5, 0.975))
  data.frame(
    parameter = parameter,
    mean = mean(draws),
    sd = stats::sd(draws),
    q2.5 = unname(quantiles[1]),
    median = unname(quantiles[2]),
    q97.5 = unname(quantiles[3])
  )
}

summarise_weekly_draws <- function(draw_matrix, week_start) {
  quantiles <- t(vapply(seq_len(ncol(draw_matrix)), function(t) {
    stats::quantile(draw_matrix[, t], probs = c(0.025, 0.5, 0.975))
  }, numeric(3)))
  data.frame(
    week_start = week_start,
    mean = colMeans(draw_matrix),
    sd = apply(draw_matrix, 2, stats::sd),
    q2.5 = quantiles[, 1],
    median = quantiles[, 2],
    q97.5 = quantiles[, 3]
  )
}

if (sys.nframe() == 0) {
  # Analysis window and renewal kernel.
  UF_CODE <- "23"
  DATE_START <- as.Date("2015-01-04")
  DATE_END <- as.Date("2019-12-29")
  G <- 8L
  SEED_WEEKS <- 8L
  GI_MEAN <- 2
  GI_SD <- 1

  # External and transportability assumptions.
  REPORTING_REGION_SD <- 0.75
  # Weekly model dates are Sundays. These endpoints represent all full model
  # weeks beginning within the reported calendar collection windows.
  SERO_SITE <- c("Juazeiro do Norte", "Quixada")
  SERO_WINDOW_START <- as.Date(c("2018-06-03", "2018-06-03"))
  SERO_WINDOW_END <- as.Date(c("2018-12-30", "2019-12-29"))
  SERO_POSITIVE <- c(103L, 289L)
  SERO_N <- c(404L, 409L)

  # Sampling defaults; environment variables support non-destructive tests.
  N_ITER <- env_integer("RENEWAL_V2_1_ITER", 2000L)
  N_WARMUP <- env_integer("RENEWAL_V2_1_WARMUP", 1000L)
  N_CHAINS <- env_integer("RENEWAL_V2_1_CHAINS", 4L)
  FIT_SEED <- 24052018L
  OUTPUT_TAG <- Sys.getenv("RENEWAL_V2_1_TAG", unset = "")

  if (N_WARMUP >= N_ITER) stop("RENEWAL_V2_1_WARMUP must be smaller than RENEWAL_V2_1_ITER")

  panel <- readRDS(here::here("01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
  ceara_weekly <- build_state_weekly(panel, UF_CODE) |>
    dplyr::filter(week_start >= DATE_START, week_start <= DATE_END)

  expected_dates <- seq(DATE_START, DATE_END, by = "week")
  if (!identical(ceara_weekly$week_start, expected_dates)) {
    stop("Ceara data do not form the expected complete, ordered weekly sequence")
  }
  if (anyNA(ceara_weekly$cases) || any(ceara_weekly$cases < 0)) {
    stop("Weekly case counts must be complete and non-negative")
  }
  if (any(abs(ceara_weekly$cases - round(ceara_weekly$cases)) > 1e-8)) {
    stop("Weekly case counts must be integer-valued")
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

  w <- discretize_gamma_generation_interval(GI_MEAN, GI_SD, G)
  N_pop_t <- as.numeric(ceara_weekly$population)
  if (length(N_pop_t) != nrow(ceara_weekly) ||
    any(!is.finite(N_pop_t)) || any(N_pop_t < 1)) {
    stop("Weekly state population must be finite and at least one in every week")
  }
  time_scaled <- seq(-0.5, 0.5, length.out = nrow(ceara_weekly))

  # With zero initial reports, a weak log-scale prior centered on 10 true
  # infections/week keeps the renewal recursion alive without fixing scale.
  log_seed_prior_mean <- rep(log(10), SEED_WEEKS)

  stan_data_v2_1 <- list(
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

  # Initialize the latent path from a five-week moving average rescaled to a
  # 25% cumulative attack rate. This is only a starting point for adaptation;
  # the posterior scale is inferred jointly from cases and serology.
  smooth_cases <- vapply(seq_len(stan_data_v2_1$N), function(t) {
    idx <- max(1L, t - 2L):min(stan_data_v2_1$N, t + 2L)
    mean(stan_data_v2_1$C[idx] + 0.5)
  }, numeric(1))
  rough_X <- smooth_cases * (0.25 * N_pop_t[1] / sum(smooth_cases))
  rough_log_hazard <- numeric(stan_data_v2_1$N)
  rough_susceptible <- N_pop_t[1]
  for (t in seq_len(stan_data_v2_1$N)) {
    rough_X[t] <- min(rough_X[t], 0.2 * rough_susceptible)
    rough_log_hazard[t] <- log(-log1p(-rough_X[t] / rough_susceptible))
    rough_susceptible <- rough_susceptible - rough_X[t]
    if (t < stan_data_v2_1$N) {
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
      rho_sym_brazil = 0.25,
      logit_rho_sym_ceara_mid = stats::qlogis(0.116),
      reporting_trend = 0.3,
      p_symp = 0.52,
      sigma_geo = 0.5,
      z_site = rep(0, length(SERO_SITE)),
      phi_obs = 20
    )
  }

  message(sprintf(
    "[data] %s to %s | N=%d | cases=%s | population range=%s to %s",
    DATE_START, DATE_END, nrow(ceara_weekly),
    format(sum(ceara_weekly$cases), big.mark = ","),
    format(min(N_pop_t), big.mark = ",", scientific = FALSE),
    format(max(N_pop_t), big.mark = ",", scientific = FALSE)
  ))
  message("[generation interval] ", paste(sprintf("%.4f", w), collapse = ", "))
  message(
    "[serology] ",
    paste(
      sprintf(
        "%s: %d/%d (%s to %s)",
        SERO_SITE, SERO_POSITIVE, SERO_N,
        SERO_WINDOW_START, SERO_WINDOW_END
      ),
      collapse = "; "
    )
  )

  stan_file <- here::here("02_Script/stan/renewal_ceara_v2_1.stan")
  fit_renewal_v2_1 <- rstan::stan(
    file = stan_file,
    data = stan_data_v2_1,
    iter = N_ITER,
    warmup = N_WARMUP,
    chains = N_CHAINS,
    seed = FIT_SEED,
    init = init_fn,
    control = list(adapt_delta = 0.97, max_treedepth = 13),
    refresh = max(1L, floor(N_ITER / 10L))
  )

  suffix <- if (nzchar(OUTPUT_TAG)) paste0("_", OUTPUT_TAG) else ""
  out_path <- here::here(
    "02_Script/stan", paste0("renewal_ceara_v2_1_fit", suffix, ".rds")
  )

  posterior_draws <- rstan::extract(
    fit_renewal_v2_1,
    pars = c("p_symp", "rho_sym_t", "overall_detection")
  )
  p_symp_posterior <- summarise_scalar_draws(
    posterior_draws$p_symp,
    "p_symp"
  )
  symptomatic_reporting_posterior <- summarise_weekly_draws(
    posterior_draws$rho_sym_t,
    ceara_weekly$week_start
  )
  overall_detection_posterior <- summarise_weekly_draws(
    posterior_draws$overall_detection,
    ceara_weekly$week_start
  )
  posterior_correlations <- data.frame(
    week_start = ceara_weekly$week_start,
    p_symp_vs_symptomatic_reporting = vapply(
      seq_len(stan_data_v2_1$N),
      function(t) stats::cor(posterior_draws$p_symp, posterior_draws$rho_sym_t[, t]),
      numeric(1)
    ),
    p_symp_vs_overall_detection = vapply(
      seq_len(stan_data_v2_1$N),
      function(t) stats::cor(posterior_draws$p_symp, posterior_draws$overall_detection[, t]),
      numeric(1)
    ),
    symptomatic_reporting_vs_overall_detection = vapply(
      seq_len(stan_data_v2_1$N),
      function(t) {
        stats::cor(
          posterior_draws$rho_sym_t[, t],
          posterior_draws$overall_detection[, t]
        )
      },
      numeric(1)
    )
  )

  fit_summary <- summary(fit_renewal_v2_1)$summary
  key_parameters <- c(
    "sigma_R", "p_symp", "rho_sym_brazil", "rho_sym_ceara_mid",
    "reporting_ceara_offset", "reporting_trend",
    "sigma_geo", "phi_obs",
    "state_attack_window[1]", "state_attack_window[2]",
    "p_site[1]", "p_site[2]"
  )

  fit_bundle <- list(
    fit = fit_renewal_v2_1,
    stan_data = stan_data_v2_1,
    weekly_data = ceara_weekly,
    posterior_summaries = list(
      p_symp = p_symp_posterior,
      symptomatic_reporting = symptomatic_reporting_posterior,
      overall_detection = overall_detection_posterior,
      correlations = posterior_correlations
    ),
    config = list(
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

  divergences <- rstan::get_num_divergent(fit_renewal_v2_1)
  treedepth_hits <- rstan::get_num_max_treedepth(fit_renewal_v2_1)
  min_bfmi <- min(rstan::get_bfmi(fit_renewal_v2_1))
  message(sprintf(
    "[diagnostics] divergences=%d | max_treedepth hits=%d | min E-BFMI=%.3f",
    divergences, treedepth_hits, min_bfmi
  ))
}
