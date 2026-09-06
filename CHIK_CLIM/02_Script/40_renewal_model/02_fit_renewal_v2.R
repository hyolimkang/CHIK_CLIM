# =============================================================================
# 02_fit_renewal_v2.R
#
# Fit the dynamic susceptible-depletion renewal model for Ceara, 2015-2018.
# This is a new scientific model version; v1 sources and results are preserved.
#
# Main assumptions
# ----------------
# * X[t] is true weekly incident infection count.
# * Closed 2015 population; infection-derived immunity does not wane.
# * R0[t] follows a smooth first-order Gaussian random walk.
# * p_symp is externally fixed (primary analysis: 0.52).
# * Symptomatic reporting changes monotonically on the logit scale.
# * Beta(20, 60) describes Brazil-wide long-run symptomatic reporting;
#   Ceara can differ through a region-level log-odds offset.
# * The 2018 Juazeiro survey (103/404) is linked to Ceara attack rate with a
#   fixed geographic heterogeneity SD. Its midpoint is represented by the
#   week ending 2018-09-30.
# * Assay sensitivity and specificity default to 1 until validated values are
#   supplied from the assay documentation.
#
# Environment overrides for a computational smoke test:
#   RENEWAL_V2_ITER, RENEWAL_V2_WARMUP, RENEWAL_V2_CHAINS,
#   RENEWAL_V2_TAG, RENEWAL_V2_P_SYMP
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

source(here::here("02_Script/40_renewal_model/01_build_ceara_state_weekly.R"))

discretize_gamma_generation_interval <- function(mean_weeks, sd_weeks, G) {
  shape <- (mean_weeks / sd_weeks)^2
  rate <- mean_weeks / sd_weeks^2
  interval_mass <- diff(stats::pgamma(0:G, shape = shape, rate = rate))
  interval_mass / sum(interval_mass)
}

env_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1) stop(name, " must be a positive integer")
  parsed
}

env_numeric <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.numeric(value))
  if (!is.finite(parsed)) stop(name, " must be numeric")
  parsed
}

if (sys.nframe() == 0) {
  # Analysis window and renewal kernel.
  UF_CODE <- "23"
  DATE_START <- as.Date("2015-01-04")
  DATE_END <- as.Date("2018-12-30")
  G <- 8L
  SEED_WEEKS <- 8L
  GI_MEAN <- 2
  GI_SD <- 1

  # External and transportability assumptions.
  P_SYMP <- env_numeric("RENEWAL_V2_P_SYMP", 0.52)
  REPORTING_REGION_SD <- 0.75
  SERO_DATE <- as.Date("2018-09-30")
  SERO_POSITIVE <- 103L
  SERO_TESTED <- 404L
  SERO_SENSITIVITY <- 1
  SERO_SPECIFICITY <- 1
  SERO_GEOGRAPHIC_SD <- 0.35

  # Sampling defaults; environment variables support non-destructive tests.
  N_ITER <- env_integer("RENEWAL_V2_ITER", 2000L)
  N_WARMUP <- env_integer("RENEWAL_V2_WARMUP", 1000L)
  N_CHAINS <- env_integer("RENEWAL_V2_CHAINS", 4L)
  FIT_SEED <- 24052018L
  OUTPUT_TAG <- Sys.getenv("RENEWAL_V2_TAG", unset = "")

  if (N_WARMUP >= N_ITER) stop("RENEWAL_V2_WARMUP must be smaller than RENEWAL_V2_ITER")
  if (P_SYMP <= 0 || P_SYMP >= 1) stop("RENEWAL_V2_P_SYMP must lie strictly between 0 and 1")

  panel <- readRDS(here::here("01_Data/chik_dlnm_panel_muni_week_2015_2024.rds"))
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

  sero_t <- match(SERO_DATE, ceara_weekly$week_start)
  if (is.na(sero_t)) stop("SERO_DATE is absent from the weekly model grid")

  w <- discretize_gamma_generation_interval(GI_MEAN, GI_SD, G)
  N_pop <- as.numeric(ceara_weekly$population[1])
  time_scaled <- seq(-0.5, 0.5, length.out = nrow(ceara_weekly))

  # With zero initial reports, a weak log-scale prior centered on 10 true
  # infections/week keeps the renewal recursion alive without fixing scale.
  log_seed_prior_mean <- rep(log(10), SEED_WEEKS)

  stan_data_v2 <- list(
    N = nrow(ceara_weekly),
    G = G,
    seed_weeks = SEED_WEEKS,
    C = as.integer(ceara_weekly$cases),
    w = as.vector(w),
    N_pop = N_pop,
    time_scaled = as.vector(time_scaled),
    p_symp = P_SYMP,
    log_seed_prior_mean = log_seed_prior_mean,
    seed_prior_sd = 1.5,
    reporting_region_sd = REPORTING_REGION_SD,
    sero_t = sero_t,
    sero_positive = SERO_POSITIVE,
    sero_tested = SERO_TESTED,
    sero_sensitivity = SERO_SENSITIVITY,
    sero_specificity = SERO_SPECIFICITY,
    sero_geographic_sd = SERO_GEOGRAPHIC_SD
  )

  # Initialize the latent path from a five-week moving average rescaled to a
  # 25% cumulative attack rate. This is only a starting point for adaptation;
  # the posterior scale is inferred jointly from cases and serology.
  smooth_cases <- vapply(seq_len(stan_data_v2$N), function(t) {
    idx <- max(1L, t - 2L):min(stan_data_v2$N, t + 2L)
    mean(stan_data_v2$C[idx] + 0.5)
  }, numeric(1))
  rough_X <- smooth_cases * (0.25 * N_pop / sum(smooth_cases))
  rough_log_hazard <- numeric(stan_data_v2$N)
  rough_susceptible <- N_pop
  for (t in seq_len(stan_data_v2$N)) {
    rough_X[t] <- min(rough_X[t], 0.2 * rough_susceptible)
    rough_log_hazard[t] <- log(-log1p(-rough_X[t] / rough_susceptible))
    rough_susceptible <- rough_susceptible - rough_X[t]
  }

  init_fn <- function() {
    list(
      sigma_R = 0.05,
      log_infection_scale = rough_log_hazard[1],
      log_hazard_relative = rough_log_hazard[-1] - rough_log_hazard[1],
      rho_sym_brazil = 0.25,
      logit_rho_sym_ceara_mid = stats::qlogis(0.116),
      reporting_trend = 0.3,
      logit_site_attack_at_serosurvey = stats::qlogis(SERO_POSITIVE / SERO_TESTED),
      phi_obs = 20
    )
  }

  message(sprintf(
    "[data] %s to %s | N=%d | cases=%s | population=%s",
    DATE_START, DATE_END, nrow(ceara_weekly),
    format(sum(ceara_weekly$cases), big.mark = ","),
    format(N_pop, big.mark = ",", scientific = FALSE)
  ))
  message("[generation interval] ", paste(sprintf("%.4f", w), collapse = ", "))
  message(sprintf(
    "[serology] %d/%d at %s | geographic log-odds SD=%.2f",
    SERO_POSITIVE, SERO_TESTED, SERO_DATE, SERO_GEOGRAPHIC_SD
  ))

  stan_file <- here::here("02_Script/stan/renewal_ceara_v2.stan")
  fit_renewal_v2 <- rstan::stan(
    file = stan_file,
    data = stan_data_v2,
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
    "02_Script/stan", paste0("renewal_ceara_v2_fit", suffix, ".rds")
  )
  fit_bundle <- list(
    fit = fit_renewal_v2,
    stan_data = stan_data_v2,
    weekly_data = ceara_weekly,
    config = list(
      uf_code = UF_CODE,
      date_start = DATE_START,
      date_end = DATE_END,
      generation_interval_mean_weeks = GI_MEAN,
      generation_interval_sd_weeks = GI_SD,
      p_symp = P_SYMP,
      reporting_region_sd = REPORTING_REGION_SD,
      sero_date = SERO_DATE,
      sero_geographic_sd = SERO_GEOGRAPHIC_SD,
      closed_population = TRUE,
      iterations = N_ITER,
      warmup = N_WARMUP,
      chains = N_CHAINS,
      seed = FIT_SEED,
      stan_source = normalizePath(stan_file, winslash = "/", mustWork = TRUE)
    )
  )
  saveRDS(fit_bundle, out_path)
  message("[save] ", out_path)

  fit_summary <- summary(fit_renewal_v2)$summary
  key_parameters <- c(
    "sigma_R", "rho_sym_brazil", "rho_sym_ceara_mid",
    "reporting_ceara_offset", "reporting_trend",
    "sero_geographic_offset", "phi_obs",
    "state_attack_at_serosurvey", "site_attack_at_serosurvey"
  )
  print(round(fit_summary[key_parameters, c("mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat")], 3))

  divergences <- rstan::get_num_divergent(fit_renewal_v2)
  treedepth_hits <- rstan::get_num_max_treedepth(fit_renewal_v2)
  min_bfmi <- min(rstan::get_bfmi(fit_renewal_v2))
  message(sprintf(
    "[diagnostics] divergences=%d | max_treedepth hits=%d | min E-BFMI=%.3f",
    divergences, treedepth_hits, min_bfmi
  ))
}
