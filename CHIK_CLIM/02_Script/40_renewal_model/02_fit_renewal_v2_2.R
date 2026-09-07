# =============================================================================
# 02_fit_renewal_v2_2.R
#
# Fit the explicit-demography renewal model for Ceara, 2015-2019 (same short
# calibration window as the v2.1 baseline fit, so the comparison isolates
# the demographic-mechanism change described in
# docs/model_development_roadmap.md -- "Run B"). v2.1 sources and results
# are preserved; this is a new file, not a modification of v2.1.
#
# What changed versus v2.1
# -------------------------
# * State population is no longer a single N_pop_t series driving a
#   net-change rule. S[t] (susceptible) and U[t] (no longer susceptible)
#   are tracked separately and updated from real weekly births (SINASC),
#   all-cause deaths (IBGE 2024-revision annual total, day-weighted), and
#   the exact reconciliation residual N_end[t]-N_start[t]-births[t]+deaths[t]
#   -- see 00_data_prep/09_fetch_ce_sinasc_weekly_births.R and
#   00_data_prep/10_build_ceara_weekly_demography.R.
# * Deaths and the reconciliation residual are split across S/U in
#   proportion to the post-infection S/U composition (a documented
#   approximation -- see the .stan file header).
# * Transmission, reporting and serology-linkage structure, priors and
#   likelihood are UNCHANGED from v2.1 -- only the demographic mechanism
#   differs, by design.
#
# Environment overrides for a computational smoke test:
#   RENEWAL_V2_2_ITER, RENEWAL_V2_2_WARMUP, RENEWAL_V2_2_CHAINS,
#   RENEWAL_V2_2_TAG
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

summarise_scalar_draws <- function(draws, parameter) {
  quantiles <- stats::quantile(draws, probs = c(0.025, 0.5, 0.975))
  data.frame(
    parameter = parameter, mean = mean(draws), sd = stats::sd(draws),
    q2.5 = unname(quantiles[1]), median = unname(quantiles[2]), q97.5 = unname(quantiles[3])
  )
}

if (sys.nframe() == 0) {
  # Analysis window and renewal kernel -- identical to the v2.1 baseline.
  UF_CODE <- "23"
  DATE_START <- as.Date("2015-01-04")
  DATE_END <- as.Date("2019-12-29")
  G <- 8L
  SEED_WEEKS <- 8L
  GI_MEAN <- 2
  GI_SD <- 1

  REPORTING_REGION_SD <- 0.75
  SERO_SITE <- c("Juazeiro do Norte", "Quixada")
  SERO_WINDOW_START <- as.Date(c("2018-06-03", "2018-06-03"))
  SERO_WINDOW_END <- as.Date(c("2018-12-30", "2019-12-29"))
  SERO_POSITIVE <- c(103L, 289L)
  SERO_N <- c(404L, 409L)

  N_ITER <- env_integer("RENEWAL_V2_2_ITER", 2000L)
  N_WARMUP <- env_integer("RENEWAL_V2_2_WARMUP", 1000L)
  N_CHAINS <- env_integer("RENEWAL_V2_2_CHAINS", 4L)
  FIT_SEED <- 24052018L
  OUTPUT_TAG <- Sys.getenv("RENEWAL_V2_2_TAG", unset = "")

  if (N_WARMUP >= N_ITER) stop("RENEWAL_V2_2_WARMUP must be smaller than RENEWAL_V2_2_ITER")

  # Cases still come from the muni-week panel (unchanged from v2.1); only
  # the demographic inputs are swapped for the explicit weekly contract.
  panel <- readRDS(here::here("01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
  ceara_weekly_cases <- build_state_weekly(panel, UF_CODE) |>
    dplyr::filter(week_start >= DATE_START, week_start <= DATE_END)

  demography <- readRDS(here::here("01_Data/ceara_weekly_demography.rds")) |>
    dplyr::filter(week_start >= DATE_START, week_start <= DATE_END)

  ceara_weekly <- ceara_weekly_cases |>
    dplyr::select(-population) |>
    dplyr::inner_join(demography, by = "week_start")

  expected_dates <- seq(DATE_START, DATE_END, by = "week")
  if (!identical(ceara_weekly$week_start, expected_dates)) {
    stop("Ceara data do not form the expected complete, ordered weekly sequence")
  }
  if (nrow(ceara_weekly) != nrow(ceara_weekly_cases) || nrow(ceara_weekly) != nrow(demography)) {
    stop("Case panel and weekly demography table have mismatched week coverage")
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
  time_scaled <- seq(-0.5, 0.5, length.out = nrow(ceara_weekly))
  log_seed_prior_mean <- rep(log(10), SEED_WEEKS)

  stan_data_v2_2 <- list(
    N = nrow(ceara_weekly),
    G = G,
    seed_weeks = SEED_WEEKS,
    C = as.integer(ceara_weekly$cases),
    w = as.vector(w),
    N_start = as.vector(ceara_weekly$N_start),
    N_end = as.vector(ceara_weekly$N_end),
    births = as.vector(ceara_weekly$births),
    deaths = as.vector(ceara_weekly$all_cause_deaths),
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

  # Same initialization strategy as v2.1: a smoothed-case-based rough
  # latent-infection path rescaled to a 25% cumulative attack rate,
  # propagated through the SAME S/U demographic mechanism the model uses,
  # so the chain starts near a self-consistent state.
  smooth_cases <- vapply(seq_len(stan_data_v2_2$N), function(t) {
    idx <- max(1L, t - 2L):min(stan_data_v2_2$N, t + 2L)
    mean(stan_data_v2_2$C[idx] + 0.5)
  }, numeric(1))
  N_start <- stan_data_v2_2$N_start
  N_end <- stan_data_v2_2$N_end
  births <- stan_data_v2_2$births
  deaths <- stan_data_v2_2$deaths
  net_migration <- N_end - N_start - births + deaths
  rough_X <- smooth_cases * (0.25 * N_start[1] / sum(smooth_cases))
  rough_log_hazard <- numeric(stan_data_v2_2$N)
  rough_S <- N_start[1]
  rough_U <- 0
  for (t in seq_len(stan_data_v2_2$N)) {
    rough_X[t] <- min(rough_X[t], 0.2 * rough_S)
    rough_log_hazard[t] <- log(-log1p(-rough_X[t] / rough_S))
    s_after <- rough_S - rough_X[t]
    u_after <- rough_U + rough_X[t]
    frac_s <- s_after / (s_after + u_after)
    rough_S <- s_after + births[t] - deaths[t] * frac_s + net_migration[t] * frac_s
    rough_U <- u_after - deaths[t] * (1 - frac_s) + net_migration[t] * (1 - frac_s)
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
    "[data] %s to %s | N=%d | cases=%s | N_start range=%s to %s | mean weekly births=%.1f | mean weekly deaths=%.1f",
    DATE_START, DATE_END, nrow(ceara_weekly),
    format(sum(ceara_weekly$cases), big.mark = ","),
    format(min(N_start), big.mark = ",", scientific = FALSE),
    format(max(N_end), big.mark = ",", scientific = FALSE),
    mean(births), mean(deaths)
  ))
  message("[generation interval] ", paste(sprintf("%.4f", w), collapse = ", "))
  message(
    "[serology] ",
    paste(
      sprintf("%s: %d/%d (%s to %s)", SERO_SITE, SERO_POSITIVE, SERO_N, SERO_WINDOW_START, SERO_WINDOW_END),
      collapse = "; "
    )
  )

  stan_file <- here::here("02_Script/stan/renewal_ceara_v2_2.stan")
  fit_renewal_v2_2 <- rstan::stan(
    file = stan_file,
    data = stan_data_v2_2,
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
    "02_Script/stan", paste0("renewal_ceara_v2_2_fit", suffix, ".rds")
  )

  posterior_draws <- rstan::extract(fit_renewal_v2_2, pars = c("p_symp"))
  p_symp_posterior <- summarise_scalar_draws(posterior_draws$p_symp, "p_symp")

  fit_summary <- summary(fit_renewal_v2_2)$summary
  key_parameters <- c(
    "sigma_R", "p_symp", "rho_sym_brazil", "rho_sym_ceara_mid",
    "reporting_ceara_offset", "reporting_trend", "sigma_geo", "phi_obs",
    "state_attack_window[1]", "state_attack_window[2]", "p_site[1]", "p_site[2]"
  )

  fit_bundle <- list(
    fit = fit_renewal_v2_2,
    stan_data = stan_data_v2_2,
    weekly_data = ceara_weekly,
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
      demographic_mechanism = "explicit S/U accounting from real weekly births (SINASC), day-weighted annual deaths (IBGE), and exact reconciliation residual",
      demography_source = "01_Data/ceara_weekly_demography.rds",
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

  divergences <- rstan::get_num_divergent(fit_renewal_v2_2)
  treedepth_hits <- rstan::get_num_max_treedepth(fit_renewal_v2_2)
  min_bfmi <- min(rstan::get_bfmi(fit_renewal_v2_2))
  message(sprintf(
    "[diagnostics] divergences=%d | max_treedepth hits=%d | min E-BFMI=%.3f",
    divergences, treedepth_hits, min_bfmi
  ))
}
