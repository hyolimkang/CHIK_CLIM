# ===========================================================================
# 02_fit_renewal_v1.R
#
# Purpose
# -------
# Fit Stage 1 of the renewal-equation transmission model (stan/renewal_ceara_v1.stan)
# to Ceará's state-weekly case series: constant R, no susceptible depletion,
# lognormal process noise on latent infections, NegBin observation model.
# This stage only checks that the renewal machinery (generation interval
# convolution + seeding + process noise + NegBin observation) fits cleanly
# before adding time-varying R, susceptible depletion, or a climate link.
#
# Generation interval (w_tau)
# ----------------------------
# Chikungunya's serial interval (extrinsic incubation in the mosquito +
# human incubation/infectious period) is commonly cited in the literature
# around 10-20 days; a discretized Gamma with mean 2 weeks and SD 1 week is
# used here as a starting assumption (NOT estimated from data at this
# stage) — discretized via the standard renewal-model convention
# w_tau = pgamma(tau) - pgamma(tau-1), truncated at G weeks and renormalized
# to sum to 1. Revisit this once the basic model fits cleanly.
#
# Why a single-outbreak window, not the full 2014-2024 series
# --------------------------------------------------------------
# First attempt fit the full 523-week series and it was numerically
# unstable (neg_binomial_2_rng exceptions during warmup from I_t exploding
# to astronomical values). This isn't a bug — it's the expected consequence
# of a *constant* R with *no* susceptible depletion applied over a decade
# containing multiple outbreak/quiet cycles: any warmup draw with R > 1
# grows the renewal recursion exponentially for the rest of the series with
# nothing to stop it (in reality, susceptible depletion is exactly what
# stops this, but that's deferred to Stage 3). So Stage 1 is fit to one
# single-outbreak window (2015-06 to 2018-06, the 2016-17 wave) where a
# constant R is at least a plausible local approximation. Stage 2 (R_t as a
# random walk) is what gets applied to the full multi-year series, since
# real dynamics need time-varying R over that horizon regardless of
# susceptible depletion.
#
# Output
# ------
#   02_Script/stan/renewal_ceara_v1.rds — the fitted stanfit object
# ===========================================================================

for (p in c("here", "rstan", "posterior", "dplyr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr) })

rstan_options(auto_write = TRUE)
options(mc.cores = min(4, parallel::detectCores()))

source(here::here("02_Script/40_renewal_model/01_build_ceara_state_weekly.R"))

discretize_gamma_generation_interval <- function(mean_weeks, sd_weeks, G) {
  shape <- (mean_weeks / sd_weeks)^2
  rate  <- mean_weeks / sd_weeks^2
  cdf <- pgamma(0:G, shape = shape, rate = rate)
  w <- diff(cdf)
  w / sum(w)
}

if (sys.nframe() == 0) {
  # ============================================================
  UF_CODE     <- "23"
  DATE_START  <- as.Date("2015-06-01")  # single-outbreak window (2016-17 wave)
  DATE_END    <- as.Date("2018-06-30")  # see header comment for why not the full series
  G           <- 8    # max generation-interval lag, weeks
  SEED_WEEKS  <- 8    # must be >= G
  GI_MEAN     <- 2    # weeks
  GI_SD       <- 1    # weeks
  N_ITER      <- 2000
  N_WARMUP    <- 1000
  N_CHAINS    <- 4
  # ============================================================

  panel <- readRDS(here::here("01_Data/chik_dlnm_panel_muni_week_2015_2024.rds"))
  ceara_weekly <- build_state_weekly(panel, UF_CODE) |>
    dplyr::filter(week_start >= DATE_START, week_start <= DATE_END)

  w <- discretize_gamma_generation_interval(GI_MEAN, GI_SD, G)
  message("[gen interval] weights: ", paste(round(w, 3), collapse = ", "))

  stan_data <- list(
    N = nrow(ceara_weekly),
    G = G,
    seed_weeks = SEED_WEEKS,
    C = ceara_weekly$cases,
    w = w
  )

  # Explicit near-stable inits: rstan's default random init (uniform on the
  # unconstrained scale) can start R noticeably above 1 with no depletion
  # to check it, which explodes the multiplicative renewal recursion within
  # the first few warmup iterations (see header comment). Starting R near 1
  # and process noise small avoids that.
  log_I_seed_init <- log(mean(stan_data$C[1:SEED_WEEKS]) + 1)
  init_fn <- function() {
    list(
      R = 1,
      log_I_seed = rep(log_I_seed_init, SEED_WEEKS),
      sigma_process = 0.1,
      phi_obs = 5,
      eta = rep(0, stan_data$N - SEED_WEEKS)
    )
  }

  stan_file <- here::here("02_Script/stan/renewal_ceara_v1.stan")
  fit_renewal_v1 <- stan(
    file = stan_file,
    data = stan_data,
    iter = N_ITER,
    warmup = N_WARMUP,
    chains = N_CHAINS,
    seed = 42,
    init = init_fn,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )

  out_path <- here::here("02_Script/stan/renewal_ceara_v1.rds")
  saveRDS(fit_renewal_v1, out_path)
  message("[save] ", out_path)

  # ---- quick diagnostics ----
  summ <- posterior::summarise_draws(fit_renewal_v1, "rhat", "ess_bulk")
  key_pars <- summ |> dplyr::filter(variable %in% c("R", "sigma_process", "phi_obs"))
  print(key_pars)

  divergences <- rstan::get_num_divergent(fit_renewal_v1)
  max_treedepth_hits <- rstan::get_num_max_treedepth(fit_renewal_v1)
  message(sprintf("[diagnostics] divergences: %d | max_treedepth hits: %d", divergences, max_treedepth_hits))
}
