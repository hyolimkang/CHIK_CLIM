# ---------------------------------------------------------------------------
# fit_chik_ceara_stan_weekly.R
#
# Weekly chikungunya transmission model for Brazilian states (UF).
# Aggregates municipality-week SINAN data to state-week, then fits Stan models.
#
# Workflow (run top to bottom; do not source chunks in isolation)
# ---------------------------------------------------------------------------
#   Step 0  Load packages
#   Step 1  USER CONFIG (UF_CODE, UF_NAME, years)
#   Step 2  UF-specific priors (UF_PRIORS)
#   Step 3  Build state-week panel  -> ce_weekly  (saved to 01_Data/)
#   Step 4  Filter to fit window    -> ce_fit, week_id (built once)
#   Step 5  Build Stan data lists   -> stan_data_v12e / v12c / v12d (separate)
#   Step 6  Fit Stan models         -> fit_v12e (main), fit_v12d (CE pulses only)
#   Step 7  Diagnostics and plots
#
# Stan data rule: each model has its own list. Nothing is overwritten.
#   - stan_data_v12e  <-> ceara_weekly_v12e.stan  (main fit for all UFs)
#   - stan_data_v12c  <-> ceara_weekly_v12.stan    (optional; not auto-fitted)
#   - stan_data_v12d  <-> ceara_weekly_v12d.stan  (CE pulse model only)
#
# Inputs
#   01_Data/chik_brazil_muni_week_2015_2024.rds
#   01_Data/ibge_pop_muni_year_*.csv
#   01_Data/{uf_name}_births_annual_ibge.csv
#
# Outputs
#   01_Data/{uf_name}_weekly_{YEAR_START}_{YEAR_END}.rds
#   fit_v12e, fit_v12d (in memory)
# ---------------------------------------------------------------------------
# ---- Step 0. Packages ------------------------------------------------------
cran_pkgs <- c(
  "here", "dplyr", "readr", "lubridate", "tibble", "tidyr",
  "ggplot2", "rstan", "posterior"
)
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages(
    "cmdstanr",
    repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
  )
}
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(lubridate)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(rstan)
  library(posterior)
  library(cmdstanr)
})

# ---- Step 1. USER CONFIG ---------------------------------------------------
UF_CODE    <- "31"    # CE=23, MG=31, SP=35, RJ=33, BA=29, ...
UF_NAME    <- "mg"    # lowercase prefix for 01_Data file names
YEAR_START <- 2014
YEAR_END   <- 2025
CASE_VAR   <- "cases_confirmed"
GAMMA_WEEK <- 1.0

# ---- Step 2. UF-specific priors --------------------------------------------
# foi_ext_mean/sdlog : serology-based annual FOI (NA -> wide defaults in Stan data)
# s0_alpha/beta      : initial susceptible fraction prior (Beta);
#                      high values reflect near-fully susceptible introduction era.
# rho_alpha/beta     : reporting probability prior (Beta)
# pulse_windows      : list(start, end) date pairs; NULL disables pulse model (v12d)
UF_PRIORS <- list(
  "23" = list(                          # Ceará
    foi_ext_mean   = 0.0120,
    foi_ext_sdlog  = 0.25,
    s0_alpha = 99, s0_beta  = 1,
    rho_alpha = 4, rho_beta = 12,
    pulse_windows = list(
      list(start = "2017-02-01", end = "2017-05-31"),
      list(start = "2022-02-01", end = "2022-05-31")
    )
  ),
  "31" = list(                          # Minas Gerais (serology TBD — uninformative)
    foi_ext_mean   = 0.012,
    foi_ext_sdlog  = 0.5,
    use_cum_attack_prior = TRUE,
    s0_alpha = 99, s0_beta  = 1,
    rho_alpha = 2, rho_beta = 8,
    pulse_windows = NULL
  ),
  "35" = list(                          # São Paulo (serology TBD — uninformative)
    foi_ext_mean   = NA,
    foi_ext_sdlog  = NA,
    s0_alpha = 99, s0_beta  = 1,
    rho_alpha = 2, rho_beta = 8,
    pulse_windows = NULL
  ),
  "33" = list(                          # Rio de Janeiro (serology TBD — uninformative)
    foi_ext_mean   = NA,
    foi_ext_sdlog  = NA,
    s0_alpha = 99, s0_beta  = 1,
    rho_alpha = 2, rho_beta = 8,
    pulse_windows = NULL
  ),
  "29" = list(                          # Bahia (serology TBD — uninformative)
    foi_ext_mean   = NA,
    foi_ext_sdlog  = NA,
    s0_alpha = 99, s0_beta  = 1,
    rho_alpha = 2, rho_beta = 8,
    pulse_windows = NULL
  )
)

uf_prior <- UF_PRIORS[[UF_CODE]]
if (is.null(uf_prior)) {
  stop(sprintf(
    "UF_CODE '%s' is not in UF_PRIORS. Add an entry before running.",
    UF_CODE
  ))
}

source(here::here("02_Script/10_sir_transmission_model/01_build_uf_weekly_panel.R"))

# ---- Step 3. State-week panel (muni -> UF aggregate) -----------------------
panel_week <- load_weekly_panel()
ce_weekly  <- build_uf_weekly_panel(
  panel_week = panel_week,
  uf_code    = UF_CODE,
  uf_name    = UF_NAME,
  year_start = YEAR_START,
  year_end   = YEAR_END,
  case_var   = CASE_VAR
)
ce_weekly  <- standardize_uf_weekly(ce_weekly)

saveRDS(ce_weekly, sprintf("01_Data/%s_weekly_%d_%d.rds", UF_NAME, YEAR_START, YEAR_END))

# ---- Step 4. Fit window + shared indices (built once, reused below) --------
FIT_START <- as.Date("2015-01-01")
FIT_END   <- as.Date("2025-12-31")

ce_fit <- ce_weekly |>
  dplyr::filter(week_start >= FIT_START, week_start <= FIT_END) |>
  dplyr::arrange(week_start) |>
  dplyr::mutate(t = dplyr::row_number())
ce_fit <- standardize_uf_weekly(ce_fit)

N <- nrow(ce_fit)
if (N < 2L) {
  stop("ce_fit is empty — check UF_CODE, births file, and FIT_START/FIT_END.", call. = FALSE)
}

years   <- sort(unique(ce_fit$year))
Y       <- length(years)
year_id <- match(ce_fit$year, years)
W       <- 52L
week_id <- make_week_id(ce_fit, W = W)

message(sprintf(
  "ce_fit: N = %d weeks, %d years, week_id length = %d",
  N, Y, length(week_id)
))

# ---- Step 5. Stan data lists (one list per model; no overwrites) ------------
STAN_CHAINS <- 2L
STAN_CORES  <- min(STAN_CHAINS, parallel::detectCores())
options(mc.cores = STAN_CORES)
rstan_options(auto_write = TRUE)

STAN_FILE_V12E_ORIGINAL <- here::here("02_Script/stan/ceara_weekly_v12e.stan")
# v12f = v12e + explicit demographic flow (births in, deaths out) so that the
# modelled S + I + R exactly tracks the observed population pop[t].
# Set USE_DEMOGRAPHIC_FLOW <- TRUE before sourcing to fit v12f instead of v12e.
# The Stan data contract is identical, so stan_data_v12e is reused unchanged.
STAN_FILE_V12F <- here::here("02_Script/stan/ceara_weekly_v12f.stan")
if (!exists("USE_DEMOGRAPHIC_FLOW")) USE_DEMOGRAPHIC_FLOW <- TRUE
STAN_FILE_V12E <- if (USE_DEMOGRAPHIC_FLOW) STAN_FILE_V12F else STAN_FILE_V12E_ORIGINAL
STAN_FILE_V12C <- here::here("02_Script/stan/ceara_weekly_v12.stan")
STAN_FILE_V12D <- here::here("02_Script/stan/ceara_weekly_v12d.stan")

# Rough natural mortality rate for v12f demographic flow (per capita per year).
# Default ~ Brazil crude death rate 6.6 / 1000 / year. v12e ignores this field.
if (!exists("BRAZIL_ANNUAL_DEATH_RATE")) BRAZIL_ANNUAL_DEATH_RATE <- 0.0066

# v12e: cumulative attack prior.
# Recovery setting: use the original infection-scale anchor.
# For UFs without FOI, this falls back to a moderate 20% cumulative attack prior.
use_cum_attack_prior <- TRUE
cum_attack_horizon_years = 11

if (!is.na(uf_prior$foi_ext_mean)) {
  cum_attack_ext_mean <- 1 - exp(-uf_prior$foi_ext_mean * cum_attack_horizon_years)
  cum_attack_ext_sdlogit <- 0.5
} else {
  cum_attack_ext_mean <- 0.20
  cum_attack_ext_sdlogit <- 0.50
}

stan_data_v12e <- list(
  N = N,
  cases = as.integer(ce_fit$cases),
  Y = Y,
  year_id = as.integer(year_id),
  W = W,
  week_id = week_id,
  pop = as.vector(ce_fit$population),
  births_weekly = as.vector(ce_fit$births_weekly),
  gamma_fixed = GAMMA_WEEK,
  mu_death_annual = BRAZIL_ANNUAL_DEATH_RATE,  # used by v12f only; v12e ignores it
  fit_start = 20L,
  
  s0_alpha = uf_prior$s0_alpha,
  s0_beta  = uf_prior$s0_beta,
  
  rho_alpha = uf_prior$rho_alpha,
  rho_beta  = uf_prior$rho_beta,
  
  use_cum_attack_prior = as.integer(use_cum_attack_prior),
  cum_attack_ext_mean = cum_attack_ext_mean,
  cum_attack_ext_sdlogit = cum_attack_ext_sdlogit,
  
  import_prior_logmedian = log(3),
  import_prior_sd = 0.5
)

check_stan_panel(stan_data_v12e, "stan_data_v12e")

# v12c: annual FOI prior (Ceará-style; optional second model)
stan_data_v12c <- list(
  N = N,
  cases = as.integer(ce_fit$cases),
  Y = Y,
  year_id = as.integer(year_id),
  W = W,
  week_id = week_id,
  pop = as.vector(ce_fit$population),
  births_weekly = as.vector(ce_fit$births_weekly),
  gamma_fixed = GAMMA_WEEK,
  fit_start = 20L,
  s0_alpha = uf_prior$s0_alpha,
  s0_beta  = uf_prior$s0_beta,
  rho_alpha = uf_prior$rho_alpha,
  rho_beta  = uf_prior$rho_beta,
  foi_ext_mean  = if (!is.na(uf_prior$foi_ext_mean)) uf_prior$foi_ext_mean else 0.01,
  foi_ext_sdlog = if (!is.na(uf_prior$foi_ext_sdlog)) uf_prior$foi_ext_sdlog else 1.0,
  import_prior_logmedian = log(3),
  import_prior_sd = 0.5
)
check_stan_panel(stan_data_v12c, "stan_data_v12c")

# v12d: pulse outbreaks (built only when pulse_windows is non-NULL for this UF)
pulse_id <- rep(0L, N)
if (!is.null(uf_prior$pulse_windows)) {
  for (k in seq_along(uf_prior$pulse_windows)) {
    pw <- uf_prior$pulse_windows[[k]]
    pulse_id[
      ce_fit$week_start >= as.Date(pw$start) &
        ce_fit$week_start <= as.Date(pw$end)
    ] <- as.integer(k)
  }
}
K_pulse <- if (is.null(uf_prior$pulse_windows)) 0L else length(uf_prior$pulse_windows)

if (K_pulse > 0L) {
  stan_data_v12d <- list(
    N = N,
    cases = as.integer(ce_fit$cases),
    Y = Y,
    year_id = as.integer(year_id),
    W = W,
    week_id = week_id,
    pop = as.vector(ce_fit$population),
    births_weekly = as.vector(ce_fit$births_weekly),
    gamma_fixed = GAMMA_WEEK,
    fit_start = 20L,
    s0_alpha = uf_prior$s0_alpha,
    s0_beta  = uf_prior$s0_beta,
    rho_alpha = uf_prior$rho_alpha,
    rho_beta  = uf_prior$rho_beta,
    foi_ext_mean  = if (!is.na(uf_prior$foi_ext_mean)) uf_prior$foi_ext_mean else 0.01,
    foi_ext_sdlog = if (!is.na(uf_prior$foi_ext_sdlog)) uf_prior$foi_ext_sdlog else 1.0,
    import_prior_logmedian = log(3),
    import_prior_sd = 0.5,
    K_pulse = K_pulse,
    pulse_id = as.integer(pulse_id),
    pulse_prior_sd = 0.6
  )
  check_stan_panel(stan_data_v12d, "stan_data_v12d")

  init_v12d <- list(
    list(
      rho = 0.15,
      mu_log_beta = -0.44,
      log_beta_season_raw = rep(0, W),
      sigma_season = 0.07,
      delta_year_raw = rep(0, Y),
      import_count_year = rep(3, Y),
      pulse_amp = rep(0.2, K_pulse),
      I0 = 5,
      S0_frac = 0.86,
      phi_cases = 2.35
    ),
    list(
      rho = 0.16,
      mu_log_beta = -0.42,
      log_beta_season_raw = rep(0, W),
      sigma_season = 0.07,
      delta_year_raw = rep(0, Y),
      import_count_year = rep(3, Y),
      pulse_amp = rep(0.3, K_pulse),
      I0 = 5,
      S0_frac = 0.86,
      phi_cases = 2.40
    )
  )
}

# ---- Step 6. Stan fits -----------------------------------------------------
# Main model: v12e + stan_data_v12e (never pass stan_data_v12c here)
fit_v12e <- stan(
  file    = STAN_FILE_V12F,
  data    = stan_data_v12e,
  chains  = STAN_CHAINS,
  cores   = STAN_CORES,
  iter    = 1000,
  warmup  = 500,
  seed    = 123,
  init = function() {
    list(
      rho = stan_data_v12e$rho_alpha /
        (stan_data_v12e$rho_alpha + stan_data_v12e$rho_beta),
      S0_frac = stan_data_v12e$s0_alpha /
        (stan_data_v12e$s0_alpha + stan_data_v12e$s0_beta),
      mu_log_beta = 0,
      log_beta_season_raw = rep(0, stan_data_v12e$W),
      sigma_season = 0.07,
      delta_year_raw = rep(0, stan_data_v12e$Y),
      import_count_year = rep(exp(stan_data_v12e$import_prior_logmedian), stan_data_v12e$Y),
      I0 = 5,
      phi_cases = 5
    )
  },
  control = list(
    adapt_delta = 0.97,
    max_treedepth = 13
  )
)

fit_v12c <- fit_v12e  # alias for legacy downstream code

if (K_pulse > 0L) {
  fit_v12d <- stan(
    file    = STAN_FILE_V12D,
    data    = stan_data_v12d,
    chains  = STAN_CHAINS,
    cores   = STAN_CORES,
    iter    = 1000,
    warmup  = 500,
    seed    = 123,
    init    = init_v12d,
    control = list(
      adapt_delta = 0.90,
      max_treedepth = 10
    ),
    refresh = 50
  )
} else {
  message("Skipping fit_v12d: no pulse_windows for UF ", UF_CODE)
}

# ---- Step 7. Diagnostics and posterior predictive plot ---------------------
print(fit_v12e, pars = c("mu_log_beta", "rho", "S0_frac", "I0", "phi_cases"))

check_hmc_diagnostics(fit_v12e)

print(
  fit_v12e,
  pars = c(
    "phi_cases", "mu_log_beta", "S0_frac", "I0",
    "import_count_year", "delta_year_raw"
  ),
  probs = c(0.025, 0.5, 0.975)
)
# Convert to long format
# Extract posterior draws
if ("CmdStanMCMC" %in% class(fit_v12e)) {
  draws_df <- fit_v12e$draws(format = "draws_df")
} else {
  draws_df <- as.data.frame(fit_v12e)
}

# Find posterior predictive columns
cases_rep_cols <- grep("^cases_rep\\[", names(draws_df), value = TRUE)


cases_rep_long <- draws_df %>%
  dplyr::select(all_of(cases_rep_cols)) %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = all_of(cases_rep_cols),
    names_to = "time",
    values_to = "cases_rep"
  ) %>%
  mutate(
    t = as.integer(gsub("cases_rep\\[|\\]", "", time))
  )

# Summarise posterior predictive intervals
cases_rep_summary <- cases_rep_long %>%
  group_by(t) %>%
  summarise(
    pred_med = median(cases_rep),
    pred_lwr = quantile(cases_rep, 0.025),
    pred_upr = quantile(cases_rep, 0.975),
    pred_lwr50 = quantile(cases_rep, 0.25),
    pred_upr50 = quantile(cases_rep, 0.75),
    .groups = "drop"
  )

# Merge with observed data
plot_df <- ce_fit %>%
  mutate(t = row_number()) %>%
  left_join(cases_rep_summary, by = "t")

# Plot
ggplot(plot_df, aes(x = t)) +
  geom_ribbon(
    aes(ymin = pred_lwr, ymax = pred_upr),
    fill = "skyblue",
    alpha = 0.18
  ) +
  geom_ribbon(
    aes(ymin = pred_lwr50, ymax = pred_upr50),
    fill = "skyblue",
    alpha = 0.35
  ) +
  geom_line(
    aes(y = pred_med, colour = "Posterior predictive median"),
    linewidth = 0.8
  ) +
  geom_point(
    aes(y = cases, colour = "Observed cases"),
    size = 0.7,
    alpha = 0.75
  ) +
  scale_colour_manual(
    values = c(
      "Observed cases" = "black",
      "Posterior predictive median" = "red"
    )
  ) +
  labs(
    x = "Week",
    y = "Reported weekly cases",
    colour = NULL,
    title = "Observed vs posterior predictive fitted cases",
    subtitle = "Blue ribbons = 50% and 95% posterior predictive intervals"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  ) +
  scale_y_log10()


