# STEP 2 -- Reusable future-renewal simulation functions.
#
# Pure functions, no file I/O, no global state beyond what is passed in.
# Mirrors (does not modify) the accepted v4.9 climate-forced Stan model's
# mechanics: same generation-interval convolution, same hazard form
# 1-exp(-FOI), same demographic S/U bookkeeping identity, same treatment
# of imports as a hazard addition (never double-counted as cases).
#
# Source this file AFTER 00_config_future_pilot.R (needs generation_weights_future(),
# WEEKS_PER_YEAR).

suppressPackageStartupMessages({ library(dplyr) })

# ================================================================
# 2.2 Future seasonal regressors -- EXACT continuation of the historical
# phase (t_index continues 574, 575, ... after the historical 1..573;
# verified against stan_data$seasonal_sin/cos in the unit tests).
# ================================================================
future_seasonal_regressors <- function(t_index) {
  list(
    sin1 = sin(2 * pi * t_index / WEEKS_PER_YEAR), cos1 = cos(2 * pi * t_index / WEEKS_PER_YEAR),
    sin2 = sin(4 * pi * t_index / WEEKS_PER_YEAR), cos2 = cos(4 * pi * t_index / WEEKS_PER_YEAR)
  )
}

# ================================================================
# 2.2 Future year effects -- Y0 (reference/neutral) vs Y1 (resampled from
# the model's OWN documented distributional assumptions: z ~ std_normal(),
# year_effect = year_effect_prior_sd*(z - mean(z)), SAME construction as
# the historical fit, applied within the future window only). The same
# treatment is applied to the year-varying seasonal AMPLITUDE hierarchy
# (A_year = exp(sigma_season_year * z_season_year)), since it uses an
# identical per-year std_normal() hierarchy in the accepted model and
# leaving it silently fixed while only forecasting the level would be an
# inconsistent half-measure. This is a declared, disclosed choice.
# ================================================================
sample_future_year_structure <- function(arm, n_future_years, year_effect_prior_sd, sigma_season_year_draw, rng) {
  stopifnot(arm %in% YEAR_EFFECT_ARMS)
  if (arm == "Y0_reference") {
    return(list(year_effect = rep(0, n_future_years), A_year = rep(1, n_future_years)))
  }
  z_year <- rng(n_future_years) # std_normal draws, caller supplies the RNG (reproducible per posterior draw/year_effect_id)
  z_season <- rng(n_future_years)
  list(
    year_effect = year_effect_prior_sd * (z_year - mean(z_year)),
    A_year = exp(sigma_season_year_draw * z_season)
  )
}

# ================================================================
# 2.3 Future climate arms.
#   C0: stationary/neutral (z_T = z_P = 0 for all future weeks).
#   C1: historical weather variability resampled in COMPLETE, seasonally
#       aligned calendar-year blocks (whole-year bootstrap; never shuffles
#       individual weeks, preserving the joint T-P temporal structure
#       since both series come from the same original historical week).
#   C2: the SAME resampled sequences as C1 (same block draw), but the
#       caller must pass beta_T_eff = beta_P_eff = 0 into compute_R0_future()
#       -- implemented there, not here, so this function only builds the
#       (possibly reused) z_T/z_P sequence.
# `hist_by_year` = list, one element per historical calendar year, each a
# data.frame(within_year_index, z_T_anom, z_P_anom) ordered by week.
# ================================================================
sample_future_climate_block <- function(arm, future_year_labels, weeks_per_future_year, hist_by_year, rng_discrete) {
  stopifnot(arm %in% CLIMATE_ARMS)
  n_future_years <- length(future_year_labels)
  if (arm == "C0_stationary_neutral") {
    total_weeks <- sum(weeks_per_future_year)
    return(list(z_T = rep(0, total_weeks), z_P = rep(0, total_weeks), source_years = rep(NA_integer_, n_future_years)))
  }
  hist_years <- as.integer(names(hist_by_year))
  chosen <- rng_discrete(n_future_years, hist_years) # which historical calendar year each future year borrows
  z_T <- vector("list", n_future_years); z_P <- vector("list", n_future_years)
  for (i in seq_len(n_future_years)) {
    block <- hist_by_year[[as.character(chosen[i])]]
    need <- weeks_per_future_year[i]
    have <- nrow(block)
    idx <- if (have >= need) seq_len(need) else c(seq_len(have), rep(have, need - have)) # pad by repeating the last week if a 52-week future year draws a 53-week historical block's complement (rare, disclosed)
    z_T[[i]] <- block$z_T_anom[idx]; z_P[[i]] <- block$z_P_anom[idx]
  }
  list(z_T = unlist(z_T), z_P = unlist(z_P), source_years = chosen)
}

# ================================================================
# log_R0_future(t) = alpha_R + future_year_effect(t) + future_seasonality(t)
#                   + beta_T_eff*z_T_future(t) + beta_P_eff*z_P_future(t)
# `beta_T_eff`/`beta_P_eff` are passed separately from the posterior draw's
# own beta_T/beta_P so that arm C2 can force them to exactly 0 while still
# using genuinely resampled (non-zero) climate input -- the forcing-
# identity unit test (2.10.F) checks this equals arm C0 for the same draw.
# ================================================================
compute_R0_future <- function(t_index, year_id_future, alpha_R, year_effect_future, A_year_future,
                               beta_sin1, beta_cos1, beta_sin2, beta_cos2, beta_T_eff, beta_P_eff, z_T_future, z_P_future) {
  reg <- future_seasonal_regressors(t_index)
  log_R0 <- alpha_R + year_effect_future[year_id_future] +
    A_year_future[year_id_future] * (beta_sin1 * reg$sin1 + beta_cos1 * reg$cos1) +
    beta_sin2 * reg$sin2 + beta_cos2 * reg$cos2 +
    beta_T_eff * z_T_future + beta_P_eff * z_P_future
  list(log_R0 = log_R0, R0 = exp(log_R0))
}

# ================================================================
# 2.1 Bridge the ONE missing demographic step between the posterior's
# terminal state and the start of the future horizon.
#
# The accepted Stan model's S[t]/U[t] represent the state AT THE START of
# week t, BEFORE that week's births/deaths/reconciliation are applied (the
# `if (t < N)` guard means the final historical week's demographic update
# is never executed, since there is no week N+1 within the fit). Handing
# S[draw, N_hist]/U[draw, N_hist] directly to the future simulator would
# silently skip that week's demographic flow. This function completes it,
# so the future horizon starts from a state consistent with N_end[N_hist]
# (verified to match exactly -- see unit tests).
# ================================================================
advance_terminal_state <- function(S_last, U_last, X_last, births_last, deaths_last, reconciliation_last) {
  S_after <- S_last - X_last; U_after <- U_last + X_last
  sfrac <- S_after / (S_after + U_after)
  list(S = S_after + births_last - deaths_last * sfrac + reconciliation_last * sfrac,
       U = U_after - deaths_last * (1 - sfrac) + reconciliation_last * (1 - sfrac))
}

# ================================================================
# 2.6 Future demography arms.
#
# D1_ibge_projection: uses the ALREADY-EXISTING, documented IBGE UF annual
# population-projection file (2024 revision) -- the only future
# demographic source available in this repo. Only its RATES (annual
# growth rate, crude birth/death rate) are borrowed, applied multiplicatively
# to this model's OWN continuing population level (anchored at N_end of the
# last historical week) -- never IBGE's absolute population level, which
# uses a different base/definition than this model's own demographic
# bookkeeping and would otherwise create an artificial discontinuity at
# the 2025/2026 boundary.
#
# D0_stationary_2025_levels: sensitivity alternative. Freezes the FINAL
# HISTORICAL YEAR's own observed mean weekly growth rate and birth/death
# rates and repeats them unchanged for the whole future horizon -- no
# future acceleration or deceleration of demographic trends. Used whenever
# IBGE-projection data are judged unavailable/unsuitable, and in any case
# always run alongside D1 as the required sensitivity comparison.
# ================================================================
build_future_demography <- function(arm, future_dates, N_end_anchor, ibge_uf_annual = NULL, hist_weekly = NULL) {
  stopifnot(arm %in% DEMOGRAPHY_ARMS)
  n <- length(future_dates)
  years <- as.integer(format(future_dates, "%Y"))
  N_start <- numeric(n); N_end <- numeric(n); births <- numeric(n); deaths <- numeric(n)

  if (arm == "D1_ibge_projection") {
    stopifnot(!is.null(ibge_uf_annual))
    ibge_uf_annual <- ibge_uf_annual |> arrange(year)
    growth_rate <- setNames(c(NA, diff(log(ibge_uf_annual$pop_total))), ibge_uf_annual$year) # log-growth per year, weekly factor = exp(this/WEEKS_PER_YEAR)
    cbr <- setNames(ibge_uf_annual$births_total / ibge_uf_annual$pop_total, ibge_uf_annual$year)
    cdr <- setNames(ibge_uf_annual$deaths_total / ibge_uf_annual$pop_total, ibge_uf_annual$year)
    N_start[1] <- N_end_anchor
    for (t in seq_len(n)) {
      y <- as.character(years[t])
      weekly_growth <- exp(unname(growth_rate[y]) / WEEKS_PER_YEAR)
      N_end[t] <- N_start[t] * weekly_growth
      births[t] <- N_start[t] * unname(cbr[y]) / WEEKS_PER_YEAR
      deaths[t] <- N_start[t] * unname(cdr[y]) / WEEKS_PER_YEAR
      if (t < n) N_start[t + 1] <- N_end[t]
    }
  } else { # D0_stationary_2025_levels
    stopifnot(!is.null(hist_weekly))
    mean_weekly_growth <- exp(mean(diff(log(hist_weekly$N_end))))
    mean_birth_rate <- mean(hist_weekly$births / hist_weekly$N_start)
    mean_death_rate <- mean(hist_weekly$all_cause_deaths / hist_weekly$N_start)
    N_start[1] <- N_end_anchor
    for (t in seq_len(n)) {
      N_end[t] <- N_start[t] * mean_weekly_growth
      births[t] <- N_start[t] * mean_birth_rate
      deaths[t] <- N_start[t] * mean_death_rate
      if (t < n) N_start[t + 1] <- N_end[t]
    }
  }
  reconciliation <- N_end - N_start - births + deaths
  tibble(week_start = future_dates, N_start = N_start, N_end = N_end, births = births, deaths = deaths, reconciliation = reconciliation)
}

# ================================================================
# 2.5-2.6 Core renewal + demography recursion.
#
# Works for BOTH the historical-replay unit test (pass is_seed/X_seed,
# deterministic=TRUE, historical R0/demography) and genuine future
# simulation (is_seed=NULL, M_t stochastic importation instead).
#
# State convention: S0/U0 are the values AT week 1 (before that week's
# hazard is applied) -- i.e. S0 = S(t=1). `X_recent_history` supplies
# X(t=0), X(t=-1), ..., X(t=-(G-1)) (most recent LAST) so infectiousness
# at t=1 can be built from genuine prior history instead of zero.
# ================================================================
simulate_renewal_core <- function(R0_vec, w, S0, U0, N_start_vec, N_end_vec, births_vec, deaths_vec,
                                   X_recent_history, lambda_import = 0, stochastic = FALSE, rng_pois = NULL, rng_binom = NULL,
                                   is_seed = NULL, X_seed = NULL) {
  N_weeks <- length(R0_vec)
  G <- length(w)
  if (is.null(is_seed)) is_seed <- integer(N_weeks)
  if (is.null(X_seed)) X_seed <- numeric(N_weeks)
  stopifnot(length(X_recent_history) == G)

  S <- numeric(N_weeks); U <- numeric(N_weeks); X <- numeric(N_weeks)
  lambda_t_vec <- numeric(N_weeks); R_eff_vec <- numeric(N_weeks); M_vec <- numeric(N_weeks)
  S[1] <- S0; U[1] <- U0
  X_hist_full <- c(X_recent_history, X) # length G + N_weeks; X_hist_full[G+t] = X[t]; X_hist_full[G+t-g] for g=1..G is history

  for (t in seq_len(N_weeks)) {
    total_t <- S[t] + U[t]
    infectiousness <- 0
    for (g in seq_len(G)) infectiousness <- infectiousness + w[g] * X_hist_full[G + t - g]

    M_t <- if (lambda_import <= 0) 0 else if (stochastic) rng_pois(1, lambda_import) else lambda_import
    M_vec[t] <- M_t
    lambda_t <- R0_vec[t] * infectiousness / total_t + M_t / total_t
    lambda_t_vec[t] <- lambda_t
    R_eff_vec[t] <- R0_vec[t] * S[t] / total_t
    p_t <- -expm1(-lambda_t)

    if (is_seed[t] == 1L) {
      X[t] <- min(X_seed[t], S[t])
    } else if (stochastic) {
      X[t] <- rng_binom(1, round(S[t]), p_t) # S(t) rounded to the nearest integer for the Binomial trial count (negligible bias at population scale; verified in unit tests)
    } else {
      X[t] <- S[t] * p_t
    }
    X_hist_full[G + t] <- X[t]

    if (t < N_weeks) {
      S_after <- S[t] - X[t]; U_after <- U[t] + X[t]
      sfrac <- S_after / (S_after + U_after)
      reconciliation_t <- N_end_vec[t] - N_start_vec[t] - births_vec[t] + deaths_vec[t]
      S[t + 1] <- S_after + births_vec[t] - deaths_vec[t] * sfrac + reconciliation_t * sfrac
      U[t + 1] <- U_after - deaths_vec[t] * (1 - sfrac) + reconciliation_t * (1 - sfrac)
    }
  }
  N_vec <- S + U
  list(S = S, U = U, N = N_vec, X = X, R0 = R0_vec, R_eff = R_eff_vec, lambda = lambda_t_vec, M = M_vec)
}

# ================================================================
# 2.7 Observation layer -- latent X kept separate from reported Y.
# ================================================================
simulate_reported_cases <- function(X, q, phi_obs, stochastic = FALSE, rng_nb = NULL) {
  mu <- q * X + 1e-9
  if (!stochastic) return(mu) # expected reported cases (deterministic mode)
  rng_nb(length(X), mu, phi_obs)
}

# ================================================================
# 2.9 Outbreak/wave detection on SIMULATED REPORTED CASES -- adapted,
# simplified single-series version of the accepted national wave-census
# "primary" peak-trough algorithm (07_national_pipeline/05_national_wave_census_v2).
# Applied to Y_t (reported), never to latent X_t directly.
# ================================================================
centered_mean_future <- function(x, width = 3L) {
  half <- width %/% 2L
  vapply(seq_along(x), function(i) mean(x[max(1L, i - half):min(length(x), i + half)]), numeric(1))
}
local_peaks_future <- function(x, neighborhood = 3L) {
  if (length(x) < 3L) return(integer())
  out <- integer()
  for (i in 2:(length(x) - 1L)) {
    lo <- max(1L, i - neighborhood); hi <- min(length(x), i + neighborhood)
    if (x[i] > 0 && x[i] >= max(x[lo:hi]) && x[i] >= x[i - 1L] && x[i] >= x[i + 1L] && (x[i] > x[i - 1L] || x[i] > x[i + 1L])) out <- c(out, i)
  }
  out
}
detect_simulated_waves <- function(week_start, reported_cases, population, settings = WAVE_DETECTION) {
  incidence <- reported_cases / population * 1e5
  smooth <- centered_mean_future(incidence, settings$smooth_weeks)
  peaks <- local_peaks_future(smooth, settings$peak_neighborhood_weeks)
  if (!length(peaks)) return(tibble(onset_week = as.Date(character()), peak_week = as.Date(character()), end_week = as.Date(character()),
                                     total_cases = numeric(), cumulative_incidence_per_100k = numeric(), major_epidemic = logical()))
  rows <- vector("list", length(peaks))
  prev_end <- 1L
  for (i in seq_along(peaks)) {
    pk <- peaks[i]
    trough_before <- prev_end - 1L + which.min(smooth[prev_end:pk])
    baseline <- smooth[trough_before]; threshold <- baseline + settings$onset_fraction * (smooth[pk] - baseline)
    onset <- pk
    if (pk > trough_before + 1L) {
      candidates <- seq.int(trough_before + 1L, pk - 1L)
      sustained <- candidates[smooth[candidates] >= threshold]
      if (length(sustained)) onset <- sustained[1L]
    }
    after <- if (pk < length(smooth)) (pk + 1L):length(smooth) else integer()
    eligible <- after[smooth[after] <= settings$trough_fraction * smooth[pk]]
    end_idx <- if (length(eligible)) eligible[1L] else length(smooth)
    rows[[i]] <- tibble(onset_week = week_start[onset], peak_week = week_start[pk], end_week = week_start[end_idx],
                         total_cases = sum(reported_cases[onset:end_idx]), cumulative_incidence_per_100k = sum(incidence[onset:end_idx]))
    prev_end <- end_idx
  }
  bind_rows(rows) |> distinct(onset_week, .keep_all = TRUE) |>
    mutate(major_epidemic = total_cases >= settings$major_abs_case_minimum & cumulative_incidence_per_100k >= settings$major_incidence_minimum_per_100k)
}
