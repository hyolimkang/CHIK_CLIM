// Bahia-specific CONSERVATIVE extension of
// renewal_ceara_v4_9_hierarchical_seasonality_optional_serology.stan.
// NOT a new transmission model: the renewal equation, R0(t)
// (harmonic + hierarchical A_year), generation interval, demographic S/U
// recursion, lifelong immunity, imports, NB2 case likelihood, phi_obs
// prior, all transmission priors, fixed-q strategy, and the conditional
// seed mechanism are BYTE-IDENTICAL to that file. S[t]/U[t] retain their
// original epidemiological meaning (actual counts of susceptible/immune
// individuals in Bahia) -- there is no weighted/effective/FOI-weighted S,
// no spatial compartments, no risk strata.
//
// The ONLY new statistical component: the single-survey serology block
// (t_sero/sero_pos/sero_n/use_serology) is replaced by a J_sero=6-survey
// block that lets each survey site's true local seroprevalence differ
// from the Bahia STATE-level infection-derived immune proportion
// (p_state[t] = immune_prop[t] = U[t]/(S[t]+U[t]), unchanged) via one
// non-centred geographic log-odds offset per survey:
//   eta_geo[j] = sero_geographic_sd * z_geo[j],  z_geo[j] ~ Normal(0,1)
// with sero_geographic_sd FIXED (=1.0, supplied as data, NOT estimated).
// For survey j, p_site[j,t] = inv_logit(logit(p_state[t]) + eta_geo[j])
// is averaged (uniform weights, built and audited in R --
// 30_bahia_v4_9_replication/scripts/06_build_bahia_serology_windows.R)
// over that survey's collection-window weeks
// (sero_window_start_idx[j]:sero_window_start_idx[j]+sero_window_n_weeks[j]-1,
// a contiguous run by construction) to give p_site_window[j], which enters
// EXACTLY the same beta-binomial likelihood family and kappa_sero
// definition/value already used by the successful v4.9 model (kappa_sero
// fixed, not tuned here).
//
// No allfoi/bottom-up-FOI information of any kind is used anywhere in
// this file.

data {
  int<lower=2> N;
  int<lower=1> G;
  array[N] int<lower=0> C;
  simplex[G] w;
  vector<lower=0>[N] N_start;
  vector<lower=0>[N] N_end;
  vector<lower=0>[N] births;
  vector<lower=0>[N] deaths;
  int<lower=1> Y;
  array[N] int<lower=1, upper=Y> year_id;
  vector[N] seasonal_sin;
  vector[N] seasonal_cos;
  vector[N] seasonal_sin2;
  vector[N] seasonal_cos2;
  real<lower=0> imports_per_week;
  real<lower=0> year_effect_prior_sd;
  real<lower=0> kappa_sero; // fixed, same value/convention as v4.9 (50)
  real<lower=0, upper=1> q;

  array[N] int<lower=0, upper=1> is_seed;
  vector<lower=0>[N] X_seed;
  int<lower=0, upper=N> N_fit;
  array[N_fit] int<lower=1, upper=N> fit_index;

  // Multi-site Bahia serology (replaces the single Ceara/Juazeiro anchor).
  int<lower=1> J_sero;
  array[J_sero] int<lower=0> sero_n_positive;
  array[J_sero] int<lower=1> sero_n_tested;
  array[J_sero] int<lower=1, upper=N> sero_window_start_idx;
  array[J_sero] int<lower=1> sero_window_n_weeks;
  real<lower=0> sero_geographic_sd; // FIXED (=1.0), not estimated
}

transformed data {
  vector[N] reconciliation;
  for (t in 1:N) {
    reconciliation[t] = N_end[t] - N_start[t] - births[t] + deaths[t];
  }
}

parameters {
  real alpha_R;
  vector[Y] z_year;
  real beta_sin1;
  real beta_cos1;
  real beta_sin2;
  real beta_cos2;
  real<lower=0> sigma_season_year;
  vector[Y] z_season_year;
  real<lower=0> phi_obs;

  vector[J_sero] z_geo; // non-centred geographic offsets
}

transformed parameters {
  vector[Y] year_effect;
  vector[Y] A_year;
  vector[N] log_R0;
  vector[N] R0_t;
  vector[N] R_eff_t;
  vector[N] force_of_infection;
  vector[N] X;
  vector[N] S;
  vector[N] U;
  vector[N] S_prop;
  vector[N] immune_prop; // = p_state[t], unchanged epidemiological meaning
  vector[N] expected_reported_cases;

  vector[J_sero] eta_geo;
  vector[J_sero] p_site_window;
  vector[J_sero] alpha_sero_site;
  vector[J_sero] beta_sero_site;

  year_effect = year_effect_prior_sd * (z_year - mean(z_year));
  A_year = exp(sigma_season_year * z_season_year);
  eta_geo = sero_geographic_sd * z_geo;

  S[1] = N_start[1];
  U[1] = 0;

  for (t in 1:N) {
    real total = S[t] + U[t];
    real infectiousness = 0;

    if (S[t] < 0 || U[t] < 0 || fabs(total - N_start[t]) > 1e-6 * fmax(1, N_start[t])) {
      reject("Invalid demographic state at week ", t);
    }
    for (g in 1:G) {
      if (t > g) infectiousness += w[g] * X[t - g];
    }

    log_R0[t] = alpha_R + year_effect[year_id[t]] +
      A_year[year_id[t]] * (beta_sin1 * seasonal_sin[t] + beta_cos1 * seasonal_cos[t]) +
      beta_sin2 * seasonal_sin2[t] + beta_cos2 * seasonal_cos2[t];
    R0_t[t] = exp(log_R0[t]);
    force_of_infection[t] = R0_t[t] * infectiousness / total + imports_per_week / total;

    if (is_seed[t] == 1) {
      if (X_seed[t] > S[t]) {
        reject("Seed infections exceed available susceptibles at week ", t);
      }
      X[t] = X_seed[t];
    } else {
      X[t] = S[t] * (-expm1(-force_of_infection[t]));
    }

    R_eff_t[t] = R0_t[t] * S[t] / total;
    S_prop[t] = S[t] / total;
    immune_prop[t] = U[t] / total;
    expected_reported_cases[t] = q * X[t] + 1e-9;

    if (t < N) {
      real susceptible_after_infection = S[t] - X[t];
      real immune_after_infection = U[t] + X[t];
      real susceptible_fraction = susceptible_after_infection /
        (susceptible_after_infection + immune_after_infection);
      S[t + 1] = susceptible_after_infection + births[t] -
        deaths[t] * susceptible_fraction + reconciliation[t] * susceptible_fraction;
      U[t + 1] = immune_after_infection - deaths[t] * (1 - susceptible_fraction) +
        reconciliation[t] * (1 - susceptible_fraction);
    }
  }

  // Geographically-adjusted site-level seroprevalence, averaged (uniform
  // weights = simple mean over a contiguous window, per the R-built
  // weight audit) over each survey's collection window.
  for (j in 1:J_sero) {
    int start = sero_window_start_idx[j];
    int nwk = sero_window_n_weeks[j];
    real acc = 0;
    for (k in 1:nwk) {
      int t = start + k - 1;
      real p_state_safe = fmin(1 - 1e-9, fmax(1e-9, immune_prop[t]));
      acc += inv_logit(logit(p_state_safe) + eta_geo[j]);
    }
    p_site_window[j] = acc / nwk;
    alpha_sero_site[j] = p_site_window[j] * kappa_sero;
    beta_sero_site[j] = (1 - p_site_window[j]) * kappa_sero;
  }
}

model {
  alpha_R ~ normal(log(1.2), 0.5);
  z_year ~ std_normal();
  beta_sin1 ~ normal(0, 0.25);
  beta_cos1 ~ normal(0, 0.25);
  beta_sin2 ~ normal(0, 0.10);
  beta_cos2 ~ normal(0, 0.10);
  sigma_season_year ~ normal(0, 0.20);
  z_season_year ~ std_normal();
  phi_obs ~ gamma(2, 0.1);

  z_geo ~ std_normal(); // non-centred; sero_geographic_sd fixed, not estimated

  C[fit_index] ~ neg_binomial_2(expected_reported_cases[fit_index], phi_obs);
  for (j in 1:J_sero) {
    sero_n_positive[j] ~ beta_binomial(sero_n_tested[j], alpha_sero_site[j], beta_sero_site[j]);
  }
}

generated quantities {
  array[N] int C_pred;
  vector[N] log_lik_cases;
  vector[J_sero] log_lik_serology;
  array[J_sero] int sero_pred;
  vector[J_sero] p_state_window; // same window average, WITHOUT the geographic offset (comparison only)

  for (t in 1:N) {
    C_pred[t] = neg_binomial_2_rng(expected_reported_cases[t], phi_obs);
    log_lik_cases[t] = neg_binomial_2_lpmf(C[t] | expected_reported_cases[t], phi_obs);
  }
  for (j in 1:J_sero) {
    int start = sero_window_start_idx[j];
    int nwk = sero_window_n_weeks[j];
    real acc = 0;
    for (k in 1:nwk) {
      acc += fmin(1 - 1e-9, fmax(1e-9, immune_prop[start + k - 1]));
    }
    p_state_window[j] = acc / nwk;
    log_lik_serology[j] = beta_binomial_lpmf(sero_n_positive[j] | sero_n_tested[j], alpha_sero_site[j], beta_sero_site[j]);
    sero_pred[j] = beta_binomial_rng(sero_n_tested[j], alpha_sero_site[j], beta_sero_site[j]);
  }
}
