// Rio de Janeiro CITY auxiliary model B: frozen v4.9 transmission
// structure (renewal equation, R0(t) harmonic + hierarchical A_year,
// generation interval, demographic S/U recursion, lifelong immunity,
// imports, NB2 case likelihood, phi_obs prior, global-q strategy) --
// BYTE-IDENTICAL to renewal_bahia_v4_9_global_q_multisite_serology.stan
// in every respect EXCEPT the serology observation model.
//
// Because U10 (Perisse 2020) used a stratified household sample with
// design/nonresponse weights, and the reported weighted prevalence (0.180,
// 95% CI [0.148,0.212]) differs from the raw count (371/2120=0.175), a
// design-aware logit-normal observation model is used on the WEIGHTED
// prevalence directly, rather than treating 371/2120 as a naive simple
// binomial sample via beta-binomial:
//
//   logit(p_hat_U10) ~ Normal(logit(p_city_window), SE_logit_U10)
//
// SE_logit_U10 is derived (in R, passed as data) from the reported 95% CI
// via the standard delta-method approximation:
//   SE_logit = (logit(CI_hi) - logit(CI_lo)) / (2 * 1.96)
//
// Because geography is matched (city cases, city serology, city latent
// S/U), there is NO geographic-offset term (eta_geo = 0 by construction,
// not estimated) -- p_city_window is the DIRECT mean city immune fraction
// over the survey window, exactly as immune_prop[t] already is for the
// homogeneous state models.

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

  array[N] int<lower=0, upper=1> is_seed;
  vector<lower=0>[N] X_seed;
  int<lower=0, upper=N> N_fit;
  array[N_fit] int<lower=1, upper=N> fit_index;

  // Rio de Janeiro CITY serology (U10, Perisse 2020) -- design-aware
  // logit-normal likelihood on the WEIGHTED prevalence, geography matched
  // (no eta_geo).
  int<lower=1> J_sero;
  array[J_sero] real<lower=0, upper=1> p_hat_sero; // observed weighted prevalence
  array[J_sero] real<lower=0> se_logit_sero; // delta-method SE on the logit scale
  array[J_sero] int<lower=1, upper=N> sero_window_start_idx;
  array[J_sero] int<lower=1> sero_window_n_weeks;

  real logit_q_prior_mean; // logit(0.10)
  real<lower=0> logit_q_prior_sd; // 1.0, SAME broad prior family as the state model
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

  real logit_q; // single global, time-constant q_city
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
  vector[N] immune_prop;
  vector[N] expected_reported_cases;
  real<lower=0, upper=1> q;

  vector[J_sero] p_city_window;

  q = inv_logit(logit_q);
  year_effect = year_effect_prior_sd * (z_year - mean(z_year));
  A_year = exp(sigma_season_year * z_season_year);

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

  // Direct city-level serology link -- NO geographic offset (geography matched).
  for (j in 1:J_sero) {
    int start = sero_window_start_idx[j];
    int nwk = sero_window_n_weeks[j];
    real acc = 0;
    for (k in 1:nwk) {
      int t = start + k - 1;
      acc += fmin(1 - 1e-9, fmax(1e-9, immune_prop[t]));
    }
    p_city_window[j] = acc / nwk;
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

  logit_q ~ normal(logit_q_prior_mean, logit_q_prior_sd);

  C[fit_index] ~ neg_binomial_2(expected_reported_cases[fit_index], phi_obs);

  for (j in 1:J_sero) {
    logit(p_hat_sero[j]) ~ normal(logit(p_city_window[j]), se_logit_sero[j]);
  }
}

generated quantities {
  array[N] int C_pred;
  vector[N] log_lik_cases;
  vector[J_sero] log_lik_serology;
  array[J_sero] real sero_pred; // predictive draw of the weighted prevalence (logit-normal)
  real logit_q_prior_draw = normal_rng(logit_q_prior_mean, logit_q_prior_sd);
  real q_prior_draw = inv_logit(logit_q_prior_draw);

  for (t in 1:N) {
    C_pred[t] = neg_binomial_2_rng(expected_reported_cases[t], phi_obs);
    log_lik_cases[t] = neg_binomial_2_lpmf(C[t] | expected_reported_cases[t], phi_obs);
  }
  for (j in 1:J_sero) {
    log_lik_serology[j] = normal_lpdf(logit(p_hat_sero[j]) | logit(p_city_window[j]), se_logit_sero[j]);
    sero_pred[j] = inv_logit(normal_rng(logit(p_city_window[j]), se_logit_sero[j]));
  }
}
