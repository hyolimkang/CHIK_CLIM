// PROPOSED Phase-2 canary extension of the frozen v4.9 model
// (renewal_ceara_v4_9_hierarchical_seasonality.stan). NOT YET FITTED --
// created for review only, per the climate-forced v4.9 extension Phase 2
// specification.
//
// Every line below is byte-identical to the frozen file EXCEPT the four
// blocks marked "NEW" / "CHANGED". No spatial compartments, no dynamic AR
// process, no interaction term, no change to q treatment, observation
// model, generation interval, demographic accounting, seed mechanism, or
// serology likelihood.
//
// Diff summary:
//   data:                + z_T_anom[N], z_P_anom[N]   (standardised climate anomalies)
//   parameters:           + beta_T, beta_P
//   transformed parameters: log_R0[t] gains "+ beta_T*z_T_anom[t] + beta_P*z_P_anom[t]"
//   model:                + beta_T ~ normal(0,0.15); beta_P ~ normal(0,0.15)
//   generated quantities: + climate_multiplier[t] = exp(beta_T*z_T_anom[t] + beta_P*z_P_anom[t])

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
  int<lower=1, upper=N> t_sero;
  int<lower=0> sero_pos;
  int<lower=1> sero_n;
  real<lower=0> kappa_sero;
  real<lower=0, upper=1> q; // FIXED, unchanged from v4.9

  array[N] int<lower=0, upper=1> is_seed;
  vector<lower=0>[N] X_seed;
  int<lower=0, upper=N> N_fit;
  array[N_fit] int<lower=1, upper=N> fit_index;

  vector[N] z_T_anom; // NEW: standardised temperature anomaly (primary lag spec)
  vector[N] z_P_anom; // NEW: standardised precipitation anomaly (primary lag spec)
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

  real beta_T; // NEW: shared climate-temperature coefficient on log_R0
  real beta_P; // NEW: shared climate-precipitation coefficient on log_R0
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
  real p_state_sero_at_anchor;
  real p_sero_safe;
  real alpha_sero;
  real beta_sero;

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

    // CHANGED: log_R0[t] gains the shared climate-anomaly forcing term.
    // Everything else in this line is unchanged from v4.9.
    log_R0[t] = alpha_R + year_effect[year_id[t]] +
      A_year[year_id[t]] * (beta_sin1 * seasonal_sin[t] + beta_cos1 * seasonal_cos[t]) +
      beta_sin2 * seasonal_sin2[t] + beta_cos2 * seasonal_cos2[t] +
      beta_T * z_T_anom[t] + beta_P * z_P_anom[t];
    R0_t[t] = exp(log_R0[t]);
    force_of_infection[t] = R0_t[t] * infectiousness / total + imports_per_week / total;

    if (is_seed[t] == 1) {
      if (X_seed[t] > S[t]) {
        reject("2022 seed infections exceed available susceptibles at week ", t);
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

  p_state_sero_at_anchor = immune_prop[t_sero];
  p_sero_safe = fmin(1 - 1e-9, fmax(1e-9, p_state_sero_at_anchor));
  alpha_sero = p_sero_safe * kappa_sero;
  beta_sero = (1 - p_sero_safe) * kappa_sero;
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

  beta_T ~ normal(0, 0.15); // NEW: regularising, centred at zero (no climate effect)
  beta_P ~ normal(0, 0.15); // NEW: regularising, centred at zero (no climate effect)

  C[fit_index] ~ neg_binomial_2(expected_reported_cases[fit_index], phi_obs);
  sero_pos ~ beta_binomial(sero_n, alpha_sero, beta_sero);
}

generated quantities {
  array[N] int C_pred;
  vector[N] log_lik_cases;
  vector[N] climate_multiplier; // NEW: must be declared here
  real log_lik_serology;
  int sero_pred;

  for (t in 1:N) {
    climate_multiplier[t] =
      exp(beta_T * z_T_anom[t] +
          beta_P * z_P_anom[t]);

    C_pred[t] = neg_binomial_2_rng(expected_reported_cases[t], phi_obs);
    log_lik_cases[t] = neg_binomial_2_lpmf(C[t] | expected_reported_cases[t], phi_obs);
  }
  log_lik_serology = beta_binomial_lpmf(sero_pos | sero_n, alpha_sero, beta_sero);
  sero_pred = beta_binomial_rng(sero_n, alpha_sero, beta_sero);
}
