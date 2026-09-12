// Ceará v4.0: parsimonious all-age, no-vaccination renewal benchmark.
//
// This is intentionally a low-dimensional retrospective transmission model.
// It does not estimate reporting, introductions, climate effects, serology
// offsets, or a weekly latent R process. Those quantities are either fixed
// scenario inputs or reserved for later, separately validated extensions.

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
  real<lower=0, upper=1> q_fixed;
  real<lower=0> imports_per_week;
  real<lower=0> year_effect_prior_sd;
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
  real beta_sin;
  real beta_cos;
  real<lower=0> phi_obs;
}

transformed parameters {
  vector[Y] year_effect;
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

  year_effect = year_effect_prior_sd * (z_year - mean(z_year));
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
      beta_sin * seasonal_sin[t] + beta_cos * seasonal_cos[t];
    R0_t[t] = exp(log_R0[t]);
    force_of_infection[t] = R0_t[t] * infectiousness / total + imports_per_week / total;
    X[t] = S[t] * (-expm1(-force_of_infection[t]));
    R_eff_t[t] = R0_t[t] * S[t] / total;
    S_prop[t] = S[t] / total;
    immune_prop[t] = U[t] / total;
    expected_reported_cases[t] = q_fixed * X[t] + 1e-9;

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
}

model {
  alpha_R ~ normal(log(1.2), 0.5);
  z_year ~ std_normal();
  beta_sin ~ normal(0, 0.25);
  beta_cos ~ normal(0, 0.25);
  phi_obs ~ gamma(2, 0.1);
  C ~ neg_binomial_2(expected_reported_cases, phi_obs);
}

generated quantities {
  array[N] int C_pred;
  vector[N] log_lik_cases;
  for (t in 1:N) {
    C_pred[t] = neg_binomial_2_rng(expected_reported_cases[t], phi_obs);
    log_lik_cases[t] = neg_binomial_2_lpmf(C[t] | expected_reported_cases[t], phi_obs);
  }
}
