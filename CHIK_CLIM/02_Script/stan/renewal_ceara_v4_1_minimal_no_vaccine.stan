// Ceará v4.1: parsimonious all-age, no-vaccination renewal benchmark.
//
// Same scientific backbone as v4.0 (renewal_ceara_v4_0_minimal_no_vaccine.stan
// -- HMC-passing at max_treedepth=14), with exactly two substantive changes:
//
// A. The annual transmission-deviation term: v4.0 sampled Y independent
//    z_year and recentred them (z_year - mean(z_year)), a sum-to-zero effect
//    that still retains one likelihood-redundant direction. v4.1 samples only
//    the Y-1 identified degrees of freedom via an explicit orthonormal
//    sum-to-zero basis B_year (built in R, passed as data), and learns a
//    single annual-heterogeneity scale sigma_year instead of fixing it.
//
// B. The reporting fraction: v4.0 fixed q_fixed = 0.10 as a data input with
//    zero uncertainty. v4.1 estimates ONE constant q with an informative
//    prior approximating the induced distribution of
//    Beta(30,28) [p_symptomatic] * Beta(20,60) [p_detect_given_symptomatic]
//    -- moment-matched and Kolmogorov-Smirnov-checked in
//    00_check_q_prior_approximation.R (Beta(16.2, 108.7) confirmed adequate,
//    KS = 0.0096). p_symptomatic and p_detect_given_symptomatic are NOT
//    estimated separately: the case likelihood identifies mainly their
//    product, so separating them would add a non-identifiable direction.
//
// No other latent process is added: still no climate effects, serology
// likelihood, age structure, vaccination, weekly latent R process, AR(1),
// fitted importation, or additional seasonal harmonics.

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
  real<lower=0> imports_per_week;
  matrix[Y, Y - 1] B_year;
  real<lower=0> q_prior_a;
  real<lower=0> q_prior_b;
}

transformed data {
  vector[N] reconciliation;
  for (t in 1:N) {
    reconciliation[t] = N_end[t] - N_start[t] - births[t] + deaths[t];
  }
}

parameters {
  real alpha_R;
  vector[Y - 1] z_year_free;
  real<lower=0> sigma_year;
  real beta_sin;
  real beta_cos;
  real<lower=0> phi_obs;
  real<lower=0, upper=1> q;
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

  year_effect = sigma_year * (B_year * z_year_free);
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
}

model {
  alpha_R ~ normal(log(1.2), 0.5);
  z_year_free ~ std_normal();
  sigma_year ~ normal(0, 0.35); // truncated at 0 by the <lower=0> declaration: half-normal
  beta_sin ~ normal(0, 0.25);
  beta_cos ~ normal(0, 0.25);
  phi_obs ~ gamma(2, 0.1);
  q ~ beta(q_prior_a, q_prior_b);
  C ~ neg_binomial_2(expected_reported_cases, phi_obs);
}

generated quantities {
  array[N] int C_pred;
  vector[N] log_lik_cases;
  real sum_year_effect = sum(year_effect); // should be ~0 for every draw; QC only
  for (t in 1:N) {
    C_pred[t] = neg_binomial_2_rng(expected_reported_cases[t], phi_obs);
    log_lik_cases[t] = neg_binomial_2_lpmf(C[t] | expected_reported_cases[t], phi_obs);
  }
}
