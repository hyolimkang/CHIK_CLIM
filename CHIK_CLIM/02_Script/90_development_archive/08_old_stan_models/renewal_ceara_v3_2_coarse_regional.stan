// Coarse regional direct-renewal model: 8 official-mapping units, 2015--2019.
data {
  int<lower=2> K;
  int<lower=2> N;
  int<lower=1> G;
  int<lower=1> seed_weeks;
  array[K, N] int<lower=0> C;
  simplex[G] w;
  simplex[seed_weeks] seed_profile;
  matrix<lower=0>[K, N] N_start;
  matrix<lower=0>[K, N] N_end;
  matrix<lower=0>[K, N] births;
  matrix<lower=0>[K, N] deaths;
  real<lower=0> q_prior_a;
  real<lower=0> q_prior_b;
  int<lower=1> J;
  array[J] int<lower=1, upper=K> sero_unit;
  array[J] int<lower=1, upper=N> sero_window_start;
  array[J] int<lower=1, upper=N> sero_window_end;
  array[J] int<lower=0> sero_positive;
  array[J] int<lower=1> sero_n;
}
transformed data {
  matrix[K, N] reconciliation;
  if (seed_weeks < G || seed_weeks >= N) reject("Require G <= seed_weeks < N");
  for (j in 1:J)
    if (sero_window_start[j] > sero_window_end[j]) reject("Invalid serology window");
  for (k in 1:K)
    for (t in 1:N)
      reconciliation[k, t] = N_end[k, t] - N_start[k, t] - births[k, t] + deaths[k, t];
}
parameters {
  real alpha_R;
  real<lower=0> sigma_R;
  vector[N - 1] z_R;
  real<lower=0> sigma_region;
  vector[K] z_region;
  real mu_seed;
  real<lower=0> sigma_seed;
  vector[K] z_seed;
  real<lower=0, upper=1> q;
  real<lower=0> phi_obs;
}
transformed parameters {
  vector[N] mu_R;
  vector[K] region_offset = sigma_region * z_region;
  vector[K] log_seed = mu_seed + sigma_seed * z_seed;
  matrix[K, N] R0;
  matrix[K, N] R_eff;
  matrix[K, N] X;
  matrix[K, N] S;
  matrix[K, N] U;
  matrix[K, N] S_prop;
  matrix[K, N] immune_prop;
  matrix[K, N] expected_reported_cases;
  vector<lower=0, upper=1>[J] p_sero;

  mu_R[1] = alpha_R;
  for (t in 2:N) mu_R[t] = mu_R[t - 1] + sigma_R * z_R[t - 1];

  for (k in 1:K) {
    S[k, 1] = N_start[k, 1];
    U[k, 1] = 0;
    for (t in 1:N) {
      real population_t = S[k, t] + U[k, t];
      if (S[k, t] < 0 || U[k, t] < 0 ||
          fabs(population_t - N_start[k, t]) > 1e-6 * fmax(1.0, N_start[k, t]))
        reject("Invalid S/U accounting for region ", k, " week ", t);
      R0[k, t] = exp(mu_R[t] + region_offset[k]);
      if (t <= seed_weeks) {
        X[k, t] = S[k, t] * (-expm1(-seed_profile[t] * exp(log_seed[k]) / population_t));
      } else {
        real infectiousness = 0;
        for (g in 1:G) infectiousness += w[g] * X[k, t - g];
        X[k, t] = S[k, t] * (-expm1(-R0[k, t] * infectiousness / population_t));
      }
      R_eff[k, t] = R0[k, t] * S[k, t] / population_t;
      S_prop[k, t] = S[k, t] / population_t;
      immune_prop[k, t] = U[k, t] / population_t;
      expected_reported_cases[k, t] = q * X[k, t] + 1e-9;
      if (t < N) {
        real s_after = S[k, t] - X[k, t];
        real u_after = U[k, t] + X[k, t];
        real frac_s = s_after / (s_after + u_after);
        S[k, t + 1] = s_after + births[k, t] - deaths[k, t] * frac_s
                      + reconciliation[k, t] * frac_s;
        U[k, t + 1] = u_after - deaths[k, t] * (1 - frac_s)
                      + reconciliation[k, t] * (1 - frac_s);
        if (S[k, t + 1] < 0 || U[k, t + 1] < 0)
          reject("Negative demographic state for region ", k, " week ", t + 1);
      }
    }
  }
  for (j in 1:J) {
    real total = 0;
    for (t in sero_window_start[j]:sero_window_end[j]) total += immune_prop[sero_unit[j], t];
    p_sero[j] = total / (sero_window_end[j] - sero_window_start[j] + 1);
  }
}
model {
  alpha_R ~ normal(log(1), 0.5);
  sigma_R ~ normal(0, 0.10);
  z_R ~ normal(0, 1);
  sigma_region ~ normal(0, 0.25);
  z_region ~ normal(0, 1);
  mu_seed ~ normal(log(10), 1.5);
  sigma_seed ~ normal(0, 1);
  z_seed ~ normal(0, 1);
  q ~ beta(q_prior_a, q_prior_b);
  phi_obs ~ gamma(2, 0.1);
  for (k in 1:K)
    for (t in 1:N)
      C[k, t] ~ neg_binomial_2(expected_reported_cases[k, t], phi_obs);
  for (j in 1:J) sero_positive[j] ~ binomial(sero_n[j], p_sero[j]);
}
generated quantities {
  array[K, N] int C_pred;
  array[N] int C_ceara;
  array[N] int C_ceara_pred;
  vector[N] X_ceara;
  vector[N] S_ceara_prop;
  vector[N] immune_ceara_prop;
  vector[N] R_eff_ceara;
  array[J] int sero_positive_pred;
  for (t in 1:N) {
    int c_obs = 0;
    int c_pred = 0;
    real s_total = 0;
    real u_total = 0;
    real x_total = 0;
    real reff_total = 0;
    for (k in 1:K) {
      C_pred[k, t] = neg_binomial_2_rng(expected_reported_cases[k, t], phi_obs);
      c_obs += C[k, t]; c_pred += C_pred[k, t];
      s_total += S[k, t]; u_total += U[k, t]; x_total += X[k, t];
      reff_total += R_eff[k, t] * (S[k, t] + U[k, t]);
    }
    C_ceara[t] = c_obs; C_ceara_pred[t] = c_pred; X_ceara[t] = x_total;
    S_ceara_prop[t] = s_total / (s_total + u_total);
    immune_ceara_prop[t] = u_total / (s_total + u_total);
    R_eff_ceara[t] = reff_total / (s_total + u_total);
  }
  for (j in 1:J) sero_positive_pred[j] = binomial_rng(sero_n[j], p_sero[j]);
}
