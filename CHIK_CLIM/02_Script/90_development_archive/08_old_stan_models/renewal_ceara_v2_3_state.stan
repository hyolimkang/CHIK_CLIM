// v2.3 state renewal model: explicit demography, direct detection, and a
// stationary seasonal AR(1) R0 process. No recurrent importation is included.
data {
  int<lower=2> N; int<lower=1> G; int<lower=1> seed_weeks;
  array[N] int<lower=0> C; simplex[G] w;
  vector<lower=0>[N] N_start; vector<lower=0>[N] N_end;
  vector<lower=0>[N] births; vector<lower=0>[N] deaths;
  vector[N] seasonal_sin; vector[N] seasonal_cos;
  vector[seed_weeks] log_seed_prior_mean; real<lower=0> seed_prior_sd;
  real<lower=0> q_prior_a; real<lower=0> q_prior_b;
  int<lower=1> J; array[J] int<lower=1,upper=N> sero_window_start; array[J] int<lower=1,upper=N> sero_window_end;
  array[J] int<lower=0> sero_positive; array[J] int<lower=1> sero_n;
  real<lower=0,upper=1> sero_power;
}
transformed data {
  vector[N] reconciliation;
  if (seed_weeks < G || seed_weeks >= N) reject("Require G <= seed_weeks < N");
  for (t in 1:N) reconciliation[t] = N_end[t] - N_start[t] - births[t] + deaths[t];
}
parameters {
  real alpha_R; real beta_sin; real beta_cos;
  real<lower=0,upper=.995> phi_R; real<lower=0> sigma_R; vector[N] z_R;
  vector[seed_weeks] log_seed_hazard;
  real<lower=0,upper=1> q;
  real<lower=0> sigma_geo; vector[J] z_site;
  real<lower=0> phi_obs;
}
transformed parameters {
  vector[N] r_deviation; vector[N] log_R0; vector[N] R0_t; vector[N] R_eff_t;
  vector[N] X; vector[N] S; vector[N] U; vector[N] immune_prop; vector[N] S_prop;
  vector[N] expected_reported_cases; vector[J] state_attack_window; vector[J] p_site;
  r_deviation[1] = sigma_R / sqrt(1 - square(phi_R)) * z_R[1];
  for (t in 2:N) r_deviation[t] = phi_R * r_deviation[t - 1] + sigma_R * z_R[t];
  for (t in 1:N) { log_R0[t] = alpha_R + beta_sin * seasonal_sin[t] + beta_cos * seasonal_cos[t] + r_deviation[t]; R0_t[t] = exp(log_R0[t]); }
  S[1] = N_start[1]; U[1] = 0;
  for (t in 1:N) {
    real total = S[t] + U[t];
    if (S[t] < 0 || U[t] < 0 || fabs(total - N_start[t]) > 1e-6 * fmax(1, N_start[t])) reject("Invalid demographic state at week ", t);
    if (t <= seed_weeks) X[t] = S[t] * (-expm1(-exp(log_seed_hazard[t])));
    else {
      real infectiousness = 0;
      for (g in 1:G) infectiousness += w[g] * X[t - g];
      X[t] = S[t] * (-expm1(-R0_t[t] * infectiousness / total));
    }
    R_eff_t[t] = R0_t[t] * S[t] / total; S_prop[t] = S[t] / total; immune_prop[t] = U[t] / total;
    expected_reported_cases[t] = q * X[t] + 1e-9;
    if (t < N) {
      real s_after = S[t] - X[t]; real u_after = U[t] + X[t]; real frac_s = s_after / (s_after + u_after);
      S[t + 1] = s_after + births[t] - deaths[t] * frac_s + reconciliation[t] * frac_s;
      U[t + 1] = u_after - deaths[t] * (1 - frac_s) + reconciliation[t] * (1 - frac_s);
      if (S[t + 1] < 0 || U[t + 1] < 0) reject("Negative demographic state at week ", t + 1);
    }
  }
  for (j in 1:J) {
    real a = 0;
    for (t in sero_window_start[j]:sero_window_end[j]) a += immune_prop[t];
    state_attack_window[j] = a / (sero_window_end[j] - sero_window_start[j] + 1);
    p_site[j] = inv_logit(logit(state_attack_window[j]) + sigma_geo * z_site[j]);
  }
}
model {
  alpha_R ~ normal(log(1.2), .5); beta_sin ~ normal(0, .25); beta_cos ~ normal(0, .25);
  phi_R ~ beta(2, 2); sigma_R ~ normal(0, .10); z_R ~ normal(0, 1);
  log_seed_hazard ~ normal(log_seed_prior_mean - log(segment(N_start, 1, seed_weeks)), seed_prior_sd);
  q ~ beta(q_prior_a, q_prior_b); sigma_geo ~ normal(0, 1); z_site ~ normal(0, 1); phi_obs ~ gamma(2, .1);
  C ~ neg_binomial_2(expected_reported_cases, phi_obs);
  for (j in 1:J) target += sero_power * binomial_lpmf(sero_positive[j] | sero_n[j], p_site[j]);
}
generated quantities {
  array[N] int C_pred; vector[N] overall_detection; array[J] int sero_positive_pred;
  for (t in 1:N) { C_pred[t] = neg_binomial_2_rng(expected_reported_cases[t], phi_obs); overall_detection[t] = q; }
  for (j in 1:J) sero_positive_pred[j] = binomial_rng(sero_n[j], p_site[j]);
}
