// DIAGNOSTIC 4-WEEK OBSERVATION SENSITIVITY (mode audit Part 4). NOT a
// production model. The latent WEEKLY renewal/S/U dynamics, transmission
// model, q definition, Juazeiro likelihood, q prior, imports, generation
// interval, seasonality, and annual effects are ALL unchanged from
// renewal_ceara_v4_3_weak_q.stan. The ONLY change: the case observation
// likelihood is aggregated into non-overlapping 4-week blocks (with its own
// estimated block-level phi_obs_block, not the weekly phi_obs) instead of
// 313 separate weekly NB likelihood terms. Tests whether the high-q mode's
// dominance depends on treating temporally autocorrelated weekly residuals
// (ACF~0.85 at lag 1, confirmed in 07_weekly_residual_acf.R) as if they were
// independent evidence.

data {
  int<lower=2> N;
  int<lower=1> G;
  array[N] int<lower=0> C; // still read in (used for C_pred/log_lik_cases diagnostics only, not the primary likelihood)
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
  real<lower=0> year_effect_prior_sd;
  int<lower=1, upper=N> t_sero;
  int<lower=0> sero_pos;
  int<lower=1> sero_n;
  real<lower=0> kappa_sero;
  // 4-week block aggregation (Section 2): non-overlapping consecutive
  // 4-week blocks. The final incomplete block (313 weeks = 78 full blocks
  // of 4 + 1 leftover week) is EXCLUDED from the observation likelihood in
  // BOTH fits (documented rule) via include_in_block=0 for that week; the
  // latent weekly renewal recursion still runs for all N weeks regardless.
  int<lower=1> n_blocks;
  array[N] int<lower=1, upper=n_blocks> block_id;
  array[N] int<lower=0, upper=1> include_in_block;
  int<lower=0> use_serology;
  array[n_blocks] int<lower=0> observed_block;
}

transformed data {
  vector[N] reconciliation;
  real eta_q_prior_mean = logit(0.05);
  for (t in 1:N) {
    reconciliation[t] = N_end[t] - N_start[t] - births[t] + deaths[t];
  }
}

parameters {
  real alpha_R;
  vector[Y] z_year;
  real beta_sin;
  real beta_cos;
  real<lower=0> phi_obs_block; // block-level dispersion, estimated fresh -- NOT copied from any weekly phi_obs
  real eta_q;
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
  vector[n_blocks] expected_block = rep_vector(0, n_blocks);
  real p_state_sero_at_anchor;
  real p_sero_safe;
  real alpha_sero;
  real beta_sero;
  real<lower=0, upper=1> q = inv_logit(eta_q);

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
    expected_reported_cases[t] = q * X[t] + 1e-9;
    expected_block[block_id[t]] += include_in_block[t] * expected_reported_cases[t];

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
  beta_sin ~ normal(0, 0.25);
  beta_cos ~ normal(0, 0.25);
  phi_obs_block ~ gamma(2, 0.1);
  eta_q ~ normal(eta_q_prior_mean, 1.25);
  observed_block ~ neg_binomial_2(expected_block, phi_obs_block); // 4-week aggregated likelihood, NOT weekly
  if (use_serology) {
    sero_pos ~ beta_binomial(sero_n, alpha_sero, beta_sero);
  }
}

generated quantities {
  array[n_blocks] int block_pred;
  vector[n_blocks] log_lik_block;
  real log_lik_serology;
  int sero_pred;
  real log_prior_q;
  for (b in 1:n_blocks) {
    block_pred[b] = neg_binomial_2_rng(expected_block[b], phi_obs_block);
    log_lik_block[b] = neg_binomial_2_lpmf(observed_block[b] | expected_block[b], phi_obs_block);
  }
  log_lik_serology = beta_binomial_lpmf(sero_pos | sero_n, alpha_sero, beta_sero);
  sero_pred = beta_binomial_rng(sero_n, alpha_sero, beta_sero);
  log_prior_q = normal_lpdf(eta_q | eta_q_prior_mean, 1.25);
}
