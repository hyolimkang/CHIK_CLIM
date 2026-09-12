// Retrospective annual FOI / susceptibility reconstruction.
// This is intentionally an annual, state-level model: no renewal process,
// R0/Rt, climate covariates, weekly transmission, or forecasting component.

data {
  int<lower=2> Y;
  array[Y] int<lower=0> C;
  vector<lower=0>[Y] N_start;
  vector<lower=0>[Y] N_end;
  vector<lower=0>[Y] births;
  vector<lower=0>[Y] deaths;
  vector[Y] reconciliation;

  real log_foi_prior_mean;
  real<lower=0> log_foi_prior_sd;
  real q_logit_prior_mean;
  real<lower=0> q_logit_prior_sd;
  int<lower=0, upper=1> use_dynamic;
}

parameters {
  real log_mu_lambda;
  real logit_q;
  real log_phi_obs;
  real<lower=-0.95, upper=0.95> rho_delta;
  real log_sigma_delta;
  vector[Y] z_delta;
}

transformed parameters {
  real<lower=0> mu_lambda = exp(log_mu_lambda);
  real<lower=0, upper=1> q = inv_logit(logit_q);
  real<lower=0> phi_obs = exp(log_phi_obs);
  real<lower=0> sigma_delta = exp(log_sigma_delta);
  vector[Y] delta;
  vector<lower=0>[Y] lambda;
  vector<lower=0, upper=1>[Y] attack_prob;
  vector<lower=0>[Y] latent_infections;
  vector<lower=0>[Y] S_start;
  vector<lower=0>[Y] S_after;
  vector<lower=0>[Y] S_end;
  vector<lower=0>[Y] U_start;
  vector<lower=0>[Y] U_end;
  vector<lower=0, upper=1>[Y] S_start_prop;
  vector<lower=0, upper=1>[Y] S_end_prop;
  vector<lower=0>[Y] expected_cases;

  if (use_dynamic == 1) {
    delta[1] = sigma_delta * z_delta[1] / sqrt(1 - square(rho_delta));
    for (y in 2:Y)
      delta[y] = rho_delta * delta[y - 1] + sigma_delta * z_delta[y];
  } else {
    delta = rep_vector(0, Y);
  }

  S_start[1] = N_start[1];
  U_start[1] = 0;
  for (y in 1:Y) {
    lambda[y] = exp(log_mu_lambda + delta[y]);
    attack_prob[y] = -expm1(-lambda[y]);
    latent_infections[y] = S_start[y] * attack_prob[y];
    S_after[y] = S_start[y] - latent_infections[y];
    U_end[y] = U_start[y] + latent_infections[y];
    S_start_prop[y] = S_start[y] / N_start[y];
    expected_cases[y] = q * latent_infections[y] + 1e-9;

    if (y < Y) {
      real frac_s = S_after[y] / N_start[y];
      real deaths_s = deaths[y] * frac_s;
      real reconciliation_s = reconciliation[y] * frac_s;
      S_end[y] = S_after[y] + births[y] - deaths_s + reconciliation_s;
      U_end[y] = U_end[y] - (deaths[y] - deaths_s)
                 + (reconciliation[y] - reconciliation_s);
      S_start[y + 1] = S_end[y];
      U_start[y + 1] = U_end[y];
    } else {
      S_end[y] = S_after[y] + births[y] - deaths[y] * (S_after[y] / N_start[y])
                 + reconciliation[y] * (S_after[y] / N_start[y]);
      U_end[y] = U_end[y] - deaths[y] * (U_end[y] / N_start[y])
                 + reconciliation[y] * (U_end[y] / N_start[y]);
    }
    S_end_prop[y] = S_end[y] / N_end[y];

    if (fabs((S_end[y] + U_end[y]) - N_end[y]) > 1e-5 * N_end[y])
      reject("Demographic accounting does not conserve population in year ", y);
  }
}

model {
  log_mu_lambda ~ normal(log_foi_prior_mean, log_foi_prior_sd);
  logit_q ~ normal(q_logit_prior_mean, q_logit_prior_sd);
  // These two priors are algebraically the same as their original-scale
  // forms, with the Jacobian terms included for log-scale sampling.
  target += gamma_lpdf(phi_obs | 2, 0.1) + log_phi_obs;
  rho_delta ~ normal(0, 0.45);
  target += normal_lpdf(sigma_delta | 0, 0.75) + log_sigma_delta;
  z_delta ~ normal(0, 1);

  C ~ neg_binomial_2(expected_cases, phi_obs);
}

generated quantities {
  array[Y] int C_pred;
  vector[Y] log_lik_cases;
  real realized_mean_lambda_2015_2025 = mean(lambda);
  real realized_geometric_mean_lambda_2015_2025 = exp(mean(log(lambda)));
  real cumulative_latent_infections = sum(latent_infections);
  real realized_mean_to_mu_ratio = realized_mean_lambda_2015_2025 / mu_lambda;
  real implied_long_run_arithmetic_mean_lambda = use_dynamic == 1
    ? mu_lambda * exp(0.5 * square(sigma_delta) / (1 - square(rho_delta)))
    : mu_lambda;

  for (y in 1:Y) {
    C_pred[y] = neg_binomial_2_rng(expected_cases[y], phi_obs);
    log_lik_cases[y] = neg_binomial_2_lpmf(C[y] | expected_cases[y], phi_obs);
  }
}
