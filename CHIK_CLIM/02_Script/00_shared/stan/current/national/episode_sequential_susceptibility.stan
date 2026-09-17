// Period-level episode-sequential susceptibility reconstruction.
//
// There is no weekly transmission process in this model. Each period has one
// cumulative force of infection (lambda), and S/U are propagated exactly over
// the observed period-level demographic accounting quantities.

data {
  int<lower=1> P;
  array[P] int<lower=0> C;
  array[P] int<lower=0, upper=1> is_epidemic;
  vector[P] year_centered;

  vector<lower=0>[P] N_start;
  vector<lower=0>[P] N_end;
  vector<lower=0>[P] births;
  vector<lower=0>[P] deaths;
  vector[P] reconciliation;

  real q_alpha_prior_mean;
  real<lower=0> q_alpha_prior_sd;
  real q_beta_prior_mean;
  real<lower=0> q_beta_prior_sd;
  real log_lambda_epi_prior_mean;
  real<lower=0> log_lambda_epi_prior_sd;
  real log_lambda_inter_prior_mean;
  real<lower=0> log_lambda_inter_prior_sd;
  real log_phi_prior_mean;
  real<lower=0> log_phi_prior_sd;
}

parameters {
  vector[P] log_lambda;
  real alpha_q;
  real beta_q;
  // This is an exact lognormal prior parameterization of phi, not a prior
  // change. Sampling the unconstrained log scale avoids a positive-scale
  // geometry that produced divergences in the pilot's first run.
  real log_phi;
}

transformed parameters {
  vector<lower=0>[P] lambda = exp(log_lambda);
  vector<lower=0, upper=1>[P] attack_prob;
  vector<lower=0, upper=1>[P] q;
  vector<lower=0>[P] S_start;
  vector<lower=0>[P] U_start;
  vector<lower=0>[P] latent_infections;
  vector<lower=0>[P] S_end;
  vector<lower=0>[P] U_end;
  vector<lower=0>[P] expected_cases;
  real<lower=0> phi = exp(log_phi);

  for (p in 1:P) {
    if (p == 1) {
      S_start[p] = N_start[p];
      U_start[p] = 0;
    } else {
      S_start[p] = S_end[p - 1];
      U_start[p] = U_end[p - 1];
    }

    attack_prob[p] = 1 - exp(-lambda[p]);
    latent_infections[p] = S_start[p] * attack_prob[p];
    q[p] = inv_logit(alpha_q + beta_q * year_centered[p]);
    expected_cases[p] = q[p] * latent_infections[p];

    {
      real S_after = S_start[p] - latent_infections[p];
      real U_after = U_start[p] + latent_infections[p];
      real susceptible_share = S_after / (S_after + U_after);

      // Births enter S. Deaths and population reconciliation are allocated
      // in proportion to the post-infection S/U composition, as in v2.2.
      S_end[p] = S_after + births[p] - deaths[p] * susceptible_share
                 + reconciliation[p] * susceptible_share;
      U_end[p] = U_after - deaths[p] * (1 - susceptible_share)
                 + reconciliation[p] * (1 - susceptible_share);
    }
  }
}

model {
  alpha_q ~ normal(q_alpha_prior_mean, q_alpha_prior_sd);
  beta_q ~ normal(q_beta_prior_mean, q_beta_prior_sd);
  log_phi ~ normal(log_phi_prior_mean, log_phi_prior_sd);

  for (p in 1:P) {
    if (is_epidemic[p] == 1) {
      log_lambda[p] ~ normal(log_lambda_epi_prior_mean, log_lambda_epi_prior_sd);
    } else {
      log_lambda[p] ~ normal(log_lambda_inter_prior_mean, log_lambda_inter_prior_sd);
    }
    C[p] ~ neg_binomial_2(expected_cases[p], phi);
  }
}

generated quantities {
  array[P] int C_rep;
  vector[P] S_prop_start;
  vector[P] S_prop_end;
  vector[P] accounting_error;

  for (p in 1:P) {
    C_rep[p] = neg_binomial_2_rng(expected_cases[p], phi);
    S_prop_start[p] = S_start[p] / N_start[p];
    S_prop_end[p] = S_end[p] / N_end[p];
    accounting_error[p] = S_end[p] + U_end[p] - N_end[p];
  }
}
