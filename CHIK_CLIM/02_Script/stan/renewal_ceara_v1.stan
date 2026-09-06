// renewal_ceara_v1.stan
//
// Stage 1 of the renewal-equation transmission model (simplest version):
//   - constant R (no time variation, no climate link yet)
//   - no susceptible depletion (S_t = 1 always)
//   - reporting fraction rho fixed at 1 (I_t is defined as "expected
//     reported infections", not true infections — separating the two
//     needs external calibration data like seroprevalence, deferred to a
//     later stage)
//
// Latent infections I_t cannot be a literal integer NegBin-distributed
// Stan *parameter* (HMC needs continuous parameters). Following the same
// approach EpiNow2 uses: I_t is continuous, with lognormal multiplicative
// process noise around the renewal-equation mean:
//   I_t = R * sum_tau w[tau] * I_{t-tau} * exp(eta_t),  eta_t ~ Normal(0, sigma)
// The observation layer C_t ~ NegBinomial2(I_t, phi) is a proper discrete
// likelihood, which is fine since C_t is data, not a parameter.

data {
  int<lower=1> N;                    // number of weeks
  int<lower=1> G;                    // max generation-interval lag (weeks)
  int<lower=1> seed_weeks;           // initial weeks seeded freely (>= G)
  array[N] int<lower=0> C;           // observed weekly confirmed cases
  vector[G] w;                       // generation-interval weights, sums to 1
}

transformed data {
  real log_I_seed_prior_mean;
  {
    real acc = 0;
    for (t in 1:seed_weeks) acc += C[t];
    log_I_seed_prior_mean = log(acc / seed_weeks + 1);
  }
}

parameters {
  real<lower=0> R;                   // constant reproduction number
  vector[seed_weeks] log_I_seed;     // free log-infections for the seed period
  real<lower=0> sigma_process;       // SD of lognormal process noise
  real<lower=0> phi_obs;             // NegBin overdispersion (observation)
  vector[N - seed_weeks] eta;        // process-noise innovations post-seeding
}

transformed parameters {
  vector[N] I;                       // latent expected-reported-infections series

  for (t in 1:seed_weeks) {
    I[t] = exp(log_I_seed[t]);
  }

  for (t in (seed_weeks + 1):N) {
    real conv = 0;
    for (tau in 1:G) {
      conv += w[tau] * I[t - tau];
    }
    I[t] = R * conv * exp(sigma_process * eta[t - seed_weeks]);
  }
}

model {
  // ---- priors ----
  R ~ lognormal(log(1.0), 0.5);
  sigma_process ~ normal(0, 0.5) T[0, ];
  phi_obs ~ gamma(2, 0.1);
  log_I_seed ~ normal(log_I_seed_prior_mean, 1);
  eta ~ std_normal();

  // ---- observation model ----
  for (t in 1:N) {
    C[t] ~ neg_binomial_2(I[t] + 1e-6, phi_obs);
  }
}

generated quantities {
  array[N] int C_pred;
  for (t in 1:N) {
    C_pred[t] = neg_binomial_2_rng(I[t] + 1e-6, phi_obs);
  }
}
