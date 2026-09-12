// Ceara episode-level early-phase renewal model.
//
// Each row is one pre-defined outbreak episode and each column is a week in
// its initial growth window. Lambda is calculated in R from the observed
// incidence history and the fixed generation-interval weights. The model
// estimates one effective reproduction number per episode, not a weekly R(t).

data {
  int<lower=1> K;                       // number of outbreak episodes
  int<lower=1> W;                       // weeks retained per episode
  array[K, W] int<lower=0> I_obs;       // reported weekly incidence
  matrix<lower=0>[K, W] Lambda;         // fixed renewal infectiousness
}

parameters {
  vector[K] log_Re;
  real<lower=0> phi;                    // shared NB2 overdispersion
}

transformed parameters {
  vector<lower=0>[K] Re = exp(log_Re);
  matrix<lower=0>[K, W] mu;

  for (k in 1:K) {
    for (w in 1:W) {
      mu[k, w] = Re[k] * Lambda[k, w];
    }
  }
}

model {
  // Weakly regularising early-epidemic prior: median 1.5, broad 95% range.
  log_Re ~ normal(log(1.5), 0.75);
  phi ~ lognormal(log(20), 1);

  for (k in 1:K) {
    for (w in 1:W) {
      I_obs[k, w] ~ neg_binomial_2(mu[k, w], phi);
    }
  }
}

generated quantities {
  array[K, W] int I_rep;

  for (k in 1:K) {
    for (w in 1:W) {
      I_rep[k, w] = neg_binomial_2_rng(mu[k, w], phi);
    }
  }
}
