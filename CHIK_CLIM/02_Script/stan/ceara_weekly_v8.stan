data {
  int<lower=2> N;                         // number of weeks
  array[N] int<lower=0> cases;            // weekly reported cases

  int<lower=1> Y;                         // number of years
  array[N] int<lower=1, upper=Y> year_id; // year index: 1,...,Y

  int<lower=2> W;                         // number of seasonal weeks, usually 52
  array[N] int<lower=1, upper=W> week_id; // week of year: 1,...,52

  vector<lower=0>[N] pop;                 // weekly population
  vector<lower=0>[N] births_weekly;       // weekly births

  real<lower=0> gamma_fixed;              // weekly recovery rate / hazard

  // Fixed reporting rate.
  // Sensitivity: rho_fixed = 0.05, 0.10, 0.20
  real<lower=0.001, upper=0.8> rho_fixed;

  // Start likelihood after initial transient.
  int<lower=2, upper=N> fit_start;

  // Prior hyperparameters for initial susceptible fraction.
  // Example: beta(8, 2) = mean 0.8, flexible but still high.
  real<lower=0> s0_alpha;
  real<lower=0> s0_beta;
}

parameters {
  // Baseline transmission level
  real mu_log_beta;

  // Cyclic 52-week seasonal beta on log scale
  vector[W] log_beta_season_raw;
  real<lower=0.01, upper=0.50> sigma_season;

  // Year-specific transmission deviation
  vector<lower=-1.5, upper=1.5>[Y] delta_year_raw;

  // Year-specific external infectious pressure.
  // Effective forcing, not observed imported cases.
  vector<lower=0.1, upper=2000>[Y] import_count_year;

  // Initial conditions
  real<lower=1, upper=100> I0;

  // More flexible than previous lower=0.6.
  // Still avoids exactly 0, which can create numerical/pathological behaviour.
  real<lower=0.05, upper=0.999> S0_frac;

  // Observation overdispersion
  real<lower=1> phi;
}

transformed parameters {
  vector[Y] delta_year;
  vector[W] log_beta_season;
  vector[N] beta_t;

  vector[N] S_frac;
  vector[N] I_frac;
  vector[N] R_frac;
  vector[N] new_inf_frac;
  vector[N] mu_cases;

  // Centre annual effects to reduce confounding with mu_log_beta.
  delta_year = delta_year_raw - mean(delta_year_raw);

  // Centre seasonal beta to reduce confounding with mu_log_beta.
  log_beta_season = log_beta_season_raw - mean(log_beta_season_raw);

  // Time-varying beta:
  // baseline + cyclic seasonal week effect + annual multiplier
  for (t in 1:N) {
    int y = year_id[t];
    int w = week_id[t];

    beta_t[t] = exp(
      mu_log_beta
      + delta_year[y]
      + log_beta_season[w]
    );
  }

  // Initial states as population fractions
  S_frac[1] = S0_frac;
  I_frac[1] = I0 / pop[1];
  R_frac[1] = fmax(1.0 - S_frac[1] - I_frac[1], 0.0);

  new_inf_frac[1] = 1e-12;
  mu_cases[1] = 1e-6;

  for (t in 2:N) {
    int y = year_id[t];

    real birth_frac;
    real import_frac;
    real infectious_pressure;
    real lambda_t;
    real p_inf;
    real p_rec;
    real recov_frac;

    birth_frac = births_weekly[t] / pop[t];

    // Year-specific external infectious pressure as population fraction
    import_frac = import_count_year[y] / pop[t];

    // Local infectious pressure + effective external pressure
    infectious_pressure = I_frac[t - 1] + import_frac;

    // Probability-based incidence
    lambda_t = beta_t[t] * infectious_pressure;
    p_inf = 1.0 - exp(-lambda_t);

    new_inf_frac[t] = S_frac[t - 1] * p_inf;

    // Recovery over one weekly interval
    p_rec = 1.0 - exp(-gamma_fixed);
    recov_frac = I_frac[t - 1] * p_rec;

    // SIR update
    S_frac[t] = S_frac[t - 1] - new_inf_frac[t] + birth_frac;
    I_frac[t] = I_frac[t - 1] + new_inf_frac[t] - recov_frac;
    R_frac[t] = R_frac[t - 1] + recov_frac;

    // Expected reported cases
    mu_cases[t] = rho_fixed * new_inf_frac[t] * pop[t] + 1e-6;
  }
}

model {
  // Baseline beta
  mu_log_beta ~ normal(log(1.0), 0.5);

  // Important anchor:
  // The random-walk prior only constrains differences between weeks.
  // This weakly anchors the overall raw seasonal level and helps geometry.
  log_beta_season_raw ~ normal(0, 1);

  // Cyclic random-walk prior for 52-week seasonal beta.
  // Adjacent weeks should have similar log seasonal beta.
  for (w in 2:W) {
    log_beta_season_raw[w] - log_beta_season_raw[w - 1] ~
      normal(0, sigma_season);
  }

  // Cyclic closure: week 1 should connect smoothly to week W.
  log_beta_season_raw[1] - log_beta_season_raw[W] ~
    normal(0, sigma_season);

  // Smoothness prior.
  // Smaller sigma_season = smoother seasonal beta.
  sigma_season ~ normal(0, 0.15);

  // Year effects: annual effective transmission multiplier.
  delta_year_raw ~ normal(0, 0.7);

  // Year-specific effective importation / outbreak pressure.
  import_count_year ~ lognormal(log(30), 1.2);

  // Initial conditions
  I0 ~ lognormal(log(5), 0.6);

  // Flexible S0 prior supplied from R.
  // Example: beta(8, 2), beta(5, 2), beta(20, 2) sensitivity.
  S0_frac ~ beta(s0_alpha, s0_beta);

  // Observation overdispersion
  phi ~ lognormal(log(8), 0.8);

  // Likelihood
  for (t in fit_start:N) {
    cases[t] ~ neg_binomial_2(mu_cases[t], phi);
  }
}

generated quantities {
  vector[N] R0_t;
  vector[N] Reff_t;
  vector[N] S_out;
  vector[N] I_out;
  vector[N] R_out;
  vector[N] new_inf;
  vector[Y] annual_infections;
  vector[Y] annual_reported_mean;
  array[N] int cases_rep;

  for (y in 1:Y) {
    annual_infections[y] = 0;
    annual_reported_mean[y] = 0;
  }

  cases_rep[1] = 0;

  for (t in 1:N) {
    R0_t[t] = beta_t[t] / gamma_fixed;
    Reff_t[t] = beta_t[t] / gamma_fixed * S_frac[t];

    S_out[t] = S_frac[t] * pop[t];
    I_out[t] = I_frac[t] * pop[t];
    R_out[t] = R_frac[t] * pop[t];
    new_inf[t] = new_inf_frac[t] * pop[t];

    if (t >= fit_start) {
      annual_infections[year_id[t]] += new_inf[t];
      annual_reported_mean[year_id[t]] += mu_cases[t];

      cases_rep[t] = neg_binomial_2_rng(mu_cases[t], phi);
    } else {
      cases_rep[t] = 0;
    }
  }
}
