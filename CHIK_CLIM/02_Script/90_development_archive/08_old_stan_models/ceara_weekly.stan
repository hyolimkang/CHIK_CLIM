data {
  int<lower=2> N;                         // number of weeks
  array[N] int<lower=0> cases;            // weekly reported cases

  int<lower=1> Y;                         // number of years
  array[N] int<lower=1, upper=Y> year_id; // year index: 1,...,Y

  vector<lower=0>[N] pop;                 // weekly population
  vector<lower=0>[N] births_weekly;       // weekly births

  real<lower=0> gamma_fixed;              // weekly recovery rate / hazard

  // Fixed reporting rate.
  // Run sensitivity: rho_fixed = 0.05, 0.10, 0.20
  real<lower=0.001, upper=0.8> rho_fixed;

  // Start likelihood after initial transient.
  // Recommended first try: 20
  int<lower=2, upper=N> fit_start;
}

parameters {
  // Baseline seasonal transmission
  real mu_log_beta;
  real a1;
  real b1;
  real a2;
  real b2;

  // Bounded year-specific transmission deviation
  vector<lower=-1.5, upper=1.5>[Y] delta_year_raw;

  // Year-specific seasonal amplitude multiplier
  vector<lower=0.3, upper=3.0>[Y] amp_year;

  // Year-specific external infectious pressure.
  // This is effective infectious pressure, not observed importation count.
  vector<lower=0.1, upper=500>[Y] import_count_year;

  // Initial conditions
  real<lower=1, upper=100> I0;
  real<lower=0.6, upper=0.999> S0_frac;

  // Observation overdispersion
  real<lower=1> phi;
}

transformed parameters {
  vector[Y] delta_year;
  vector[N] beta_t;

  vector[N] S_frac;
  vector[N] I_frac;
  vector[N] R_frac;
  vector[N] new_inf_frac;
  vector[N] mu_cases;

  // Centre annual effects to reduce confounding with mu_log_beta.
  // Scale fixed at 1.0 to avoid unstable sigma_year estimation.
  delta_year = delta_year_raw - mean(delta_year_raw);

  // Time-varying beta: weak shared seasonality + year-specific amplitude + year effect
  for (t in 1:N) {
    real angle = 2.0 * pi() * (t - 1) / 52.18;
    real seasonal;
    int y = year_id[t];

    seasonal =
      a1 * sin(angle) + b1 * cos(angle)
      + a2 * sin(2.0 * angle) + b2 * cos(2.0 * angle);

    beta_t[t] = exp(
      mu_log_beta
      + delta_year[y]
      + amp_year[y] * seasonal
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
  mu_log_beta ~ normal(log(1.0), 0.35);

  // Shared seasonal shape.
  // Still weak, but year-specific amp_year can amplify it in outbreak years.
  a1 ~ normal(0, 0.10);
  b1 ~ normal(0, 0.10);
  a2 ~ normal(0, 0.05);
  b2 ~ normal(0, 0.05);

  // Year effects: allow outbreak years, but bounded.
  delta_year_raw ~ normal(0, 0.7);

  // Year-specific seasonal amplitude.
  amp_year ~ lognormal(log(1.0), 0.35);

  // Year-specific effective importation / outbreak pressure.
  // Most years low, outbreak years can be higher.
  import_count_year ~ lognormal(log(20), 1.0);

  // Initial conditions
  I0 ~ lognormal(log(5), 0.6);
  S0_frac ~ beta(60, 3);

  // Observation overdispersion.
  // Relaxed from lower=5 because previous fit hit the boundary.
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

    // Only accumulate diagnostics over the likelihood-fitted period
    if (t >= fit_start) {
      annual_infections[year_id[t]] += new_inf[t];
      annual_reported_mean[year_id[t]] += mu_cases[t];

      cases_rep[t] = neg_binomial_2_rng(mu_cases[t], phi);
    } else {
      cases_rep[t] = 0;
    }
  }
}

