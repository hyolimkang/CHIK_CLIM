data {
  int<lower=2> N;
  array[N] int<lower=0> cases;

  int<lower=1> Y;
  array[N] int<lower=1, upper=Y> year_id;

  int<lower=2> W;
  array[N] int<lower=1, upper=W> week_id;

  vector<lower=0>[N] pop;
  vector<lower=0>[N] births_weekly;

  real<lower=0> gamma_fixed;

  int<lower=2, upper=N> fit_start;

  // Prior hyperparameters for reporting probability rho
  real<lower=0> rho_alpha;
  real<lower=0> rho_beta;

  // Prior hyperparameters for initial susceptible fraction
  real<lower=0> s0_alpha;
  real<lower=0> s0_beta;
}

parameters {
  // Reporting probability to be inferred
  real<lower=0.001, upper=0.95> rho;

  // Baseline transmission level
  real mu_log_beta;

  // Cyclic 52-week seasonal beta on log scale
  vector[W] log_beta_season_raw;
  real<lower=0.01, upper=0.50> sigma_season;

  // Year-specific annual effective transmission multiplier
  vector<lower=-1.5, upper=1.5>[Y] delta_year_raw;

  // Year-specific effective external infectious pressure
  vector<lower=0.1, upper=2000>[Y] import_count_year;

  // Initial conditions
  real<lower=1, upper=100> I0;
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

  vector[Y] annual_infections_model;
  real mean_annual_infections_model;

  delta_year = delta_year_raw - mean(delta_year_raw);
  log_beta_season = log_beta_season_raw - mean(log_beta_season_raw);

  for (t in 1:N) {
    int y = year_id[t];
    int w = week_id[t];

    beta_t[t] = exp(
      mu_log_beta
      + delta_year[y]
      + log_beta_season[w]
    );
  }

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
    import_frac = import_count_year[y] / pop[t];

    infectious_pressure = I_frac[t - 1] + import_frac;

    lambda_t = beta_t[t] * infectious_pressure;
    p_inf = 1.0 - exp(-lambda_t);

    new_inf_frac[t] = S_frac[t - 1] * p_inf;

    p_rec = 1.0 - exp(-gamma_fixed);
    recov_frac = I_frac[t - 1] * p_rec;

    S_frac[t] = S_frac[t - 1] - new_inf_frac[t] + birth_frac;
    I_frac[t] = I_frac[t - 1] + new_inf_frac[t] - recov_frac;
    R_frac[t] = R_frac[t - 1] + recov_frac;

    // Expected confirmed cases using inferred rho
    mu_cases[t] = rho * new_inf_frac[t] * pop[t] + 1e-6;
  }

  annual_infections_model = rep_vector(0.0, Y);

  for (t in fit_start:N) {
    annual_infections_model[year_id[t]] += new_inf_frac[t] * pop[t];
  }

  mean_annual_infections_model = mean(annual_infections_model);
}

model {
  // Confirmed-case reporting prior informed by external infection estimate.
  // Suggested first try: rho_alpha = 8, rho_beta = 23, mean approx 0.26.
  rho ~ beta(rho_alpha, rho_beta);

  mu_log_beta ~ normal(log(1.0), 0.5);

  // Anchor seasonal raw levels
  log_beta_season_raw ~ normal(0, 1);

  // Cyclic random-walk prior for seasonal beta
  for (w in 2:W) {
    log_beta_season_raw[w] - log_beta_season_raw[w - 1] ~
      normal(0, sigma_season);
  }

  log_beta_season_raw[1] - log_beta_season_raw[W] ~
    normal(0, sigma_season);

  sigma_season ~ normal(0, 0.15);

  delta_year_raw ~ normal(0, 0.7);

  import_count_year ~ lognormal(log(30), 1.2);

  I0 ~ lognormal(log(5), 0.6);

  S0_frac ~ beta(s0_alpha, s0_beta);

  phi ~ lognormal(log(8), 0.8);

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
  vector[Y] annual_reported_mean;
  array[N] int cases_rep;

  for (t in 1:N) {
    R0_t[t] = beta_t[t] / gamma_fixed;
    Reff_t[t] = beta_t[t] / gamma_fixed * S_frac[t];

    S_out[t] = S_frac[t] * pop[t];
    I_out[t] = I_frac[t] * pop[t];
    R_out[t] = R_frac[t] * pop[t];
    new_inf[t] = new_inf_frac[t] * pop[t];
  }

  for (y in 1:Y) {
    annual_reported_mean[y] = 0;
  }

  cases_rep[1] = 0;

  for (t in 2:N) {
    if (t >= fit_start) {
      annual_reported_mean[year_id[t]] += mu_cases[t];
      cases_rep[t] = neg_binomial_2_rng(mu_cases[t], phi);
    } else {
      cases_rep[t] = 0;
    }
  }
}
