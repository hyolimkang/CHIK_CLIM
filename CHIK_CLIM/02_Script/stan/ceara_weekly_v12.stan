data {
  int<lower=2> N;                         // number of weeks
  array[N] int<lower=0> cases;            // weekly confirmed/reported cases

  int<lower=1> Y;                         // number of years
  array[N] int<lower=1, upper=Y> year_id; // year index: 1,...,Y

  int<lower=2> W;                         // number of seasonal weeks, usually 52
  array[N] int<lower=1, upper=W> week_id; // week of year: 1,...,52

  vector<lower=0>[N] pop;                 // weekly population
  vector<lower=0>[N] births_weekly;       // weekly births

  real<lower=0> gamma_fixed;              // weekly recovery hazard

  // Start likelihood after initial transient
  int<lower=2, upper=N> fit_start;

  // Initial susceptible fraction prior
  real<lower=0> s0_alpha;
  real<lower=0> s0_beta;

  // Reporting probability prior
  real<lower=0> rho_alpha;
  real<lower=0> rho_beta;

  // Ceará-specific external annual FOI
  real<lower=0> foi_ext_mean;
  real<lower=0> foi_ext_sdlog;

  // Prior for annual importation / external pressure
  real import_prior_logmedian;
  real<lower=0> import_prior_sd;
}

parameters {
  // Reporting probability: confirmed cases / true infections
  real<lower=0.01, upper=0.6> rho;

  // Baseline transmission level
  // Regularised to prevent high-beta explosive modes.
  real<lower=-2, upper=1.5> mu_log_beta;

  // Cyclic 52-week seasonal beta on log scale
  vector[W] log_beta_season_raw;
  real<lower=0.01, upper=0.35> sigma_season;

  // Year-specific transmission deviation
  vector<lower=-1.2, upper=1.2>[Y] delta_year_raw;

  // Year-specific effective external infectious pressure
  // Strongly regularised so it does not become a hidden outbreak generator.
  vector<lower=0.01, upper=75>[Y] import_count_year;

  // Initial infectious individuals
  real<lower=1, upper=100> I0;

  // Initial susceptible fraction
  real<lower=0.05, upper=0.999> S0_frac;

  // Observation overdispersion for reported cases
  real<lower=1> phi_cases;
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

  vector[N] lambda_out;
  vector[Y] annual_foi;
  real mean_annual_foi_model;

  // Centre annual effects to reduce confounding with mu_log_beta.
  delta_year = delta_year_raw - mean(delta_year_raw);

  // Centre seasonal beta to reduce confounding with mu_log_beta.
  log_beta_season = log_beta_season_raw - mean(log_beta_season_raw);

  // Time-varying beta.
  for (t in 1:N) {
    int y = year_id[t];
    int w = week_id[t];

    beta_t[t] = exp(
      mu_log_beta
      + delta_year[y]
      + log_beta_season[w]
    );
  }

  // Initial states as population fractions.
  S_frac[1] = S0_frac;
  I_frac[1] = I0 / pop[1];
  R_frac[1] = fmax(1.0 - S_frac[1] - I_frac[1], 0.0);

  new_inf_frac[1] = 1e-12;
  mu_cases[1] = 1e-6;
  lambda_out[1] = 0.0;

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

    // Year-specific external infectious pressure as population fraction.
    import_frac = import_count_year[y] / pop[t];

    // Local infectious pressure + effective external pressure.
    infectious_pressure = I_frac[t - 1] + import_frac;

    // Weekly FOI hazard.
    lambda_t = beta_t[t] * infectious_pressure;
    lambda_out[t] = lambda_t;

    // Probability-based incidence.
    p_inf = 1.0 - exp(-lambda_t);
    new_inf_frac[t] = S_frac[t - 1] * p_inf;

    // Recovery over one weekly interval.
    p_rec = 1.0 - exp(-gamma_fixed);
    recov_frac = I_frac[t - 1] * p_rec;

    // SIR update.
    S_frac[t] = S_frac[t - 1] - new_inf_frac[t] + birth_frac;
    I_frac[t] = I_frac[t - 1] + new_inf_frac[t] - recov_frac;
    R_frac[t] = R_frac[t - 1] + recov_frac;

    // Expected confirmed/reported cases.
    mu_cases[t] = rho * new_inf_frac[t] * pop[t] + 1e-6;
  }

  // Annual cumulative FOI: sum of weekly FOI hazards within each year.
  annual_foi = rep_vector(0.0, Y);

  for (t in fit_start:N) {
    annual_foi[year_id[t]] += lambda_out[t];
  }

  // Simple prototype: average across all model years.
  mean_annual_foi_model = mean(annual_foi);
}

model {
  // Reporting probability.
  rho ~ beta(rho_alpha, rho_beta);

  // Baseline beta.
  mu_log_beta ~ normal(0, 0.4);

  // Anchor raw seasonal beta to improve geometry.
  log_beta_season_raw ~ normal(0, 1);

  // Cyclic random-walk prior for seasonal beta.
  for (w in 2:W) {
    log_beta_season_raw[w] - log_beta_season_raw[w - 1] ~
      normal(0, sigma_season);
  }

  // Cyclic closure.
  log_beta_season_raw[1] - log_beta_season_raw[W] ~
    normal(0, sigma_season);

  // Smoother seasonal beta than v12 simple.
  sigma_season ~ normal(0, 0.12);

  // Year effects: slightly tightened.
  delta_year_raw ~ normal(0, 0.5);

  // Effective importation / external pressure.
  import_count_year ~ lognormal(import_prior_logmedian, import_prior_sd);

  // Initial conditions.
  I0 ~ lognormal(log(5), 0.6);
  S0_frac ~ beta(s0_alpha, s0_beta);

  // Observation overdispersion.
  phi_cases ~ lognormal(log(10), 0.5);

  // Reported case likelihood.
  for (t in fit_start:N) {
    cases[t] ~ neg_binomial_2(mu_cases[t], phi_cases);
  }

  // External Ceará serology-informed annual FOI constraint.
  target += normal_lpdf(
    log(foi_ext_mean) |
    log(mean_annual_foi_model + 1e-9),
    foi_ext_sdlog
  );
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
  vector[Y] annual_foi_out;

  real mean_annual_infections_model;
  real mean_annual_reported_model;
  real cumulative_attack_model;

  real foi_ext_log_lik;
  vector[N] log_lik_cases;

  array[N] int cases_rep;

  for (y in 1:Y) {
    annual_infections[y] = 0.0;
    annual_reported_mean[y] = 0.0;
    annual_foi_out[y] = annual_foi[y];
  }

  cumulative_attack_model = 0.0;

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
      cumulative_attack_model += new_inf_frac[t];

      log_lik_cases[t] = neg_binomial_2_lpmf(
        cases[t] | mu_cases[t], phi_cases
      );

      cases_rep[t] = neg_binomial_2_rng(mu_cases[t], phi_cases);
    } else {
      log_lik_cases[t] = 0.0;
      cases_rep[t] = 0;
    }
  }

  mean_annual_infections_model = mean(annual_infections);
  mean_annual_reported_model = mean(annual_reported_mean);

  foi_ext_log_lik = normal_lpdf(
    log(foi_ext_mean) |
    log(mean_annual_foi_model + 1e-9),
    foi_ext_sdlog
  );
}
