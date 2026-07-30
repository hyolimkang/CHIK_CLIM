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

  real<lower=0> s0_alpha;
  real<lower=0> s0_beta;

  real<lower=0> rho_alpha;
  real<lower=0> rho_beta;

  real import_prior_logmedian;
  real<lower=0> import_prior_sd;

  // Raw serology data
  int<lower=0> J_sero;
  array[J_sero] int<lower=1, upper=N> sero_t;
  array[J_sero] int<lower=0> sero_pos;
  array[J_sero] int<lower=0> sero_n;

  // Test sensitivity / specificity
  real<lower=0, upper=1> sero_sens;
  real<lower=0, upper=1> sero_spec;

  // Outbreak pulse structure
  int<lower=1> K_pulse;
  array[N] int<lower=0, upper=K_pulse> pulse_id;
  real<lower=0> pulse_prior_sd;
}

parameters {
  // Reporting probability
  real<lower=0.01, upper=0.45> rho;

  // Baseline transmission
  real<lower=-2, upper=1.5> mu_log_beta;

  // Cyclic seasonal beta
  vector[W] log_beta_season_raw;
  real<lower=0.01, upper=0.35> sigma_season;

  // Annual transmission deviation
  vector<lower=-1.2, upper=1.2>[Y] delta_year_raw;

  // v12c-lite import structure
  vector<lower=0.01, upper=75>[Y] import_count_year;

  // Moderate outbreak pulse on log beta
  vector<lower=0, upper=1.2>[K_pulse] pulse_amp;

  // Initial conditions
  real<lower=1, upper=100> I0;
  real<lower=0.05, upper=0.999> S0_frac;

  // Observation overdispersion
  real<lower=1> phi_cases;
}

transformed parameters {
  vector[Y] delta_year;
  vector[W] log_beta_season;
  vector[N] beta_t;
  vector[N] pulse_effect_t;

  vector[N] S_frac;
  vector[N] I_frac;
  vector[N] R_frac;
  vector[N] new_inf_frac;
  vector[N] mu_cases;

  vector[N] lambda_out;
  vector[Y] annual_foi;
  real mean_annual_foi_model;

  real cumulative_attack_model;

  // Centre annual and seasonal effects
  delta_year = delta_year_raw - mean(delta_year_raw);
  log_beta_season = log_beta_season_raw - mean(log_beta_season_raw);

  for (t in 1:N) {
    int y = year_id[t];
    int w = week_id[t];

    if (pulse_id[t] == 0) {
      pulse_effect_t[t] = 0.0;
    } else {
      pulse_effect_t[t] = pulse_amp[pulse_id[t]];
    }

    beta_t[t] = exp(
      mu_log_beta
      + delta_year[y]
      + log_beta_season[w]
      + pulse_effect_t[t]
    );
  }

  // Initial states
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

    import_frac = import_count_year[y] / pop[t];

    infectious_pressure = I_frac[t - 1] + import_frac;

    lambda_t = beta_t[t] * infectious_pressure;
    lambda_out[t] = lambda_t;

    p_inf = 1.0 - exp(-lambda_t);
    new_inf_frac[t] = S_frac[t - 1] * p_inf;

    p_rec = 1.0 - exp(-gamma_fixed);
    recov_frac = I_frac[t - 1] * p_rec;

    S_frac[t] = S_frac[t - 1] - new_inf_frac[t] + birth_frac;
    I_frac[t] = I_frac[t - 1] + new_inf_frac[t] - recov_frac;
    R_frac[t] = R_frac[t - 1] + recov_frac;

    mu_cases[t] = rho * new_inf_frac[t] * pop[t] + 1e-6;
  }

  annual_foi = rep_vector(0.0, Y);

  for (t in fit_start:N) {
    annual_foi[year_id[t]] += lambda_out[t];
  }

  mean_annual_foi_model = mean(annual_foi);

  cumulative_attack_model = 0.0;
  for (t in fit_start:N) {
    cumulative_attack_model += new_inf_frac[t];
  }
}

model {
  rho ~ beta(rho_alpha, rho_beta);

  mu_log_beta ~ normal(0, 0.4);

  log_beta_season_raw ~ normal(0, 1);

  for (w in 2:W) {
    log_beta_season_raw[w] - log_beta_season_raw[w - 1] ~
      normal(0, sigma_season);
  }

  log_beta_season_raw[1] - log_beta_season_raw[W] ~
    normal(0, sigma_season);

  sigma_season ~ normal(0, 0.12);

  delta_year_raw ~ normal(0, 0.5);

  import_count_year ~ lognormal(import_prior_logmedian, import_prior_sd);

  pulse_amp ~ normal(0, pulse_prior_sd);

  I0 ~ lognormal(log(5), 0.6);
  S0_frac ~ beta(s0_alpha, s0_beta);

  phi_cases ~ lognormal(log(10), 0.5);

  // Reported case likelihood
  for (t in fit_start:N) {
    cases[t] ~ neg_binomial_2(mu_cases[t], phi_cases);
  }

  // Raw serology likelihood
  // v13a simplification:
  // serology observes total immune / seropositive fraction at survey time.
  for (j in 1:J_sero) {
    real p_true;
    real p_obs;

    p_true = 1.0 - S_frac[sero_t[j]];
    p_true = fmin(fmax(p_true, 1e-6), 1.0 - 1e-6);

    p_obs =
      sero_sens * p_true
      + (1.0 - sero_spec) * (1.0 - p_true);

    p_obs = fmin(fmax(p_obs, 1e-6), 1.0 - 1e-6);

    sero_pos[j] ~ binomial(sero_n[j], p_obs);
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
  vector[Y] annual_foi_out;

  vector[K_pulse] pulse_multiplier;

  real mean_annual_infections_model;
  real mean_annual_reported_model;
  real cumulative_attack_out;

  vector[N] log_lik_cases;
  vector[J_sero] log_lik_sero;

  array[N] int cases_rep;

  for (k in 1:K_pulse) {
    pulse_multiplier[k] = exp(pulse_amp[k]);
  }

  for (y in 1:Y) {
    annual_infections[y] = 0.0;
    annual_reported_mean[y] = 0.0;
    annual_foi_out[y] = annual_foi[y];
  }

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

      log_lik_cases[t] = neg_binomial_2_lpmf(
        cases[t] | mu_cases[t], phi_cases
      );

      cases_rep[t] = neg_binomial_2_rng(mu_cases[t], phi_cases);
    } else {
      log_lik_cases[t] = 0.0;
      cases_rep[t] = 0;
    }
  }

  for (j in 1:J_sero) {
    real p_true;
    real p_obs;

    p_true = 1.0 - S_frac[sero_t[j]];
    p_true = fmin(fmax(p_true, 1e-6), 1.0 - 1e-6);

    p_obs =
      sero_sens * p_true
      + (1.0 - sero_spec) * (1.0 - p_true);

    p_obs = fmin(fmax(p_obs, 1e-6), 1.0 - 1e-6);

    log_lik_sero[j] = binomial_lpmf(sero_pos[j] | sero_n[j], p_obs);
  }

  mean_annual_infections_model = mean(annual_infections);
  mean_annual_reported_model = mean(annual_reported_mean);

  cumulative_attack_out = cumulative_attack_model;
}
