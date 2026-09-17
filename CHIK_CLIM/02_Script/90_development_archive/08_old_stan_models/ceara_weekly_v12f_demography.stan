data {
  int<lower=2> N;
  int<lower=0> cases[N];

  int<lower=1> Y;
  int<lower=1, upper=Y> year_id[N];

  int<lower=2> W;
  int<lower=1, upper=W> week_id[N];

  vector<lower=0>[N] pop;
  vector<lower=0>[N] births_weekly;

  real<lower=0> gamma_fixed;

  int<lower=2, upper=N> fit_start;

  real<lower=0> s0_alpha;
  real<lower=0> s0_beta;

  real<lower=0> rho_alpha;
  real<lower=0> rho_beta;

  // External cumulative attack constraint derived from serology-informed FOI.
  real<lower=0, upper=1> cum_attack_ext_mean;
  real<lower=0> cum_attack_ext_sdlogit;

  real import_prior_logmedian;
  real<lower=0> import_prior_sd;

  int<lower=0, upper=1> use_cum_attack_prior;
}

parameters {
  // Reporting probability.
  // Upper bound reflects external infection-scale plausibility.
  real<lower=0.01, upper=0.45> rho;

  real<lower=-2, upper=1.5> mu_log_beta;

  vector[W] log_beta_season_raw;
  real<lower=0.01, upper=0.35> sigma_season;

  vector<lower=-1.2, upper=1.2>[Y] delta_year_raw;

  // Annual imported infectious equivalents.
  vector<lower=0.01, upper=75>[Y] import_count_year;

  real<lower=1, upper=100> I0;
  real<lower=0.05, upper=0.999> S0_frac;

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
  vector[N] new_inf_count;
  vector[N] mu_cases;

  vector[N] lambda_out;
  vector[Y] annual_foi;
  real mean_annual_foi_model;

  real cumulative_attack_model;
  real cumulative_attack_model_clamped;

  delta_year = delta_year_raw - mean(delta_year_raw);
  log_beta_season = log_beta_season_raw - mean(log_beta_season_raw);

  for (t in 1:N) {
    int y;
    int w;

    y = year_id[t];
    w = week_id[t];

    beta_t[t] = exp(
      mu_log_beta
      + delta_year[y]
      + log_beta_season[w]
    );
  }

  S_frac[1] = S0_frac;
  I_frac[1] = I0 / pop[1];
  R_frac[1] = fmax(1.0 - S_frac[1] - I_frac[1], 0.0);

  new_inf_count[1] = 1e-12;
  new_inf_frac[1] = new_inf_count[1] / pop[1];
  mu_cases[1] = 1e-6;
  lambda_out[1] = 0.0;

  for (t in 2:N) {
    int y;
    real import_count_week;
    real import_frac;
    real infectious_pressure;
    real lambda_t;
    real p_inf;
    real p_rec;
    real S_prev_count;
    real I_prev_count;
    real R_prev_count;
    real recov_count;
    real non_birth_demo_count;
    real S_next_count;
    real I_next_count;
    real R_next_count;
    real total_next_count;

    y = year_id[t];

    import_count_week = import_count_year[y] / 52.0;
    import_frac = import_count_week / pop[t - 1];

    S_prev_count = S_frac[t - 1] * pop[t - 1];
    I_prev_count = I_frac[t - 1] * pop[t - 1];
    R_prev_count = R_frac[t - 1] * pop[t - 1];

    infectious_pressure = I_frac[t - 1] + import_frac;

    lambda_t = beta_t[t] * infectious_pressure;
    lambda_out[t] = lambda_t;

    p_inf = 1.0 - exp(-lambda_t);
    new_inf_count[t] = S_prev_count * p_inf;
    new_inf_frac[t] = new_inf_count[t] / pop[t - 1];

    p_rec = 1.0 - exp(-gamma_fixed);
    recov_count = I_prev_count * p_rec;

    S_next_count = S_prev_count - new_inf_count[t] + births_weekly[t];
    I_next_count = I_prev_count + new_inf_count[t] - recov_count;
    R_next_count = R_prev_count + recov_count;

    // IBGE population mixes annual estimates, Census 2022, and interpolation.
    // Treat the residual population change as non-birth demographic adjustment
    // and distribute it proportionally across compartments, instead of letting
    // the denominator alone create artificial jumps in S(t) / N(t).
    non_birth_demo_count = pop[t] - pop[t - 1] - births_weekly[t];

    S_next_count += non_birth_demo_count * S_frac[t - 1];
    I_next_count += non_birth_demo_count * I_frac[t - 1];
    R_next_count += non_birth_demo_count * R_frac[t - 1];

    S_next_count = fmax(S_next_count, 1e-9);
    I_next_count = fmax(I_next_count, 1e-9);
    R_next_count = fmax(R_next_count, 1e-9);
    total_next_count = S_next_count + I_next_count + R_next_count;

    S_frac[t] = S_next_count / total_next_count;
    I_frac[t] = I_next_count / total_next_count;
    R_frac[t] = R_next_count / total_next_count;

    mu_cases[t] = rho * new_inf_count[t] + 1e-6;
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

  cumulative_attack_model_clamped =
    fmin(fmax(cumulative_attack_model, 1e-6), 1.0 - 1e-6);
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

  I0 ~ lognormal(log(5), 0.6);
  S0_frac ~ beta(s0_alpha, s0_beta);

  phi_cases ~ lognormal(log(10), 0.5);

  for (t in fit_start:N) {
    cases[t] ~ neg_binomial_2(mu_cases[t], phi_cases);
  }

  if (use_cum_attack_prior == 1) {
    target += normal_lpdf(
      logit(cumulative_attack_model_clamped) |
      logit(cum_attack_ext_mean),
      cum_attack_ext_sdlogit
    );
  }
}

generated quantities {
  vector[N] R0_t;
  vector[N] Reff_t;

  vector[N] S_out;
  vector[N] I_out;
  vector[N] R_out;
  vector[N] new_inf;

  vector[N] local_infectious_frac_out;
  vector[N] import_frac_out;
  vector[N] infectious_pressure_out;

  vector[Y] annual_infections;
  vector[Y] annual_reported_mean;
  vector[Y] annual_foi_out;

  real mean_annual_infections_model;
  real mean_annual_reported_model;

  real cumulative_attack_out;
  real cum_attack_ext_log_lik;

  vector[N] log_lik_cases;

  int cases_rep[N];

  for (y in 1:Y) {
    annual_infections[y] = 0.0;
    annual_reported_mean[y] = 0.0;
    annual_foi_out[y] = annual_foi[y];
  }

  for (t in 1:N) {
    R0_t[t] = beta_t[t] / (1.0 - exp(-gamma_fixed));
    Reff_t[t] = beta_t[t] / (1.0 - exp(-gamma_fixed)) * S_frac[t];

    S_out[t] = S_frac[t] * pop[t];
    I_out[t] = I_frac[t] * pop[t];
    R_out[t] = R_frac[t] * pop[t];
    new_inf[t] = new_inf_count[t];

    local_infectious_frac_out[t] = I_frac[t];
    import_frac_out[t] = (import_count_year[year_id[t]] / 52.0) / pop[t];
    infectious_pressure_out[t] =
      local_infectious_frac_out[t] + import_frac_out[t];

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

  mean_annual_infections_model = mean(annual_infections);
  mean_annual_reported_model = mean(annual_reported_mean);

  cumulative_attack_out = cumulative_attack_model;

  if (use_cum_attack_prior == 1) {
    cum_attack_ext_log_lik = normal_lpdf(
      logit(cumulative_attack_model_clamped) |
      logit(cum_attack_ext_mean),
      cum_attack_ext_sdlogit
    );
  } else {
    cum_attack_ext_log_lik = 0.0;
  }
}
