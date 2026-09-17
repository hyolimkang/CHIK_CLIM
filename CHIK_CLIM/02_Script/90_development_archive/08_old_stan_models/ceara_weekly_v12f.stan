// ---------------------------------------------------------------------------
// ceara_weekly_v12f.stan
//
// v12e + explicit demographic flow (births in, NATURAL deaths out).
//
// What changed vs v12e
// --------------------
// v12e added births to S but had NO deaths, so the compartment counts summed
// to pop[t-1] + births_weekly[t] while being divided by pop[t]. That broke
// population conservation.
//
// v12f adds mechanistic demography:
//   1. epidemic transitions (infection, recovery) conserve the total,
//   2. births_weekly[t] enter the susceptible compartment,
//   3. a NATURAL per-capita mortality rate (mu_death_annual, supplied as data)
//      removes people from S, I and R every week.
//
// Because deaths are driven by a real (rough) Brazilian annual death rate and
// births by the observed birth series, the modelled population evolves on its
// own: pop_model[t] = (pop_model[t-1] + births[t]) * (1 - p_death). It is NOT
// pinned to the observed pop[t]; the gap (pop_balance_resid) is reported so the
// demographic drift can be inspected. Frequency-dependent terms (prevalence and
// imports) use the modelled population pop_model so the accounting is internally
// consistent.
//
// mu_death_annual is a rough natural mortality rate, e.g. Brazil crude death
// rate ~ 6.6 / 1000 / year (0.0066) or 1 / life_expectancy (~1/76 ~ 0.0132).
// It is passed as data so it can be changed without recompiling.
//
// Data contract adds ONE field vs v12e: mu_death_annual. The existing R data
// list only needs that single extra entry; v12e ignores it.
// ---------------------------------------------------------------------------
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

  // Rough natural mortality rate (per capita, per YEAR). Converted to a weekly
  // survival probability internally. e.g. 0.0066 (crude death rate ~6.6/1000).
  real<lower=0> mu_death_annual;

  int<lower=2, upper=N> fit_start;

  real<lower=0> s0_alpha;
  real<lower=0> s0_beta;

  real<lower=0> rho_alpha;
  real<lower=0> rho_beta;

  // External cumulative attack constraint derived from serology-informed FOI
  real<lower=0, upper=1> cum_attack_ext_mean;
  real<lower=0> cum_attack_ext_sdlogit;

  real import_prior_logmedian;
  real<lower=0> import_prior_sd;

  int<lower=0, upper=1> use_cum_attack_prior;
}

transformed data {
  // Weekly natural mortality probability from the annual rate.
  real p_death = 1.0 - exp(-mu_death_annual / 52.0);
}

parameters {
  // Reporting probability.
  // Upper bound now reflects external infection-scale plausibility.
  real<lower=0.01, upper=0.45> rho; // 10-year constant reporting is a strong assumption

  real<lower=-2, upper=1.5> mu_log_beta; // beta = exp(mu_log_beta)

  vector[W] log_beta_season_raw;
  real<lower=0.01, upper=0.35> sigma_season;

  vector<lower=-2.0, upper=2.0>[Y] delta_year_raw;

  // v12c-lite import structure
  vector<lower=0.01, upper=75>[Y] import_count_year;

  // beta_t = exp(mu_log_beta + delta_year[y] + log_beta_season[w])
  // weekly transmission rate = baseline transmission * year-specific factor * seasonal factor

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

  // Demographic flow (counts)
  vector[N] pop_model;          // modelled population, evolves by births - deaths
  vector[N] births_count;
  vector[N] deaths_count;       // total natural deaths this week
  vector[N] deaths_S;
  vector[N] deaths_I;
  vector[N] deaths_R;

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

  // Initial state: anchor the modelled population to the first observed value.
  pop_model[1] = pop[1];
  S_frac[1] = S0_frac;
  I_frac[1] = I0 / pop[1];
  R_frac[1] = fmax(1.0 - S_frac[1] - I_frac[1], 0.0);

  new_inf_count[1] = 1e-12;
  new_inf_frac[1] = new_inf_count[1] / pop[1];
  mu_cases[1] = 1e-6;
  lambda_out[1] = 0.0;

  births_count[1] = 0.0;
  deaths_count[1] = 0.0;
  deaths_S[1] = 0.0;
  deaths_I[1] = 0.0;
  deaths_R[1] = 0.0;

  for (t in 2:N) {
    int y;
    real N_prev;
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

    // Post-epidemic, post-birth counts (before natural deaths)
    real S_post;
    real I_post;
    real R_post;

    // Post-death counts
    real S_new;
    real I_new;
    real R_new;
    real N_new;

    y = year_id[t];
    N_prev = pop_model[t - 1];

    import_count_week = import_count_year[y] / 52.0;
    import_frac       = import_count_week / N_prev;

    S_prev_count = S_frac[t - 1] * N_prev;
    I_prev_count = I_frac[t - 1] * N_prev;
    R_prev_count = R_frac[t - 1] * N_prev;

    infectious_pressure = I_frac[t - 1] + import_frac;

    lambda_t = beta_t[t] * infectious_pressure;
    lambda_out[t] = lambda_t;

    p_inf = 1.0 - exp(-lambda_t);
    new_inf_count[t] = S_prev_count * p_inf;

    p_rec = 1.0 - exp(-gamma_fixed);
    recov_count = I_prev_count * p_rec;

    // ----- Epidemic transitions + births -----------------------------------
    // Infection and recovery only move people between S, I, R (total conserved).
    // Newborns enter S.
    births_count[t] = births_weekly[t];

    S_post = S_prev_count - new_inf_count[t] + births_count[t];
    I_post = I_prev_count + new_inf_count[t] - recov_count;
    R_post = R_prev_count + recov_count;

    // ----- Natural mortality -----------------------------------------------
    // Remove a constant per-capita fraction p_death from every compartment.
    deaths_S[t] = S_post * p_death;
    deaths_I[t] = I_post * p_death;
    deaths_R[t] = R_post * p_death;
    deaths_count[t] = deaths_S[t] + deaths_I[t] + deaths_R[t];

    S_new = S_post - deaths_S[t];
    I_new = I_post - deaths_I[t];
    R_new = R_post - deaths_R[t];

    N_new = S_new + I_new + R_new;   // = (N_prev + births) * (1 - p_death)
    pop_model[t] = N_new;

    // Fractions of the MODELLED population (sum to 1 by construction).
    S_frac[t] = S_new / N_new;
    I_frac[t] = I_new / N_new;
    R_frac[t] = R_new / N_new;

    new_inf_frac[t] = new_inf_count[t] / N_new;

    // Reported cases track incident infections during the week (pre-death),
    // consistent with v12e.
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
  rho ~ beta(rho_alpha, rho_beta); // reported confirmed cases / true infections

  mu_log_beta ~ normal(0, 0.4);

  log_beta_season_raw ~ normal(0, 1);

  for (w in 2:W) {
    log_beta_season_raw[w] - log_beta_season_raw[w - 1] ~
      normal(0, sigma_season);
  }

  log_beta_season_raw[1] - log_beta_season_raw[W] ~
    normal(0, sigma_season);

  sigma_season ~ normal(0, 0.12);

  delta_year_raw ~ normal(0, 0.8);

  import_count_year ~ lognormal(import_prior_logmedian, import_prior_sd);

  I0 ~ lognormal(log(5), 0.6);
  S0_frac ~ beta(s0_alpha, s0_beta);

  phi_cases ~ lognormal(log(10), 0.5);

  for (t in fit_start:N) {
    cases[t] ~ neg_binomial_2(mu_cases[t], phi_cases);
  }

  // Optional external FOI-derived cumulative attack constraint.
  // Used only when use_cum_attack_prior == 1.
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

  // Demographic flow outputs
  vector[N] births_out;
  vector[N] deaths_out;
  vector[N] deaths_S_out;
  vector[N] deaths_I_out;
  vector[N] deaths_R_out;
  vector[N] pop_model_out;      // modelled population (births - natural deaths)
  vector[N] pop_balance_resid;  // pop_model - observed pop, demographic drift check

  vector[Y] annual_births;
  vector[Y] annual_deaths;

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
    annual_births[y] = 0.0;
    annual_deaths[y] = 0.0;
  }

  for (t in 1:N) {
    R0_t[t] = beta_t[t] / (1.0 - exp(-gamma_fixed));
    Reff_t[t] = beta_t[t] / (1.0 - exp(-gamma_fixed)) * S_frac[t];

    // Counts are expressed on the modelled population.
    S_out[t] = S_frac[t] * pop_model[t];
    I_out[t] = I_frac[t] * pop_model[t];
    R_out[t] = R_frac[t] * pop_model[t];
    new_inf[t] = new_inf_count[t];

    births_out[t] = births_count[t];
    deaths_out[t] = deaths_count[t];
    deaths_S_out[t] = deaths_S[t];
    deaths_I_out[t] = deaths_I[t];
    deaths_R_out[t] = deaths_R[t];

    pop_model_out[t] = pop_model[t];
    pop_balance_resid[t] = pop_model[t] - pop[t];

    annual_births[year_id[t]] += births_count[t];
    annual_deaths[year_id[t]] += deaths_count[t];

    local_infectious_frac_out[t] = I_frac[t];
    import_frac_out[t] = (import_count_year[year_id[t]] / 52.0) / pop_model[t];
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
