// renewal_ceara_v2.stan
//
// Dynamic susceptible-depletion renewal model for weekly chikungunya cases
// in Ceara. X[t] is true incident infection count, not reported incidence.
// A smooth weekly R0[t] random walk provides process variation. For efficient
// HMC, the model is parameterized by weekly log infection hazards and R0[t] is
// recovered from the renewal equation; this is the same deterministic renewal
// process expressed in coordinates directly informed by observed incidence.
// Reporting is separated into a fixed externally supplied symptomatic
// probability and a monotonically increasing symptomatic-case reporting
// probability.
//
// The population is closed over this first implementation. Immunity is
// lifelong over the fitting period and S_prop[1] starts at one immediately
// before the first seeded infections. A 2018 serosurvey is linked to state
// attack rate through a Juazeiro-versus-Ceara log-odds offset.

data {
  int<lower=2> N;
  int<lower=1> G;
  int<lower=1> seed_weeks;
  array[N] int<lower=0> C;
  simplex[G] w;

  real<lower=1> N_pop;
  vector[N] time_scaled;
  real<lower=0, upper=1> p_symp;

  vector[seed_weeks] log_seed_prior_mean;
  real<lower=0> seed_prior_sd;
  real<lower=0> reporting_region_sd;

  int<lower=1, upper=N> sero_t;
  int<lower=0> sero_positive;
  int<lower=1> sero_tested;
  real<lower=0, upper=1> sero_sensitivity;
  real<lower=0, upper=1> sero_specificity;
  real<lower=0> sero_geographic_sd;
}

parameters {
  real<lower=0> sigma_R;
  real log_infection_scale;
  vector[N - 1] log_hazard_relative;

  real<lower=0, upper=1> rho_sym_brazil;
  real logit_rho_sym_ceara_mid;
  real<lower=0> reporting_trend;

  real logit_site_attack_at_serosurvey;
  real<lower=0> phi_obs;
}

transformed parameters {
  vector[N] log_infection_hazard;
  vector[N] log_R0;
  vector[N] R0_t;
  vector[N] R_eff_t;
  vector[N] X;
  vector[N] S_prop;
  vector[N] immune_prop;
  vector[N] rho_sym_t;
  vector[N] rho_total_t;
  vector[N] expected_reported_cases;
  real<lower=0, upper=1> rho_sym_ceara_mid;
  real reporting_ceara_offset;
  real state_attack_at_serosurvey;
  real site_attack_at_serosurvey;
  real sero_geographic_offset;
  real serosurvey_positive_probability;

  {
    real susceptible_count = N_pop;

    log_infection_hazard[1] = log_infection_scale;
    for (t in 2:N) {
      log_infection_hazard[t] = log_infection_scale
                                + log_hazard_relative[t - 1];
    }

    rho_sym_ceara_mid = inv_logit(logit_rho_sym_ceara_mid);
    reporting_ceara_offset = logit_rho_sym_ceara_mid
                             - logit(rho_sym_brazil);

    for (t in 1:N) {
      rho_sym_t[t] = inv_logit(
        logit_rho_sym_ceara_mid + reporting_trend * time_scaled[t]
      );
      rho_total_t[t] = p_symp * rho_sym_t[t];

      X[t] = susceptible_count * (-expm1(-exp(log_infection_hazard[t])));

      if (t > seed_weeks) {
        real infectiousness = 0;
        for (g in 1:G) {
          infectiousness += w[g] * X[t - g];
        }
        log_R0[t] = log_infection_hazard[t] + log(N_pop)
                    - log(infectiousness);
        R0_t[t] = exp(log_R0[t]);
      } else {
        // R0 is not identified before a complete generation-interval history
        // exists. Fill these output positions after the recursion below.
        log_R0[t] = 0;
        R0_t[t] = 1;
      }

      R_eff_t[t] = R0_t[t] * susceptible_count / N_pop;
      susceptible_count -= X[t];
      S_prop[t] = susceptible_count / N_pop;
      immune_prop[t] = 1 - S_prop[t];
      expected_reported_cases[t] = rho_total_t[t] * X[t] + 1e-9;
    }

    // The seed-period values are display placeholders, explicitly copied from
    // the first estimable R0. Downstream plots omit these first seed weeks.
    for (t in 1:seed_weeks) {
      log_R0[t] = log_R0[seed_weeks + 1];
      R0_t[t] = R0_t[seed_weeks + 1];
      R_eff_t[t] = R0_t[t];
    }
  }

  state_attack_at_serosurvey = immune_prop[sero_t];
  site_attack_at_serosurvey = inv_logit(logit_site_attack_at_serosurvey);
  sero_geographic_offset = logit_site_attack_at_serosurvey
                           - logit(state_attack_at_serosurvey);
  serosurvey_positive_probability =
    sero_sensitivity * site_attack_at_serosurvey
    + (1 - sero_specificity) * (1 - site_attack_at_serosurvey);
}

model {
  // Smooth time-varying basic reproduction number.
  sigma_R ~ normal(0, 0.10);
  log_R0[seed_weeks + 1] ~ normal(log(1.0), 0.5);
  for (t in (seed_weeks + 2):N) {
    log_R0[t] ~ normal(log_R0[t - 1], sigma_R);
  }

  // Weak seed information is expressed on the corresponding low-incidence
  // count scale: log(hazard) is approximately log(X) - log(N_pop).
  segment(log_infection_hazard, 1, seed_weeks) ~ normal(
    log_seed_prior_mean - rep_vector(log(N_pop), seed_weeks),
    seed_prior_sd
  );

  // Brazil-wide long-run symptomatic reporting, with a Ceara offset and a
  // non-negative linear change on the logit scale during 2015-2018.
  rho_sym_brazil ~ beta(20, 60);
  logit_rho_sym_ceara_mid ~ normal(
    logit(rho_sym_brazil), reporting_region_sd
  );
  reporting_trend ~ normal(0, 0.75);

  // Geographic transport from state attack rate to the Juazeiro survey.
  logit_site_attack_at_serosurvey ~ normal(
    logit(state_attack_at_serosurvey), sero_geographic_sd
  );

  phi_obs ~ gamma(2, 0.1);

  C ~ neg_binomial_2(expected_reported_cases, phi_obs);
  sero_positive ~ binomial(sero_tested, serosurvey_positive_probability);
}

generated quantities {
  array[N] int C_pred;
  vector[N] log_lik_cases;
  int sero_positive_pred;
  real log_lik_serology;

  for (t in 1:N) {
    C_pred[t] = neg_binomial_2_rng(expected_reported_cases[t], phi_obs);
    log_lik_cases[t] = neg_binomial_2_lpmf(
      C[t] | expected_reported_cases[t], phi_obs
    );
  }
  sero_positive_pred = binomial_rng(
    sero_tested, serosurvey_positive_probability
  );
  log_lik_serology = binomial_lpmf(
    sero_positive | sero_tested, serosurvey_positive_probability
  );
}
