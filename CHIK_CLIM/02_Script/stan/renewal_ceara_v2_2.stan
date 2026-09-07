// renewal_ceara_v2_2.stan
//
// Explicit demographic accounting for the Ceara renewal model, per
// docs/model_development_roadmap.md sections 3-4 ("Run B"). This changes
// ONLY the demographic mechanism versus v2.1: transmission, reporting and
// serology-linkage structure, priors and likelihood are unchanged, so a
// v2.1-vs-v2.2 comparison isolates the effect of that one change.
//
// v2.1 updated S[t] from the net population change alone (all of a
// positive delta entering S, all of a negative delta leaving S
// proportionally) -- an approximation the roadmap shows misses real
// demographic turnover even when total population is unchanged (its
// section 3.3 worked example). v2.2 instead tracks S[t] (susceptible) and
// U[t] (no longer susceptible under the model's infection-derived
// protection assumption) as two conserved compartments:
//   S[t+1] = S[t] - X[t] + births[t] - deaths_S[t] + migration_S[t]
//   U[t+1] = U[t] + X[t]              - deaths_U[t] + migration_U[t]
// Births enter S in full (an aggregate approximation -- no neonatal
// protection is modelled). All-cause deaths and the net population
// reconciliation residual (net_migration_reconciliation[t] = N_end[t] -
// N_start[t] - births[t] + deaths[t], an exact identity given externally
// fixed N, not an estimate of true migration) are both split across S/U in
// proportion to the *post-infection* S/U composition -- a documented
// approximation (roadmap section 4.4), since the age/immune composition of
// deaths and unmeasured net migration is not separately observed here.
//
// N[t] = S[t] + U[t] is therefore conserved by construction and is not
// supplied as a separate external series (unlike v2.1's N_pop_t) --
// N_start/N_end only seed S[1] and define the reconciliation residual.

data {
  int<lower=2> N;
  int<lower=1> G;
  int<lower=1> seed_weeks;
  array[N] int<lower=0> C;
  simplex[G] w;

  vector<lower=0>[N] N_start;
  vector<lower=0>[N] N_end;
  vector<lower=0>[N] births;
  vector<lower=0>[N] deaths;
  vector[N] time_scaled;

  vector[seed_weeks] log_seed_prior_mean;
  real<lower=0> seed_prior_sd;
  real<lower=0> reporting_region_sd;

  int<lower=1> J;
  array[J] int<lower=1, upper=N> sero_window_start;
  array[J] int<lower=1, upper=N> sero_window_end;
  array[J] int<lower=0> sero_positive;
  array[J] int<lower=1> sero_n;
}

transformed data {
  // Exact residual, not an estimated migration flow: whatever the fixed
  // N series does that births and deaths alone do not explain.
  vector[N] net_migration_reconciliation;
  for (t in 1:N) {
    net_migration_reconciliation[t] = N_end[t] - N_start[t] - births[t] + deaths[t];
  }
}

parameters {
  real<lower=0> sigma_R;
  real log_infection_scale;
  vector[N - 1] log_hazard_relative;

  real<lower=0, upper=1> p_symp;
  real<lower=0, upper=1> rho_sym_brazil;
  real logit_rho_sym_ceara_mid;
  real<lower=0> reporting_trend;

  real<lower=0> sigma_geo;
  vector[J] z_site;
  real<lower=0> phi_obs;
}

transformed parameters {
  vector[N] log_infection_hazard;
  vector[N] log_R0;
  vector[N] R0_t;
  vector[N] R_eff_t;
  vector[N] X;
  vector[N] S;
  vector[N] U;
  vector[N] immune_prop;
  vector[N] rho_sym_t;
  vector[N] rho_total_t;
  vector[N] expected_reported_cases;
  real<lower=0, upper=1> rho_sym_ceara_mid;
  real reporting_ceara_offset;
  vector[J] state_attack_window;
  vector[J] sero_geographic_offset;
  vector[J] p_site;

  {
    S[1] = N_start[1];
    U[1] = 0;

    log_infection_hazard[1] = log_infection_scale;
    for (t in 2:N) {
      log_infection_hazard[t] = log_infection_scale
                                + log_hazard_relative[t - 1];
    }

    rho_sym_ceara_mid = inv_logit(logit_rho_sym_ceara_mid);
    reporting_ceara_offset = logit_rho_sym_ceara_mid
                             - logit(rho_sym_brazil);

    for (t in 1:N) {
      if (S[t] < 0 || U[t] < 0) {
        reject("Negative demographic state at week ", t,
               ": S=", S[t], " U=", U[t]);
      }

      rho_sym_t[t] = inv_logit(
        logit_rho_sym_ceara_mid + reporting_trend * time_scaled[t]
      );
      rho_total_t[t] = p_symp * rho_sym_t[t];

      X[t] = S[t] * (-expm1(-exp(log_infection_hazard[t])));

      if (t > seed_weeks) {
        real infectiousness = 0;
        for (g in 1:G) {
          infectiousness += w[g] * X[t - g];
        }
        log_R0[t] = log_infection_hazard[t] + log(S[t] + U[t])
                    - log(infectiousness);
        R0_t[t] = exp(log_R0[t]);
      } else {
        // R0 is not identified before a complete generation-interval history
        // exists. Fill these output positions after the recursion below.
        log_R0[t] = 0;
        R0_t[t] = 1;
      }

      R_eff_t[t] = R0_t[t] * S[t] / (S[t] + U[t]);
      immune_prop[t] = U[t] / (S[t] + U[t]);
      expected_reported_cases[t] = rho_total_t[t] * X[t] + 1e-9;

      if (t < N) {
        real s_after = S[t] - X[t];
        real u_after = U[t] + X[t];
        // Post-infection composition (pre-vital-event population = N_start[t]
        // exactly, since infection only moves people between S and U).
        real frac_s = s_after / (s_after + u_after);

        real deaths_s = deaths[t] * frac_s;
        real deaths_u = deaths[t] - deaths_s;
        real migration_s = net_migration_reconciliation[t] * frac_s;
        real migration_u = net_migration_reconciliation[t] - migration_s;

        S[t + 1] = s_after + births[t] - deaths_s + migration_s;
        U[t + 1] = u_after - deaths_u + migration_u;

        if (S[t + 1] < 0 || U[t + 1] < 0) {
          reject("Demographic accounting produced a negative state at week ",
                 t + 1, ": S=", S[t + 1], " U=", U[t + 1],
                 " (check net_migration_reconciliation at week ", t, ")");
        }
      }
    }

    // The seed-period values are display placeholders, explicitly copied from
    // the first estimable R0. Downstream plots omit these first seed weeks.
    for (t in 1:seed_weeks) {
      log_R0[t] = log_R0[seed_weeks + 1];
      R0_t[t] = R0_t[seed_weeks + 1];
      R_eff_t[t] = R0_t[t];
    }
  }

  for (j in 1:J) {
    state_attack_window[j] = mean(segment(
      immune_prop,
      sero_window_start[j],
      sero_window_end[j] - sero_window_start[j] + 1
    ));
    sero_geographic_offset[j] = sigma_geo * z_site[j];
    p_site[j] = inv_logit(
      logit(state_attack_window[j]) + sero_geographic_offset[j]
    );
  }
}

model {
  // Smooth time-varying basic reproduction number.
  sigma_R ~ normal(0, 0.10);
  log_R0[seed_weeks + 1] ~ normal(log(1.0), 0.5);
  for (t in (seed_weeks + 2):N) {
    log_R0[t] ~ normal(log_R0[t - 1], sigma_R);
  }

  // Weak seed information is expressed on the corresponding low-incidence
  // count scale: log(hazard) is approximately log(X) - log(N_start).
  segment(log_infection_hazard, 1, seed_weeks) ~ normal(
    log_seed_prior_mean - log(segment(N_start, 1, seed_weeks)),
    seed_prior_sd
  );

  // Brazil-wide long-run symptomatic reporting, with a Ceara offset and a
  // non-negative linear change on the logit scale during the fitting window.
  p_symp ~ beta(30, 28);
  rho_sym_brazil ~ beta(20, 60);
  logit_rho_sym_ceara_mid ~ normal(
    logit(rho_sym_brazil), reporting_region_sd
  );
  reporting_trend ~ normal(0, 0.75);

  // Non-centered geographic heterogeneity around each window-mean state
  // attack rate.
  sigma_geo ~ normal(0, 1);
  z_site ~ normal(0, 1);

  phi_obs ~ gamma(2, 0.1);

  C ~ neg_binomial_2(expected_reported_cases, phi_obs);
  for (j in 1:J) {
    sero_positive[j] ~ binomial(sero_n[j], p_site[j]);
  }
}

generated quantities {
  array[N] int C_pred;
  vector[N] log_lik_cases;
  vector[N] overall_detection;
  array[J] int sero_positive_pred;
  vector[J] log_lik_serology;

  for (t in 1:N) {
    C_pred[t] = neg_binomial_2_rng(expected_reported_cases[t], phi_obs);
    log_lik_cases[t] = neg_binomial_2_lpmf(
      C[t] | expected_reported_cases[t], phi_obs
    );
    overall_detection[t] = p_symp * rho_sym_t[t];
  }
  for (j in 1:J) {
    sero_positive_pred[j] = binomial_rng(sero_n[j], p_site[j]);
    log_lik_serology[j] = binomial_lpmf(
      sero_positive[j] | sero_n[j], p_site[j]
    );
  }
}
