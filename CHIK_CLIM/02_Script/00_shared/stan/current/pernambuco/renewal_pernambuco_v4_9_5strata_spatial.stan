// Pernambuco Spatial Model 1: minimal 5-stratum extension of the frozen
// v4.9 global-q architecture (renewal_bahia_v4_9_global_q_multisite_serology.stan).
//
// Purpose (per the pre-registered spatial-turnover-diagnostic follow-up
// spec): test whether allowing geographically distinct susceptible pools
// removes the pathological low-q / high-depletion / high-R0 geometry found
// in the homogeneous-state PE model, WITHOUT optimising fit flexibility.
//
// What is NEW relative to the frozen v4.9 architecture:
//   - R (=5) independent renewal/S-U recursions, one per stratum, using the
//     SAME demographic bookkeeping and lifelong-immunity assumptions.
//   - log R0[r,t] = alpha_global + delta_region[r] + year_effect[year] +
//     A_year[year]*(seasonal harmonics) + (fixed second harmonic), with
//     delta_region a non-centred, sum-to-zero-constrained regional
//     intercept (region_contrast_basis, built in R exactly like the
//     existing Y-1 orthonormal year-effect basis). Harmonics, A_year, and
//     year effects remain COMMON across strata (not region-specific).
//   - ONE common q_PE (logit_q) across all 5 strata -- same prior as the
//     homogeneous PE global-q model.
//   - ONE shared phi_obs across all 5 regional case series.
//   - Statewide importation is allocated across strata proportional to
//     population (imports_r, precomputed in R) -- no new import
//     parameters, no inter-region transmission coupling, no mobility
//     matrix.
//   - Recife U14 serology (Aug2018-Feb2019, n=2070/770) links DIRECTLY to
//     the Recife stratum's own immune_prop window average -- NO
//     state-to-Recife geographic-offset (eta_geo) transport step, because
//     Recife is now its own stratum (the whole point of this model).
//
// What is UNCHANGED: renewal equation, generation interval, NB2 case
// likelihood, phi_obs prior, beta-binomial serology convention and
// kappa_sero, lifelong-immunity S/U bookkeeping, logit_q prior. No
// region-specific q, harmonics, seasonal phase, A_year, or year effects.
// No mobility matrix, no hand-selected regional seeds (is_seed/X_seed
// retained but zero everywhere, exactly as in the validated homogeneous PE
// model, which found no qualifying seed gap).

data {
  int<lower=2> N; // weeks
  int<lower=2> R; // strata (5)
  int<lower=1> G; // generation-interval length

  array[N, R] int<lower=0> C;
  simplex[G] w;

  matrix<lower=0>[N, R] N_start;
  matrix<lower=0>[N, R] N_end;
  matrix<lower=0>[N, R] births;
  matrix<lower=0>[N, R] deaths;
  matrix<lower=0>[N, R] imports_r; // population-scaled statewide imports

  int<lower=2> Y;
  array[N] int<lower=1, upper=Y> year_id;
  vector[N] seasonal_sin;
  vector[N] seasonal_cos;
  vector[N] seasonal_sin2;
  vector[N] seasonal_cos2;
  real<lower=0> year_effect_prior_sd;
  matrix[Y, Y - 1] year_contrast_basis;
  matrix[R, R - 1] region_contrast_basis;

  array[N, R] int<lower=0, upper=1> is_seed;
  matrix<lower=0>[N, R] X_seed;

  int<lower=1, upper=R> RECIFE_INDEX;
  int<lower=0> J_sero;
  array[J_sero] int<lower=0> sero_n_positive;
  array[J_sero] int<lower=1> sero_n_tested;
  array[J_sero] int<lower=1, upper=N> sero_window_start_idx;
  array[J_sero] int<lower=1> sero_window_n_weeks;
  real<lower=0> kappa_sero;

  real logit_q_prior_mean;
  real<lower=0> logit_q_prior_sd;
}

transformed data {
  matrix[N, R] reconciliation;
  for (t in 1:N) {
    for (r in 1:R) {
      reconciliation[t, r] = N_end[t, r] - N_start[t, r] - births[t, r] + deaths[t, r];
    }
  }
}

parameters {
  real alpha_global;
  vector[R - 1] z_region_contrast;
  real<lower=0> sigma_region;

  vector[Y - 1] z_year_contrast;
  real beta_sin1;
  real beta_cos1;
  real beta_sin2;
  real beta_cos2;
  real<lower=0> sigma_season_year;
  vector[Y] z_season_year;
  real<lower=0> phi_obs;

  real logit_q;
}

transformed parameters {
  vector[R] delta_region;
  vector[Y] year_effect;
  vector[Y] A_year;
  matrix[N, R] log_R0;
  matrix[N, R] R0_t;
  matrix[N, R] R_eff_t;
  matrix[N, R] force_of_infection;
  matrix[N, R] X;
  matrix[N, R] S;
  matrix[N, R] U;
  matrix[N, R] S_prop;
  matrix[N, R] immune_prop;
  matrix[N, R] expected_reported_cases;
  real<lower=0, upper=1> q_PE;

  vector[J_sero] p_site_window;
  vector[J_sero] alpha_sero_site;
  vector[J_sero] beta_sero_site;

  q_PE = inv_logit(logit_q);
  delta_region = sigma_region * (region_contrast_basis * z_region_contrast);
  year_effect = year_effect_prior_sd * (year_contrast_basis * z_year_contrast);
  A_year = exp(sigma_season_year * z_season_year);

  for (r in 1:R) {
    S[1, r] = N_start[1, r];
    U[1, r] = 0;
  }

  for (t in 1:N) {
    for (r in 1:R) {
      real total = S[t, r] + U[t, r];
      real infectiousness = 0;

      if (S[t, r] < 0 || U[t, r] < 0 || fabs(total - N_start[t, r]) > 1e-6 * fmax(1, N_start[t, r])) {
        reject("Invalid demographic state at week ", t, " stratum ", r);
      }
      for (g in 1:G) {
        if (t > g) infectiousness += w[g] * X[t - g, r];
      }

      log_R0[t, r] = alpha_global + delta_region[r] + year_effect[year_id[t]] +
        A_year[year_id[t]] * (beta_sin1 * seasonal_sin[t] + beta_cos1 * seasonal_cos[t]) +
        beta_sin2 * seasonal_sin2[t] + beta_cos2 * seasonal_cos2[t];
      R0_t[t, r] = exp(log_R0[t, r]);
      force_of_infection[t, r] = R0_t[t, r] * infectiousness / total + imports_r[t, r] / total;

      if (is_seed[t, r] == 1) {
        if (X_seed[t, r] > S[t, r]) {
          reject("Seed infections exceed available susceptibles at week ", t, " stratum ", r);
        }
        X[t, r] = X_seed[t, r];
      } else {
        X[t, r] = S[t, r] * (-expm1(-force_of_infection[t, r]));
      }

      R_eff_t[t, r] = R0_t[t, r] * S[t, r] / total;
      S_prop[t, r] = S[t, r] / total;
      immune_prop[t, r] = U[t, r] / total;
      expected_reported_cases[t, r] = q_PE * X[t, r] + 1e-9;

      if (t < N) {
        real susceptible_after_infection = S[t, r] - X[t, r];
        real immune_after_infection = U[t, r] + X[t, r];
        real susceptible_fraction = susceptible_after_infection /
          (susceptible_after_infection + immune_after_infection);
        S[t + 1, r] = susceptible_after_infection + births[t, r] -
          deaths[t, r] * susceptible_fraction + reconciliation[t, r] * susceptible_fraction;
        U[t + 1, r] = immune_after_infection - deaths[t, r] * (1 - susceptible_fraction) +
          reconciliation[t, r] * (1 - susceptible_fraction);
      }
    }
  }

  // Recife U14 serology: DIRECT link to the Recife stratum's own
  // immune_prop window average. No geographic-offset transport (Section 9).
  for (j in 1:J_sero) {
    int start = sero_window_start_idx[j];
    int nwk = sero_window_n_weeks[j];
    real acc = 0;
    for (k in 1:nwk) {
      int t = start + k - 1;
      acc += fmin(1 - 1e-9, fmax(1e-9, immune_prop[t, RECIFE_INDEX]));
    }
    p_site_window[j] = acc / nwk;
    alpha_sero_site[j] = p_site_window[j] * kappa_sero;
    beta_sero_site[j] = (1 - p_site_window[j]) * kappa_sero;
  }
}

model {
  alpha_global ~ normal(log(1.2), 0.5);
  z_region_contrast ~ std_normal();
  sigma_region ~ normal(0, 0.3); // conservative weakly informative prior (half-normal via <lower=0>)

  z_year_contrast ~ std_normal();
  beta_sin1 ~ normal(0, 0.25);
  beta_cos1 ~ normal(0, 0.25);
  beta_sin2 ~ normal(0, 0.10);
  beta_cos2 ~ normal(0, 0.10);
  sigma_season_year ~ normal(0, 0.20);
  z_season_year ~ std_normal();
  phi_obs ~ gamma(2, 0.1);

  logit_q ~ normal(logit_q_prior_mean, logit_q_prior_sd);

  for (r in 1:R) {
    C[, r] ~ neg_binomial_2(expected_reported_cases[, r], phi_obs);
  }
  if (J_sero > 0) {
    for (j in 1:J_sero) {
      sero_n_positive[j] ~ beta_binomial(sero_n_tested[j], alpha_sero_site[j], beta_sero_site[j]);
    }
  }
}

generated quantities {
  array[N, R] int C_pred;
  matrix[N, R] log_lik_cases;
  vector[J_sero] log_lik_serology;
  array[J_sero] int sero_pred;
  vector[N] S_PE_prop;
  vector[N] U_PE_prop;
  array[N] int C_pred_state_sum;
  real logit_q_prior_draw = normal_rng(logit_q_prior_mean, logit_q_prior_sd);
  real q_prior_draw = inv_logit(logit_q_prior_draw);

  for (t in 1:N) {
    real S_state = 0;
    real U_state = 0;
    real total_state = 0;
    int csum = 0;
    for (r in 1:R) {
      C_pred[t, r] = neg_binomial_2_rng(expected_reported_cases[t, r], phi_obs);
      log_lik_cases[t, r] = neg_binomial_2_lpmf(C[t, r] | expected_reported_cases[t, r], phi_obs);
      csum += C_pred[t, r];
      S_state += S[t, r];
      U_state += U[t, r];
      total_state += S[t, r] + U[t, r];
    }
    C_pred_state_sum[t] = csum;
    S_PE_prop[t] = S_state / total_state;
    U_PE_prop[t] = U_state / total_state;
  }
  for (j in 1:J_sero) {
    log_lik_serology[j] = beta_binomial_lpmf(sero_n_positive[j] | sero_n_tested[j], alpha_sero_site[j], beta_sero_site[j]);
    sero_pred[j] = beta_binomial_rng(sero_n_tested[j], alpha_sero_site[j], beta_sero_site[j]);
  }
}
