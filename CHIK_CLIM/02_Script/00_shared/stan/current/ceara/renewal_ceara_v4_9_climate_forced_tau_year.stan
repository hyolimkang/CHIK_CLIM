// Derived from the accepted renewal_ceara_v4_9_climate_forced_canary.stan.
// ONLY the year-effect variance specification changes -- everything else
// (renewal mechanism, seasonal harmonics, seasonal-amplitude hierarchy,
// climate coefficients/priors, generation interval, demography,
// observation likelihood, serology likelihood, fixed q) is byte-identical.
//
// Diff summary vs. the accepted file:
//   data:    - year_effect_prior_sd (now a free parameter, not fixed data)
//             + tau_prior_scale (the HalfNormal prior's scale -- varies
//               across P1/P2/P3 sensitivity priors without recompiling)
//   parameters: + tau_year (lower=0 -> HalfNormal via the declared bound)
//   transformed parameters: year_effect = tau_year * (z_year - mean(z_year))
//                            [was: year_effect_prior_sd * (z_year - mean(z_year))]
//   model:   + tau_year ~ normal(0, tau_prior_scale)  [HalfNormal via lower=0]
//
// Non-centred parameterisation: z_year ~ std_normal() (unit-scale raw
// parameter, as in the accepted model) multiplied by the free scale
// tau_year -- this is the standard non-centred hierarchical form.
//
// This file is also reused, UNCHANGED, for the training-only holdout
// fits (shorter N, is_seed all 0, no 2022 seed) -- no Stan-level change
// is needed for that; only the R-side data preparation differs.

data {
  int<lower=2> N;
  int<lower=1> G;
  array[N] int<lower=0> C;
  simplex[G] w;
  vector<lower=0>[N] N_start;
  vector<lower=0>[N] N_end;
  vector<lower=0>[N] births;
  vector<lower=0>[N] deaths;
  int<lower=1> Y;
  array[N] int<lower=1, upper=Y> year_id;
  vector[N] seasonal_sin;
  vector[N] seasonal_cos;
  vector[N] seasonal_sin2;
  vector[N] seasonal_cos2;
  real<lower=0> imports_per_week;
  real<lower=0> tau_prior_scale; // CHANGED: was year_effect_prior_sd (fixed data); now the HalfNormal prior scale for the ESTIMATED tau_year
  int<lower=1, upper=N> t_sero;
  int<lower=0> sero_pos;
  int<lower=1> sero_n;
  real<lower=0> kappa_sero;
  real<lower=0, upper=1> q; // FIXED, unchanged from v4.9

  array[N] int<lower=0, upper=1> is_seed;
  vector<lower=0>[N] X_seed;
  int<lower=0, upper=N> N_fit;
  array[N_fit] int<lower=1, upper=N> fit_index;

  vector[N] z_T_anom;
  vector[N] z_P_anom;
}

transformed data {
  vector[N] reconciliation;
  for (t in 1:N) {
    reconciliation[t] = N_end[t] - N_start[t] - births[t] + deaths[t];
  }
}

parameters {
  real alpha_R;
  vector[Y] z_year;
  real<lower=0> tau_year; // NEW: estimated year-effect scale (was fixed data year_effect_prior_sd)
  real beta_sin1;
  real beta_cos1;
  real beta_sin2;
  real beta_cos2;
  real<lower=0> sigma_season_year;
  vector[Y] z_season_year;
  real<lower=0> phi_obs;

  real beta_T;
  real beta_P;
}

transformed parameters {
  vector[Y] year_effect;
  vector[Y] A_year;
  vector[N] log_R0;
  vector[N] R0_t;
  vector[N] R_eff_t;
  vector[N] force_of_infection;
  vector[N] X;
  vector[N] S;
  vector[N] U;
  vector[N] S_prop;
  vector[N] immune_prop;
  vector[N] expected_reported_cases;
  real p_state_sero_at_anchor;
  real p_sero_safe;
  real alpha_sero;
  real beta_sero;

  year_effect = tau_year * (z_year - mean(z_year)); // CHANGED: tau_year (estimated) replaces year_effect_prior_sd (fixed)
  A_year = exp(sigma_season_year * z_season_year);

  S[1] = N_start[1];
  U[1] = 0;

  for (t in 1:N) {
    real total = S[t] + U[t];
    real infectiousness = 0;

    if (S[t] < 0 || U[t] < 0 || fabs(total - N_start[t]) > 1e-6 * fmax(1, N_start[t])) {
      reject("Invalid demographic state at week ", t);
    }
    for (g in 1:G) {
      if (t > g) infectiousness += w[g] * X[t - g];
    }

    log_R0[t] = alpha_R + year_effect[year_id[t]] +
      A_year[year_id[t]] * (beta_sin1 * seasonal_sin[t] + beta_cos1 * seasonal_cos[t]) +
      beta_sin2 * seasonal_sin2[t] + beta_cos2 * seasonal_cos2[t] +
      beta_T * z_T_anom[t] + beta_P * z_P_anom[t];
    R0_t[t] = exp(log_R0[t]);
    force_of_infection[t] = R0_t[t] * infectiousness / total + imports_per_week / total;

    if (is_seed[t] == 1) {
      if (X_seed[t] > S[t]) {
        reject("2022 seed infections exceed available susceptibles at week ", t);
      }
      X[t] = X_seed[t];
    } else {
      X[t] = S[t] * (-expm1(-force_of_infection[t]));
    }

    R_eff_t[t] = R0_t[t] * S[t] / total;
    S_prop[t] = S[t] / total;
    immune_prop[t] = U[t] / total;
    expected_reported_cases[t] = q * X[t] + 1e-9;

    if (t < N) {
      real susceptible_after_infection = S[t] - X[t];
      real immune_after_infection = U[t] + X[t];
      real susceptible_fraction = susceptible_after_infection /
        (susceptible_after_infection + immune_after_infection);
      S[t + 1] = susceptible_after_infection + births[t] -
        deaths[t] * susceptible_fraction + reconciliation[t] * susceptible_fraction;
      U[t + 1] = immune_after_infection - deaths[t] * (1 - susceptible_fraction) +
        reconciliation[t] * (1 - susceptible_fraction);
    }
  }

  p_state_sero_at_anchor = immune_prop[t_sero];
  p_sero_safe = fmin(1 - 1e-9, fmax(1e-9, p_state_sero_at_anchor));
  alpha_sero = p_sero_safe * kappa_sero;
  beta_sero = (1 - p_sero_safe) * kappa_sero;
}

model {
  alpha_R ~ normal(log(1.2), 0.5);
  z_year ~ std_normal();
  tau_year ~ normal(0, tau_prior_scale); // NEW: HalfNormal(0, tau_prior_scale) via the lower=0 declaration
  beta_sin1 ~ normal(0, 0.25);
  beta_cos1 ~ normal(0, 0.25);
  beta_sin2 ~ normal(0, 0.10);
  beta_cos2 ~ normal(0, 0.10);
  sigma_season_year ~ normal(0, 0.20);
  z_season_year ~ std_normal();
  phi_obs ~ gamma(2, 0.1);

  beta_T ~ normal(0, 0.15);
  beta_P ~ normal(0, 0.15);

  C[fit_index] ~ neg_binomial_2(expected_reported_cases[fit_index], phi_obs);
  sero_pos ~ beta_binomial(sero_n, alpha_sero, beta_sero);
}

generated quantities {
  array[N] int C_pred;
  vector[N] log_lik_cases;
  vector[N] climate_multiplier;
  real log_lik_serology;
  int sero_pred;

  for (t in 1:N) {
    climate_multiplier[t] =
      exp(beta_T * z_T_anom[t] +
          beta_P * z_P_anom[t]);

    C_pred[t] = neg_binomial_2_rng(expected_reported_cases[t], phi_obs);
    log_lik_cases[t] = neg_binomial_2_lpmf(C[t] | expected_reported_cases[t], phi_obs);
  }
  log_lik_serology = beta_binomial_lpmf(sero_pos | sero_n, alpha_sero, beta_sero);
  sero_pred = beta_binomial_rng(sero_n, alpha_sero, beta_sero);
}
