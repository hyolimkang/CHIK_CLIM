// Climate-forced v4.9 CANARY (Phase 2 of the climate-forcing extension,
// see 02_Script/40_renewal_model/36_climate_forced_v4_9/PHASE0_AUDIT.md).
//
// This file is BYTE-IDENTICAL to the frozen
// renewal_ceara_v4_9_hierarchical_seasonality.stan EXCEPT for three
// additions, each marked "// CLIMATE:" below:
//   1. two new data vectors: z_T_anom, z_P_anom (standardised climate
//      anomalies, primary lag spec, from 01_build_climate_anomaly_covariates.R)
//   2. two new parameters: beta_T, beta_P, each with a regularising
//      normal(0, 0.15) prior (prior-predictive-checked in
//      02_prior_predictive_climate_multiplier.R: central 50% of
//      climate_multiplier = exp(beta_T*zT + beta_P*zP) falls in
//      [0.93, 1.08], 90% in [0.77, 1.30], extreme multipliers rare)
//   3. one existing line extended: log_R0[t] gains the additive term
//      "+ beta_T * z_T_anom[t] + beta_P * z_P_anom[t]"
//
// NOTHING else changes: the S/U recursion, generation-interval convolution,
// demographic reconciliation, NB2 observation likelihood, Juazeiro serology
// treatment, fixed q, 2022 seed mechanism, and all other priors are
// untouched. climate_multiplier[t] itself is NOT computed inside Stan (to
// keep this diff minimal); it is reconstructed post-hoc in R from the
// posterior draws of beta_T/beta_P and the (fixed, known) z_T_anom/z_P_anom
// data vectors -- exp(a+b) = exp(a)*exp(b), so R0_t[t] = exp(log_R0[t])
// already equals baseline_R0[t] * climate_multiplier[t] with no other
// change required to the R0_t/R_eff_t computation below.
//
// ---- original v4.9 header, unchanged, for provenance ----
// v4.9 = v4.8 + a parsimonious hierarchical seasonal-SHAPE model, replacing
// ONLY the temporal R0(t) parameterisation. v4.8 showed that conditioning
// on an externally seeded 2022 reintroduction leaves the post-seed 2022
// burden essentially unchanged from v4.7 (still ~half the observed
// total), AND that even the already-fitted 2017 epidemic's WEEKLY peak is
// underpredicted by the posterior median (~half the observed ~10,000/week
// peak) despite an adequately calibrated ANNUAL burden -- i.e. the broad
// NB2 predictive interval covers the peak, but the mean trajectory is too
// flat. v4.8 (and v4.0-v4.8 more generally) shares one first-harmonic
// seasonal curve across every year, with only an additive annual
// intercept (year_effect). v4.9 tests whether allowing the SEASONAL
// AMPLITUDE itself to vary (weakly, hierarchically) between years can
// sharpen the mean epidemic shape, without adding a high-dimensional
// weekly R process, region effects, or outbreak-specific pulses.
//
// Everything except the log_R0(t) formula is BYTE-IDENTICAL to v4.8:
// fixed q (data), demographic S/U recursion, lifelong immunity, weekly NB
// observation likelihood (seed weeks excluded via fit_index), phi_obs
// prior, Juazeiro serology treatment, the 2022 conditional seed mechanism
// (is_seed/X_seed), imports_per_week, generation-interval kernel.

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
  vector[N] seasonal_sin;   // first annual harmonic (unchanged from v4.0-v4.8)
  vector[N] seasonal_cos;
  vector[N] seasonal_sin2;  // second annual harmonic (period ~26.1 weeks)
  vector[N] seasonal_cos2;
  real<lower=0> imports_per_week;
  real<lower=0> year_effect_prior_sd;
  int<lower=1, upper=N> t_sero;
  int<lower=0> sero_pos;
  int<lower=1> sero_n;
  real<lower=0> kappa_sero;
  real<lower=0, upper=1> q; // FIXED, same grid as v4.7/v4.8

  // 2022 seed window -- identical mechanism to v4.8, unchanged.
  array[N] int<lower=0, upper=1> is_seed;
  vector<lower=0>[N] X_seed;
  int<lower=0, upper=N> N_fit;
  array[N_fit] int<lower=1, upper=N> fit_index;

  // CLIMATE: standardised anomalies, primary lag spec (Phase 1). Aligned
  // 1:1 to this state's own N weekly rows (same week_start ordering as
  // every other per-week vector in this data block).
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
  real beta_sin1;             // first-harmonic coefficients (year-modulated by A_year)
  real beta_cos1;
  real beta_sin2;              // second-harmonic coefficients (shared across years, NOT year-modulated)
  real beta_cos2;
  real<lower=0> sigma_season_year; // hierarchical SD of log-amplitude deviations across years
  vector[Y] z_season_year;         // non-centred raw year-specific log-amplitude deviations
  real<lower=0> phi_obs;
  real beta_T; // CLIMATE: temperature-anomaly coefficient on log_R0
  real beta_P; // CLIMATE: precipitation-anomaly coefficient on log_R0
}

transformed parameters {
  vector[Y] year_effect;
  vector[Y] A_year; // year-specific first-harmonic amplitude multiplier, centred at 1
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

  year_effect = year_effect_prior_sd * (z_year - mean(z_year));
  A_year = exp(sigma_season_year * z_season_year); // non-centred: A_year = 1 when z_season_year = 0

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
      beta_T * z_T_anom[t] + beta_P * z_P_anom[t]; // CLIMATE: only change to this line
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
  beta_sin1 ~ normal(0, 0.25); // identical scale to v4.0-v4.8's shared beta_sin/beta_cos prior
  beta_cos1 ~ normal(0, 0.25);
  beta_sin2 ~ normal(0, 0.10); // conservative shrinkage: second harmonic is new, higher-order flexibility
  beta_cos2 ~ normal(0, 0.10);
  sigma_season_year ~ normal(0, 0.20); // half-normal (sigma>0): typical year stays close to the shared curve
  z_season_year ~ std_normal();        // non-centred
  phi_obs ~ gamma(2, 0.1); // UNCHANGED from v4.8 (diagnostic-only audit target, not modified)
  beta_T ~ normal(0, 0.15); // CLIMATE: regularising, prior-predictive-checked
  beta_P ~ normal(0, 0.15); // CLIMATE: regularising, prior-predictive-checked

  C[fit_index] ~ neg_binomial_2(expected_reported_cases[fit_index], phi_obs);
  sero_pos ~ beta_binomial(sero_n, alpha_sero, beta_sero);
}

generated quantities {
  array[N] int C_pred;
  vector[N] log_lik_cases;
  real log_lik_serology;
  int sero_pred;
  for (t in 1:N) {
    C_pred[t] = neg_binomial_2_rng(expected_reported_cases[t], phi_obs);
    log_lik_cases[t] = neg_binomial_2_lpmf(C[t] | expected_reported_cases[t], phi_obs);
  }
  log_lik_serology = beta_binomial_lpmf(sero_pos | sero_n, alpha_sero, beta_sero);
  sero_pred = beta_binomial_rng(sero_n, alpha_sero, beta_sero);
}
