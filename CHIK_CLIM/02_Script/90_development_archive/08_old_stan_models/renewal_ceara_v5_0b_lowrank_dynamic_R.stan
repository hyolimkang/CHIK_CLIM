// v5.0b = v4.9 + a LOW-RANK (K=8 cyclic basis coefficients per year, NOT a
// weekly latent state) deviation from its R0 backbone. Supersedes v5.0a
// (weekly AR(1), N=573 sequentially-correlated latent parameters), which
// was archived as computationally infeasible (~85% of transitions hit
// max_treedepth=14, real divergences during warmup, projected many-hour
// runtime) -- see 29_v5_0_dynamic_R/outputs/FAILED_v5_0a_.../FAILURE_SUMMARY.md.
// v5.0b tests the SAME scientific hypothesis (can a limited short-timescale
// R0 deviation reproduce the sharp 2017 epidemic while retaining the good
// 2022 fit?) with only Y*K = 11*8 = 88 latent deviation parameters
// (z_dynamic) plus one global shrinkage scale (tau_dynamic), instead of a
// weekly AR(1) state.
//
// mu_R[t] is EXACTLY v4.9's log_R0(t) formula (unchanged). The dynamic
// deviation is:
//   delta_R[t] = sum_k B_dynamic[t,k] * theta[year_id[t], k]
//   theta[y,k] = tau_dynamic * z_dynamic[y,k]      (non-centred)
// where B_dynamic (data, built in R by
// 29_v5_0_dynamic_R/scripts/04_build_lowrank_basis_and_prior_predictive.R)
// is a periodic piecewise-linear ("hat") basis over 8 equally spaced knots
// across the annual cycle, PRE-CENTRED WITHIN EACH CALENDAR YEAR so that
// ANY theta[y,] choice yields a within-year zero-mean deviation by
// construction -- this is what prevents delta_R from duplicating
// year_effect's or A_year's role (no additional Stan-side centring needed).
//
// Everything else (fixed q, demographic S/U recursion, lifelong immunity,
// GI kernel, 2022 seed mechanism + likelihood exclusion, NB2 observation
// model, phi_obs prior, harmonic/A_year structure, imports, direct
// Juazeiro serology treatment, HMC configuration) is UNCHANGED from v4.9.

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
  real<lower=0> year_effect_prior_sd;
  int<lower=1, upper=N> t_sero;
  int<lower=0> sero_pos;
  int<lower=1> sero_n;
  real<lower=0> kappa_sero;
  real<lower=0, upper=1> q;

  array[N] int<lower=0, upper=1> is_seed;
  vector<lower=0>[N] X_seed;
  int<lower=0, upper=N> N_fit;
  array[N_fit] int<lower=1, upper=N> fit_index;

  int<lower=1> K;              // number of cyclic basis knots (=8)
  matrix[N, K] B_dynamic;      // pre-centred-within-year cyclic basis (fixed data)
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
  real beta_sin1;
  real beta_cos1;
  real beta_sin2;
  real beta_cos2;
  real<lower=0> sigma_season_year;
  vector[Y] z_season_year;
  real<lower=0> phi_obs;

  // Low-rank dynamic R0 deviation.
  real<lower=0> tau_dynamic;
  matrix[Y, K] z_dynamic;
}

transformed parameters {
  vector[Y] year_effect;
  vector[Y] A_year;
  matrix[Y, K] theta;
  vector[N] mu_R;
  vector[N] delta_R;
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
  A_year = exp(sigma_season_year * z_season_year);
  theta = tau_dynamic * z_dynamic; // non-centred

  for (t in 1:N) {
    mu_R[t] = alpha_R + year_effect[year_id[t]] +
      A_year[year_id[t]] * (beta_sin1 * seasonal_sin[t] + beta_cos1 * seasonal_cos[t]) +
      beta_sin2 * seasonal_sin2[t] + beta_cos2 * seasonal_cos2[t];
    delta_R[t] = 0;
    for (k in 1:K) {
      delta_R[t] += B_dynamic[t, k] * theta[year_id[t], k];
    }
  }

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

    log_R0[t] = mu_R[t] + delta_R[t];
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
  beta_sin1 ~ normal(0, 0.25);
  beta_cos1 ~ normal(0, 0.25);
  beta_sin2 ~ normal(0, 0.10);
  beta_cos2 ~ normal(0, 0.10);
  sigma_season_year ~ normal(0, 0.20);
  z_season_year ~ std_normal();
  phi_obs ~ gamma(2, 0.1); // UNCHANGED from v4.9

  // Low-rank dynamic-R priors -- pre-registered via prior-predictive
  // simulation (04_build_lowrank_basis_and_prior_predictive.R) BEFORE this
  // model was fit; NOT tuned after inspecting the 2017 fit.
  tau_dynamic ~ normal(0, 0.15); // half-normal (tau_dynamic > 0)
  to_vector(z_dynamic) ~ std_normal();

  C[fit_index] ~ neg_binomial_2(expected_reported_cases[fit_index], phi_obs);
  sero_pos ~ beta_binomial(sero_n, alpha_sero, beta_sero); // direct anchor, unchanged from v4.9
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
