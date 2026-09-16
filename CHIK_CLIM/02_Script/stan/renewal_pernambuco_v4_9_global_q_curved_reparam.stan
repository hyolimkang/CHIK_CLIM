// Pernambuco v4.9 global-q: FREE quadratic ridge reparameterisation
// (Model D). Supersedes the affine-only reparameterisation, which left 26
// post-warmup divergences concentrated at the low-q/high-depletion tail
// where the true alpha_R vs logit_q relationship is visibly curved, not
// linear.
//
// This preserves the frozen v4.9 epidemiological model exactly: renewal
// transmission, demographic S/U recursion, lifelong immunity, NB2 cases,
// harmonic/year-specific seasonality, one global q, and the geographic-
// offset serology observation layer are unchanged. The Y-1 orthonormal
// sum-to-zero year-effect basis (from the affine-reparam model) is
// retained unchanged.
//
// Computational change only: alpha_R is now represented as
//   d = logit_q - u0
//   alpha_R = c0 + gamma_curve + b1*d + b2*d^2
// where c0, u0, b1, b2 are FIXED DATA (estimated externally, via OLS on
// non-divergent draws from the affine-reparam full run, re-centred at
// u0=median(logit_q) -- NOT the old affine reference/slope, which the
// user has explicitly noted is only a computational coordinate constant).
// The map (gamma_curve, logit_q) -> (alpha_R, logit_q) is one-to-one with
// unit Jacobian (d alpha_R / d gamma_curve = 1, d logit_q/d logit_q = 1,
// d logit_q/d gamma_curve = 0), so no Jacobian correction is required.
// The ORIGINAL independent priors are applied to alpha_R (post-transform)
// and logit_q; gamma_curve itself receives NO separate prior (same
// convention as the affine model's gamma_ridge).

data {
  int<lower=2> N;
  int<lower=1> G;
  array[N] int<lower=0> C;
  simplex[G] w;
  vector<lower=0>[N] N_start;
  vector<lower=0>[N] N_end;
  vector<lower=0>[N] births;
  vector<lower=0>[N] deaths;
  int<lower=2> Y;
  array[N] int<lower=1, upper=Y> year_id;
  vector[N] seasonal_sin;
  vector[N] seasonal_cos;
  vector[N] seasonal_sin2;
  vector[N] seasonal_cos2;
  real<lower=0> imports_per_week;
  real<lower=0> year_effect_prior_sd;
  matrix[Y, Y - 1] year_contrast_basis;

  array[N] int<lower=0, upper=1> is_seed;
  vector<lower=0>[N] X_seed;
  int<lower=0, upper=N> N_fit;
  array[N_fit] int<lower=1, upper=N> fit_index;

  int<lower=0> J_sero;
  array[J_sero] int<lower=0> sero_n_positive;
  array[J_sero] int<lower=1> sero_n_tested;
  array[J_sero] int<lower=1, upper=N> sero_window_start_idx;
  array[J_sero] int<lower=1> sero_window_n_weeks;
  real<lower=0> sero_geographic_sd;

  real logit_q_prior_mean;
  real<lower=0> logit_q_prior_sd;

  // Fixed quadratic-ridge coordinate constants (NOT estimated in Stan).
  real c0;
  real u0;
  real b1;
  real b2;
}

transformed data {
  vector[N] reconciliation;
  for (t in 1:N) {
    reconciliation[t] = N_end[t] - N_start[t] - births[t] + deaths[t];
  }
}

parameters {
  real gamma_curve;
  vector[Y - 1] z_year_contrast;
  real beta_sin1;
  real beta_cos1;
  real beta_sin2;
  real beta_cos2;
  real<lower=0> sigma_season_year;
  vector[Y] z_season_year;
  real<lower=0> phi_obs;
  vector[J_sero] z_geo;
  real logit_q;
}

transformed parameters {
  real d_ridge;
  real alpha_R;
  real<lower=0, upper=1> q;
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
  vector[J_sero] eta_geo;
  vector[J_sero] p_site_window;
  vector[J_sero] alpha_sero_site;
  vector[J_sero] beta_sero_site;

  // Unit-Jacobian coordinate transform (quadratic in logit_q, linear in
  // gamma_curve): d(alpha_R)/d(gamma_curve) = 1 exactly. The alpha_R prior
  // below is therefore unchanged in its effect on the joint posterior.
  d_ridge = logit_q - u0;
  alpha_R = c0 + gamma_curve + b1 * d_ridge + b2 * square(d_ridge);
  q = inv_logit(logit_q);
  year_effect = year_effect_prior_sd * (year_contrast_basis * z_year_contrast);
  A_year = exp(sigma_season_year * z_season_year);
  eta_geo = sero_geographic_sd * z_geo;

  S[1] = N_start[1];
  U[1] = 0;
  for (t in 1:N) {
    real total = S[t] + U[t];
    real infectiousness = 0;
    if (S[t] < 0 || U[t] < 0 || fabs(total - N_start[t]) > 1e-6 * fmax(1, N_start[t])) {
      reject("Invalid demographic state at week ", t);
    }
    for (g in 1:G) if (t > g) infectiousness += w[g] * X[t - g];

    log_R0[t] = alpha_R + year_effect[year_id[t]] +
      A_year[year_id[t]] * (beta_sin1 * seasonal_sin[t] + beta_cos1 * seasonal_cos[t]) +
      beta_sin2 * seasonal_sin2[t] + beta_cos2 * seasonal_cos2[t];
    R0_t[t] = exp(log_R0[t]);
    force_of_infection[t] = R0_t[t] * infectiousness / total + imports_per_week / total;
    X[t] = is_seed[t] == 1 ? X_seed[t] : S[t] * (-expm1(-force_of_infection[t]));
    R_eff_t[t] = R0_t[t] * S[t] / total;
    S_prop[t] = S[t] / total;
    immune_prop[t] = U[t] / total;
    expected_reported_cases[t] = q * X[t] + 1e-9;

    if (t < N) {
      real susceptible_after = S[t] - X[t];
      real immune_after = U[t] + X[t];
      real susceptible_fraction = susceptible_after / (susceptible_after + immune_after);
      S[t + 1] = susceptible_after + births[t] - deaths[t] * susceptible_fraction + reconciliation[t] * susceptible_fraction;
      U[t + 1] = immune_after - deaths[t] * (1 - susceptible_fraction) + reconciliation[t] * (1 - susceptible_fraction);
    }
  }

  for (j in 1:J_sero) {
    int start = sero_window_start_idx[j];
    int nwk = sero_window_n_weeks[j];
    real acc = 0;
    for (k in 1:nwk) {
      real p_state_safe = fmin(1 - 1e-9, fmax(1e-9, immune_prop[start + k - 1]));
      acc += inv_logit(logit(p_state_safe) + eta_geo[j]);
    }
    p_site_window[j] = acc / nwk;
    alpha_sero_site[j] = p_site_window[j] * 50;
    beta_sero_site[j] = (1 - p_site_window[j]) * 50;
  }
}

model {
  // These are the original v4.9 priors, deliberately applied to alpha_R/q.
  // gamma_curve itself receives NO separate prior (same convention as the
  // affine model's gamma_ridge) -- the alpha_R prior below, evaluated on
  // the derived alpha_R, plays that role via the unit-Jacobian transform.
  alpha_R ~ normal(log(1.2), 0.5);
  z_year_contrast ~ std_normal();
  beta_sin1 ~ normal(0, 0.25);
  beta_cos1 ~ normal(0, 0.25);
  beta_sin2 ~ normal(0, 0.10);
  beta_cos2 ~ normal(0, 0.10);
  sigma_season_year ~ normal(0, 0.20);
  z_season_year ~ std_normal();
  phi_obs ~ gamma(2, 0.1);
  z_geo ~ std_normal();
  logit_q ~ normal(logit_q_prior_mean, logit_q_prior_sd);

  C[fit_index] ~ neg_binomial_2(expected_reported_cases[fit_index], phi_obs);
  if (J_sero > 0) {
    for (j in 1:J_sero) {
      sero_n_positive[j] ~ beta_binomial(sero_n_tested[j], alpha_sero_site[j], beta_sero_site[j]);
    }
  }
}

generated quantities {
  array[N] int C_pred;
  vector[N] log_lik_cases;
  vector[J_sero] log_lik_serology;
  array[J_sero] int sero_pred;
  vector[J_sero] p_state_window;
  for (t in 1:N) {
    C_pred[t] = neg_binomial_2_rng(expected_reported_cases[t], phi_obs);
    log_lik_cases[t] = neg_binomial_2_lpmf(C[t] | expected_reported_cases[t], phi_obs);
  }
  for (j in 1:J_sero) {
    int start = sero_window_start_idx[j];
    int nwk = sero_window_n_weeks[j];
    real acc = 0;
    for (k in 1:nwk) acc += fmin(1 - 1e-9, fmax(1e-9, immune_prop[start + k - 1]));
    p_state_window[j] = acc / nwk;
    log_lik_serology[j] = beta_binomial_lpmf(sero_n_positive[j] | sero_n_tested[j], alpha_sero_site[j], beta_sero_site[j]);
    sero_pred[j] = beta_binomial_rng(sero_n_tested[j], alpha_sero_site[j], beta_sero_site[j]);
  }
}
