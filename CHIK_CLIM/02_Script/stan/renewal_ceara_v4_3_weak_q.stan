// Ceará v4.3: short-period (2015-2019) q-identifiability calibration.
//
// EXACT scientific backbone of the v4.2/Stage-B4 model (verified against
// 22_v4_3_short_2015_2019_q_calibration/archive/
// renewal_ceara_v4_2_q_juazeiro_serology_ARCHIVED.stan): latent infection
// renewal, susceptible/immune bookkeeping, fixed generation interval,
// fixed imports_per_week, fixed year_effect_prior_sd annual effects, same
// seasonal terms, same NB case likelihood, same beta-binomial serology
// observation model. No AR(1), no fitted seed, no climate, no spatial
// effects, no sigma_year, no time-varying q -- the only scientific change
// is (a) the observation window (truncated to 2015-2019 by the CALLER via
// N/weekly data, not by this file), (b) the reporting-fraction prior, and
// (c) an optional on/off switch for the serology likelihood so ONE file
// serves both Fit A (cases only) and Fit B (cases + Juazeiro).
//
// (b) The legacy q ~ Beta(16.2,108.7) prior is REMOVED here: the repo audit
// found no documented external evidence tying that prior to the actual
// observational endpoint of C[t] (see 21_v4_2_.../DEVELOPMENT_LOG.md and
// the conversation record). Replaced by a deliberately WEAK calibration
// prior on the logit scale: eta_q ~ normal(logit(0.05), 1.25), q =
// inv_logit(eta_q). This is wide enough to give real prior mass to q
// anywhere from ~0.02 to >0.20 (quantified in
// scripts/01_plot_weak_q_prior.R) and is NOT tuned after seeing the
// posterior.
//
// (c) use_serology (0/1, data) turns the ONE Juazeiro serology term on or
// off without duplicating the whole file. p_state_sero_at_anchor and its
// generated-quantities diagnostics are always computed either way (harmless
// -- a pure function of the deterministic S/U trajectory), so Fit A's
// implied 2018 seroprevalence remains directly comparable to Fit B's, even
// though only Fit B's likelihood is informed by the observed 103/404.

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
  real<lower=0> imports_per_week;
  real<lower=0> year_effect_prior_sd;
  int<lower=1, upper=N> t_sero;
  int<lower=0> sero_pos;
  int<lower=1> sero_n;
  real<lower=0> kappa_sero;
  int<lower=0, upper=1> use_serology; // Fit A = 0, Fit B = 1
}

transformed data {
  vector[N] reconciliation;
  real eta_q_prior_mean = logit(0.05);
  for (t in 1:N) {
    reconciliation[t] = N_end[t] - N_start[t] - births[t] + deaths[t];
  }
}

parameters {
  real alpha_R;
  vector[Y] z_year;
  real beta_sin;
  real beta_cos;
  real<lower=0> phi_obs;
  real eta_q;
}

transformed parameters {
  vector[Y] year_effect;
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
  real<lower=0, upper=1> q = inv_logit(eta_q);

  year_effect = year_effect_prior_sd * (z_year - mean(z_year));
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
      beta_sin * seasonal_sin[t] + beta_cos * seasonal_cos[t];
    R0_t[t] = exp(log_R0[t]);
    force_of_infection[t] = R0_t[t] * infectiousness / total + imports_per_week / total;
    X[t] = S[t] * (-expm1(-force_of_infection[t]));
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
  beta_sin ~ normal(0, 0.25);
  beta_cos ~ normal(0, 0.25);
  phi_obs ~ gamma(2, 0.1);
  eta_q ~ normal(eta_q_prior_mean, 1.25); // weak calibration prior -- Section 2
  C ~ neg_binomial_2(expected_reported_cases, phi_obs);
  if (use_serology) {
    sero_pos ~ beta_binomial(sero_n, alpha_sero, beta_sero);
  }
}

generated quantities {
  array[N] int C_pred;
  vector[N] log_lik_cases;
  real log_lik_serology;
  int sero_pred;
  real log_prior_q;
  for (t in 1:N) {
    C_pred[t] = neg_binomial_2_rng(expected_reported_cases[t], phi_obs);
    log_lik_cases[t] = neg_binomial_2_lpmf(C[t] | expected_reported_cases[t], phi_obs);
  }
  log_lik_serology = beta_binomial_lpmf(sero_pos | sero_n, alpha_sero, beta_sero);
  sero_pred = beta_binomial_rng(sero_n, alpha_sero, beta_sero);
  log_prior_q = normal_lpdf(eta_q | eta_q_prior_mean, 1.25);
}
