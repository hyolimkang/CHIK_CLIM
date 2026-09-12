// v4.5: FINAL 2015-2019 identifiability test (mode audit decision-gate
// Outcome B+D follow-up). The latent WEEKLY renewal process (X[t], S[t],
// U[t], R0_t, R_eff_t, demographics, generation interval, seasonality,
// annual effects, imports) is IDENTICAL to renewal_ceara_v4_3_weak_q.stan.
// The q definition and its weak logit-normal calibration prior are also
// unchanged, and q remains ONE freely estimated parameter (not fixed).
//
// The ONLY change is the case observation model, split into two parts to
// remove the weekly-independence pseudo-replication identified in the
// mode-validity audit (weekly residual ACF ~0.81-0.86 at lag 1; 4-week
// block aggregation collapsed the pathological bimodality but left block
// residual ACF ~0.8, i.e. still severe):
//
//   1. Absolute burden: for each of the Y (=5) calendar years,
//        year_total_C[y] ~ NegBinomial(q * X_total_year[y], phi_burden)
//      -- only Y=5 nominal observations, matching how much independent
//      information the data can plausibly offer about ascertainment.
//
//   2. Weekly temporal shape CONDITIONAL ON THE OBSERVED yearly total:
//        weekly counts within year y ~ DirichletMultinomial(year_total_C[y],
//                                                             pi_y, kappa_shape)
//      where pi_y = X[t] / sum(X[t] over year y). This uses the full
//      within-year temporal pattern of X (informative about R0_t/seasonality)
//      without pretending each week is an independent NB draw -- the DM
//      contributes exactly one joint log-density term per year (5 total),
//      correctly modeling the negative covariance among weeks that must sum
//      to a fixed total.
//
// Juazeiro 2018 serology (103/404, beta-binomial, kappa_sero=50) is
// included unconditionally with the same conservative treatment as v4.2-v4.4.

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
  array[N] int<lower=1, upper=Y> year_id; // must be non-decreasing in t (chronological data)
  vector[N] seasonal_sin;
  vector[N] seasonal_cos;
  real<lower=0> imports_per_week;
  real<lower=0> year_effect_prior_sd;
  int<lower=1, upper=N> t_sero;
  int<lower=0> sero_pos;
  int<lower=1> sero_n;
  real<lower=0> kappa_sero;
}

transformed data {
  vector[N] reconciliation;
  real eta_q_prior_mean = logit(0.05);
  array[Y] int<lower=1, upper=N> year_start;
  array[Y] int<lower=1, upper=N> year_end;
  array[Y] int<lower=1> year_n_weeks;
  array[Y] int<lower=0> year_total_C;

  for (t in 1:N) {
    reconciliation[t] = N_end[t] - N_start[t] - births[t] + deaths[t];
  }

  {
    int y_cur = 1;
    year_start[1] = 1;
    for (t in 2:N) {
      if (year_id[t] != y_cur) {
        year_end[y_cur] = t - 1;
        y_cur += 1;
        year_start[y_cur] = t;
      }
    }
    year_end[Y] = N;
    for (y in 1:Y) {
      year_n_weeks[y] = year_end[y] - year_start[y] + 1;
      year_total_C[y] = sum(C[year_start[y]:year_end[y]]);
    }
  }
}

parameters {
  real alpha_R;
  vector[Y] z_year;
  real beta_sin;
  real beta_cos;
  real eta_q;
  real<lower=0> phi_burden;   // yearly burden NB dispersion, freshly estimated
  real<lower=0> kappa_shape;  // within-year Dirichlet-multinomial concentration, freshly estimated
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
  real p_state_sero_at_anchor;
  real p_sero_safe;
  real alpha_sero;
  real beta_sero;
  real<lower=0, upper=1> q = inv_logit(eta_q);
  vector[Y] X_total_year;
  vector[Y] expected_burden;

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

  for (y in 1:Y) {
    X_total_year[y] = sum(X[year_start[y]:year_end[y]]);
  }
  expected_burden = q * X_total_year + 1e-9;
}

model {
  alpha_R ~ normal(log(1.2), 0.5);
  z_year ~ std_normal();
  beta_sin ~ normal(0, 0.25);
  beta_cos ~ normal(0, 0.25);
  eta_q ~ normal(eta_q_prior_mean, 1.25);
  phi_burden ~ gamma(2, 0.1);
  kappa_shape ~ gamma(2, 0.1);

  // 1. Absolute yearly burden (Y=5 terms only).
  year_total_C ~ neg_binomial_2(expected_burden, phi_burden);

  // 2. Weekly shape within each year, conditional on the OBSERVED yearly
  // total (manual Dirichlet-multinomial log density: one joint term per year).
  for (y in 1:Y) {
    int Ty = year_n_weeks[y];
    vector[Ty] Xy = segment(X, year_start[y], Ty);
    vector[Ty] pi_y = Xy / sum(Xy);
    vector[Ty] alpha_y = kappa_shape * pi_y;
    real sum_y = year_total_C[y];
    real log_dm = lgamma(kappa_shape) - lgamma(sum_y + kappa_shape);
    for (k in 1:Ty) {
      int idx = year_start[y] + k - 1;
      log_dm += lgamma(C[idx] + alpha_y[k]) - lgamma(alpha_y[k]);
    }
    target += log_dm;
  }

  sero_pos ~ beta_binomial(sero_n, alpha_sero, beta_sero);
}

generated quantities {
  array[Y] int C_total_pred;
  vector[Y] log_lik_burden;
  vector[Y] log_lik_shape;
  real log_lik_serology;
  int sero_pred;
  real log_prior_q;
  array[N] int weekly_shape_pred = rep_array(0, N);

  for (y in 1:Y) {
    int Ty = year_n_weeks[y];
    vector[Ty] Xy = segment(X, year_start[y], Ty);
    vector[Ty] pi_y = Xy / sum(Xy);
    vector[Ty] alpha_y = kappa_shape * pi_y;
    real sum_y = year_total_C[y];
    real log_dm = lgamma(kappa_shape) - lgamma(sum_y + kappa_shape);

    C_total_pred[y] = neg_binomial_2_rng(expected_burden[y], phi_burden);
    log_lik_burden[y] = neg_binomial_2_lpmf(year_total_C[y] | expected_burden[y], phi_burden);

    for (k in 1:Ty) {
      int idx = year_start[y] + k - 1;
      log_dm += lgamma(C[idx] + alpha_y[k]) - lgamma(alpha_y[k]);
    }
    log_lik_shape[y] = log_dm;

    {
      vector[Ty] theta_y = dirichlet_rng(alpha_y);
      array[Ty] int counts_y = multinomial_rng(theta_y, year_total_C[y]);
      for (k in 1:Ty) {
        weekly_shape_pred[year_start[y] + k - 1] = counts_y[k];
      }
    }
  }

  log_lik_serology = beta_binomial_lpmf(sero_pos | sero_n, alpha_sero, beta_sero);
  sero_pred = beta_binomial_rng(sero_n, alpha_sero, beta_sero);
  log_prior_q = normal_lpdf(eta_q | eta_q_prior_mean, 1.25);
}
