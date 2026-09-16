// v5.0a = v4.9 + a strongly regularized weekly AR(1) residual on log R0(t),
// keeping the direct Juazeiro serology anchor unchanged ("direct_serology"
// -- see v5.0b for the site-to-state transport-offset variant, built only
// if the Phase-1 serology ablation shows serology materially suppresses
// the 2017 burden).
//
// v4.9 showed that a smooth two-harmonic seasonal curve, even with
// year-specific amplitude (A_year), cannot sharpen the 2017 epidemic peak
// (posterior median ~53% of the observed ~10,000/week peak) while it DOES
// reproduce the 2022 recurrence and 2018 serology almost exactly at
// q=0.05. v5.0a tests whether a regularized short-timescale departure
// from that smooth backbone -- NOT a free weekly random walk -- can close
// the 2017 gap without degrading 2022 or requiring implausible phi_obs.
//
// mu_R[t] is EXACTLY v4.9's log_R0(t) formula (renamed, unchanged):
//   mu_R[t] = alpha_R + year_effect[year_id[t]]
//             + A_year[year_id[t]] * (beta_sin1*sin1[t] + beta_cos1*cos1[t])
//             + beta_sin2*sin2[t] + beta_cos2*cos2[t]
// log_R0[t] = mu_R[t] + delta_R[t], where delta_R is a non-centred
// stationary AR(1) process, RE-CENTRED TO ZERO MEAN WITHIN EACH CALENDAR
// YEAR so it cannot reproduce year-level intercepts (year_effect's job)
// or the broad seasonal amplitude (A_year's job) -- it can only represent
// shorter-timescale within-year departures from that backbone.
//
// Everything else (fixed q, demographic S/U recursion, lifelong immunity,
// 2022 seed mechanism + likelihood exclusion, NB2 observation model,
// phi_obs prior, harmonic/A_year structure, imports, GI kernel, Juazeiro
// serology treatment) is UNCHANGED from v4.9.

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
  array[N] int<lower=1, upper=Y> year_id; // must be non-decreasing in t
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
}

transformed data {
  vector[N] reconciliation;
  array[Y] int<lower=1, upper=N> year_start;
  array[Y] int<lower=1, upper=N> year_end;

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

  // Dynamic short-timescale R0 residual (non-centred stationary AR(1)).
  real<lower=0, upper=1> rho_R;
  real<lower=0> sigma_R_dynamic;
  vector[N] z_delta;
}

transformed parameters {
  vector[Y] year_effect;
  vector[Y] A_year;
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

  for (t in 1:N) {
    mu_R[t] = alpha_R + year_effect[year_id[t]] +
      A_year[year_id[t]] * (beta_sin1 * seasonal_sin[t] + beta_cos1 * seasonal_cos[t]) +
      beta_sin2 * seasonal_sin2[t] + beta_cos2 * seasonal_cos2[t];
  }

  {
    vector[N] delta_raw;
    real stationary_sd = sigma_R_dynamic / sqrt(1 - rho_R * rho_R);
    delta_raw[1] = z_delta[1] * stationary_sd;
    for (t in 2:N) {
      delta_raw[t] = rho_R * delta_raw[t - 1] + sigma_R_dynamic * z_delta[t];
    }
    // Re-centre to zero mean WITHIN each calendar year so delta_R cannot
    // absorb year_effect's or A_year's role (identifiability, Section
    // "IDENTIFIABILITY BETWEEN A_year AND delta_R").
    for (y in 1:Y) {
      real year_mean = mean(delta_raw[year_start[y]:year_end[y]]);
      for (t in year_start[y]:year_end[y]) {
        delta_R[t] = delta_raw[t] - year_mean;
      }
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

  // Dynamic-R priors -- pre-registered here, NOT tuned after inspecting
  // the 2017 fit (Section "PRIORS"/"Do NOT tune sigma_R_dynamic").
  rho_R ~ beta(8, 2);              // mean 0.8: favours substantial positive autocorrelation
  sigma_R_dynamic ~ normal(0, 0.15); // conservative half-normal (sigma_R_dynamic > 0)
  z_delta ~ std_normal();

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
