// v4.8 = v4.7 + externally conditioned 2022 epidemic re-seeding.
//
// v4.7 (renewal_ceara_v4_6_fixed_q_sweep.stan run over 2015-2025) showed
// that the continuous renewal process, with q fixed anywhere in 0.05-0.30,
// systematically underpredicts the 2022 recurrence: by late 2021 the
// convolution term (sum_g w[g]*X[t-g]) has decayed to near zero after
// ~4 years of very low incidence, so no fixed R0_t magnitude can reignite
// growth without an external seed -- imports_per_week is a tiny constant
// (1/week) unrelated to any real 2022-specific reintroduction event.
//
// v4.8 tests ONE structural hypothesis: given an externally conditioned
// 2022 reintroduction, can the EXISTING susceptibility (continuous, never
// reset) and the EXISTING renewal/transmission model reproduce the
// subsequent 2022 wave? Everything else is byte-for-byte identical to
// v4.6/v4.7: fixed q, same priors, same GI kernel, same R0(t)
// parameterisation, same demographic S/U bookkeeping, same lifelong
// infection-derived immunity, same Juazeiro serology treatment.
//
// SEED MECHANISM (Option A from the design spec): for a short, externally
// specified window of weeks at the start of the 2022 recurrence
// (is_seed[t]=1), X[t] is NOT computed from the renewal equation -- it is
// set directly from the observed case count under the SAME q mapping used
// everywhere else (X_seed[t] = C[t]/q, computed in R). Those weeks are
// EXCLUDED from the case likelihood (via fit_index) to avoid double use of
// the same observations both to define the latent state and to score it.
// The demographic S/U recursion is completely unaffected: it consumes
// whatever X[t] is available (renewal-derived or seed-conditioned)
// exactly as before, so susceptibility never resets and remains
// continuous through 2022. After the seed window, X[t] is again generated
// entirely by the renewal equation, using the seeded infections (now part
// of the convolution history) to restart local transmission.

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
  real<lower=0, upper=1> q; // FIXED, swept across the same grid as v4.7

  // 2022 seed window (Section: "2022 SEED IMPLEMENTATION"). Built entirely
  // in R (see 01_fit_v4_8.R / make_v4_8_data()) -- never hard-coded here.
  array[N] int<lower=0, upper=1> is_seed;   // 1 for the ~8-week 2022 seed window
  vector<lower=0>[N] X_seed;                // = C[t]/q for seed weeks; 0 elsewhere (unused when is_seed=0)
  int<lower=0, upper=N> N_fit;              // N - (length of seed window)
  array[N_fit] int<lower=1, upper=N> fit_index; // non-seed week indices; seed weeks excluded from the likelihood
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
  real beta_sin;
  real beta_cos;
  real<lower=0> phi_obs;
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

    if (is_seed[t] == 1) {
      // Externally conditioned 2022 seed: X[t] comes from data, not the
      // renewal equation. force_of_infection[t] above is still computed
      // (diagnostic continuity / unused here) but does not drive X[t].
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
  beta_sin ~ normal(0, 0.25);
  beta_cos ~ normal(0, 0.25);
  phi_obs ~ gamma(2, 0.1);
  // Seed weeks are excluded here (fit_index omits them) -- their C[t] was
  // already used to define X_seed, so they cannot also count as evidence
  // for or against the model.
  C[fit_index] ~ neg_binomial_2(expected_reported_cases[fit_index], phi_obs);
  sero_pos ~ beta_binomial(sero_n, alpha_sero, beta_sero);
}

generated quantities {
  array[N] int C_pred;
  vector[N] log_lik_cases;
  real log_lik_serology;
  int sero_pred;
  for (t in 1:N) {
    // Computed for ALL weeks (including seed weeks) for plotting
    // continuity only. Seed-week C_pred/log_lik_cases are NOT part of the
    // model's target and must not be read as predictive performance
    // (expected_reported_cases there is, by construction, ~= C[t]).
    C_pred[t] = neg_binomial_2_rng(expected_reported_cases[t], phi_obs);
    log_lik_cases[t] = neg_binomial_2_lpmf(C[t] | expected_reported_cases[t], phi_obs);
  }
  log_lik_serology = beta_binomial_lpmf(sero_pos | sero_n, alpha_sero, beta_sero);
  sero_pred = beta_binomial_rng(sero_n, alpha_sero, beta_sero);
}
