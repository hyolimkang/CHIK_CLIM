// renewal_ceara_v3_0_spatial.stan
//
// Municipality-level renewal model for Ceara, 2015--2019. This is a new
// spatial version of v2.2: each municipality has its own S/U/X trajectory,
// while transmission and reporting are strongly pooled. The demographic
// recursion is the v2.2 recursion applied independently by municipality.

data {
  int<lower=2> M;
  int<lower=2> N;
  int<lower=1> G;
  int<lower=1> seed_weeks;
  array[M, N] int<lower=0> C;
  simplex[G] w;

  matrix<lower=0>[M, N] N_start;
  matrix<lower=0>[M, N] N_end;
  matrix<lower=0>[M, N] births;
  matrix<lower=0>[M, N] deaths;

  vector[seed_weeks] log_seed_prior_mean;
  real<lower=0> seed_prior_sd;

  int<lower=1> J;
  array[J] int<lower=1, upper=M> sero_muni;
  array[J] int<lower=1, upper=N> sero_window_start;
  array[J] int<lower=1, upper=N> sero_window_end;
  array[J] int<lower=0> sero_positive;
  array[J] int<lower=1> sero_n;
}

transformed data {
  matrix[M, N] net_migration_reconciliation;

  if (seed_weeks < G || seed_weeks >= N) {
    reject("Require G <= seed_weeks < N");
  }
  for (j in 1:J) {
    if (sero_window_start[j] > sero_window_end[j]) {
      reject("Serology window start occurs after its end for survey ", j);
    }
  }
  for (m in 1:M) {
    for (t in 1:N) {
      net_migration_reconciliation[m, t] = N_end[m, t] - N_start[m, t]
                                               - births[m, t] + deaths[m, t];
    }
  }
}

parameters {
  // Shared temporal R0 process and strongly pooled municipality offsets.
  real<lower=0> sigma_R;
  vector[N - seed_weeks] mu_R;
  real<lower=0> sigma_muni_R;
  vector[M] z_muni_R;

  // The only municipality-time-specific transmission parameters are the
  // regularised initial seed hazards. Later infections follow renewal.
  matrix[M, seed_weeks] log_seed_hazard;

  real<lower=0, upper=1> p_symp;
  real<lower=0, upper=1> rho_sym_global;
  real<lower=0> phi_obs;
}

transformed parameters {
  vector[M] municipality_R_offset;
  matrix[M, N] log_R0;
  matrix[M, N] R0;
  matrix[M, N] R_eff;
  matrix[M, N] X;
  matrix[M, N] S;
  matrix[M, N] U;
  matrix[M, N] S_prop;
  matrix[M, N] immune_prop;
  matrix[M, N] expected_reported_cases;
  vector<lower=0, upper=1>[J] p_sero;
  real<lower=0, upper=1> overall_detection;

  municipality_R_offset = sigma_muni_R * z_muni_R;
  overall_detection = p_symp * rho_sym_global;

  for (m in 1:M) {
    S[m, 1] = N_start[m, 1];
    U[m, 1] = 0;

    for (t in 1:N) {
      real population_t = S[m, t] + U[m, t];

      if (S[m, t] < 0 || U[m, t] < 0 ||
          fabs(population_t - N_start[m, t]) > 1e-6 * fmax(1.0, N_start[m, t])) {
        reject("Invalid demographic state for municipality ", m, " week ", t,
               ": S=", S[m, t], " U=", U[m, t],
               " N_start=", N_start[m, t]);
      }

      if (t <= seed_weeks) {
        // R0 is unidentified until a complete local generation interval is
        // observed. These display values are filled from the first
        // estimable week below.
        log_R0[m, t] = 0;
        R0[m, t] = 1;
        X[m, t] = S[m, t] * (-expm1(-exp(log_seed_hazard[m, t])));
      } else {
        real infectiousness = 0;
        real force_of_infection;
        for (g in 1:G) {
          infectiousness += w[g] * X[m, t - g];
        }
        log_R0[m, t] = mu_R[t - seed_weeks] + municipality_R_offset[m];
        R0[m, t] = exp(log_R0[m, t]);
        force_of_infection = R0[m, t] * infectiousness / population_t;
        X[m, t] = S[m, t] * (-expm1(-force_of_infection));
      }

      expected_reported_cases[m, t] = overall_detection * X[m, t] + 1e-9;

      if (t < N) {
        real s_after = S[m, t] - X[m, t];
        real u_after = U[m, t] + X[m, t];
        real frac_s = s_after / (s_after + u_after);
        real deaths_s = deaths[m, t] * frac_s;
        real deaths_u = deaths[m, t] - deaths_s;
        real migration_s = net_migration_reconciliation[m, t] * frac_s;
        real migration_u = net_migration_reconciliation[m, t] - migration_s;

        S[m, t + 1] = s_after + births[m, t] - deaths_s + migration_s;
        U[m, t + 1] = u_after - deaths_u + migration_u;

        if (S[m, t + 1] < 0 || U[m, t + 1] < 0) {
          reject("Demographic accounting produced a negative state for municipality ",
                 m, " week ", t + 1);
        }
      }
    }

    for (t in 1:seed_weeks) {
      log_R0[m, t] = log_R0[m, seed_weeks + 1];
      R0[m, t] = R0[m, seed_weeks + 1];
    }
    for (t in 1:N) {
      R_eff[m, t] = R0[m, t] * S[m, t] / (S[m, t] + U[m, t]);
      S_prop[m, t] = S[m, t] / (S[m, t] + U[m, t]);
      immune_prop[m, t] = U[m, t] / (S[m, t] + U[m, t]);
    }
  }

  for (j in 1:J) {
    real immune_sum = 0;
    for (t in sero_window_start[j]:sero_window_end[j]) {
      immune_sum += immune_prop[sero_muni[j], t];
    }
    p_sero[j] = immune_sum / (sero_window_end[j] - sero_window_start[j] + 1);
  }
}

model {
  // Same temporal smoothness scale as v2.2, now shared across municipalities.
  sigma_R ~ normal(0, 0.10);
  mu_R[1] ~ normal(log(1.0), 0.5);
  if (N - seed_weeks > 1) {
    for (t in 2:(N - seed_weeks)) {
      mu_R[t] ~ normal(mu_R[t - 1], sigma_R);
    }
  }

  // New, intentionally tight scale prior for the new municipality intercepts.
  sigma_muni_R ~ normal(0, 0.25);
  z_muni_R ~ normal(0, 1);

  for (m in 1:M) {
    for (t in 1:seed_weeks) {
      log_seed_hazard[m, t] ~ normal(
        log_seed_prior_mean[t] - log(N_start[m, t]), seed_prior_sd
      );
    }
  }

  p_symp ~ beta(30, 28);
  rho_sym_global ~ beta(20, 60);
  phi_obs ~ gamma(2, 0.1);

  for (m in 1:M) {
    for (t in 1:N) {
      C[m, t] ~ neg_binomial_2(expected_reported_cases[m, t], phi_obs);
    }
  }
  for (j in 1:J) {
    sero_positive[j] ~ binomial(sero_n[j], p_sero[j]);
  }
}

generated quantities {
  array[M, N] int C_pred;
  matrix[M, N] log_lik_cases;
  array[N] int C_ceara;
  array[N] int C_ceara_pred;
  vector[N] X_ceara_total;
  vector[N] S_ceara_prop;
  vector[N] immune_ceara_prop;
  array[J] int sero_positive_pred;
  vector[J] log_lik_serology;

  for (t in 1:N) {
    int observed_total = 0;
    int predicted_total = 0;
    real s_total = 0;
    real u_total = 0;
    real x_total = 0;
    for (m in 1:M) {
      C_pred[m, t] = neg_binomial_2_rng(expected_reported_cases[m, t], phi_obs);
      log_lik_cases[m, t] = neg_binomial_2_lpmf(
        C[m, t] | expected_reported_cases[m, t], phi_obs
      );
      observed_total += C[m, t];
      predicted_total += C_pred[m, t];
      s_total += S[m, t];
      u_total += U[m, t];
      x_total += X[m, t];
    }
    C_ceara[t] = observed_total;
    C_ceara_pred[t] = predicted_total;
    X_ceara_total[t] = x_total;
    S_ceara_prop[t] = s_total / (s_total + u_total);
    immune_ceara_prop[t] = u_total / (s_total + u_total);
  }
  for (j in 1:J) {
    sero_positive_pred[j] = binomial_rng(sero_n[j], p_sero[j]);
    log_lik_serology[j] = binomial_lpmf(sero_positive[j] | sero_n[j], p_sero[j]);
  }
}

