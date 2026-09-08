// renewal_ceara_v3_1_minimal_spatial.stan
//
// Six-unit stratified renewal model for Ceara, 2015--2019. It retains the
// explicit v2.2 S/U demographic accounting but replaces 184 weakly identified
// municipality trajectories with three named municipalities and three fixed
// burden strata. The statewide susceptibility trajectory is derived only by
// aggregating these six unit states.

data {
  int<lower=6, upper=6> K;
  int<lower=2> N;
  int<lower=1> G;
  int<lower=1> seed_weeks;
  array[K, N] int<lower=0> C;
  simplex[G] w;

  matrix<lower=0>[K, N] N_start;
  matrix<lower=0>[K, N] N_end;
  matrix<lower=0>[K, N] births;
  matrix<lower=0>[K, N] deaths;

  vector[seed_weeks] log_seed_prior_mean;
  real<lower=0> seed_prior_sd;

  int<lower=1> J;
  array[J] int<lower=1, upper=K> sero_unit;
  array[J] int<lower=1, upper=N> sero_window_start;
  array[J] int<lower=1, upper=N> sero_window_end;
  array[J] int<lower=0> sero_positive;
  array[J] int<lower=1> sero_n;
}

transformed data {
  matrix[K, N] reconciliation;

  if (seed_weeks < G || seed_weeks >= N) {
    reject("Require G <= seed_weeks < N");
  }
  for (j in 1:J) {
    if (sero_window_start[j] > sero_window_end[j]) {
      reject("Serology window start occurs after its end for survey ", j);
    }
  }
  for (k in 1:K) {
    for (t in 1:N) {
      reconciliation[k, t] = N_end[k, t] - N_start[k, t]
                              - births[k, t] + deaths[k, t];
    }
  }
}

parameters {
  real<lower=0> sigma_R;
  vector[N - seed_weeks] mu_R;
  real<lower=0> sigma_unit;
  vector[K] z_unit;

  matrix[K, seed_weeks] log_seed_hazard;

  real<lower=0, upper=1> p_symp;
  real<lower=0, upper=1> rho_sym;
  real<lower=0> phi_obs;
}

transformed parameters {
  vector[K] unit_offset;
  matrix[K, N] log_R0;
  matrix[K, N] R0;
  matrix[K, N] R_eff;
  matrix[K, N] X;
  matrix[K, N] S;
  matrix[K, N] U;
  matrix[K, N] S_prop;
  matrix[K, N] immune_prop;
  matrix[K, N] expected_reported_cases;
  vector<lower=0, upper=1>[J] p_sero;
  real<lower=0, upper=1> overall_detection;

  unit_offset = sigma_unit * z_unit;
  overall_detection = p_symp * rho_sym;

  for (k in 1:K) {
    S[k, 1] = N_start[k, 1];
    U[k, 1] = 0;

    for (t in 1:N) {
      real population_t = S[k, t] + U[k, t];

      if (S[k, t] < 0 || U[k, t] < 0 ||
          fabs(population_t - N_start[k, t]) > 1e-6 * fmax(1.0, N_start[k, t])) {
        reject("Invalid demographic state for unit ", k, " week ", t);
      }

      if (t <= seed_weeks) {
        log_R0[k, t] = 0;
        R0[k, t] = 1;
        X[k, t] = S[k, t] * (-expm1(-exp(log_seed_hazard[k, t])));
      } else {
        real infectiousness = 0;
        real force_of_infection;
        for (g in 1:G) {
          infectiousness += w[g] * X[k, t - g];
        }
        log_R0[k, t] = mu_R[t - seed_weeks] + unit_offset[k];
        R0[k, t] = exp(log_R0[k, t]);
        force_of_infection = R0[k, t] * infectiousness / population_t;
        X[k, t] = S[k, t] * (-expm1(-force_of_infection));
      }

      expected_reported_cases[k, t] = overall_detection * X[k, t] + 1e-9;

      if (t < N) {
        real s_after = S[k, t] - X[k, t];
        real u_after = U[k, t] + X[k, t];
        real frac_s = s_after / (s_after + u_after);
        real deaths_s = deaths[k, t] * frac_s;
        real deaths_u = deaths[k, t] - deaths_s;
        real reconciliation_s = reconciliation[k, t] * frac_s;
        real reconciliation_u = reconciliation[k, t] - reconciliation_s;

        S[k, t + 1] = s_after + births[k, t] - deaths_s + reconciliation_s;
        U[k, t + 1] = u_after - deaths_u + reconciliation_u;

        if (S[k, t + 1] < 0 || U[k, t + 1] < 0) {
          reject("Demographic accounting produced a negative state for unit ",
                 k, " week ", t + 1);
        }
      }
    }

    for (t in 1:seed_weeks) {
      log_R0[k, t] = log_R0[k, seed_weeks + 1];
      R0[k, t] = R0[k, seed_weeks + 1];
    }
    for (t in 1:N) {
      R_eff[k, t] = R0[k, t] * S[k, t] / (S[k, t] + U[k, t]);
      S_prop[k, t] = S[k, t] / (S[k, t] + U[k, t]);
      immune_prop[k, t] = U[k, t] / (S[k, t] + U[k, t]);
    }
  }

  for (j in 1:J) {
    real immune_sum = 0;
    for (t in sero_window_start[j]:sero_window_end[j]) {
      immune_sum += immune_prop[sero_unit[j], t];
    }
    p_sero[j] = immune_sum / (sero_window_end[j] - sero_window_start[j] + 1);
  }
}

model {
  sigma_R ~ normal(0, 0.10);
  mu_R[1] ~ normal(log(1.0), 0.5);
  if (N - seed_weeks > 1) {
    for (t in 2:(N - seed_weeks)) {
      mu_R[t] ~ normal(mu_R[t - 1], sigma_R);
    }
  }

  sigma_unit ~ normal(0, 0.25);
  z_unit ~ normal(0, 1);

  for (k in 1:K) {
    for (t in 1:seed_weeks) {
      log_seed_hazard[k, t] ~ normal(
        log_seed_prior_mean[t] - log(N_start[k, t]), seed_prior_sd
      );
    }
  }

  p_symp ~ beta(30, 28);
  rho_sym ~ beta(20, 60);
  phi_obs ~ gamma(2, 0.1);

  for (k in 1:K) {
    for (t in 1:N) {
      C[k, t] ~ neg_binomial_2(expected_reported_cases[k, t], phi_obs);
    }
  }
  for (j in 1:J) {
    sero_positive[j] ~ binomial(sero_n[j], p_sero[j]);
  }
}

generated quantities {
  array[K, N] int C_pred;
  matrix[K, N] log_lik_cases;
  array[N] int C_ceara;
  array[N] int C_ceara_pred;
  vector[N] X_ceara_total;
  vector[N] S_ceara_prop;
  vector[N] immune_ceara_prop;
  vector[N] R0_shared;
  vector[N] R_eff_ceara;
  array[J] int sero_positive_pred;
  vector[J] log_lik_serology;

  for (t in 1:N) {
    int observed_total = 0;
    int predicted_total = 0;
    real s_total = 0;
    real u_total = 0;
    real x_total = 0;
    real reff_numerator = 0;
    for (k in 1:K) {
      C_pred[k, t] = neg_binomial_2_rng(expected_reported_cases[k, t], phi_obs);
      log_lik_cases[k, t] = neg_binomial_2_lpmf(
        C[k, t] | expected_reported_cases[k, t], phi_obs
      );
      observed_total += C[k, t];
      predicted_total += C_pred[k, t];
      s_total += S[k, t];
      u_total += U[k, t];
      x_total += X[k, t];
      reff_numerator += R_eff[k, t] * (S[k, t] + U[k, t]);
    }
    C_ceara[t] = observed_total;
    C_ceara_pred[t] = predicted_total;
    X_ceara_total[t] = x_total;
    S_ceara_prop[t] = s_total / (s_total + u_total);
    immune_ceara_prop[t] = u_total / (s_total + u_total);
    R_eff_ceara[t] = reff_numerator / (s_total + u_total);
    if (t <= seed_weeks) {
      R0_shared[t] = 1;
    } else {
      R0_shared[t] = exp(mu_R[t - seed_weeks]);
    }
  }
  for (t in 1:seed_weeks) {
    R0_shared[t] = R0_shared[seed_weeks + 1];
  }
  for (j in 1:J) {
    sero_positive_pred[j] = binomial_rng(sero_n[j], p_sero[j]);
    log_lik_serology[j] = binomial_lpmf(sero_positive[j] | sero_n[j], p_sero[j]);
  }
}

