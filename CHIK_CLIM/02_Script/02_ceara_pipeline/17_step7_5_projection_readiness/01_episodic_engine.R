# STEP 7.5 -- episodic/bursty introduction-process engine. ADDITIVE
# ONLY: does not edit 16_step7_importation_persistence/01_simulation_engine.R.
# Reuses that file's simulate_renewal_core()/compute_R0_future()/
# advance_terminal_state()/generation_weights_future() (sourced
# transitively) and its make_rng_pois_step7() SHARED-UNIFORM-COUPLING
# convention, extended with one new closure type.
#
# MODEL B (episodic): B_t ~ Bernoulli(p); M_t | B_t=1 ~ Poisson(lambda/p);
# M_t | B_t=0 = 0. E[M_t] = lambda exactly, matching Model A's mean.
# Draws TWO uniforms per week (u1 for the burst indicator, u2 for the
# conditional magnitude, u2 ALWAYS drawn even when no burst occurs, to
# keep weekly RNG consumption constant within this model family) -- so
# different p values sharing the same seed see the SAME u1 sequence
# (common random numbers WITHIN Model B's own p-variants), but Model A
# (1 uniform/week) and Model B (2 uniforms/week) are NOT claimed to be
# paired against each other beyond sharing posterior draws, initial
# states, climate, year effects and seasonal amplitude -- disclosed
# explicitly, consistent with the task's "do not claim equal expected
# pressure implies identical/paired outcomes."
make_rng_pois_episodic <- function(p, lambda_mean) {
  week_counter <- 0L
  lambda_burst <- lambda_mean / p
  function(n, lam) {
    week_counter <<- week_counter + 1L
    u1 <- stats::runif(n); burst <- u1 < p
    u2 <- stats::runif(n)
    ifelse(burst, stats::qpois(u2, lambda_burst), 0)
  }
}

#' Copy of 16_step7_importation_persistence's simulate_holdout_step7(),
#' extended with lambda_mode="episodic" (p, lambda_mean). All other
#' modes/behaviour byte-identical to that file -- duplicated here rather
#' than editing it, to leave Step 7's own module completely untouched.
simulate_holdout_step75 <- function(training_bundle, holdout_spec, lambda_mode = c("constant", "episodic"),
                                     lambda_const = NULL, burst_p = NULL, burst_lambda_mean = NULL,
                                     year_effect_mode = c("centered", "uncentered"),
                                     n_draws, n_replicates, seed_offset) {
  lambda_mode <- match.arg(lambda_mode); year_effect_mode <- match.arg(year_effect_mode)
  sd_train <- training_bundle$stan_data; N_train <- sd_train$N
  draws <- rstan::extract(training_bundle$fit, pars = c("alpha_R", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                                                          "sigma_season_year", "tau_year", "beta_T", "beta_P", "phi_obs", "S", "U", "X"), permuted = TRUE)
  n_draws_available <- length(draws$alpha_R)

  holdout_weekly <- load_observed_weekly_step7(holdout_spec$holdout_start, holdout_spec$holdout_end)
  n_holdout <- nrow(holdout_weekly)
  t_index_holdout <- (N_train + 1):(N_train + n_holdout)
  holdout_years <- as.integer(format(holdout_weekly$week_start, "%Y"))
  holdout_year_labels <- sort(unique(holdout_years)); n_holdout_years <- length(holdout_year_labels)
  year_id_holdout <- match(holdout_years, holdout_year_labels)

  obs_clim <- load_observed_climate_step7(holdout_spec$holdout_start, holdout_spec$holdout_end)
  stopifnot(nrow(obs_clim) == n_holdout, all(obs_clim$week_start == holdout_weekly$week_start))
  z_T_holdout <- obs_clim$z_T_anom; z_P_holdout <- obs_clim$z_P_anom

  reconciliation_last_train <- sd_train$N_end[N_train] - sd_train$N_start[N_train] - sd_train$births[N_train] + sd_train$deaths[N_train]

  set.seed(seed_offset)
  draw_ids <- sample.int(n_draws_available, n_draws)
  traj_list <- vector("list", n_draws * n_replicates); slot <- 0L
  year_effect_rows <- vector("list", length(draw_ids))
  for (di in seq_along(draw_ids)) {
    d <- draw_ids[di]
    adv <- advance_terminal_state(draws$S[d, N_train], draws$U[d, N_train], draws$X[d, N_train], sd_train$births[N_train], sd_train$deaths[N_train], reconciliation_last_train)
    X_recent <- draws$X[d, (N_train - G_GEN + 1):N_train]

    set.seed(seed_offset + d)
    z_year_h <- rnorm(n_holdout_years); z_season_h <- rnorm(n_holdout_years)
    year_effect_h <- if (year_effect_mode == "centered") draws$tau_year[d] * (z_year_h - mean(z_year_h)) else draws$tau_year[d] * z_year_h
    A_year_h <- exp(draws$sigma_season_year[d] * z_season_h)
    year_effect_rows[[di]] <- tibble::tibble(posterior_draw_id = d, holdout_year = holdout_year_labels, year_effect = year_effect_h, A_year = A_year_h, tau_used = draws$tau_year[d])

    r0_h <- compute_R0_future(t_index_holdout, year_id_holdout, draws$alpha_R[d], year_effect_h, A_year_h,
                               draws$beta_sin1[d], draws$beta_cos1[d], draws$beta_sin2[d], draws$beta_cos2[d],
                               draws$beta_T[d], draws$beta_P[d], z_T_holdout, z_P_holdout)

    for (r in seq_len(n_replicates)) {
      slot <- slot + 1L
      set.seed(seed_offset + d * 1000L + r)
      rng_pois <- switch(lambda_mode,
        constant = make_rng_pois_step7(lambda_fn = function(wk) lambda_const),
        episodic = make_rng_pois_episodic(p = burst_p, lambda_mean = burst_lambda_mean))
      sim <- simulate_renewal_core(
        R0_vec = r0_h$R0, w = W_GEN, S0 = adv$S, U0 = adv$U,
        N_start_vec = holdout_weekly$N_start, N_end_vec = holdout_weekly$N_end, births_vec = holdout_weekly$births, deaths_vec = holdout_weekly$all_cause_deaths,
        X_recent_history = X_recent, lambda_import = 1, stochastic = TRUE,
        rng_pois = rng_pois, rng_binom = function(n, size, p) rbinom(n, size, p)
      )
      Y_rep <- simulate_reported_cases(sim$X, q = QVAL, phi_obs = draws$phi_obs[d], stochastic = TRUE, rng_nb = function(n, mu, phi) rnbinom(n, mu = mu, size = phi))
      traj_list[[slot]] <- tibble::tibble(posterior_draw_id = d, stochastic_replicate_id = r, week_start = holdout_weekly$week_start,
                                   S = sim$S, U = sim$U, N = sim$N, X_latent = sim$X, Y_reported = Y_rep,
                                   R0 = sim$R0, R_eff = sim$R_eff, lambda_t = sim$lambda, M_import = sim$M)
    }
  }
  list(trajectories = dplyr::bind_rows(traj_list), year_effects = dplyr::bind_rows(year_effect_rows),
       observed = holdout_weekly |> dplyr::select(week_start, observed_cases = cases, population = N_start),
       n_draws_available = n_draws_available, N_train = N_train)
}
