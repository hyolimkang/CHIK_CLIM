# STEP 7 shared simulation engine. Reuses simulate_renewal_core(),
# compute_R0_future(), advance_terminal_state(), generation_weights_future()
# from 10_future_projection_pilot/03_future_renewal_functions.R UNMODIFIED
# (sourced via 00_config_year_effect_validation.R). This file adds only
# NEW, NON-INVASIVE wrapper logic: it never edits the core simulator.
#
# VARIANCE-REDUCTION METHOD (documented, per the task's requirement):
# simulate_renewal_core() calls `rng_pois(1, lambda_import)` once per
# week for the introduction term. Rather than letting this draw a fresh
# Poisson(lambda) variate per arm (which would NOT keep arms coupled --
# different lambda values consume randomness differently and, worse,
# once M_t differs between arms in week t, the whole downstream X(t)
# trajectory and hence the demographic state diverges, decoupling all
# LATER weeks' random draws even if the same seed was used), every
# rng_pois closure below instead does:
#     u <- runif(1); qpois(u, lambda_for_this_week)
# i.e. a SHARED-UNIFORM (inverse-CDF) coupling: exactly one uniform is
# consumed per week regardless of lambda's value (including lambda=0,
# where qpois(u,0)=0 deterministically), so every arm's random-number
# STREAM POSITION remains aligned going into the following rbinom()
# infection draw for as long as the trajectories have not yet diverged,
# and the coupling is monotone (a larger u maps to a weakly larger M_t
# at every lambda) which is the standard antithetic/common-random-number
# construction for comparing Poisson arms. Marginally, each arm's M_t
# is still exactly Poisson(lambda_for_this_week) -- only the DEPENDENCE
# across arms is induced, never the marginal distribution.

G_GEN <- 8L
W_GEN <- generation_weights_future()

load_observed_weekly_step7 <- function(date_start, date_end) {
  panel <- readRDS(file.path(ROOT, "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |> dplyr::filter(week_start >= date_start, week_start <= date_end)
  weekly_demography <- readRDS(file.path(ROOT, "01_Data/ceara_weekly_demography_2015_2025.rds")) |>
    dplyr::filter(week_start >= date_start, week_start <= date_end)
  weekly <- weekly_cases |> dplyr::select(-population) |> dplyr::inner_join(weekly_demography, by = "week_start") |> dplyr::arrange(week_start)
  stopifnot(nrow(weekly) == nrow(weekly_cases), nrow(weekly) == nrow(weekly_demography))
  weekly
}
load_observed_climate_step7 <- function(date_start, date_end) {
  anom <- readr::read_csv(ACCEPTED_ANOMALY_CSV, show_col_types = FALSE) |>
    dplyr::filter(state == "CE", lag_spec == "primary") |> dplyr::mutate(week_start = as.Date(week_start)) |>
    dplyr::filter(week_start >= date_start, week_start <= date_end) |> dplyr::arrange(week_start)
  anom$z_T_anom[is.na(anom$z_T_anom)] <- 0; anom$z_P_anom[is.na(anom$z_P_anom)] <- 0
  anom
}

#' Build a shared-uniform-coupled rng_pois closure.
#' `lambda_fn(week_index)` returns the intended lambda for holdout week
#' `week_index` (1-based, within the holdout only). `fixed_M` (optional)
#' is a pre-realised M_t vector to REPLAY deterministically (Phase 3's
#' timing-permutation arm) -- a uniform is still drawn each week (for
#' RNG-stream alignment with the trajectory `fixed_M` was captured from)
#' but discarded in favour of the fixed value.
make_rng_pois_step7 <- function(lambda_fn = NULL, fixed_M = NULL) {
  week_counter <- 0L
  function(n, lam) {
    week_counter <<- week_counter + 1L
    u <- stats::runif(n)
    if (!is.null(fixed_M)) return(fixed_M[week_counter])
    stats::qpois(u, lambda_fn(week_counter))
  }
}

#' Core Step 7 holdout simulator. A generalisation of
#' 13_retrospective_resurgence_diagnostic's simulate_holdout_diag() /
#' 14_hmc_reliability_audit's simulate_holdout_ctrl() that (a) allows
#' lambda_import to vary by week via `lambda_fn`/`fixed_M`, instead of a
#' single scalar, and (b) allows the future year-effect centring
#' convention to be switched (Phase 5), WITHOUT modifying
#' simulate_renewal_core() itself. Reuses the SAME training-only
#' posterior and seeding conventions as Step 6/13/14.
#'
#' `lambda_mode`: "constant" (single value throughout), "split" (two
#' values either side of `split_date`), or "fixed_M" (replay a supplied
#' per-draw-per-replicate M_t matrix, Phase 3 only).
simulate_holdout_step7 <- function(training_bundle, holdout_spec, lambda_mode = c("constant", "split", "fixed_M"),
                                    lambda_const = NULL, lambda_pre = NULL, lambda_post = NULL, split_date = NULL,
                                    fixed_M_lookup = NULL, # function(posterior_draw_id, replicate_id) -> numeric vector, used iff lambda_mode=="fixed_M"
                                    year_effect_mode = c("centered", "uncentered"),
                                    n_draws, n_replicates, seed_offset, capture_M = TRUE) {
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

  if (lambda_mode == "split") {
    is_post <- holdout_weekly$week_start >= split_date
    lambda_week_vec <- ifelse(is_post, lambda_post, lambda_pre)
  }

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
    year_effect_rows[[di]] <- tibble::tibble(posterior_draw_id = d, holdout_year = holdout_year_labels, year_effect = year_effect_h, A_year = A_year_h,
                                              tau_used = draws$tau_year[d], year_effect_mode = year_effect_mode)

    r0_h <- compute_R0_future(t_index_holdout, year_id_holdout, draws$alpha_R[d], year_effect_h, A_year_h,
                               draws$beta_sin1[d], draws$beta_cos1[d], draws$beta_sin2[d], draws$beta_cos2[d],
                               draws$beta_T[d], draws$beta_P[d], z_T_holdout, z_P_holdout)

    for (r in seq_len(n_replicates)) {
      slot <- slot + 1L
      set.seed(seed_offset + d * 1000L + r)
      rng_pois <- switch(lambda_mode,
        constant = make_rng_pois_step7(lambda_fn = function(wk) lambda_const),
        split    = make_rng_pois_step7(lambda_fn = function(wk) lambda_week_vec[wk]),
        fixed_M  = make_rng_pois_step7(fixed_M = fixed_M_lookup(d, r)))
      sim <- simulate_renewal_core(
        R0_vec = r0_h$R0, w = W_GEN, S0 = adv$S, U0 = adv$U,
        N_start_vec = holdout_weekly$N_start, N_end_vec = holdout_weekly$N_end, births_vec = holdout_weekly$births, deaths_vec = holdout_weekly$all_cause_deaths,
        X_recent_history = X_recent, lambda_import = 1, stochastic = TRUE, # dummy positive placeholder -- the ACTUAL per-week lambda is applied inside rng_pois, ignoring this argument
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
