# PHASE 1 -- WHY does H1 under-predict the 2022 resurgence while H2
# performs better? NO NEW STAN FITTING: reuses the training-only
# tau_year posteriors (ce_tau_year_training_{H1,H2}.rds) already fitted
# in 12_future_year_effect_validation. Re-simulates the holdout period
# with an EXTENDED version of simulate_holdout() that additionally
# captures realised weekly imports (M), susceptible/infectious state
# (U), and per-week/per-year transmission inputs (year_effect, A_year,
# lambda hazard) -- none of which the original saved
# PHASE3_4_holdout_simulation_results.rds trajectories retained.
#
# Uses IDENTICAL seeds (seed_offset = BASE_SEED_YV + 777L, same draw/
# replicate loop order) to the original reference cells, so S/X/R0/R_eff
# reproduce BIT-FOR-BIT -- verified explicitly below -- and only the
# NEW diagnostic columns are genuinely new computation.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/13_retrospective_resurgence_diagnostic/00_config_retro_diagnostic.R"))

G <- 8L
w <- generation_weights_future()

load_observed_weekly <- function(date_start, date_end) {
  panel <- readRDS(file.path(ROOT, "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |> dplyr::filter(week_start >= date_start, week_start <= date_end)
  weekly_demography <- readRDS(file.path(ROOT, "01_Data/ceara_weekly_demography_2015_2025.rds")) |>
    dplyr::filter(week_start >= date_start, week_start <= date_end)
  weekly <- weekly_cases |> dplyr::select(-population) |> dplyr::inner_join(weekly_demography, by = "week_start") |> dplyr::arrange(week_start)
  stopifnot(nrow(weekly) == nrow(weekly_cases), nrow(weekly) == nrow(weekly_demography))
  weekly
}
load_observed_climate <- function(date_start, date_end) {
  anom <- read_csv(ACCEPTED_ANOMALY_CSV, show_col_types = FALSE) |>
    dplyr::filter(state == "CE", lag_spec == "primary") |> dplyr::mutate(week_start = as.Date(week_start)) |>
    dplyr::filter(week_start >= date_start, week_start <= date_end) |> dplyr::arrange(week_start)
  anom$z_T_anom[is.na(anom$z_T_anom)] <- 0; anom$z_P_anom[is.na(anom$z_P_anom)] <- 0
  anom
}

#' Extended simulate_holdout() -- SAME mechanics/seeds as
#' 12_future_year_effect_validation/06_phase3_4_retrospective_validation.R,
#' with added diagnostic capture (U, M realised imports, lambda hazard,
#' per-week year_effect/A_year) needed for Phase 1's mechanism decomposition.
simulate_holdout_diag <- function(training_bundle, holdout_spec, generator_arm, lambda_import, climate_arm = "observed",
                                   n_draws, n_replicates, seed_offset) {
  sd_train <- training_bundle$stan_data
  N_train <- sd_train$N
  draws <- rstan::extract(training_bundle$fit, pars = c("alpha_R", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                                                          "sigma_season_year", "tau_year", "beta_T", "beta_P", "phi_obs", "S", "U", "X"), permuted = TRUE)
  n_draws_available <- length(draws$alpha_R)

  holdout_weekly <- load_observed_weekly(holdout_spec$holdout_start, holdout_spec$holdout_end)
  n_holdout <- nrow(holdout_weekly)
  t_index_holdout <- (N_train + 1):(N_train + n_holdout)
  holdout_years <- as.integer(format(holdout_weekly$week_start, "%Y"))
  holdout_year_labels <- sort(unique(holdout_years))
  n_holdout_years <- length(holdout_year_labels)
  year_id_holdout <- match(holdout_years, holdout_year_labels)

  if (climate_arm == "observed") {
    obs_clim <- load_observed_climate(holdout_spec$holdout_start, holdout_spec$holdout_end)
    stopifnot(nrow(obs_clim) == n_holdout, all(obs_clim$week_start == holdout_weekly$week_start))
    z_T_holdout <- obs_clim$z_T_anom; z_P_holdout <- obs_clim$z_P_anom
  } else {
    train_clim <- load_observed_climate(holdout_spec$train_start, holdout_spec$train_end) |> mutate(year = as.integer(format(week_start, "%Y")))
    train_by_year <- split(train_clim, train_clim$year)
    z_T_holdout <- numeric(n_holdout); z_P_holdout <- numeric(n_holdout)
    set.seed(seed_offset + 9999L)
    for (yi in seq_along(holdout_year_labels)) {
      src_year <- sample(names(train_by_year), 1L); block <- train_by_year[[src_year]]
      idx_in_year <- which(holdout_years == holdout_year_labels[yi])
      need <- length(idx_in_year); have <- nrow(block)
      take <- if (have >= need) seq_len(need) else c(seq_len(have), rep(have, need - have))
      z_T_holdout[idx_in_year] <- block$z_T_anom[take]; z_P_holdout[idx_in_year] <- block$z_P_anom[take]
    }
  }

  reconciliation_last_train <- sd_train$N_end[N_train] - sd_train$N_start[N_train] - sd_train$births[N_train] + sd_train$deaths[N_train]

  set.seed(seed_offset)
  draw_ids <- sample.int(n_draws_available, n_draws)
  traj_list <- vector("list", n_draws * n_replicates); slot <- 0L
  year_effect_rows <- vector("list", length(draw_ids))
  init_state_rows <- vector("list", length(draw_ids))
  for (di in seq_along(draw_ids)) {
    d <- draw_ids[di]
    adv <- advance_terminal_state(draws$S[d, N_train], draws$U[d, N_train], draws$X[d, N_train], sd_train$births[N_train], sd_train$deaths[N_train], reconciliation_last_train)
    X_recent <- draws$X[d, (N_train - G + 1):N_train]
    init_state_rows[[di]] <- tibble(posterior_draw_id = d, S_init = adv$S, U_init = adv$U, N_init = adv$S + adv$U, S_frac_init = adv$S / (adv$S + adv$U))

    set.seed(seed_offset + d)
    z_year_h <- rnorm(n_holdout_years); z_season_h <- rnorm(n_holdout_years)
    tau_for_arm <- if (generator_arm == "G1_fixed_SD_040") 0.40 else draws$tau_year[d]
    year_effect_h <- tau_for_arm * (z_year_h - mean(z_year_h))
    A_year_h <- exp(draws$sigma_season_year[d] * z_season_h)
    year_effect_rows[[di]] <- tibble(posterior_draw_id = d, holdout_year = holdout_year_labels, year_effect = year_effect_h, A_year = A_year_h, tau_used = tau_for_arm)

    r0_h <- compute_R0_future(t_index_holdout, year_id_holdout, draws$alpha_R[d], year_effect_h, A_year_h,
                               draws$beta_sin1[d], draws$beta_cos1[d], draws$beta_sin2[d], draws$beta_cos2[d],
                               draws$beta_T[d], draws$beta_P[d], z_T_holdout, z_P_holdout)

    for (r in seq_len(n_replicates)) {
      slot <- slot + 1L
      set.seed(seed_offset + d * 1000L + r)
      sim <- simulate_renewal_core(
        R0_vec = r0_h$R0, w = w, S0 = adv$S, U0 = adv$U,
        N_start_vec = holdout_weekly$N_start, N_end_vec = holdout_weekly$N_end, births_vec = holdout_weekly$births, deaths_vec = holdout_weekly$all_cause_deaths,
        X_recent_history = X_recent, lambda_import = lambda_import, stochastic = TRUE,
        rng_pois = function(n, lam) rpois(n, lam), rng_binom = function(n, size, p) rbinom(n, size, p)
      )
      Y_rep <- simulate_reported_cases(sim$X, q = QVAL, phi_obs = draws$phi_obs[d], stochastic = TRUE, rng_nb = function(n, mu, phi) rnbinom(n, mu = mu, size = phi))
      traj_list[[slot]] <- tibble(posterior_draw_id = d, stochastic_replicate_id = r, week_start = holdout_weekly$week_start,
                                   S = sim$S, U = sim$U, N = sim$N, X_latent = sim$X, Y_reported = Y_rep,
                                   R0 = sim$R0, R_eff = sim$R_eff, lambda_t = sim$lambda, M_import = sim$M)
    }
  }
  list(trajectories = bind_rows(traj_list), year_effects = bind_rows(year_effect_rows), init_states = bind_rows(init_state_rows),
       observed = holdout_weekly |> select(week_start, observed_cases = cases, population = N_start), n_draws_available = n_draws_available)
}

message("[Phase1] Loading training-only fits (reused, NOT refit)...")
training_bundles <- lapply(names(HOLDOUTS), function(h) readRDS(file.path(DIR_FITS_YV, sprintf("ce_tau_year_training_%s.rds", h))))
names(training_bundles) <- names(HOLDOUTS)

CELLS <- expand.grid(holdout = c("H1", "H2"), generator = c("G1_fixed_SD_040", "G2_estimated_tau"), stringsAsFactors = FALSE)
diag_results <- list()
for (i in seq_len(nrow(CELLS))) {
  h <- CELLS$holdout[i]; g <- CELLS$generator[i]
  message(sprintf("  [%s / %s] diagnostic re-simulation (lambda=%s, observed climate)...", h, g, IMPORTATION_REFERENCE))
  diag_results[[paste(h, g, sep = "_")]] <- simulate_holdout_diag(training_bundles[[h]], HOLDOUTS[[h]], g, IMPORTATION_REFERENCE, "observed",
                                                                    DIAG_ENSEMBLE$n_draws, DIAG_ENSEMBLE$n_replicates, DIAG_SEED_OFFSET)
}
saveRDS(diag_results, file.path(DIR_RESULTS_RD, "PHASE1_diagnostic_simulation_results.rds"))

# ================================================================
# Integrity check: new diagnostic re-simulation must reproduce the
# ORIGINAL saved S/X_latent/Y_reported/R0/R_eff bit-for-bit (same
# seeds/order) -- proves this is a faithful re-derivation, not a
# different simulation, before trusting the new diagnostic columns.
# ================================================================
orig <- readRDS(file.path(DIR_RESULTS_YV_PRIOR, "PHASE3_4_holdout_simulation_results.rds"))
integrity <- bind_rows(lapply(seq_len(nrow(CELLS)), function(i) {
  h <- CELLS$holdout[i]; g <- CELLS$generator[i]
  nm <- paste(h, g, sep = "_"); orig_key <- paste(h, g, "observed", IMPORTATION_REFERENCE, sep = "_")
  if (!orig_key %in% names(orig)) return(tibble(cell = nm, orig_key = orig_key, identical_S = NA, identical_X = NA, identical_Y = NA, identical_R0 = NA, identical_Reff = NA, max_abs_diff_S = NA_real_))
  a <- diag_results[[nm]]$trajectories |> arrange(posterior_draw_id, stochastic_replicate_id, week_start)
  b <- orig[[orig_key]]$trajectories |> arrange(posterior_draw_id, stochastic_replicate_id, week_start)
  tibble(cell = nm, orig_key = orig_key, n_rows_match = nrow(a) == nrow(b),
         identical_S = isTRUE(all.equal(a$S, b$S)), identical_X = isTRUE(all.equal(a$X_latent, b$X_latent)),
         identical_Y = isTRUE(all.equal(a$Y_reported, b$Y_reported)), identical_R0 = isTRUE(all.equal(a$R0, b$R0)), identical_Reff = isTRUE(all.equal(a$R_eff, b$R_eff)),
         max_abs_diff_S = max(abs(a$S - b$S)))
}))
write_csv(integrity, file.path(DIR_TABLES_RD, "PHASE1_integrity_check_vs_original.csv"))
message("\n=== Integrity check: diagnostic re-simulation vs original saved PHASE3_4 trajectories ===")
print(as.data.frame(integrity), digits = 6)
if (!all(unlist(integrity[, c("identical_S","identical_X","identical_Y","identical_R0","identical_Reff")]), na.rm = TRUE)) {
  warning("Diagnostic re-simulation does NOT exactly reproduce the original saved trajectories -- investigate before trusting Phase 1 results.")
}

message("\n[Phase1] Diagnostic re-simulation complete -> ", DIR_RESULTS_RD)
