# PHASE 2 -- Controlled, single-mechanism sensitivity experiments. NO
# NEW STAN FITTING: reuses the training-only tau_year posteriors.
# Experiment A: importation grid (prespecified 0.25/1/4, NOT tuned to
#   holdout outcomes -- values fixed in 00_config_year_effect_validation.R
#   BEFORE any holdout comparison was computed).
# Experiment B: G1 (fixed SD=0.40) vs G2 (estimated tau_year), paired on
#   common z_year/z_season/climate/import draws per (posterior_draw,
#   replicate) -- built into simulate_holdout()'s own seeding, not a new design.
# Experiment C: impact of posterior uncertainty in the training-cutoff
#   susceptible population, isolated by fixing S/U at the cutoff to the
#   POSTERIOR MEAN terminal state while leaving every other
#   draw-specific parameter (tau_year, alpha_R, beta_T/P, phi_obs) free.
# Experiment D: joint timing of realised imports vs R_eff>1 vs
#   susceptible depletion in 2022 -- pure analysis of Phase 1's
#   diagnostic trajectories, no new simulation.

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

#' Same mechanics as Phase 1's simulate_holdout_diag(), plus an optional
#' `fixed_init = list(S=, U=)` that REPLACES every draw's own
#' advance_terminal_state() output with a single common value -- used
#' ONLY by Experiment C to isolate S/U-uncertainty's contribution.
simulate_holdout_ctrl <- function(training_bundle, holdout_spec, generator_arm, lambda_import, climate_arm = "observed",
                                   n_draws, n_replicates, seed_offset, fixed_init = NULL) {
  sd_train <- training_bundle$stan_data
  N_train <- sd_train$N
  draws <- rstan::extract(training_bundle$fit, pars = c("alpha_R", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                                                          "sigma_season_year", "tau_year", "beta_T", "beta_P", "phi_obs", "S", "U", "X"), permuted = TRUE)
  n_draws_available <- length(draws$alpha_R)

  holdout_weekly <- load_observed_weekly(holdout_spec$holdout_start, holdout_spec$holdout_end)
  n_holdout <- nrow(holdout_weekly)
  t_index_holdout <- (N_train + 1):(N_train + n_holdout)
  holdout_years <- as.integer(format(holdout_weekly$week_start, "%Y"))
  holdout_year_labels <- sort(unique(holdout_years)); n_holdout_years <- length(holdout_year_labels)
  year_id_holdout <- match(holdout_years, holdout_year_labels)

  obs_clim <- load_observed_climate(holdout_spec$holdout_start, holdout_spec$holdout_end)
  stopifnot(nrow(obs_clim) == n_holdout, all(obs_clim$week_start == holdout_weekly$week_start))
  z_T_holdout <- obs_clim$z_T_anom; z_P_holdout <- obs_clim$z_P_anom

  reconciliation_last_train <- sd_train$N_end[N_train] - sd_train$N_start[N_train] - sd_train$births[N_train] + sd_train$deaths[N_train]

  set.seed(seed_offset)
  draw_ids <- sample.int(n_draws_available, n_draws)
  traj_list <- vector("list", n_draws * n_replicates); slot <- 0L
  for (di in seq_along(draw_ids)) {
    d <- draw_ids[di]
    if (is.null(fixed_init)) {
      adv <- advance_terminal_state(draws$S[d, N_train], draws$U[d, N_train], draws$X[d, N_train], sd_train$births[N_train], sd_train$deaths[N_train], reconciliation_last_train)
    } else {
      adv <- fixed_init # SAME S/U for every draw -- isolates parameter (not initial-state) uncertainty
    }
    X_recent <- draws$X[d, (N_train - G + 1):N_train]

    set.seed(seed_offset + d)
    z_year_h <- rnorm(n_holdout_years); z_season_h <- rnorm(n_holdout_years)
    tau_for_arm <- if (generator_arm == "G1_fixed_SD_040") 0.40 else draws$tau_year[d]
    year_effect_h <- tau_for_arm * (z_year_h - mean(z_year_h))
    A_year_h <- exp(draws$sigma_season_year[d] * z_season_h)

    r0_h <- compute_R0_future(t_index_holdout, year_id_holdout, draws$alpha_R[d], year_effect_h, A_year_h,
                               draws$beta_sin1[d], draws$beta_cos1[d], draws$beta_sin2[d], draws$beta_cos2[d],
                               draws$beta_T[d], draws$beta_P[d], z_T_holdout, z_P_holdout)

    for (r in seq_len(n_replicates)) {
      slot <- slot + 1L
      set.seed(seed_offset + d * 1000L + r)
      sim <- simulate_renewal_core(R0_vec = r0_h$R0, w = w, S0 = adv$S, U0 = adv$U,
        N_start_vec = holdout_weekly$N_start, N_end_vec = holdout_weekly$N_end, births_vec = holdout_weekly$births, deaths_vec = holdout_weekly$all_cause_deaths,
        X_recent_history = X_recent, lambda_import = lambda_import, stochastic = TRUE,
        rng_pois = function(n, lam) rpois(n, lam), rng_binom = function(n, size, p) rbinom(n, size, p))
      Y_rep <- simulate_reported_cases(sim$X, q = QVAL, phi_obs = draws$phi_obs[d], stochastic = TRUE, rng_nb = function(n, mu, phi) rnbinom(n, mu = mu, size = phi))
      traj_list[[slot]] <- tibble(posterior_draw_id = d, stochastic_replicate_id = r, week_start = holdout_weekly$week_start,
                                   S = sim$S, N = sim$N, X_latent = sim$X, Y_reported = Y_rep, R_eff = sim$R_eff)
    }
  }
  list(trajectories = bind_rows(traj_list), observed = holdout_weekly |> select(week_start, observed_cases = cases, population = N_start))
}

summarise_cell <- function(res, label) {
  cum <- res$trajectories |> group_by(posterior_draw_id, stochastic_replicate_id) |> summarise(cum_Y = sum(Y_reported), cum_X = sum(X_latent), .groups = "drop")
  tibble(cell = label, median_cum_Y = median(cum$cum_Y), sd_cum_Y = sd(cum$cum_Y), p975_cum_Y = quantile(cum$cum_Y, .975),
         median_cum_X = median(cum$cum_X), n = nrow(cum))
}

message("[Phase2] Loading training-only fits (reused, NOT refit)...")
training_bundles <- lapply(names(HOLDOUTS), function(h) readRDS(file.path(DIR_FITS_YV, sprintf("ce_tau_year_training_%s.rds", h))))
names(training_bundles) <- names(HOLDOUTS)

# ================================================================
# Experiment A -- importation grid, G1 AND G2 (G2's grid already run in
# 12_future_year_effect_validation at smoke scale; G1's grid is NEW here,
# smoke scale, for completeness of the G1-vs-G2 x lambda comparison).
# ================================================================
message("\n[Experiment A] Importation sensitivity grid for G1 (new, smoke scale; G2 grid reused from Phase 12)...")
expA_new <- list()
for (h in names(HOLDOUTS)) for (lam in setdiff(IMPORTATION_SENSITIVITY_GRID, IMPORTATION_REFERENCE)) {
  message(sprintf("  [%s / G1 / lambda=%.2f] smoke-scale run...", h, lam))
  res <- simulate_holdout_ctrl(training_bundles[[h]], HOLDOUTS[[h]], "G1_fixed_SD_040", lam, "observed",
                                ENSEMBLE_RETRO$smoke$n_draws, ENSEMBLE_RETRO$smoke$n_replicates, DIAG_SEED_OFFSET)
  expA_new[[paste(h, "G1_fixed_SD_040", lam, sep = "_")]] <- res
}
# G1 reference (lambda=1) and full G2 grid -- reuse already-saved results, no recomputation.
orig_retro <- readRDS(file.path(DIR_RESULTS_YV_PRIOR, "PHASE3_4_holdout_simulation_results.rds"))
expA_rows <- bind_rows(
  bind_rows(lapply(names(expA_new), function(nm) summarise_cell(expA_new[[nm]], nm))),
  bind_rows(lapply(names(orig_retro), function(nm) if (grepl("^(H1|H2)_(G1_fixed_SD_040|G2_estimated_tau)_observed_", nm)) summarise_cell(orig_retro[[nm]], nm) else NULL))
) |> separate(cell, into = c("holdout", "generator", "climate_or_lambda", "lambda_tail"), sep = "_(?=observed|training|0\\.|[0-9])", extra = "merge", fill = "right") |>
  mutate(lambda = suppressWarnings(as.numeric(gsub("^observed_|^training_resampled_", "", paste(climate_or_lambda, lambda_tail))))) |>
  select(-climate_or_lambda, -lambda_tail)
write_csv(expA_rows, file.path(DIR_TABLES_RD, "PHASE2_experimentA_importation_grid.csv"))
message("\n=== Experiment A: importation grid, G1 vs G2 ==="); print(as.data.frame(expA_rows), digits = 4)

# ================================================================
# Experiment B -- G1 vs G2, paired on common z_year/z_season/climate/
# import draws (same seed per posterior_draw_id -> same random shocks,
# only the year_effect SCALE differs). Uses Phase 1's diagnostic
# results (already run with these exact seeds).
# ================================================================
message("\n[Experiment B] G1 vs G2 paired decomposition (from Phase 1 diagnostic trajectories)...")
diag1 <- readRDS(file.path(DIR_RESULTS_RD, "PHASE1_diagnostic_simulation_results.rds"))
expB_rows <- bind_rows(lapply(c("H1", "H2"), function(h) {
  g1 <- diag1[[paste0(h, "_G1_fixed_SD_040")]]$trajectories |> group_by(posterior_draw_id, stochastic_replicate_id) |> summarise(cum_Y_G1 = sum(Y_reported), cum_X_G1 = sum(X_latent), .groups = "drop")
  g2 <- diag1[[paste0(h, "_G2_estimated_tau")]]$trajectories |> group_by(posterior_draw_id, stochastic_replicate_id) |> summarise(cum_Y_G2 = sum(Y_reported), cum_X_G2 = sum(X_latent), .groups = "drop")
  paired <- inner_join(g1, g2, by = c("posterior_draw_id", "stochastic_replicate_id")) |> mutate(holdout = h, delta_Y = cum_Y_G2 - cum_Y_G1, delta_X = cum_X_G2 - cum_X_G1)
  tibble(holdout = h, median_delta_Y = median(paired$delta_Y), lo95_delta_Y = quantile(paired$delta_Y, .025), hi95_delta_Y = quantile(paired$delta_Y, .975),
         median_cum_Y_G1 = median(paired$cum_Y_G1), median_cum_Y_G2 = median(paired$cum_Y_G2), n_paired = nrow(paired))
}))
write_csv(expB_rows, file.path(DIR_TABLES_RD, "PHASE2_experimentB_G1_vs_G2_paired.csv"))
message("\n=== Experiment B: paired G2-G1 contrast (common random numbers) ==="); print(as.data.frame(expB_rows), digits = 4)

# ================================================================
# Experiment C -- impact of posterior uncertainty in the training-
# cutoff susceptible population. Compares "full posterior S/U" vs
# "fixed-at-posterior-mean S/U" (both re-run at SAME smoke scale/seed
# for an exact, self-contained paired comparison), for the most
# consequential case: H1 / G2 (the arm failing to reproduce 2022), and
# H2 / G2 for contrast.
# ================================================================
message("\n[Experiment C] Initial-susceptibility uncertainty decomposition (H1 and H2, G2)...")
expC_rows <- list()
for (h in names(HOLDOUTS)) {
  tb <- training_bundles[[h]]; sd_train <- tb$stan_data; N_train <- sd_train$N
  draws_full <- rstan::extract(tb$fit, pars = c("S", "U", "X"), permuted = TRUE)
  reconciliation_last_train <- sd_train$N_end[N_train] - sd_train$N_start[N_train] - sd_train$births[N_train] + sd_train$deaths[N_train]
  adv_all <- lapply(seq_len(nrow(draws_full$S)), function(d) advance_terminal_state(draws_full$S[d, N_train], draws_full$U[d, N_train], draws_full$X[d, N_train], sd_train$births[N_train], sd_train$deaths[N_train], reconciliation_last_train))
  mean_S_init <- mean(vapply(adv_all, `[[`, numeric(1), "S")); mean_U_init <- mean(vapply(adv_all, `[[`, numeric(1), "U"))
  message(sprintf("  [%s] posterior mean S_init=%.0f, U_init=%.0f (S_frac=%.4f); across-draw SD of S_frac_init=%.5f",
                   h, mean_S_init, mean_U_init, mean_S_init / (mean_S_init + mean_U_init),
                   sd(vapply(adv_all, function(a) a$S / (a$S + a$U), numeric(1)))))

  res_full  <- simulate_holdout_ctrl(tb, HOLDOUTS[[h]], "G2_estimated_tau", IMPORTATION_REFERENCE, "observed", EXPERIMENT_C_ENSEMBLE$n_draws, EXPERIMENT_C_ENSEMBLE$n_replicates, EXPERIMENT_C_SEED_OFFSET, fixed_init = NULL)
  res_fixed <- simulate_holdout_ctrl(tb, HOLDOUTS[[h]], "G2_estimated_tau", IMPORTATION_REFERENCE, "observed", EXPERIMENT_C_ENSEMBLE$n_draws, EXPERIMENT_C_ENSEMBLE$n_replicates, EXPERIMENT_C_SEED_OFFSET, fixed_init = list(S = mean_S_init, U = mean_U_init))

  cum_full <- res_full$trajectories |> group_by(posterior_draw_id, stochastic_replicate_id) |> summarise(cum_Y = sum(Y_reported), .groups = "drop")
  cum_fixed <- res_fixed$trajectories |> group_by(posterior_draw_id, stochastic_replicate_id) |> summarise(cum_Y = sum(Y_reported), .groups = "drop")
  expC_rows[[h]] <- tibble(holdout = h, arm = c("full_posterior_S_init", "fixed_mean_S_init"),
                            median_cum_Y = c(median(cum_full$cum_Y), median(cum_fixed$cum_Y)),
                            sd_cum_Y = c(sd(cum_full$cum_Y), sd(cum_fixed$cum_Y)),
                            var_ratio_fixed_over_full = sd(cum_fixed$cum_Y)^2 / sd(cum_full$cum_Y)^2)
}
expC_rows <- bind_rows(expC_rows)
write_csv(expC_rows, file.path(DIR_TABLES_RD, "PHASE2_experimentC_init_susceptibility.csv"))
message("\n=== Experiment C: S_init uncertainty contribution to outcome variance ==="); print(as.data.frame(expC_rows), digits = 4)

# ================================================================
# Experiment D -- joint timing: realised imports vs R_eff>1 vs
# susceptible depletion, in 2022, from Phase 1's diagnostic trajectories
# (pure analysis, no new simulation).
# ================================================================
message("\n[Experiment D] Joint timing decomposition for 2022 (from Phase 1 diagnostics)...")
expD_rows <- bind_rows(lapply(names(diag1), function(nm) {
  traj <- diag1[[nm]]$trajectories |> mutate(year = as.integer(format(week_start, "%Y"))) |> dplyr::filter(year == 2022)
  tibble(cell = nm, n_week_draw_rep = nrow(traj),
         frac_weeks_Reff_gt1 = mean(traj$R_eff > 1), frac_weeks_any_import = mean(traj$M_import >= 1),
         frac_weeks_import_during_Reff_gt1 = mean(traj$M_import >= 1 & traj$R_eff > 1),
         frac_weeks_import_during_Reff_le1 = mean(traj$M_import >= 1 & traj$R_eff <= 1),
         median_S_frac_start_2022 = median(traj$S[traj$week_start == min(traj$week_start)] / traj$N[traj$week_start == min(traj$week_start)]),
         median_total_imports_2022 = median(tapply(traj$M_import, paste(traj$posterior_draw_id, traj$stochastic_replicate_id), sum)))
}))
write_csv(expD_rows, file.path(DIR_TABLES_RD, "PHASE2_experimentD_timing_2022.csv"))
message("\n=== Experiment D: 2022 timing decomposition ==="); print(as.data.frame(expD_rows), digits = 4)

message("\n[Phase2] Controlled sensitivity experiments complete -> ", DIR_TABLES_RD)
