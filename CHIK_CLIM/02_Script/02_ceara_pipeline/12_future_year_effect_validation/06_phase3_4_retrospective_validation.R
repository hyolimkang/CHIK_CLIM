# PHASE 3/4 -- Retrospective (held-out) epidemic validation.
#
# THE MOST IMPORTANT VALIDATION IN THIS STUDY. For H1 and H2, simulates
# FORWARD from the TRAINING-ONLY posterior's own terminal state into the
# holdout period, WITHOUT ever reading holdout-period observed cases and
# WITHOUT the 2022 conditioned seed anywhere. Observed holdout cases are
# used ONLY afterwards, to evaluate (never to fit/calibrate/tune) the
# simulated trajectories.
#
# Primary arm: observed historical climate during the holdout (available
# because this is retrospective -- the holdout period has already
# happened), cases withheld. Explicitly labelled a RETROSPECTIVE,
# OBSERVED-CLIMATE-CONDITIONAL forecast, NOT a genuine prospective one.
# Secondary arm (smoke scale, "if feasible"): climate resampled from
# TRAINING-period blocks only (no information from the holdout period's
# actual weather at all).
#
# Importation: prespecified sensitivity grid {0.25, 1, 4}/week, NOT
# tuned to holdout outcomes (grid fixed in 00_config_year_effect_validation.R
# before any holdout comparison was computed).

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/12_future_year_effect_validation/00_config_year_effect_validation.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R")) # build_state_weekly(), generation_weights()

G <- 8L
w <- generation_weights_future()

#' Observed weekly CE case/demography series for ANY date range (used
#' both to build holdout-period structural inputs -- demography/climate,
#' never cases for conditioning -- and, separately, to evaluate against).
load_observed_weekly <- function(date_start, date_end) {
  panel <- readRDS(file.path(ROOT, "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |> dplyr::filter(week_start >= date_start, week_start <= date_end)
  weekly_demography <- readRDS(file.path(ROOT, "01_Data/ceara_weekly_demography_2015_2025.rds")) |>
    dplyr::filter(week_start >= date_start, week_start <= date_end)
  weekly <- weekly_cases |> dplyr::select(-population) |> dplyr::inner_join(weekly_demography, by = "week_start") |> dplyr::arrange(week_start)
  stopifnot(nrow(weekly) == nrow(weekly_cases), nrow(weekly) == nrow(weekly_demography))
  weekly
}

#' Observed climate anomalies (primary lag, CE) for ANY date range --
#' used as the PRIMARY retrospective arm's "observed-climate-conditional" input.
load_observed_climate <- function(date_start, date_end) {
  anom <- read_csv(ACCEPTED_ANOMALY_CSV, show_col_types = FALSE) |>
    dplyr::filter(state == "CE", lag_spec == "primary") |> dplyr::mutate(week_start = as.Date(week_start)) |>
    dplyr::filter(week_start >= date_start, week_start <= date_end) |> dplyr::arrange(week_start)
  anom$z_T_anom[is.na(anom$z_T_anom)] <- 0; anom$z_P_anom[is.na(anom$z_P_anom)] <- 0
  anom
}

GEN_ARMS <- c("G1_fixed_SD_040", "G2_estimated_tau")

#' Core: simulate the holdout period forward from a TRAINING-ONLY fit,
#' for one (holdout, generator-arm, lambda_import, climate-arm) cell.
simulate_holdout <- function(training_bundle, holdout_spec, generator_arm, lambda_import, climate_arm = "observed",
                              n_draws, n_replicates, seed_offset) {
  sd_train <- training_bundle$stan_data
  N_train <- sd_train$N
  draws <- rstan::extract(training_bundle$fit, pars = c("alpha_R", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                                                          "sigma_season_year", "tau_year", "beta_T", "beta_P", "phi_obs", "S", "U", "X"), permuted = TRUE)
  n_draws_available <- length(draws$alpha_R)

  holdout_weekly <- load_observed_weekly(holdout_spec$holdout_start, holdout_spec$holdout_end) # demography (structural input) + cases (EVALUATION ONLY, never read below for conditioning)
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
  } else { # "training_resampled": whole-year blocks resampled from TRAINING period only, no holdout-period weather used at all
    train_clim <- load_observed_climate(holdout_spec$train_start, holdout_spec$train_end) |> mutate(year = as.integer(format(week_start, "%Y")), woy = as.integer(format(week_start, "%U")))
    train_by_year <- split(train_clim, train_clim$year)
    z_T_holdout <- numeric(n_holdout); z_P_holdout <- numeric(n_holdout)
    set.seed(seed_offset + 9999L)
    for (yi in seq_along(holdout_year_labels)) {
      src_year <- sample(names(train_by_year), 1L)
      block <- train_by_year[[src_year]]
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
  for (di in seq_along(draw_ids)) {
    d <- draw_ids[di]
    adv <- advance_terminal_state(draws$S[d, N_train], draws$U[d, N_train], draws$X[d, N_train], sd_train$births[N_train], sd_train$deaths[N_train], reconciliation_last_train) # TRAINING-ONLY posterior terminal state -- never the full-period fit
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
      sim <- simulate_renewal_core(
        R0_vec = r0_h$R0, w = w, S0 = adv$S, U0 = adv$U,
        N_start_vec = holdout_weekly$N_start, N_end_vec = holdout_weekly$N_end, births_vec = holdout_weekly$births, deaths_vec = holdout_weekly$all_cause_deaths,
        X_recent_history = X_recent, lambda_import = lambda_import, stochastic = TRUE, # NO is_seed/X_seed passed -- the 2022 conditioned seed is NEVER used here, by construction
        rng_pois = function(n, lam) rpois(n, lam), rng_binom = function(n, size, p) rbinom(n, size, p)
      )
      Y_rep <- simulate_reported_cases(sim$X, q = QVAL, phi_obs = draws$phi_obs[d], stochastic = TRUE, rng_nb = function(n, mu, phi) rnbinom(n, mu = mu, size = phi))
      traj_list[[slot]] <- tibble(posterior_draw_id = d, stochastic_replicate_id = r, week_start = holdout_weekly$week_start,
                                   S = sim$S, N = sim$N, X_latent = sim$X, Y_reported = Y_rep, R0 = sim$R0, R_eff = sim$R_eff)
    }
  }
  list(trajectories = bind_rows(traj_list), observed = holdout_weekly |> select(week_start, observed_cases = cases, population = N_start), n_draws_available = n_draws_available)
}

# ================================================================
# Observed holdout "major epidemic" ground truth, SAME detection rule as simulations
# ================================================================
observed_major_epidemic_summary <- function(holdout_spec) {
  obs <- load_observed_weekly(holdout_spec$holdout_start, holdout_spec$holdout_end)
  waves <- detect_simulated_waves(obs$week_start, obs$cases, obs$N_start)
  list(observed = obs, waves = waves, has_major = if (nrow(waves)) any(waves$major_epidemic) else FALSE,
       first_major_onset = if (nrow(waves) && any(waves$major_epidemic)) min(waves$onset_week[waves$major_epidemic]) else as.Date(NA),
       largest_epidemic_cases = if (nrow(waves)) max(waves$total_cases) else 0,
       annual_cases = obs |> mutate(year = as.integer(format(week_start, "%Y"))) |> group_by(year) |> summarise(cases = sum(cases), .groups = "drop"))
}

# ================================================================
# Run all cells
# ================================================================
message("[Phase3/4] Loading training-only fits...")
training_bundles <- lapply(names(HOLDOUTS), function(h) readRDS(file.path(DIR_FITS_YV, sprintf("ce_tau_year_training_%s.rds", h))))
names(training_bundles) <- names(HOLDOUTS)

all_results <- list()
all_obs <- list()
for (h in names(HOLDOUTS)) {
  hs <- HOLDOUTS[[h]]
  message(sprintf("\n=== %s ===", hs$label))
  all_obs[[h]] <- observed_major_epidemic_summary(hs)
  message(sprintf("  Observed: %d major epidemic(s) in holdout; largest=%.0f cases; first onset=%s",
                   sum(all_obs[[h]]$waves$major_epidemic), all_obs[[h]]$largest_epidemic_cases, as.character(all_obs[[h]]$first_major_onset)))

  for (g in GEN_ARMS) {
    message(sprintf("  [%s / %s] full-scale run (observed climate, lambda_import=%s)...", h, g, IMPORTATION_REFERENCE))
    res <- simulate_holdout(training_bundles[[h]], hs, g, IMPORTATION_REFERENCE, "observed",
                             ENSEMBLE_RETRO$full$n_draws, ENSEMBLE_RETRO$full$n_replicates, BASE_SEED_YV + 777L)
    all_results[[paste(h, g, "observed", IMPORTATION_REFERENCE, sep = "_")]] <- res
  }

  # Importation sensitivity grid, smoke scale, NOT tuned to the result above (grid fixed in config before this loop ran)
  for (lam in setdiff(IMPORTATION_SENSITIVITY_GRID, IMPORTATION_REFERENCE)) {
    message(sprintf("  [%s / G2 / importation sensitivity lambda=%.2f] smoke-scale run...", h, lam))
    res <- simulate_holdout(training_bundles[[h]], hs, "G2_estimated_tau", lam, "observed",
                             ENSEMBLE_RETRO$smoke$n_draws, ENSEMBLE_RETRO$smoke$n_replicates, BASE_SEED_YV + 777L)
    all_results[[paste(h, "G2_estimated_tau", "observed", lam, sep = "_")]] <- res
  }

  # Training-resampled climate, secondary arm, smoke scale ("if feasible")
  message(sprintf("  [%s / G2 / training-resampled climate] smoke-scale run...", h))
  res <- simulate_holdout(training_bundles[[h]], hs, "G2_estimated_tau", IMPORTATION_REFERENCE, "training_resampled",
                           ENSEMBLE_RETRO$smoke$n_draws, ENSEMBLE_RETRO$smoke$n_replicates, BASE_SEED_YV + 777L)
  all_results[[paste(h, "G2_estimated_tau", "training_resampled", IMPORTATION_REFERENCE, sep = "_")]] <- res
}
saveRDS(all_results, file.path(DIR_RESULTS_YV, "PHASE3_4_holdout_simulation_results.rds"))
saveRDS(all_obs, file.path(DIR_RESULTS_YV, "PHASE3_4_observed_holdout_summaries.rds"))

# ================================================================
# Validation metrics: predictive medians, 50/80/95% PI, coverage, against observed
# ================================================================
compute_validation_metrics <- function(cell_name, res, obs_summary) {
  traj <- res$trajectories |> mutate(year = as.integer(format(week_start, "%Y")))
  annual_sim <- traj |> group_by(posterior_draw_id, stochastic_replicate_id, year) |> summarise(annual_cases = sum(Y_reported), .groups = "drop")
  annual_pi <- annual_sim |> group_by(year) |>
    summarise(median = median(annual_cases), lo25 = quantile(annual_cases, .25), hi75 = quantile(annual_cases, .75),
              lo10 = quantile(annual_cases, .10), hi90 = quantile(annual_cases, .90), lo025 = quantile(annual_cases, .025), hi975 = quantile(annual_cases, .975), .groups = "drop") |>
    left_join(obs_summary$annual_cases, by = "year") |>
    mutate(cell = cell_name, covered_50 = cases >= lo25 & cases <= hi75, covered_80 = cases >= lo10 & cases <= hi90, covered_95 = cases >= lo025 & cases <= hi975)

  cum_sim <- traj |> group_by(posterior_draw_id, stochastic_replicate_id) |> summarise(cum_cases = sum(Y_reported), .groups = "drop")
  observed_cum <- sum(obs_summary$annual_cases$cases)

  traj_split <- split(traj, list(traj$posterior_draw_id, traj$stochastic_replicate_id), sep = "||", drop = TRUE)
  major <- bind_rows(lapply(traj_split, function(g) {
    waves <- detect_simulated_waves(g$week_start, g$Y_reported, g$N)
    has_major_epi <- if (nrow(waves)) any(waves$major_epidemic) else FALSE
    tibble(posterior_draw_id = g$posterior_draw_id[1], stochastic_replicate_id = g$stochastic_replicate_id[1],
           has_major = has_major_epi, largest = if (nrow(waves)) max(waves$total_cases) else 0,
           first_onset = if (has_major_epi) min(waves$onset_week[waves$major_epidemic]) else as.Date(NA))
  }))

  fadeout <- traj |> group_by(posterior_draw_id, stochastic_replicate_id) |> summarise(weeks_zero_cases = sum(Y_reported < 1), .groups = "drop")

  list(
    annual_pi = annual_pi,
    cell_summary = tibble(cell = cell_name,
                           prob_major_epidemic_simulated = mean(major$has_major), observed_has_major = obs_summary$has_major,
                           median_cum_cases = median(cum_sim$cum_cases), lo95_cum = quantile(cum_sim$cum_cases, .025), hi95_cum = quantile(cum_sim$cum_cases, .975), observed_cum_cases = observed_cum,
                           cum_cases_covered_95 = observed_cum >= quantile(cum_sim$cum_cases, .025) && observed_cum <= quantile(cum_sim$cum_cases, .975),
                           median_largest_epidemic = median(major$largest), observed_largest_epidemic = obs_summary$largest_epidemic_cases,
                           median_weeks_zero_cases = median(fadeout$weeks_zero_cases), frac_weeks_zero_cases = median(fadeout$weeks_zero_cases) / n_distinct(traj$week_start),
                           annual_coverage_50 = mean(annual_pi$covered_50, na.rm = TRUE), annual_coverage_80 = mean(annual_pi$covered_80, na.rm = TRUE), annual_coverage_95 = mean(annual_pi$covered_95, na.rm = TRUE)),
    timing = major |> dplyr::filter(!is.na(first_onset)) |> mutate(cell = cell_name)
  )
}

metrics_list <- lapply(names(all_results), function(nm) {
  h <- strsplit(nm, "_")[[1]][1]
  compute_validation_metrics(nm, all_results[[nm]], all_obs[[h]])
})
names(metrics_list) <- names(all_results)

annual_pi_all <- bind_rows(lapply(metrics_list, `[[`, "annual_pi"))
cell_summary_all <- bind_rows(lapply(metrics_list, `[[`, "cell_summary"))
timing_all <- bind_rows(lapply(metrics_list, `[[`, "timing"))

write_csv(annual_pi_all, file.path(DIR_TABLES_YV, "PHASE4_annual_PI_vs_observed.csv"))
write_csv(cell_summary_all, file.path(DIR_TABLES_YV, "PHASE4_cell_validation_summary.csv"))
write_csv(timing_all, file.path(DIR_RESULTS_YV, "PHASE4_first_major_epidemic_timing.csv"))

message("\n=== PHASE 4 validation summary (all cells) ===")
print(as.data.frame(cell_summary_all), digits = 4)

message("\n[Phase3/4] Retrospective validation complete. Tables -> ", DIR_TABLES_YV)
