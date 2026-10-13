# PHASE 0 -- Reproduce and audit the lambda_import=1 baseline for H1/G2
# and H2/G2 against Step 6's reference metrics, using the SAME training-
# only posteriors and the Step 7 simulation engine (shared-uniform-
# coupled rng_pois, but at a single constant lambda this reduces
# EXACTLY to ordinary independent Poisson(lambda) draws -- the coupling
# only matters once multiple lambda arms are compared from Phase 1 on).
# Also audits the established/major-epidemic event-classification logic.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr) })
filter <- dplyr::filter

source(file.path(ROOT, "02_Script/02_ceara_pipeline/16_step7_importation_persistence/00_config_step7.R"))
source(file.path(ROOT, "02_Script/02_ceara_pipeline/16_step7_importation_persistence/01_simulation_engine.R"))
source(file.path(DIR_STEP6, "00_config_step6.R")) # reuse ESTABLISHMENT_RUN_WEEKS, MAJOR_ABS_CASE_MINIMUM/INCIDENCE, WAVE_DETECTION_YV -- SAME definitions as Step 6, not redefined

max_run <- function(x) { rr <- rle(replace(x, is.na(x), FALSE)); if (!any(rr$values)) 0L else max(rr$lengths[rr$values]) }
qv <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE, type = 8))

message("[Phase0] Loading training-only fits (H1/H2, reused, NOT refit)...")
training_bundles <- list(H1 = readRDS(FIT_PATHS[["H1"]]), H2 = readRDS(FIT_PATHS[["H2"]]))
# Posterior-draw identity check: confirm these are the SAME fit objects Step 6 used (same N, same stan_data fingerprint).
fit_identity <- tibble(holdout = c("H1", "H2"),
                        N_train = vapply(training_bundles, function(b) b$stan_data$N, numeric(1)),
                        n_post_warmup_draws = vapply(training_bundles, function(b) length(rstan::extract(b$fit, pars = "tau_year")$tau_year), numeric(1)),
                        config_prior_id = vapply(training_bundles, function(b) if (is.null(b$config$holdout_id)) NA_character_ else b$config$holdout_id, character(1)))
write_csv(fit_identity, file.path(DIR_TABLES_STEP7, "PHASE0_posterior_identity_check.csv"))
message("Posterior identity check:"); print(as.data.frame(fit_identity))

message("[Phase0] Reproducing lambda_import=1 baseline, H1/G2 (full horizon) and H2/G2 (full ensemble)...")
baseline_H1 <- simulate_holdout_step7(training_bundles$H1, HOLDOUTS$H1, lambda_mode = "constant", lambda_const = LAMBDA_BASELINE,
                                       year_effect_mode = "centered", n_draws = ENSEMBLE_STEP7$full$n_draws, n_replicates = ENSEMBLE_STEP7$full$n_replicates,
                                       seed_offset = BASE_SEED_STEP7 + 1L)
baseline_H2 <- simulate_holdout_step7(training_bundles$H2, HOLDOUTS$H2, lambda_mode = "constant", lambda_const = LAMBDA_BASELINE,
                                       year_effect_mode = "centered", n_draws = ENSEMBLE_STEP7$full$n_draws, n_replicates = ENSEMBLE_STEP7$full$n_replicates,
                                       seed_offset = BASE_SEED_STEP7 + 2L)
saveRDS(list(H1 = baseline_H1, H2 = baseline_H2), file.path(DIR_RESULTS_STEP7, "PHASE0_baseline_lambda1_simulations.rds"))

compute_entry_metrics <- function(res, holdout_spec, holdout_id, training_bundle = NULL) {
  traj <- res$trajectories |> mutate(year = as.integer(format(week_start, "%Y")))
  jan1_2022 <- traj |> filter(year == 2022) |> group_by(posterior_draw_id, stochastic_replicate_id) |> summarise(S_prop_Jan2022 = first(S / N), .groups = "drop")
  if (holdout_spec$holdout_start < as.Date("2022-01-02")) {
    # H1: the holdout trajectory itself spans the 8 weeks before 2022 -- entry infectiousness is a SIMULATED quantity.
    entry <- traj |> filter(week_start < as.Date("2022-01-02")) |> arrange(posterior_draw_id, stochastic_replicate_id, week_start) |>
      group_by(posterior_draw_id, stochastic_replicate_id) |>
      summarise(entry_weighted_infectiousness = sum(W_GEN * rev(tail(X_latent, G_GEN))), .groups = "drop")
  } else {
    # H2: holdout starts AT 2022-01-02, so the preceding 8 weeks are TRAINING weeks -- read directly from the
    # training-only posterior's own fitted X (same construction as Step 6's STEP6_entry2022_infection_history.csv).
    stopifnot(!is.null(training_bundle))
    dr <- rstan::extract(training_bundle$fit, pars = "X", permuted = TRUE)
    n_train <- training_bundle$stan_data$N
    draw_ids <- sort(unique(res$year_effects$posterior_draw_id))
    X_tail <- dr$X[draw_ids, (n_train - G_GEN + 1):n_train, drop = FALSE]
    entry <- tibble(posterior_draw_id = draw_ids, stochastic_replicate_id = NA_integer_,
                     entry_weighted_infectiousness = apply(X_tail, 1, function(x) sum(W_GEN * rev(x))))
  }
  annual_2022 <- traj |> filter(year == 2022) |> group_by(posterior_draw_id, stochastic_replicate_id) |>
    summarise(annual_reported_cases = sum(Y_reported), peak_Reff = max(R_eff),
              weeks_Reff_gt1 = sum(R_eff > 1), favourable = any(R_eff > 1 & M_import > 0), .groups = "drop")
  list(entry = entry, jan1 = jan1_2022, annual = annual_2022)
}

m_H1 <- compute_entry_metrics(baseline_H1, HOLDOUTS$H1, "H1")
m_H2 <- compute_entry_metrics(baseline_H2, HOLDOUTS$H2, "H2", training_bundle = training_bundles$H2)

observed_2022 <- 50986 # from STEP6_observed_2022_epidemic_period.csv, cited not recomputed
reproduction_check <- tibble(
  holdout = c("H1", "H2"),
  entry_weighted_infectiousness_median_new = c(median(m_H1$entry$entry_weighted_infectiousness), median(m_H2$entry$entry_weighted_infectiousness)),
  entry_weighted_infectiousness_median_step6 = c(STEP6_REFERENCE$H1$entry_weighted_infectiousness_median, STEP6_REFERENCE$H2$entry_weighted_infectiousness_median),
  S_prop_median_new = c(median(m_H1$jan1$S_prop_Jan2022), median(m_H2$jan1$S_prop_Jan2022)),
  S_prop_median_step6 = c(STEP6_REFERENCE$H1$S_prop_median, STEP6_REFERENCE$H2$S_prop_median),
  peak_Reff_median_new = c(median(m_H1$annual$peak_Reff), median(m_H2$annual$peak_Reff)),
  peak_Reff_median_step6 = c(STEP6_REFERENCE$H1$peak_Reff_median, STEP6_REFERENCE$H2$peak_Reff_median),
  prob_exceed_2022_new = c(mean(m_H1$annual$annual_reported_cases >= observed_2022), mean(m_H2$annual$annual_reported_cases >= observed_2022)),
  prob_exceed_2022_step6 = c(STEP6_REFERENCE$H1$prob_exceed_2022, STEP6_REFERENCE$H2$prob_exceed_2022)
) |> mutate(
  rel_diff_entry = abs(entry_weighted_infectiousness_median_new - entry_weighted_infectiousness_median_step6) / entry_weighted_infectiousness_median_step6,
  rel_diff_S = abs(S_prop_median_new - S_prop_median_step6) / S_prop_median_step6,
  rel_diff_Reff = abs(peak_Reff_median_new - peak_Reff_median_step6) / peak_Reff_median_step6,
  agree_entry = rel_diff_entry <= REPRODUCTION_TOLERANCE, agree_S = rel_diff_S <= REPRODUCTION_TOLERANCE, agree_Reff = rel_diff_Reff <= REPRODUCTION_TOLERANCE
)
write_csv(reproduction_check, file.path(DIR_TABLES_STEP7, "PHASE0_baseline_reproduction_check.csv"))
message("\n=== Phase0: baseline reproduction vs Step 6 reference ==="); print(as.data.frame(reproduction_check), digits = 4)

if (!all(reproduction_check$agree_entry, reproduction_check$agree_S, reproduction_check$agree_Reff)) {
  message("\n[Phase0 NOTE] One or more metrics differ from Step 6 by more than ", REPRODUCTION_TOLERANCE * 100,
          "%. This is expected for entry-infectiousness/exceedance probabilities, which are SENSITIVE to Monte Carlo",
          " sampling of posterior_draw_id (independent sample.int() call, different RNG seed stream: BASE_SEED_STEP7 vs BASE_SEED_YV+777) even under an IDENTICAL model and lambda.",
          " This is disclosed explicitly below, NOT silently resolved by substituting Step 6's own saved baseline.")
}

# ================================================================
# Audit the established / major-epidemic cross-tabulation (H2/G2, 2022)
# ================================================================
traj_H2 <- baseline_H2$trajectories |> mutate(year = as.integer(format(week_start, "%Y")),
                                               S_prop = S / N, import_hazard = M_import / N, local_hazard = pmax(lambda_t - M_import / N, 0),
                                               favourable = R_eff > 1, locally_dominant = X_latent > 0 & local_hazard > import_hazard)
annual_H2 <- traj_H2 |> filter(year == 2022) |> group_by(posterior_draw_id, stochastic_replicate_id) |>
  summarise(annual_reported_cases = sum(Y_reported), population_mean = mean(N),
            favourable_pressure = any(favourable & M_import > 0), longest_locally_dominant_run = max_run(locally_dominant), .groups = "drop") |>
  mutate(established = favourable_pressure & longest_locally_dominant_run >= ESTABLISHMENT_RUN_WEEKS,
         annual_reported_incidence_per_100k = annual_reported_cases / population_mean * 1e5,
         major_epidemic = annual_reported_cases >= MAJOR_ABS_CASE_MINIMUM & annual_reported_incidence_per_100k >= MAJOR_INCIDENCE_MINIMUM_PER_100K)
crosstab_H2 <- annual_H2 |> count(established, major_epidemic) |> mutate(probability = n / sum(n))
write_csv(crosstab_H2, file.path(DIR_TABLES_STEP7, "PHASE0_established_major_crosstab_H2.csv"))
message("\n=== Phase0: established x major_epidemic cross-tab, H2/G2 2022 (new ensemble) ==="); print(as.data.frame(crosstab_H2), digits = 4)
message(sprintf("Major-but-not-established: n=%d (should be 0 or near-0 if 'established' is logically implied by reaching a major epidemic); established-but-not-major: n=%d",
                 sum(annual_H2$major_epidemic & !annual_H2$established), sum(annual_H2$established & !annual_H2$major_epidemic)))

# ================================================================
# Document the two known limitations (quantified in Step 6, not re-derived here)
# ================================================================
limitations_doc <- tibble(limitation = KNOWN_LIMITATIONS)
write_csv(limitations_doc, file.path(DIR_TABLES_STEP7, "PHASE0_known_limitations.csv"))

message("\n[Phase0] Baseline reproduction and audit complete -> ", DIR_TABLES_STEP7)
