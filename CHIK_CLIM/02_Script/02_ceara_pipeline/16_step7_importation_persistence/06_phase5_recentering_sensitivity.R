# PHASE 5 -- generator-definition sensitivity (OPTIONAL diagnostic arm).
# Audits the known future-year-effect recentering issue (Step 6 FAIL
# item): the accepted generator draws z_year_h ~ std_normal() per
# holdout year and CENTRES it (z_year_h - mean(z_year_h)) within the
# finite holdout horizon, shrinking the marginal SD of each year's
# realised effect from tau_year to tau_year*sqrt((K-1)/K). This arm
# instead uses UNCENTRED year effects (year_effect = tau_year*z_year_h,
# no within-horizon mean subtraction) -- genuine iid marginal draws
# around alpha_R, respecting alpha_R's own historical-centring
# convention (alpha_R itself is untouched; only the per-year deviation
# is no longer forced to sum to zero over the finite future horizon).
# lambda_import=1 throughout (ONLY the year-effect convention changes --
# never combined with an importation change in the same arm, per the
# task's instruction not to combine multiple structural changes).

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr) })
filter <- dplyr::filter

source(file.path(ROOT, "02_Script/02_ceara_pipeline/16_step7_importation_persistence/00_config_step7.R"))
source(file.path(ROOT, "02_Script/02_ceara_pipeline/16_step7_importation_persistence/01_simulation_engine.R"))
source(file.path(DIR_STEP6, "00_config_step6.R"))

max_run <- function(x) { rr <- rle(replace(x, is.na(x), FALSE)); if (!any(rr$values)) 0L else max(rr$lengths[rr$values]) }
qv <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE, type = 8))

training_H1 <- readRDS(FIT_PATHS[["H1"]])
phase0 <- readRDS(file.path(DIR_RESULTS_STEP7, "PHASE0_baseline_lambda1_simulations.rds"))
SEED_H1 <- BASE_SEED_STEP7 + 1L

message("[Phase5] Simulating H1/G2 with UNCENTRED future year effects (lambda_import=1, isolated single change)...")
uncentered <- simulate_holdout_step7(training_H1, HOLDOUTS$H1, lambda_mode = "constant", lambda_const = LAMBDA_BASELINE,
                                      year_effect_mode = "uncentered", n_draws = ENSEMBLE_STEP7$full$n_draws, n_replicates = ENSEMBLE_STEP7$full$n_replicates,
                                      seed_offset = SEED_H1)
saveRDS(uncentered, file.path(DIR_RESULTS_STEP7, "PHASE5_H1_uncentered_year_effect_trajectories.rds"))

summarise_recenter_arm <- function(res, label) {
  traj <- res$trajectories |> mutate(year = as.integer(format(week_start, "%Y")), import_hazard = M_import / N, local_hazard = pmax(lambda_t - M_import / N, 0),
                                      favourable = R_eff > 1, locally_dominant = X_latent > 0 & local_hazard > import_hazard)
  pre2022 <- traj |> filter(week_start < as.Date("2022-01-02"))
  entry <- pre2022 |> arrange(posterior_draw_id, stochastic_replicate_id, week_start) |> group_by(posterior_draw_id, stochastic_replicate_id) |>
    summarise(entry_weighted_infectiousness = sum(W_GEN * rev(tail(X_latent, G_GEN))), .groups = "drop")
  annual_2022 <- traj |> filter(year == 2022) |> group_by(posterior_draw_id, stochastic_replicate_id) |>
    summarise(peak_R0 = max(R0), peak_Reff = max(R_eff), annual_reported_cases = sum(Y_reported), cumulative_latent_infections = sum(X_latent),
              favourable_pressure = any(favourable & M_import > 0), longest_locally_dominant_run = max_run(locally_dominant), population_mean = mean(N), .groups = "drop") |>
    mutate(established = favourable_pressure & longest_locally_dominant_run >= ESTABLISHMENT_RUN_WEEKS,
           annual_reported_incidence_per_100k = annual_reported_cases / population_mean * 1e5,
           major_epidemic = annual_reported_cases >= MAJOR_ABS_CASE_MINIMUM & annual_reported_incidence_per_100k >= MAJOR_INCIDENCE_MINIMUM_PER_100K)
  ye <- res$year_effects |> group_by(posterior_draw_id) |> summarise(sd_year_effect_within_horizon = sd(year_effect), mean_year_effect = mean(year_effect), .groups = "drop")
  tibble(arm = label, median_sd_year_effect_realised = median(ye$sd_year_effect_within_horizon),
         median_entry_weighted_infectiousness = median(entry$entry_weighted_infectiousness),
         median_peak_R0 = median(annual_2022$peak_R0), median_peak_Reff = median(annual_2022$peak_Reff),
         median_cumulative_latent_2022 = median(annual_2022$cumulative_latent_infections),
         median_annual_reported_cases_2022 = median(annual_2022$annual_reported_cases),
         prob_established = mean(annual_2022$established), prob_major_epidemic = mean(annual_2022$major_epidemic),
         prob_exceed_observed_2022 = mean(annual_2022$annual_reported_cases >= 50986), n = nrow(annual_2022))
}
recenter_summary <- bind_rows(summarise_recenter_arm(phase0$H1, "centered_accepted_convention"), summarise_recenter_arm(uncentered, "uncentered_diagnostic_arm"))
write_csv(recenter_summary, file.path(DIR_TABLES_STEP7, "PHASE5_recentering_sensitivity_summary.csv"))
message("\n=== PHASE 5: centred (accepted) vs uncentred future year-effect convention, H1/G2, lambda=1 ===")
print(as.data.frame(recenter_summary), digits = 4)

message("\n[Phase5] Recentering-sensitivity diagnostic arm complete -> ", DIR_TABLES_STEP7)
