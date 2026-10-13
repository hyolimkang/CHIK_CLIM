# PHASE 2 -- separate pre-2022 persistence from 2022-onward
# reintroduction. H1/G2 only (H2's holdout starts exactly at 2022, so
# the pre/post split is not meaningful there -- by construction H2 has
# no "pre-2022" holdout period to vary).
#
# Experiment A: vary lambda_import for 2020-2021, FIX lambda_import=1
# from 2022 onward. Isolates the effect of pre-2022 introduction
# pressure/standing infections on 2022 outcomes.
# Experiment B: FIX lambda_import=1 for 2020-2021, vary lambda_import
# from 2022 onward. Isolates the effect of introduction pressure DURING
# the resurgence period itself.
#
# Both experiments' lambda=1/lambda=1 arm is IDENTICAL to Phase 1's
# H1 lambda=1 baseline (same seed, same split point degenerates to a
# constant) -- reused, not resimulated.

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
SEED_H1 <- BASE_SEED_STEP7 + 1L # SAME seed as Phase 0/1's H1 arms -- common random numbers throughout
SPLIT_DATE <- as.Date("2022-01-02")

summarise_split_arm <- function(res, exp_id, arm_lambda) {
  traj <- res$trajectories |> mutate(year = as.integer(format(week_start, "%Y")),
                                      import_hazard = M_import / N, local_hazard = pmax(lambda_t - M_import / N, 0),
                                      favourable = R_eff > 1, locally_dominant = X_latent > 0 & local_hazard > import_hazard)
  pre2022 <- traj |> filter(week_start < SPLIT_DATE)
  entry <- pre2022 |> arrange(posterior_draw_id, stochastic_replicate_id, week_start) |> group_by(posterior_draw_id, stochastic_replicate_id) |>
    summarise(entry_weighted_infectiousness = sum(W_GEN * rev(tail(X_latent, G_GEN))), .groups = "drop")
  annual_2022 <- traj |> filter(year == 2022) |> group_by(posterior_draw_id, stochastic_replicate_id) |>
    summarise(peak_Reff = max(R_eff), cumulative_latent_infections = sum(X_latent), annual_reported_cases = sum(Y_reported),
              favourable_pressure = any(favourable & M_import > 0), longest_locally_dominant_run = max_run(locally_dominant),
              population_mean = mean(N),
              first_week_sustained_dominance = { idx <- which(locally_dominant); if (!length(idx)) NA_integer_ else idx[1] },
              .groups = "drop") |>
    mutate(established = favourable_pressure & longest_locally_dominant_run >= ESTABLISHMENT_RUN_WEEKS,
           annual_reported_incidence_per_100k = annual_reported_cases / population_mean * 1e5,
           major_epidemic = annual_reported_cases >= MAJOR_ABS_CASE_MINIMUM & annual_reported_incidence_per_100k >= MAJOR_INCIDENCE_MINIMUM_PER_100K)
  observed_2022 <- 50986
  tibble(experiment = exp_id, arm = arm_lambda,
         median_entry_weighted_infectiousness = median(entry$entry_weighted_infectiousness),
         lo95_entry = qv(entry$entry_weighted_infectiousness, .025), hi95_entry = qv(entry$entry_weighted_infectiousness, .975),
         median_time_to_sustained_dominance = median(annual_2022$first_week_sustained_dominance, na.rm = TRUE),
         prob_no_sustained_dominance = mean(is.na(annual_2022$first_week_sustained_dominance)),
         median_peak_Reff = median(annual_2022$peak_Reff),
         median_cumulative_latent_2022 = median(annual_2022$cumulative_latent_infections),
         median_annual_reported_cases_2022 = median(annual_2022$annual_reported_cases),
         prob_established = mean(annual_2022$established), prob_major_epidemic = mean(annual_2022$major_epidemic),
         prob_exceed_observed_2022 = mean(annual_2022$annual_reported_cases >= observed_2022), n = nrow(annual_2022))
}

run_experiment <- function(exp_id, lambda_pre_values, lambda_post_values) {
  rows <- list()
  for (i in seq_along(lambda_pre_values)) {
    lp <- lambda_pre_values[i]; lo <- lambda_post_values[i]
    arm_label <- sprintf("pre=%.2f_post=%.2f", lp, lo)
    if (lp == LAMBDA_BASELINE && lo == LAMBDA_BASELINE) {
      message(sprintf("  [%s / %s] reusing Phase0/1 lambda=1 baseline (same seed)", exp_id, arm_label))
      res <- phase0$H1
    } else {
      message(sprintf("  [%s / %s] simulating...", exp_id, arm_label))
      res <- simulate_holdout_step7(training_H1, HOLDOUTS$H1, lambda_mode = "split", lambda_pre = lp, lambda_post = lo, split_date = SPLIT_DATE,
                                     year_effect_mode = "centered", n_draws = ENSEMBLE_STEP7$full$n_draws, n_replicates = ENSEMBLE_STEP7$full$n_replicates,
                                     seed_offset = SEED_H1)
    }
    rows[[arm_label]] <- summarise_split_arm(res, exp_id, arm_label)
    saveRDS(res, file.path(DIR_RESULTS_STEP7, sprintf("PHASE2_%s_%s_trajectories.rds", exp_id, gsub("[.=]", "", arm_label))))
  }
  bind_rows(rows)
}

message("[Phase2] Experiment A: vary PRE-2022 lambda_import, fix POST-2022 at 1...")
expA <- run_experiment("ExperimentA_pre2022_varied", LAMBDA_GRID, rep(LAMBDA_BASELINE, length(LAMBDA_GRID))) |> mutate(varied_period = "pre-2022 (2020-2021)", lambda_pre = LAMBDA_GRID, lambda_post = LAMBDA_BASELINE)

message("[Phase2] Experiment B: fix PRE-2022 lambda_import at 1, vary POST-2022...")
expB <- run_experiment("ExperimentB_post2022_varied", rep(LAMBDA_BASELINE, length(LAMBDA_GRID)), LAMBDA_GRID) |> mutate(varied_period = "2022-onward", lambda_pre = LAMBDA_BASELINE, lambda_post = LAMBDA_GRID)

phase2_summary <- bind_rows(expA, expB)
write_csv(phase2_summary, file.path(DIR_TABLES_STEP7, "PHASE2_pre_post_2022_split_summary.csv"))
message("\n=== PHASE 2: pre-2022 vs 2022-onward importation experiments (H1/G2) ===")
print(as.data.frame(phase2_summary |> select(experiment, varied_period, lambda_pre, lambda_post, median_entry_weighted_infectiousness,
                                               prob_established, prob_major_epidemic, prob_exceed_observed_2022, median_annual_reported_cases_2022)), digits = 4)

message("\n[Phase2] Pre/post-2022 split experiments complete -> ", DIR_TABLES_STEP7)
