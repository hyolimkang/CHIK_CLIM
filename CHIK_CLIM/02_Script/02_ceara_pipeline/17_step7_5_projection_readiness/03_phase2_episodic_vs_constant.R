# PHASE 2 -- Model A (existing homogeneous Poisson, lambda=1) vs Model
# B (episodic/bursty, SAME mean pressure lambda=1, prespecified p=0.10
# and p=0.25). H1/G2, full horizon (2020-2025), full ensemble.
# Model A reuses Step 7's already-run, already-validated baseline
# (IDENTICAL seed -- same posterior draws/z_year/z_season/climate/
# initial state) rather than resimulating it.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr) })
filter <- dplyr::filter

source(file.path(ROOT, "02_Script/02_ceara_pipeline/17_step7_5_projection_readiness/00_config_step7_5.R"))
source(file.path(ROOT, "02_Script/02_ceara_pipeline/17_step7_5_projection_readiness/01_episodic_engine.R"))

training_H1 <- readRDS(FIT_PATHS[["H1"]])
modelA <- readRDS(file.path(DIR_RESULTS_STEP7, "PHASE0_baseline_lambda1_simulations.rds"))$H1 # REUSED, not resimulated

message("[Phase2] Simulating Model B (episodic), p=0.10, H1, full ensemble...")
modelB_p10 <- simulate_holdout_step75(training_H1, HOLDOUTS$H1, lambda_mode = "episodic", burst_p = BURST_PROBS[["p10"]], burst_lambda_mean = EPISODIC_LAMBDA,
                                       year_effect_mode = "centered", n_draws = ENSEMBLE_FULL$n_draws, n_replicates = ENSEMBLE_FULL$n_replicates,
                                       seed_offset = SEED_H1_STEP75)
message("[Phase2] Simulating Model B (episodic), p=0.25, H1, full ensemble...")
modelB_p25 <- simulate_holdout_step75(training_H1, HOLDOUTS$H1, lambda_mode = "episodic", burst_p = BURST_PROBS[["p25"]], burst_lambda_mean = EPISODIC_LAMBDA,
                                       year_effect_mode = "centered", n_draws = ENSEMBLE_FULL$n_draws, n_replicates = ENSEMBLE_FULL$n_replicates,
                                       seed_offset = SEED_H1_STEP75)
saveRDS(list(ModelA_constant = modelA, ModelB_p10 = modelB_p10, ModelB_p25 = modelB_p25), file.path(DIR_RESULTS_STEP75, "PHASE2_structural_arms.rds"))

# ================================================================
# Realised-process numerical check: mean M_t and burst frequency match targets.
# ================================================================
process_check <- bind_rows(
  modelA$trajectories |> summarise(arm = "ModelA_constant", mean_M = mean(M_import), frac_weeks_M_gt_0 = mean(M_import > 0)),
  modelB_p10$trajectories |> summarise(arm = "ModelB_p10", mean_M = mean(M_import), frac_weeks_M_gt_0 = mean(M_import > 0)),
  modelB_p25$trajectories |> summarise(arm = "ModelB_p25", mean_M = mean(M_import), frac_weeks_M_gt_0 = mean(M_import > 0))
) |> mutate(target_mean_M = EPISODIC_LAMBDA)
write_csv(process_check, file.path(DIR_TABLES_STEP75, "PHASE2_process_realisation_check.csv"))
message("\n=== Phase2: realised M_t process check (mean should be ~1 for all three arms) ==="); print(as.data.frame(process_check), digits = 4)

message("\n[Phase2] Model A vs Model B simulations complete -> ", DIR_RESULTS_STEP75)
