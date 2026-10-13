# STEP 7.5 -- structural discrimination and projection-readiness:
# isolated, read-only configuration. Reuses H1/H2 training-only
# posteriors and the Step 7 simulation engine UNMODIFIED (only a new,
# additive rng_pois closure for the episodic/bursty process is added in
# 01_episodic_engine.R, never editing 16_step7_importation_persistence's
# own files). NO Stan refitting. PROVISIONAL per the HMC audit.

if (!exists("ROOT")) stop("ROOT not found -- source CHIK_CLIM/00_project_setup.R first.")

DIR_STEP6 <- file.path(ROOT, "02_Script/02_ceara_pipeline/15_step6_retrospective_failure_diagnosis")
DIR_STEP7 <- file.path(ROOT, "02_Script/02_ceara_pipeline/16_step7_importation_persistence")

DIR_FITS_YV       <- file.path(ROOT, "03_Output/02_ceara_pipeline/model_fits/future_year_effect_validation")
DIR_TABLES_STEP6  <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/step6_retrospective_failure_diagnosis")
DIR_TABLES_STEP7  <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/step7_importation_persistence")
DIR_RESULTS_STEP7 <- file.path(ROOT, "03_Output/02_ceara_pipeline/results/step7_importation_persistence")
DIR_TABLES_HMC    <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/hmc_reliability_audit")

DIR_SCRIPTS_STEP75 <- file.path(ROOT, "02_Script/02_ceara_pipeline/17_step7_5_projection_readiness")
DIR_REPORTS_STEP75 <- file.path(ROOT, "03_Output/02_ceara_pipeline/reports/step7_5_projection_readiness")
DIR_TABLES_STEP75  <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/step7_5_projection_readiness")
DIR_FIGURES_STEP75 <- file.path(ROOT, "03_Output/02_ceara_pipeline/figures/step7_5_projection_readiness")
DIR_RESULTS_STEP75 <- file.path(ROOT, "03_Output/02_ceara_pipeline/results/step7_5_projection_readiness")
for (d in c(DIR_REPORTS_STEP75, DIR_TABLES_STEP75, DIR_FIGURES_STEP75, DIR_RESULTS_STEP75)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

FIT_PATHS <- c(H1 = file.path(DIR_FITS_YV, "ce_tau_year_training_H1.rds"), H2 = file.path(DIR_FITS_YV, "ce_tau_year_training_H2.rds"))
required_inputs <- c(FIT_PATHS,
  file.path(DIR_RESULTS_STEP7, "PHASE0_baseline_lambda1_simulations.rds"),
  file.path(DIR_RESULTS_STEP7, "PHASE1_H1_importation_grid.rds"),
  file.path(DIR_TABLES_STEP7, "PHASE1_importation_grid_summary.csv"),
  file.path(DIR_TABLES_STEP7, "PHASE2_pre_post_2022_split_summary.csv"),
  file.path(DIR_TABLES_HMC, "HMC_AUDIT_model_level.csv"))
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Missing required Step 7.5 input(s) -- cannot reliably reproduce the baseline, stopping per instruction: ", paste(missing_inputs, collapse = ", "))

source(file.path(DIR_STEP7, "00_config_step7.R"))     # HOLDOUTS, LAMBDA_GRID, ENSEMBLE_STEP7, BASE_SEED_STEP7, PROVISIONAL_NOTE, KNOWN_LIMITATIONS
source(file.path(DIR_STEP7, "01_simulation_engine.R")) # simulate_holdout_step7(), load_observed_weekly_step7(), make_rng_pois_step7(), G_GEN, W_GEN
source(file.path(DIR_STEP6, "00_config_step6.R"))      # ESTABLISHMENT_RUN_WEEKS, MAJOR_ABS_CASE_MINIMUM/INCIDENCE

# ---- Episodic/bursty introduction process (Phase 2), prespecified BEFORE any outcome is inspected ----
# E[M_t] = lambda is preserved exactly for every (p, lambda) pair tested, matching Model A's mean.
BURST_PROBS <- c(p10 = 0.10, p25 = 0.25)
EPISODIC_LAMBDA <- LAMBDA_BASELINE # = 1, the SAME mean pressure as Step 7's primary comparison

SEED_H1_STEP75 <- BASE_SEED_STEP7 + 1L # IDENTICAL seed to Step 7's H1 baseline -- same sampled posterior draws/z_year/z_season/climate
ENSEMBLE_SMOKE <- list(n_draws = 10L, n_replicates = 5L)
ENSEMBLE_FULL  <- ENSEMBLE_STEP7$full # 100 x 20, reused definition

message(sprintf("[step7_5 config] ready. Burst probs: %s | episodic lambda (mean)=%.2f",
                 paste(names(BURST_PROBS), BURST_PROBS, sep = "=", collapse = ", "), EPISODIC_LAMBDA))
