# Future year-effect validation -- standalone runner.
#
# Historical inference + retrospective holdout validation only. No SSP/
# ISIMIP forcing, no vaccination. Phase 1 fits tau_year (3 priors) on
# the full 2015-2025 period (2022 seed retained, disclosed). Phase 3
# fits tau_year on TRAINING-ONLY data for 2 holdout experiments (seed-
# free builder). This is computationally expensive (5 new Stan fits,
# ~10-14 min each) -- run Phases in order; each is independently
# re-runnable once its upstream fits exist.
#
# Usage:
#   source("02_Script/02_ceara_pipeline/12_future_year_effect_validation/00_run_year_effect_validation.R")

if (!exists("ROOT")) stop("ROOT not found -- open CHIK_CLIM.Rproj or source CHIK_CLIM/00_project_setup.R first.")
DIR_YV <- file.path(ROOT, "02_Script/02_ceara_pipeline/12_future_year_effect_validation")

message("=== PHASE 0: audit + new-model smoke test (GATE) ===")
source(file.path(DIR_YV, "01_phase0_audit.R"))

message("\n=== PHASE 1: fit tau_year (P1-P3), full period, 2022 seed retained ===")
source(file.path(DIR_YV, "02_phase1_fit_tau_year.R"))
source(file.path(DIR_YV, "07_phase1_tau_analysis.R"))

message("\n=== PHASE 2: posterior-predictive generator comparison (G0/G1/G2) ===")
source(file.path(DIR_YV, "05_phase2_generator_comparison.R"))

message("\n=== PHASE 3: fit tau_year on TRAINING-ONLY data (H1, H2; seed-free) ===")
source(file.path(DIR_YV, "04_phase3_fit_holdout_training.R"))

message("\n=== PHASE 3/4: retrospective holdout simulation + validation ===")
source(file.path(DIR_YV, "06_phase3_4_retrospective_validation.R"))

message("\n=== PHASE 5: figures ===")
source(file.path(DIR_YV, "08_make_figures.R"))

message("\n=== Complete. See 03_Output/02_ceara_pipeline/reports/future_year_effect_validation/CE_FUTURE_YEAR_EFFECT_VALIDATION.md ===")
