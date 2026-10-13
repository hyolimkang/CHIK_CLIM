# Climate robustness assessment -- standalone runner.
#
# Historical inference only (2015-2025). Phase 1 is a GATE: if it fails,
# this stops before any new Stan fitting. Phase 2 fits M1-M4 (~1 hour on
# this machine, 4 chains each); M0 is read-only (never refit).
#
# Usage:
#   source("02_Script/02_ceara_pipeline/11_climate_robustness/00_run_climate_robustness.R")

if (!exists("ROOT")) stop("ROOT not found -- open CHIK_CLIM.Rproj or source CHIK_CLIM/00_project_setup.R first.")
DIR_CR <- file.path(ROOT, "02_Script/02_ceara_pipeline/11_climate_robustness")

message("=== PHASE 1: climate input consistency audit (GATE) ===")
source(file.path(DIR_CR, "01_phase1_climate_input_audit.R")) # stops with an error if the gate fails

message("\n=== PHASE 2: fit M1-M4 (M0 reused, read-only) ===")
source(file.path(DIR_CR, "02_phase2_fit_M1_M4.R"))

message("\n=== PHASE 3: scientific comparison ===")
source(file.path(DIR_CR, "03_phase3_scientific_comparison.R"))

message("\n=== PHASE 4: figures ===")
source(file.path(DIR_CR, "04_phase4_figures.R"))

message("\n=== Climate robustness assessment complete. See 03_Output/02_ceara_pipeline/reports/climate_robustness/CE_CLIMATE_ROBUSTNESS_ASSESSMENT.md ===")
