# ============================================================
# STAGE 04 -- AGE RECONSTRUCTION (orchestrator only)
# ============================================================
#
# Execution order is 01, 03, 04, 05. Script 02 (02_age_cohort_simulator.R)
# is a HELPER (defines simulate_age_cohort()) sourced by 03/04 where
# needed -- it is never run as its own standalone analysis step.
# Read-only against the accepted climate fit; does not fit anything.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}

STAGE_DIR <- file.path(project_path("02_Script/40_renewal_model/38_ceara_current_pipeline"), "04_age_reconstruction")

climate_fit_path <- project_path("02_Script/40_renewal_model/36_climate_forced_v4_9/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds")
age_pop_source_path <- project_path("01_Data/ibge_pop_uf_single_age_expanded_2015_2024.rds")
missing <- c(climate_fit_path, age_pop_source_path)[!file.exists(c(climate_fit_path, age_pop_source_path))]
if (length(missing) > 0) {
  stop("Missing prerequisite(s):\n  ", paste(missing, collapse = "\n  "),
       "\n\nIf the climate fit is missing, run Stage 02 first:\n  02_historical_transmission_fitting (00_run_stage_02.R)", call. = FALSE)
}
message("[Stage 04] Prerequisites OK: accepted climate fit and age population source both found.")

message("\n[Stage 04] 01/04: build_age_population_shares ...")
source(file.path(STAGE_DIR, "01_build_age_population_shares.R"), local = .GlobalEnv)

message("\n[Stage 04] 02/04: validate_age_aggregation_identity ...")
source(file.path(STAGE_DIR, "03_validate_age_aggregation_identity.R"), local = .GlobalEnv) # sources 02_age_cohort_simulator.R internally as a helper

message("\n[Stage 04] 03/04: reconstruct_age_specific_susceptibility ...")
source(file.path(STAGE_DIR, "04_reconstruct_age_specific_susceptibility.R"), local = .GlobalEnv) # sources 02_age_cohort_simulator.R internally as a helper

message("\n[Stage 04] 04/04: plot_age_specific_susceptibility ...")
source(file.path(STAGE_DIR, "05_plot_age_specific_susceptibility.R"), local = .GlobalEnv)

message("\n[Stage 04] Complete.")
