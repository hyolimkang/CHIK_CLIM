# ============================================================
# STAGE 06 -- VACCINE COUNTERFACTUAL (orchestrator only)
# ============================================================
#
# Replaces the old, ambiguous 37/.../scripts/99_run_all.R. Sources 01-07 in
# dependency order and runs the three-arm posterior simulation (paired
# no-vaccine / full-feedback-vaccine / direct-only-diagnostic arms across
# 200 posterior draws). This is SIMULATION using the accepted fit, NOT a
# Stan refit -- no scientific change from the original script sequence.
#
# Output location is intentionally unchanged (Section 19): results are
# still written under 37_ceara_historical_age12_vaccine_counterfactual/results/
# (see 01_config.R's ANALYSIS_DIR), not under this pipeline directory.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}

STAGE_DIR <- file.path(project_path("02_Script/40_renewal_model/38_ceara_current_pipeline"), "06_vaccine_counterfactual")

climate_fit_path <- project_path("02_Script/40_renewal_model/36_climate_forced_v4_9/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds")
age_shares_path <- project_path("03_Output/tables/climate_forced_v4_9/CE_age_population_shares.csv")
missing <- c(climate_fit_path, age_shares_path)[!file.exists(c(climate_fit_path, age_shares_path))]
if (length(missing) > 0) {
  stop("Missing prerequisite(s):\n  ", paste(missing, collapse = "\n  "),
       "\n\nIf the age population shares are missing, run Stage 04 first:\n  04_age_reconstruction (00_run_stage_04.R)\n",
       "If the climate fit is missing, run Stage 02 first:\n  02_historical_transmission_fitting (00_run_stage_02.R)", call. = FALSE)
}
message("[Stage 06] Prerequisites OK: accepted climate fit and age population shares both found.")

setwd(ROOT)
source(file.path(STAGE_DIR, "01_config.R"), local = .GlobalEnv)
source(file.path(STAGE_DIR, "03_generic_age_vaccine_simulator.R"), local = .GlobalEnv)
source(file.path(STAGE_DIR, "02_load_accepted_inputs.R"), local = .GlobalEnv)
source(file.path(STAGE_DIR, "04_run_no_vaccine_arm.R"), local = .GlobalEnv)
source(file.path(STAGE_DIR, "05_run_full_feedback_vaccine_arm.R"), local = .GlobalEnv)
source(file.path(STAGE_DIR, "06_run_direct_only_diagnostic_arm.R"), local = .GlobalEnv)
source(file.path(STAGE_DIR, "07_execute_three_arms_all_draws.R"), local = .GlobalEnv)

inputs <- load_accepted_inputs()
execute_three_arms_all_draws(inputs)

message("\n[Stage 06] Complete. Three-arm execution finished; results saved under the existing 37_.../results/ directory.")
