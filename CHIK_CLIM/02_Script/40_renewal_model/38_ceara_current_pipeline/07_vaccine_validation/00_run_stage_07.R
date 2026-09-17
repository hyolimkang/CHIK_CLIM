# ============================================================
# STAGE 07 -- VACCINE VALIDATION (orchestrator only)
# ============================================================
#
# Order: 01 vaccine engine unit tests -> 02 routine age-12 eligibility QA
# -> 03 actual historical A/B/C counterfactual QA. Each script fails
# loudly (stop()) on a substantive QA failure; this runner does NOT catch
# or suppress those errors -- a failure here must halt the pipeline, not
# be silently passed through to Stage 08 impact/NNV summaries.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}

STAGE_DIR <- file.path(project_path("02_Script/40_renewal_model/38_ceara_current_pipeline"), "07_vaccine_validation")

climate_fit_path <- project_path("02_Script/40_renewal_model/36_climate_forced_v4_9/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds")
if (!file.exists(climate_fit_path)) {
  stop("Missing prerequisite:\n  ", climate_fit_path,
       "\n\nRun Stage 02 first:\n  02_historical_transmission_fitting (00_run_stage_02.R)", call. = FALSE)
}

message("[Stage 07] 01/03: vaccine_engine_unit_tests ...")
source(file.path(STAGE_DIR, "01_vaccine_engine_unit_tests.R"), local = .GlobalEnv)

message("\n[Stage 07] 02/03: routine_age12_eligibility_qa ...")
source(file.path(STAGE_DIR, "02_routine_age12_eligibility_qa.R"), local = .GlobalEnv)

three_arm_results_path <- project_path("02_Script/40_renewal_model/37_ceara_historical_age12_vaccine_counterfactual/results/three_arm_results.rds")
if (!file.exists(three_arm_results_path)) {
  stop("Missing prerequisite:\n  ", three_arm_results_path,
       "\n\nRun Stage 06 first:\n  06_vaccine_counterfactual (00_run_stage_06.R)", call. = FALSE)
}
message("\n[Stage 07] 03/03: validate_historical_counterfactual (fails loudly on any QA failure) ...")
source(file.path(STAGE_DIR, "03_validate_historical_counterfactual.R"), local = .GlobalEnv)

message("\n[Stage 07] Complete. All vaccine validation checks passed.")
