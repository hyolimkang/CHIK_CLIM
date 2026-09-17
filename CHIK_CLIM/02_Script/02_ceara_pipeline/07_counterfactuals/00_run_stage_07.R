# ============================================================
# STAGE 07 -- COUNTERFACTUALS (orchestrator only)
# ============================================================
#
# Merges what were previously two separate stages in 38_ceara_current_pipeline
# (06_vaccine_counterfactual, 07_vaccine_validation) into one, since the
# validation scripts only make sense run against this stage's own output.
# Execution order: 01-07 (config, simulator, load inputs, three arms,
# execute across posterior draws -> produces three_arm_results.rds), then
# 08-10 (vaccine engine unit tests, routine age-12 eligibility QA, then the
# actual historical A/B/C counterfactual QA -- fails loudly on any
# substantive QA failure; this runner does NOT catch or suppress that
# error, since a failure here must halt the pipeline, not be silently
# passed through to Stage 08 impact/NNV summaries).
#
# This is SIMULATION using the accepted fit, NOT a Stan refit. Results are
# written under 03_Output/results/ceara_historical_age12_vaccine_counterfactual/
# (see 01_config.R's DIR_RESULTS), not under this pipeline directory.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}

STAGE_DIR <- project_path("02_Script/02_ceara_pipeline/07_counterfactuals")

climate_fit_path <- project_path("03_Output/model_fits/ceara/climate_forced_v4_9/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds")
age_shares_path <- project_path("03_Output/tables/climate_forced_v4_9/CE_age_population_shares.csv")
missing <- c(climate_fit_path, age_shares_path)[!file.exists(c(climate_fit_path, age_shares_path))]
if (length(missing) > 0) {
  stop("Missing prerequisite(s):\n  ", paste(missing, collapse = "\n  "),
       "\n\nIf the age population shares are missing, run Stage 06 first:\n  06_age_demographic_extension (00_run_stage_06.R)\n",
       "If the climate fit is missing, run Stage 02 first:\n  02_transmission_fitting (00_run_stage_02.R)", call. = FALSE)
}
message("[Stage 07] Prerequisites OK: accepted climate fit and age population shares both found.")

setwd(ROOT)
message("\n[Stage 07] 01/10: config ...")
source(file.path(STAGE_DIR, "01_config.R"), local = .GlobalEnv)
message("[Stage 07] 02/10: generic_age_vaccine_simulator ...")
source(file.path(STAGE_DIR, "03_generic_age_vaccine_simulator.R"), local = .GlobalEnv)
message("[Stage 07] 03/10: load_accepted_inputs ...")
source(file.path(STAGE_DIR, "02_load_accepted_inputs.R"), local = .GlobalEnv)
message("[Stage 07] 04/10: run_no_vaccine_arm ...")
source(file.path(STAGE_DIR, "04_run_no_vaccine_arm.R"), local = .GlobalEnv)
message("[Stage 07] 05/10: run_full_feedback_vaccine_arm ...")
source(file.path(STAGE_DIR, "05_run_full_feedback_vaccine_arm.R"), local = .GlobalEnv)
message("[Stage 07] 06/10: run_direct_only_diagnostic_arm ...")
source(file.path(STAGE_DIR, "06_run_direct_only_diagnostic_arm.R"), local = .GlobalEnv)
message("[Stage 07] 07/10: execute_three_arms_all_draws ...")
source(file.path(STAGE_DIR, "07_execute_three_arms_all_draws.R"), local = .GlobalEnv)

inputs <- load_accepted_inputs()
execute_three_arms_all_draws(inputs)
message("[Stage 07] Three-arm execution complete.")

message("\n[Stage 07] 08/10: vaccine_engine_unit_tests ...")
source(file.path(STAGE_DIR, "08_vaccine_engine_unit_tests.R"), local = .GlobalEnv)

message("\n[Stage 07] 09/10: routine_age12_eligibility_qa ...")
source(file.path(STAGE_DIR, "09_routine_age12_eligibility_qa.R"), local = .GlobalEnv)

message("\n[Stage 07] 10/10: validate_historical_counterfactual (fails loudly on any QA failure) ...")
source(file.path(STAGE_DIR, "10_validate_historical_counterfactual.R"), local = .GlobalEnv)

message("\n[Stage 07] Complete. Results saved under 03_Output/results/ceara_historical_age12_vaccine_counterfactual/.")
