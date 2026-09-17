# ============================================================
# STAGE 05 -- CLOSED-LOOP VALIDATION (orchestrator only)
# ============================================================
#
# Answers only: does the closed-loop age simulator reproduce the accepted
# fitted historical transmission when vaccination is OFF? No vaccine unit
# tests here (those are Stage 07), and no Stan refit.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}

STAGE_DIR <- file.path(project_path("02_Script/40_renewal_model/38_ceara_current_pipeline"), "05_closed_loop_validation")

climate_fit_path <- project_path("02_Script/40_renewal_model/36_climate_forced_v4_9/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds")
age_shares_path <- project_path("03_Output/tables/climate_forced_v4_9/CE_age_population_shares.csv")
missing <- c(climate_fit_path, age_shares_path)[!file.exists(c(climate_fit_path, age_shares_path))]
if (length(missing) > 0) {
  stop("Missing prerequisite(s):\n  ", paste(missing, collapse = "\n  "),
       "\n\nRun Stage 04 first:\n  04_age_reconstruction (00_run_stage_04.R)", call. = FALSE)
}
message("[Stage 05] Prerequisites OK: accepted climate fit and age reconstruction inputs both found.")

message("\n[Stage 05] 01/01: validate_closed_loop_replay ...")
source(file.path(STAGE_DIR, "01_validate_closed_loop_replay.R"), local = .GlobalEnv)

message("\n[Stage 05] Complete.")
