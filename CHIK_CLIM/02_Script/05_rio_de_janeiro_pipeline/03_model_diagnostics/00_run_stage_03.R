# ============================================================
# STAGE 03 -- MODEL DIAGNOSTICS (Rio de Janeiro, orchestrator only)
# ============================================================
# Read-only: validates/plots the accepted state and city fits. Does not
# fit anything.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/05_rio_de_janeiro_pipeline/03_model_diagnostics")
BASE_OUT <- project_path("03_Output/05_rio_de_janeiro_pipeline/model_fits/v4_9_replication")
state_fit <- file.path(BASE_OUT, "outputs/caseonly/rj_global_q_case_only.rds")
if (!file.exists(state_fit)) {
  stop("Missing prerequisite:\n  ", state_fit, "\n\nRun Stage 02 first:\n  02_transmission_fitting (00_run_stage_02.R)", call. = FALSE)
}
message("[Stage 03] Prerequisites OK.")

message("\n[Stage 03] 01/07: rj_global_q_diagnostics (state) ...")
source(file.path(STAGE_DIR, "01_rj_global_q_diagnostics.R"), local = .GlobalEnv)

message("\n[Stage 03] 02/07: plot_rj_six_panel (state, caseonly) ...")
source(file.path(STAGE_DIR, "02_plot_rj_six_panel.R"), local = .GlobalEnv)
plot_rj_panel("caseonly")

message("\n[Stage 03] 03/07: plot_rj_fixedq_panels (state, q = 0.010-0.020) ...")
source(file.path(STAGE_DIR, "03_plot_rj_fixedq_panels.R"), local = .GlobalEnv)

message("\n[Stage 03] 04/07: rj_u10_caseonly_validation (state vs city U10 external check) ...")
source(file.path(STAGE_DIR, "04_rj_u10_caseonly_validation.R"), local = .GlobalEnv)

if (file.exists(file.path(BASE_OUT, "outputs/city_caseonly/rj_city_global_q_case_only.rds"))) {
  message("\n[Stage 03] 05/07: rj_city_caseonly_diagnostics ...")
  source(file.path(STAGE_DIR, "05_rj_city_caseonly_diagnostics.R"), local = .GlobalEnv)

  message("\n[Stage 03] 06/07: rj_city_u10_diagnostics ...")
  source(file.path(STAGE_DIR, "06_rj_city_u10_diagnostics.R"), local = .GlobalEnv)

  message("\n[Stage 03] 07/07: rj_city_ab_final_comparison ...")
  source(file.path(STAGE_DIR, "07_rj_city_ab_final_comparison.R"), local = .GlobalEnv)
} else {
  message("\n[Stage 03] 05-07/07: SKIPPED -- city case-only fit not present.")
}

message("\n[Stage 03] Complete.")
