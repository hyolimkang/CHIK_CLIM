# ============================================================
# STAGE 03 -- MODEL DIAGNOSTICS (Mato Grosso, orchestrator only)
# ============================================================
# Read-only: validates/plots the accepted fits. Does not fit anything.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/06_mato_grosso_pipeline/03_model_diagnostics")
BASE_OUT <- project_path("03_Output/06_mato_grosso_pipeline/model_fits/v4_9_replication")
primary_fit <- file.path(BASE_OUT, "outputs/caseonly/mt_global_q_case_only.rds")
if (!file.exists(primary_fit)) {
  stop("Missing prerequisite:\n  ", primary_fit, "\n\nRun Stage 02 first:\n  02_transmission_fitting (00_run_stage_02.R)", call. = FALSE)
}
message("[Stage 03] Prerequisites OK.")

message("\n[Stage 03] 01/04: mt_global_q_diagnostics ...")
source(file.path(STAGE_DIR, "01_mt_global_q_diagnostics.R"), local = .GlobalEnv)

message("\n[Stage 03] 02/04: plot_mt_hmc_diagnostics ...")
source(file.path(STAGE_DIR, "02_plot_mt_hmc_diagnostics.R"), local = .GlobalEnv)

message("\n[Stage 03] 03/04: plot_mt_six_panel (caseonly) ...")
source(file.path(STAGE_DIR, "03_plot_mt_six_panel.R"), local = .GlobalEnv)
plot_mt_panel("caseonly")

if (file.exists(file.path(BASE_OUT, "outputs/fixedq/mt_fixedq_q0.020.rds"))) {
  message("\n[Stage 03] 04/04: mt_fixed_q_analysis ...")
  source(file.path(STAGE_DIR, "04_mt_fixed_q_analysis.R"), local = .GlobalEnv)
} else {
  message("\n[Stage 03] 04/04: SKIPPED -- fixed-q profile fits not present.")
}

message("\n[Stage 03] Complete.")
