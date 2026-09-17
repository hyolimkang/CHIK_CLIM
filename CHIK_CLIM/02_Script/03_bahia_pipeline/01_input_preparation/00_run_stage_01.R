# ============================================================
# STAGE 01 -- INPUT PREPARATION (Bahia, orchestrator only)
# ============================================================
if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/03_bahia_pipeline/01_input_preparation")

panel_path <- project_path("01_Data/chik_dlnm_panel_muni_week_2015_2025.rds")
if (!file.exists(panel_path)) {
  stop("Missing prerequisite:\n  ", panel_path, "\n\nRun the data pipeline first:\n  01_data_pipeline (99_run_data_pipeline.R)", call. = FALSE)
}
message("\n[Stage 01] 01/01: bahia_data_audit ...")
source(file.path(STAGE_DIR, "01_bahia_data_audit.R"), local = .GlobalEnv)
message("\n[Stage 01] Complete.")
