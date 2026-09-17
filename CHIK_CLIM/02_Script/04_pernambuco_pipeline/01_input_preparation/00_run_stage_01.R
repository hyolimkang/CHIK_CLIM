# ============================================================
# STAGE 01 -- INPUT PREPARATION (Pernambuco, orchestrator only)
# ============================================================
# Data audit only. Does not fit or refit any Stan model.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/04_pernambuco_pipeline/01_input_preparation")

message("\n[Stage 01] 01/01: pe_data_audit ...")
source(file.path(STAGE_DIR, "01_pe_data_audit.R"), local = .GlobalEnv)

message("\n[Stage 01] Complete.")
