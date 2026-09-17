# ============================================================
# STAGE 01 -- INPUT PREPARATION (Mato Grosso, orchestrator only)
# ============================================================
# Data/episode/input audits only. Does not fit or refit any Stan model.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/06_mato_grosso_pipeline/01_input_preparation")

message("\n[Stage 01] 01/03: mt_data_audit ...")
source(file.path(STAGE_DIR, "01_mt_data_audit.R"), local = .GlobalEnv)

message("\n[Stage 01] 02/03: mt_episode_audit ...")
source(file.path(STAGE_DIR, "02_mt_episode_audit.R"), local = .GlobalEnv)

message("\n[Stage 01] 03/03: mt_v49_input_audit ...")
source(file.path(STAGE_DIR, "03_mt_v49_input_audit.R"), local = .GlobalEnv)

message("\n[Stage 01] Complete.")
