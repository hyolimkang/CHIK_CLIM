# ============================================================
# STAGE 01 -- INPUT PREPARATION (Rio de Janeiro, orchestrator only)
# ============================================================
# Data audits only (state + city). Does not fit or refit any Stan model.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/05_rio_de_janeiro_pipeline/01_input_preparation")

message("\n[Stage 01] 01/04: rj_data_audit ...")
source(file.path(STAGE_DIR, "01_rj_data_audit.R"), local = .GlobalEnv)

message("\n[Stage 01] 02/04: rj_episode_audit ...")
source(file.path(STAGE_DIR, "02_rj_episode_audit.R"), local = .GlobalEnv)

message("\n[Stage 01] 03/04: rj_v49_input_audit ...")
source(file.path(STAGE_DIR, "03_rj_v49_input_audit.R"), local = .GlobalEnv)

message("\n[Stage 01] 04/04: rj_city_data_audit ...")
source(file.path(STAGE_DIR, "04_rj_city_data_audit.R"), local = .GlobalEnv)

message("\n[Stage 01] Complete.")
