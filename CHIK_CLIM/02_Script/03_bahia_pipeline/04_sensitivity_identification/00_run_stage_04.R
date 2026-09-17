# ============================================================
# STAGE 04 -- SENSITIVITY / IDENTIFICATION (Bahia, orchestrator only)
# ============================================================
# The global-q (estimated, not fixed) identification experiment: prior
# predictive check -> identification diagnostics -> plots. Read-only
# against the global-q multisite fit; does not fit anything itself.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/03_bahia_pipeline/04_sensitivity_identification")
globalq_fit <- project_path("03_Output/03_bahia_pipeline/model_fits/v4_9_replication/outputs_global_q_multisite_serology/on_full2/renewal_bahia_global_q_fit_on_full2.rds")
if (!file.exists(globalq_fit)) {
  stop("Missing prerequisite:\n  ", globalq_fit, "\n\nRun Stage 02 first:\n  02_transmission_fitting (00_run_stage_02.R)", call. = FALSE)
}
message("[Stage 04] Prerequisites OK.")

message("\n[Stage 04] 01/03: prior_predictive_global_q ...")
source(file.path(STAGE_DIR, "01_prior_predictive_global_q.R"), local = .GlobalEnv)

message("\n[Stage 04] 02/03: global_q_identification_diagnostics ...")
source(file.path(STAGE_DIR, "02_global_q_identification_diagnostics.R"), local = .GlobalEnv)

message("\n[Stage 04] 03/03: plot_global_q_diagnostics ...")
source(file.path(STAGE_DIR, "03_plot_global_q_diagnostics.R"), local = .GlobalEnv)

message("\n[Stage 04] Complete.")
