# ============================================================
# STAGE 03 -- MODEL DIAGNOSTICS (Pernambuco, orchestrator only)
# ============================================================
# Read-only: validates/plots the accepted global-q and fixed-q-sweep fits.
# Does not fit anything.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/04_pernambuco_pipeline/03_model_diagnostics")
BASE_OUT <- project_path("03_Output/model_fits/pernambuco/v4_9_replication")
required_tags <- c("modelA_caseonly", "modelB_U14")
fit_paths <- file.path(BASE_OUT, "outputs", required_tags, paste0("renewal_pe_global_q_fit_", required_tags, ".rds"))
missing <- fit_paths[!file.exists(fit_paths)]
if (length(missing) > 0) {
  stop("Missing prerequisite(s):\n  ", paste(missing, collapse = "\n  "),
       "\n\nRun Stage 02 first:\n  02_transmission_fitting (00_run_stage_02.R)", call. = FALSE)
}
message("[Stage 03] Prerequisites OK.")

message("\n[Stage 03] 01/06: run_pe_diagnostics ...")
source(file.path(STAGE_DIR, "01_run_pe_diagnostics.R"), local = .GlobalEnv)

message("\n[Stage 03] 02/06: plot_pe_figures ...")
source(file.path(STAGE_DIR, "02_plot_pe_figures.R"), local = .GlobalEnv)

if (file.exists(project_path("03_Output/model_fits/bahia/v4_9_replication/outputs/q0.05/renewal_bahia_v4_9_fit_q0.05.rds"))) {
  message("\n[Stage 03] 03/06: cross_state_comparison (CE/BA/PE) ...")
  source(file.path(STAGE_DIR, "03_cross_state_comparison.R"), local = .GlobalEnv)
} else {
  message("\n[Stage 03] 03/06: SKIPPED -- Bahia fit not present (needed for the cross-state comparison).")
}

message("\n[Stage 03] 04/06: plot_pe_six_seven_panel (modelA_caseonly) ...")
source(file.path(STAGE_DIR, "04_plot_pe_six_seven_panel.R"), local = .GlobalEnv)
plot_pe_panel("modelA_caseonly", "pilot")

message("\n[Stage 03] 05/06: summarise_pe_fixed_q_sweep ...")
source(file.path(STAGE_DIR, "05_summarise_pe_fixed_q_sweep.R"), local = .GlobalEnv)

message("\n[Stage 03] 06/06: plot_pe_fixedq_panel (q = 0.05-0.30) ...")
source(file.path(STAGE_DIR, "06_plot_pe_fixedq_panel.R"), local = .GlobalEnv)
for (qv in c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)) plot_pe_fixedq_panel(qv, stage = "pilot")

message("\n[Stage 03] Complete.")
