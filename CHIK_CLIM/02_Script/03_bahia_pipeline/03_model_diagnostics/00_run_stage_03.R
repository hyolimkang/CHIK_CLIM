# ============================================================
# STAGE 03 -- MODEL DIAGNOSTICS (Bahia, orchestrator only)
# ============================================================
# Read-only: validates/plots the accepted fits. Does not fit anything.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/03_bahia_pipeline/03_model_diagnostics")
BASE_OUT <- project_path("03_Output/model_fits/bahia/v4_9_replication")
primary_fit <- file.path(BASE_OUT, "outputs/q0.05/renewal_bahia_v4_9_fit_q0.05.rds")
if (!file.exists(primary_fit)) {
  stop("Missing prerequisite:\n  ", primary_fit, "\n\nRun Stage 02 first:\n  02_transmission_fitting (00_run_stage_02.R)", call. = FALSE)
}
message("[Stage 03] Prerequisites OK.")

message("\n[Stage 03] 01/08: plot_bahia_six_panel (q=0.05) ...")
source(file.path(STAGE_DIR, "01_plot_bahia_six_panel.R"), local = .GlobalEnv)
plot_bahia_six_panel(0.05, "pilot")

message("\n[Stage 03] 02/08: bahia_wave_ppc_and_figures ...")
source(file.path(STAGE_DIR, "02_bahia_wave_ppc_and_figures.R"), local = .GlobalEnv)

message("\n[Stage 03] 03/08: compare_ceara_bahia_posteriors ...")
source(file.path(STAGE_DIR, "03_compare_ceara_bahia_posteriors.R"), local = .GlobalEnv)

if (file.exists(file.path(BASE_OUT, "outputs_multisite_serology/q0.05/renewal_bahia_multisite_fit_q0.05.rds"))) {
  message("\n[Stage 03] 04/08: plot_bahia_multisite_six_panel (q=0.05) ...")
  source(file.path(STAGE_DIR, "04_plot_bahia_multisite_six_panel.R"), local = .GlobalEnv)
  plot_bahia_multisite_six_panel(0.05, "pilot")

  message("\n[Stage 03] 05/08: bahia_multisite_serology_diagnostics ...")
  source(file.path(STAGE_DIR, "05_bahia_multisite_serology_diagnostics.R"), local = .GlobalEnv)

  message("\n[Stage 03] 06/08: step0_audit_multisite_serology_wiring ...")
  source(file.path(STAGE_DIR, "06_step0_audit_multisite_serology_wiring.R"), local = .GlobalEnv)
} else {
  message("\n[Stage 03] 04-06/08: SKIPPED -- multisite serology fit not present (optional sensitivity track).")
}

if (file.exists(file.path(BASE_OUT, "outputs_global_q_multisite_serology/on_full2/renewal_bahia_global_q_fit_on_full2.rds"))) {
  message("\n[Stage 03] 07/08: plot_bahia_global_q_seven_panel ...")
  source(file.path(STAGE_DIR, "07_plot_bahia_global_q_seven_panel.R"), local = .GlobalEnv)
  plot_bahia_global_q_seven_panel("on_full2", "pilot") # "on_full2" matches the actual saved fit tag (outputs_global_q_multisite_serology/on_full2/); the script's own QTAG default ("pilot") does not correspond to any saved fit
} else {
  message("\n[Stage 03] 07/08: SKIPPED -- global-q multisite fit not present (optional sensitivity track).")
}

message("\n[Stage 03] 08/08: plot_sero_ablation_comparison ...")
source(file.path(STAGE_DIR, "08_plot_sero_ablation_comparison.R"), local = .GlobalEnv)

message("\n[Stage 03] Complete.")
