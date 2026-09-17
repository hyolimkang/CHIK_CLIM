# ============================================================
# STAGE 04 -- SENSITIVITY / IDENTIFICATION (Pernambuco, orchestrator only)
# ============================================================
# Global-q and curved ridge-reparameterisation identification diagnostics
# (computationally improved coordinate systems for the frozen v4.9 model,
# not new epidemiological models -- see PE_V4_9_GLOBAL_Q_RIDGE_REPARAM.md
# and PE_CURVED_RIDGE_REPARAM_TEST.md). modelC (ridge) and modelD (curved)
# pilot/full fits are already on disk; the fit steps below reuse them and
# only refit when missing AND PERNAMBUCO_ALLOW_REFIT="TRUE".

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/04_pernambuco_pipeline/04_sensitivity_identification")
ALLOW_REFIT <- identical(Sys.getenv("PERNAMBUCO_ALLOW_REFIT", "FALSE"), "TRUE")
BASE_OUT <- project_path("03_Output/model_fits/pernambuco/v4_9_replication")

message("\n[Stage 04] 01/09: fit_pe_global_q_ridge_reparam (modelC) ...")
modelC_fit <- file.path(BASE_OUT, "outputs/modelC/pe_ridge_U14_pilot.rds")
source(file.path(STAGE_DIR, "01_fit_pe_global_q_ridge_reparam.R"), local = .GlobalEnv) # defines run_pe_v4_9_global_q_ridge_reparam(); does not auto-run under source()
if (file.exists(modelC_fit)) {
  message("  Found: ", modelC_fit, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  run_pe_v4_9_global_q_ridge_reparam()
} else {
  stop("Missing prerequisite:\n  ", modelC_fit, "\n\nSet Sys.setenv(PERNAMBUCO_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 04] 02/09: diagnose_pe_global_q_ridge_reparam ...")
source(file.path(STAGE_DIR, "02_diagnose_pe_global_q_ridge_reparam.R"), local = .GlobalEnv)

message("\n[Stage 04] 03/09: plot_pe_ridge_reparam_panel ...")
source(file.path(STAGE_DIR, "03_plot_pe_ridge_reparam_panel.R"), local = .GlobalEnv)

message("\n[Stage 04] 04/09: plot_pe_ridge_reparam_traceplots ...")
source(file.path(STAGE_DIR, "04_plot_pe_ridge_reparam_traceplots.R"), local = .GlobalEnv)

message("\n[Stage 04] 05/09: characterise_curved_ridge ...")
source(file.path(STAGE_DIR, "05_characterise_curved_ridge.R"), local = .GlobalEnv)

message("\n[Stage 04] 06/09: free_quadratic_ridge_gate ...")
source(file.path(STAGE_DIR, "06_free_quadratic_ridge_gate.R"), local = .GlobalEnv)

message("\n[Stage 04] 07/09: fit_pe_curved_ridge_reparam (modelD pilot) ...")
modelD_fit <- file.path(BASE_OUT, "outputs/modelD/pe_curved_U14_pilot.rds")
source(file.path(STAGE_DIR, "07_fit_pe_curved_ridge_reparam.R"), local = .GlobalEnv) # defines run_pe_v4_9_curved_ridge_reparam(); does not auto-run under source()
if (file.exists(modelD_fit)) {
  message("  Found: ", modelD_fit, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  run_pe_v4_9_curved_ridge_reparam()
} else {
  stop("Missing prerequisite:\n  ", modelD_fit, "\n\nSet Sys.setenv(PERNAMBUCO_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 04] 08/09: plot_curved_pilot_diagnostics ...")
source(file.path(STAGE_DIR, "08_plot_curved_pilot_diagnostics.R"), local = .GlobalEnv)

message("\n[Stage 04] 09/09: fit_pe_curved_longwarmup_diagnostic ...")
longwarmup_fit <- file.path(BASE_OUT, "outputs/modelD/pe_curved_U14_longwarmup.rds")
source(file.path(STAGE_DIR, "09_fit_pe_curved_longwarmup_diagnostic.R"), local = .GlobalEnv) # defines run_longwarmup_diagnostic(); does not auto-run under source()
if (file.exists(longwarmup_fit)) {
  message("  Found: ", longwarmup_fit, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  run_longwarmup_diagnostic()
} else {
  stop("Missing prerequisite:\n  ", longwarmup_fit, "\n\nSet Sys.setenv(PERNAMBUCO_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 04] Complete.")
