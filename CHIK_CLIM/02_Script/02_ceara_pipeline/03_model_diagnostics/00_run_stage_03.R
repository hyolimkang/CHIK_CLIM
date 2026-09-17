# ============================================================
# STAGE 03 -- MODEL DIAGNOSTICS (orchestrator only)
# ============================================================
#
# Read-only: validates the two accepted fits, then plots the historical
# reconstruction. Does NOT fit anything. Run diagnostics FIRST, then the
# reconstruction plot, per the task specification.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}

STAGE_DIR <- project_path("02_Script/02_ceara_pipeline/03_model_diagnostics")

v49_baseline_path <- project_path("03_Output/02_ceara_pipeline/model_fits/baseline/v4_9_hierarchical_seasonality/outputs/q0.05/renewal_ceara_v4_9_fit_q0.05.rds")
climate_fit_path <- project_path("03_Output/02_ceara_pipeline/model_fits/climate_forced/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds")
missing <- c(v49_baseline_path, climate_fit_path)[!file.exists(c(v49_baseline_path, climate_fit_path))]
if (length(missing) > 0) {
  stop("Missing prerequisite(s):\n  ", paste(missing, collapse = "\n  "),
       "\n\nRun Stage 02 first:\n  02_transmission_fitting (00_run_stage_02.R)", call. = FALSE)
}
message("[Stage 03] Prerequisites OK: baseline v4.9 and climate-forced v4.9 fits both found.")

message("\n[Stage 03] 01/02: validate_climate_v49_fit ...")
source(file.path(STAGE_DIR, "01_validate_climate_v49_fit.R"), local = .GlobalEnv)

message("\n[Stage 03] 02/02: plot_historical_reconstruction ...")
source(file.path(STAGE_DIR, "02_plot_historical_reconstruction.R"), local = .GlobalEnv) # defines plot_ce_climate_canary_six_panel(); does not auto-run under source()
plot_ce_climate_canary_six_panel()

message("\n[Stage 03] Complete.")
