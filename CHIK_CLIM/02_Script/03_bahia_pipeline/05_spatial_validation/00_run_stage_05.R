# ============================================================
# STAGE 05 -- SPATIAL VALIDATION (Bahia, orchestrator only)
# ============================================================
# Municipality-level spatial diagnostics: turnover diagnostic (01-06,
# formerly 31_bahia_spatial_turnover_diagnostic/) then global-q local
# consistency audit (07-13, formerly 32_bahia_global_q_local_consistency_audit/).
# All descriptive/diagnostic -- does NOT fit or refit any Stan model.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/03_bahia_pipeline/05_spatial_validation")
primary_fit <- project_path("03_Output/03_bahia_pipeline/model_fits/v4_9_replication/outputs/q0.05/renewal_bahia_v4_9_fit_q0.05.rds")
panel_path <- project_path("01_Data/chik_dlnm_panel_muni_week_2015_2025.rds")
missing <- c(primary_fit, panel_path)[!file.exists(c(primary_fit, panel_path))]
if (length(missing) > 0) {
  stop("Missing prerequisite(s):\n  ", paste(missing, collapse = "\n  "),
       "\n\nIf the fit is missing, run Stage 02 first:\n  02_transmission_fitting (00_run_stage_02.R)", call. = FALSE)
}
message("[Stage 05] Prerequisites OK.")

message("\n--- Turnover diagnostic (01-06) ---")
for (i in 1:6) {
  f <- list.files(STAGE_DIR, pattern = sprintf("^0%d_", i), full.names = TRUE)
  message(sprintf("[Stage 05] %02d/13: %s ...", i, basename(f)))
  source(f, local = .GlobalEnv)
}

message("\n--- Global-q local consistency audit (07-13) ---")
for (i in 7:13) {
  f <- list.files(STAGE_DIR, pattern = sprintf("^%02d_", i), full.names = TRUE)
  message(sprintf("[Stage 05] %02d/13: %s ...", i, basename(f)))
  source(f, local = .GlobalEnv)
}

message("\n[Stage 05] Complete.")
