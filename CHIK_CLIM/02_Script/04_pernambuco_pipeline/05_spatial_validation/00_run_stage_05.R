# ============================================================
# STAGE 05 -- SPATIAL VALIDATION (Pernambuco, orchestrator only)
# ============================================================
# Municipality-level turnover diagnostic (01-08), 5-strata / Spatial-1
# hierarchical-region model attempt (09-16), Spatial-0 common-q
# model-criticism baseline (17-20), and two deterministic feasibility
# stops (21-22). All descriptive/diagnostic.
#
# HALTED WORK, KEPT LIVE HERE (per explicit decision -- not archived):
#   - 12-15 (5-strata / "Spatial-1" fit and its downstream HMC/structural/PPC
#     diagnostics): chain 1 of the 4-chain pilot hit a severe, unrecovered
#     warmup pathology (78/1500 iterations; see
#     PE_SPATIAL1_GEOMETRY_DIAGNOSTIC.md) and the run was stopped. No
#     pe_5strata_pilot.rds was ever produced. Step 16 is the diagnostic that
#     explains this using the raw (partial) chain CSVs, which DO exist.
#     Steps 13-15 need the bundled fit and are skipped below.
#   - 21 (Q2 feasibility) and 22 (regional-R feasibility): both are
#     deterministic, no-Stan-refit feasibility stops (see
#     PE_SPATIAL0_Q2_FEASIBILITY.md / PE_REGIONAL_R_FEASIBILITY.md) that
#     concluded further HMC fitting was not warranted. They always run --
#     the stop IS their documented result, not a missing prerequisite.
#   - Spatial-0 (17-20) is a completed "model-criticism baseline" whose
#     accepted artifact is its raw per-chain CSVs (never bundled into an
#     .rds by design -- steps 19/20 read outputs/spatial0/chains/ directly).

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/04_pernambuco_pipeline/05_spatial_validation")
ALLOW_REFIT <- identical(Sys.getenv("PERNAMBUCO_ALLOW_REFIT", "FALSE"), "TRUE")
BASE_OUT <- project_path("03_Output/model_fits/pernambuco/v4_9_replication")
panel_path <- project_path("01_Data/chik_dlnm_panel_muni_week_2015_2025.rds")
modelB_fit <- file.path(BASE_OUT, "outputs/modelB_U14/renewal_pe_global_q_fit_modelB_U14.rds")
missing <- c(panel_path, modelB_fit)[!file.exists(c(panel_path, modelB_fit))]
if (length(missing) > 0) {
  stop("Missing prerequisite(s):\n  ", paste(missing, collapse = "\n  "),
       "\n\nIf the fit is missing, run Stage 02 first:\n  02_transmission_fitting (00_run_stage_02.R)", call. = FALSE)
}
message("[Stage 05] Prerequisites OK.")

message("\n--- Municipality-level turnover diagnostic (01-08) ---")
for (i in 1:8) {
  f <- list.files(STAGE_DIR, pattern = sprintf("^%02d_", i), full.names = TRUE)
  message(sprintf("[Stage 05] %02d/22: %s ...", i, basename(f)))
  source(f, local = .GlobalEnv)
}

message("\n--- 5-strata / Spatial-1 hierarchical-region attempt (09-16) ---")
for (i in 9:11) {
  f <- list.files(STAGE_DIR, pattern = sprintf("^%02d_", i), full.names = TRUE)
  message(sprintf("[Stage 05] %02d/22: %s ...", i, basename(f)))
  source(f, local = .GlobalEnv) # 11 defines build_pe_5strata_stan_data(); no auto-run under source()
}

message(sprintf("[Stage 05] %02d/22: %s ...", 12, "12_fit_pe_5strata_spatial.R"))
spatial5_fit <- file.path(BASE_OUT, "outputs/spatial5/pe_5strata_pilot.rds")
source(file.path(STAGE_DIR, "12_fit_pe_5strata_spatial.R"), local = .GlobalEnv) # defines run_pilot(); no auto-run under source()
if (file.exists(spatial5_fit)) {
  message("  Found: ", spatial5_fit, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  message("  PERNAMBUCO_ALLOW_REFIT=TRUE: attempting the 5-strata pilot fit.",
          " NOTE: this run previously stopped with a chain-1 warmup pathology (see PE_SPATIAL1_GEOMETRY_DIAGNOSTIC.md) -- a repeat attempt is a deliberate, explicit choice, not a routine refit.")
  run_pilot()
} else {
  message("  SKIPPED -- ", spatial5_fit, " does not exist.",
          " This fit was previously stopped (chain 1 warmup pathology, 78/1500 iterations; see PE_SPATIAL1_GEOMETRY_DIAGNOSTIC.md) and never completed.",
          " Steps 13-15 (which require this fit) are skipped below. Set Sys.setenv(PERNAMBUCO_ALLOW_REFIT = \"TRUE\") to deliberately re-attempt it.")
}

if (file.exists(spatial5_fit)) {
  for (i in 13:15) {
    f <- list.files(STAGE_DIR, pattern = sprintf("^%02d_", i), full.names = TRUE)
    message(sprintf("[Stage 05] %02d/22: %s ...", i, basename(f)))
    source(f, local = .GlobalEnv)
  }
} else {
  message("\n[Stage 05] 13-15/22: SKIPPED -- pe_5strata_pilot.rds not present (see note above).")
}

message(sprintf("[Stage 05] %02d/22: %s ...", 16, "16_pe_spatial1_geometry_diagnostic.R"))
source(file.path(STAGE_DIR, "16_pe_spatial1_geometry_diagnostic.R"), local = .GlobalEnv)

message("\n--- Spatial-0 common-q model-criticism baseline (17-20) ---")
for (i in c(17, 18)) {
  f <- list.files(STAGE_DIR, pattern = sprintf("^%02d_", i), full.names = TRUE)
  message(sprintf("[Stage 05] %02d/22: %s ...", i, basename(f)))
  source(f, local = .GlobalEnv) # 17 defines build_pe_spatial0_stan_data(); 18 defines run_spatial0_pilot() -- neither auto-runs under source(); existing raw chains under outputs/spatial0/chains/ are this model's accepted final artifact (no bundled .rds is produced by design)
}
for (i in 19:20) {
  f <- list.files(STAGE_DIR, pattern = sprintf("^%02d_", i), full.names = TRUE)
  message(sprintf("[Stage 05] %02d/22: %s ...", i, basename(f)))
  source(f, local = .GlobalEnv)
}

message("\n--- Deterministic feasibility stops (21-22) ---")
for (i in 21:22) {
  f <- list.files(STAGE_DIR, pattern = sprintf("^%02d_", i), full.names = TRUE)
  message(sprintf("[Stage 05] %02d/22: %s ...", i, basename(f)))
  source(f, local = .GlobalEnv)
}

message("\n[Stage 05] Complete.")
