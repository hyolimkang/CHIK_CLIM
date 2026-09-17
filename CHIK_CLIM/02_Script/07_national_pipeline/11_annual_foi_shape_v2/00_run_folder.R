# ============================================================
# 11_annual_foi_shape_v2 (National pipeline, orchestrator only)
# ============================================================
# State-specific batch extension of the existing chik_dynamic_annual_foi_shape_v2.stan
# model (SHAPE + NB2 sensitivity variants). Prepare (01) and fit (02) are
# existence-gated; diagnose (03) and the Ceara episode-driver overlay (04)
# always re-run (read-only against the fits).

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this runner.")
  source(setup_path, local = .GlobalEnv)
}
FOLDER_DIR <- project_path("02_Script/07_national_pipeline/11_annual_foi_shape_v2")
ALLOW_REFIT <- identical(Sys.getenv("NATIONAL_ALLOW_REFIT", "FALSE"), "TRUE")

message("\n[11_annual_foi_shape_v2] 01/04: prepare_annual_foi_shape_v2 ...")
fp <- project_path("03_Output/tables/annual_foi_shape_v2/annual_foi_shape_v2_input.rds")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not rebuilt).")
} else if (ALLOW_REFIT) {
  source(file.path(FOLDER_DIR, "01_prepare_annual_foi_shape_v2.R"), local = .GlobalEnv) # top-level executing script, no guard
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(NATIONAL_ALLOW_REFIT = \"TRUE\") to allow rebuilding, then re-source this runner.", call. = FALSE)
}

message("\n[11_annual_foi_shape_v2] 02/04: fit_annual_foi_shape_v2 (6 SHAPE/NB2 sensitivity variants) ...")
fp <- project_path("03_Output/tables/annual_foi_shape_v2/annual_foi_shape_v2_fit_manifest.csv")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(FOLDER_DIR, "02_fit_annual_foi_shape_v2.R"), local = .GlobalEnv) # top-level executing script, no guard -- fits all 6 variants
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(NATIONAL_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this runner.", call. = FALSE)
}

message("\n[11_annual_foi_shape_v2] 03/04: diagnose_annual_foi_shape_v2 ...")
source(file.path(FOLDER_DIR, "03_diagnose_annual_foi_shape_v2.R"), local = .GlobalEnv) # defines diagnose_annual_foi_shape_v2(); does not auto-run under source()
diagnose_annual_foi_shape_v2()

message("\n[11_annual_foi_shape_v2] 04/04: build_ceara_episode_driver_overlay ...")
source(file.path(FOLDER_DIR, "04_build_ceara_episode_driver_overlay.R"), local = .GlobalEnv) # defines build_ceara_episode_driver_overlay(); does not auto-run under source()
build_ceara_episode_driver_overlay()

message("\n[11_annual_foi_shape_v2] Complete.")
