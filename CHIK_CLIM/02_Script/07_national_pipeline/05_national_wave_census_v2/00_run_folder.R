# ============================================================
# 05_national_wave_census_v2 (National pipeline, orchestrator only)
# ============================================================
# Builds the frozen Brazil-wide chikungunya wave census (v2). Reuses the
# existing census if already built; only rebuilds when missing AND
# NATIONAL_ALLOW_REFIT="TRUE". Every downstream state/national script that
# reads major_epidemic_primary/wave_id is frozen against this exact table.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this runner.")
  source(setup_path, local = .GlobalEnv)
}
FOLDER_DIR <- project_path("02_Script/07_national_pipeline/05_national_wave_census_v2")
ALLOW_REFIT <- identical(Sys.getenv("NATIONAL_ALLOW_REFIT", "FALSE"), "TRUE")
fp <- project_path("03_Output/tables/national_wave_census_v2/brazil_chik_wave_census_v2.csv")

message("\n[05_national_wave_census_v2] build_brazil_chik_wave_census_v2 ...")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not rebuilt).")
} else if (ALLOW_REFIT) {
  source(file.path(FOLDER_DIR, "01_build_brazil_chik_wave_census_v2.R"), local = .GlobalEnv) # defines build_brazil_chik_wave_census_v2(); does not auto-run under source()
  build_brazil_chik_wave_census_v2()
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(NATIONAL_ALLOW_REFIT = \"TRUE\") to allow rebuilding, then re-source this runner.", call. = FALSE)
}

message("\n[05_national_wave_census_v2] Complete.")
