# ============================================================
# 09_recurrent_climate_susceptibility_phase (National pipeline, orchestrator only)
# ============================================================
# Recurrent-epidemic climate x susceptibility phase analysis for BA/RJ/MT
# (the primary recurrent-analysis states). All 4 scripts execute
# top-to-bottom on source() (no guarded functions); each is existence-gated.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this runner.")
  source(setup_path, local = .GlobalEnv)
}
FOLDER_DIR <- project_path("02_Script/07_national_pipeline/09_recurrent_climate_susceptibility_phase")
ALLOW_REFIT <- identical(Sys.getenv("NATIONAL_ALLOW_REFIT", "FALSE"), "TRUE")
TABLE_DIR <- project_path("03_Output/tables/national_wave_analysis")

run_gated <- function(step, fname, fp) {
  message(sprintf("\n[09_recurrent_climate_susceptibility_phase] %s: %s ...", step, fname))
  if (file.exists(fp)) {
    message("  Found: ", fp, " -- reusing (not rebuilt).")
  } else if (ALLOW_REFIT) {
    source(file.path(FOLDER_DIR, fname), local = .GlobalEnv)
  } else {
    stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(NATIONAL_ALLOW_REFIT = \"TRUE\") to allow rebuilding, then re-source this runner.", call. = FALSE)
  }
}

run_gated("01/04", "01_build_recurrent_climate_susceptibility_episodes.R", file.path(TABLE_DIR, "RECURRENT_CLIMATE_SUSCEPTIBILITY_EPISODES.csv"))
run_gated("02/04", "02_plot_recurrent_climate_susceptibility_phase.R", file.path(TABLE_DIR, "RECURRENT_RE_VALIDATION_METRICS.csv"))
run_gated("03/04", "03_build_recurrent_weekly_opportunity_panel.R", file.path(TABLE_DIR, "RECURRENT_CLIMATE_SUSCEPTIBILITY_OPPORTUNITIES_WEEKLY.csv"))
run_gated("04/04", "04_plot_recurrent_climate_susceptibility_opportunities.R", file.path(TABLE_DIR, "RECURRENT_ALL_OPPORTUNITIES_OUTCOME_SUMMARY.csv"))

message("\n[09_recurrent_climate_susceptibility_phase] Complete.")
