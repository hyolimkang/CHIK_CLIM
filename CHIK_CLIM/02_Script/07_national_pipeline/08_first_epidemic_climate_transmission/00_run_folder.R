# ============================================================
# 08_first_epidemic_climate_transmission (National pipeline, orchestrator only)
# ============================================================
# Brazil first-epidemic climate-transmission pilot (see
# FIRST_EPIDEMIC_CLIMATE_TRANSMISSION_REPORT.md). All 7 scripts execute
# top-to-bottom on source() (no guarded functions); each is existence-gated
# except the final plotting step (06), which always re-runs since it only
# reads already-built tables.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this runner.")
  source(setup_path, local = .GlobalEnv)
}
FOLDER_DIR <- project_path("02_Script/07_national_pipeline/08_first_epidemic_climate_transmission")
ALLOW_REFIT <- identical(Sys.getenv("NATIONAL_ALLOW_REFIT", "FALSE"), "TRUE")
TABLE_DIR <- project_path("03_Output/tables/national_wave_analysis")

run_gated <- function(step, fname, fp) {
  message(sprintf("\n[08_first_epidemic_climate_transmission] %s: %s ...", step, fname))
  if (file.exists(fp)) {
    message("  Found: ", fp, " -- reusing (not rebuilt).")
  } else if (ALLOW_REFIT) {
    source(file.path(FOLDER_DIR, fname), local = .GlobalEnv)
  } else {
    stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(NATIONAL_ALLOW_REFIT = \"TRUE\") to allow rebuilding, then re-source this runner.", call. = FALSE)
  }
}

run_gated("01/07", "01_build_first_episode_eligibility.R", file.path(TABLE_DIR, "FIRST_EPISODE_CLIMATE_ELIGIBILITY.csv"))
run_gated("02/07", "02_build_first_episode_weekly_re.R", file.path(TABLE_DIR, "FIRST_EPISODE_WEEKLY_RE.csv"))
run_gated("03/07", "03_build_first_episode_climate_covariates.R", file.path(TABLE_DIR, "first_epidemic_climate_training_data.csv"))
run_gated("04/07", "04_fit_first_episode_climate_gam.R", file.path(TABLE_DIR, "first_epidemic_climate_gam.rds"))
run_gated("05/07", "05_fit_first_episode_climate_gam_mc.R", file.path(TABLE_DIR, "first_epidemic_climate_gam_mc_summary.csv"))

message("\n[08_first_epidemic_climate_transmission] 06/07: plot_first_episode_climate_figures ...")
source(file.path(FOLDER_DIR, "06_plot_first_episode_climate_figures.R"), local = .GlobalEnv)

run_gated("07/07", "07_diagnose_first_episode_climate_model.R", file.path(TABLE_DIR, "first_epidemic_climate_diagnostics_overview.csv"))

message("\n[08_first_epidemic_climate_transmission] Complete.")
