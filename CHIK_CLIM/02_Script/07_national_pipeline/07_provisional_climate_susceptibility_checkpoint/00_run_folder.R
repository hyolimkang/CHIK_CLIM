# ============================================================
# 07_provisional_climate_susceptibility_checkpoint (National pipeline, orchestrator only)
# ============================================================
# Exploratory checkpoint (see the original task's "Provisional episode-level
# Re x climate x susceptibility integration checkpoint" plan): does the
# current FOI-informed susceptibility add explanatory power for early
# epidemic growth (Re) beyond climate alone? Additive/exploratory only --
# does not alter wave definitions, early-Re, long-term FOI, or the
# susceptibility model. Each step is existence-gated (reused unless missing
# AND NATIONAL_ALLOW_REFIT="TRUE") except the final plotting step.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this runner.")
  source(setup_path, local = .GlobalEnv)
}
FOLDER_DIR <- project_path("02_Script/07_national_pipeline/07_provisional_climate_susceptibility_checkpoint")
ALLOW_REFIT <- identical(Sys.getenv("NATIONAL_ALLOW_REFIT", "FALSE"), "TRUE")
TABLE_DIR <- project_path("03_Output/tables/national_wave_analysis")

run_gated <- function(step, fname, fp, fun_name) {
  message(sprintf("\n[07_provisional_climate_susceptibility_checkpoint] %s ...", fname))
  if (file.exists(fp)) {
    message("  Found: ", fp, " -- reusing (not rebuilt).")
  } else if (ALLOW_REFIT) {
    source(file.path(FOLDER_DIR, fname), local = .GlobalEnv)
    get(fun_name, envir = .GlobalEnv)()
  } else {
    stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(NATIONAL_ALLOW_REFIT = \"TRUE\") to allow rebuilding, then re-source this runner.", call. = FALSE)
  }
}

run_gated("01/06", "01_build_brazil_uf_weekly_climate.R", file.path(TABLE_DIR, "brazil_chik_uf_weekly_climate.csv"), "run_build_brazil_uf_weekly_climate")
run_gated("02/06", "02_extract_wave_onset_susceptibility_draws.R", file.path(TABLE_DIR, "brazil_chik_wave_onset_susceptibility_draws_validation.csv"), "run_extract_wave_onset_susceptibility_draws")
run_gated("03/06", "03_build_episode_climate_dataset.R", file.path(TABLE_DIR, "brazil_chik_Re_S_climate_episode_analysis_scale1_0.csv"), "run_build_episode_climate_dataset")
run_gated("04/06", "04_fit_provisional_gam_models.R", file.path(TABLE_DIR, "brazil_chik_provisional_gam_comparison.csv"), "run_provisional_gam_comparison")
run_gated("05/06", "05_fit_provisional_gam_monte_carlo.R", file.path(TABLE_DIR, "brazil_chik_provisional_gam_monte_carlo_summary.csv"), "run_provisional_gam_monte_carlo")

message("\n[07_provisional_climate_susceptibility_checkpoint] 06/06: plot_provisional_episode_figures ...")
source(file.path(FOLDER_DIR, "06_plot_provisional_episode_figures.R"), local = .GlobalEnv) # defines run_plot_provisional_episode_figures(); does not auto-run under source()
run_plot_provisional_episode_figures()

message("\n[07_provisional_climate_susceptibility_checkpoint] Complete.")
