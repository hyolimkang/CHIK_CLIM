# ============================================================
# 04_national_episode_census (National pipeline, orchestrator only)
# ============================================================
# Builds the frozen Brazil-wide chikungunya episode census. Reuses the
# existing census if already built; only rebuilds when missing AND
# NATIONAL_ALLOW_REFIT="TRUE" (this is a deterministic table build, not a
# Stan fit, but the same explicit-opt-in gate is used for consistency and
# because downstream wave/episode work is frozen against this exact table).

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this runner.")
  source(setup_path, local = .GlobalEnv)
}
FOLDER_DIR <- project_path("02_Script/07_national_pipeline/04_national_episode_census")
ALLOW_REFIT <- identical(Sys.getenv("NATIONAL_ALLOW_REFIT", "FALSE"), "TRUE")
fp <- project_path("03_Output/07_national_pipeline/tables/national_episode_census/brazil_chik_episode_census.csv")

message("\n[04_national_episode_census] build_brazil_chik_episode_census ...")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not rebuilt).")
} else if (ALLOW_REFIT) {
  source(file.path(FOLDER_DIR, "01_build_brazil_chik_episode_census.R"), local = .GlobalEnv) # defines build_brazil_chik_episode_census(); does not auto-run under source()
  build_brazil_chik_episode_census()
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(NATIONAL_ALLOW_REFIT = \"TRUE\") to allow rebuilding, then re-source this runner.", call. = FALSE)
}

message("\n[04_national_episode_census] Complete.")
