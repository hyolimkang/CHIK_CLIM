# ============================================================
# 01_episode_renewal (National pipeline, orchestrator only)
# ============================================================
# Ceara early-phase episode-level Re fit + diagnostics. Reuses the fit if it
# already exists; only refits when missing AND NATIONAL_ALLOW_REFIT="TRUE".

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this runner.")
  source(setup_path, local = .GlobalEnv)
}
FOLDER_DIR <- project_path("02_Script/07_national_pipeline/01_episode_renewal")
ALLOW_REFIT <- identical(Sys.getenv("NATIONAL_ALLOW_REFIT", "FALSE"), "TRUE")
fp <- project_path("03_Output/07_national_pipeline/model_fits/episode_renewal/renewal_ceara_episode_re_fit.rds")

message("\n[01_episode_renewal] 01/02: fit_ceara_episode_re ...")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(FOLDER_DIR, "01_fit_ceara_episode_re.R"), local = .GlobalEnv) # defines run_episode_renewal_fit(); does not auto-run under source()
  run_episode_renewal_fit()
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(NATIONAL_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this runner.", call. = FALSE)
}

message("\n[01_episode_renewal] 02/02: plot_ceara_episode_re ...")
source(file.path(FOLDER_DIR, "02_plot_ceara_episode_re.R"), local = .GlobalEnv) # defines run_episode_renewal_diagnostics(); does not auto-run under source()
run_episode_renewal_diagnostics()

message("\n[01_episode_renewal] Complete.")
