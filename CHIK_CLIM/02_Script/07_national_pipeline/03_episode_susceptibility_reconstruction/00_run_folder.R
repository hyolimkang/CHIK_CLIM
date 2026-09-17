# ============================================================
# 03_episode_susceptibility_reconstruction (National pipeline, orchestrator only)
# ============================================================
# Ceara episode-sequential susceptibility pilot: prepare periods, fit, and
# diagnose. Reuses the fit if it already exists; only refits when missing
# AND NATIONAL_ALLOW_REFIT="TRUE".
#
# PRE-EXISTING, NOT-TO-BE-FIXED FINDING: the HMC pass/fail gate in
# 03_diagnose_ceara_episode_susceptibility.R reports that only 3 of 7
# episodes pass the standard HMC gate (divergences/treedepth/Rhat/ESS/BFMI).
# This is a property of the frozen pilot-v1 fit itself, not something this
# refactor introduces or should silently fix -- it is reported as-is.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this runner.")
  source(setup_path, local = .GlobalEnv)
}
FOLDER_DIR <- project_path("02_Script/07_national_pipeline/03_episode_susceptibility_reconstruction")
ALLOW_REFIT <- identical(Sys.getenv("NATIONAL_ALLOW_REFIT", "FALSE"), "TRUE")
fp <- project_path("02_Script/stan/ceara_episode_susceptibility_pilot_v1_fit.rds")

message("\n[03_episode_susceptibility_reconstruction] 01/02: fit_ceara_episode_susceptibility ...")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(FOLDER_DIR, "02_fit_ceara_episode_susceptibility.R"), local = .GlobalEnv) # defines run_ceara_episode_susceptibility_fit(); internally sources 01_prepare_ceara_episode_periods.R; does not auto-run under source()
  run_ceara_episode_susceptibility_fit()
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(NATIONAL_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this runner.", call. = FALSE)
}

message("\n[03_episode_susceptibility_reconstruction] 02/02: diagnose_ceara_episode_susceptibility ...")
source(file.path(FOLDER_DIR, "03_diagnose_ceara_episode_susceptibility.R"), local = .GlobalEnv) # defines run_ceara_episode_susceptibility_diagnostics(); does not auto-run under source()
run_ceara_episode_susceptibility_diagnostics()

message("\n[03_episode_susceptibility_reconstruction] Complete. (Expect the HMC gate to report 3/7 episodes passing -- pre-existing, not fixed here.)")
