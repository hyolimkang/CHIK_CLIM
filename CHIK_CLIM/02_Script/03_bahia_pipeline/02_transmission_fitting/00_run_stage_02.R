# ============================================================
# STAGE 02 -- TRANSMISSION FITTING (Bahia, orchestrator only)
# ============================================================
#
# Can refit Stan (case-only v4.9, multisite serology, global-q multisite).
# Reuses each fit if it already exists; only (re)fits when missing AND
# BAHIA_ALLOW_REFIT="TRUE" -- a full refit is never triggered just because
# this file is sourced. Primary reference fit is q=0.05 (matching every
# other state's convention); the q-sweep (0.10-0.30) and the global-q
# on/off variants are sensitivity runs, not blocking prerequisites for
# downstream stages.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/03_bahia_pipeline/02_transmission_fitting")
ALLOW_REFIT <- identical(Sys.getenv("BAHIA_ALLOW_REFIT", "FALSE"), "TRUE")
BASE_OUT <- project_path("03_Output/03_bahia_pipeline/model_fits/v4_9_replication")

message("\n[Stage 02] 01/03: fit_bahia_v4_9 (case-only, q=0.05 primary) ...")
primary_fit <- file.path(BASE_OUT, "outputs/q0.05/renewal_bahia_v4_9_fit_q0.05.rds")
if (file.exists(primary_fit)) {
  message("  Found: ", primary_fit, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(STAGE_DIR, "01_fit_bahia_v4_9.R"), local = .GlobalEnv) # defines run_fit_bahia_v4_9(); does not auto-run under source()
  run_fit_bahia_v4_9()
} else {
  stop("Missing prerequisite:\n  ", primary_fit, "\n\nSet Sys.setenv(BAHIA_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] 02/03: fit_bahia_multisite_serology (q=0.05) ...")
multisite_fit <- file.path(BASE_OUT, "outputs_multisite_serology/q0.05/renewal_bahia_multisite_fit_q0.05.rds")
if (file.exists(multisite_fit)) {
  message("  Found: ", multisite_fit, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(STAGE_DIR, "02_build_bahia_serology_windows.R"), local = .GlobalEnv) # helper, no auto-run under source()
  source(file.path(STAGE_DIR, "03_fit_bahia_multisite_serology.R"), local = .GlobalEnv) # defines run_fit_bahia_multisite()
  run_fit_bahia_multisite()
} else {
  stop("Missing prerequisite:\n  ", multisite_fit, "\n\nSet Sys.setenv(BAHIA_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] 03/03: fit_bahia_global_q_multisite ...")
globalq_fit <- file.path(BASE_OUT, "outputs_global_q_multisite_serology/on_full2/renewal_bahia_global_q_fit_on_full2.rds")
if (file.exists(globalq_fit)) {
  message("  Found: ", globalq_fit, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  if (!exists("build_bahia_sero_stan_fields")) source(file.path(STAGE_DIR, "02_build_bahia_serology_windows.R"), local = .GlobalEnv)
  source(file.path(STAGE_DIR, "04_fit_bahia_global_q_multisite.R"), local = .GlobalEnv) # defines run_fit_bahia_global_q()
  run_fit_bahia_global_q()
} else {
  stop("Missing prerequisite:\n  ", globalq_fit, "\n\nSet Sys.setenv(BAHIA_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] Complete. (q-sweep 0.10-0.30 and the global-q 'off' variant are optional sensitivity runs, not checked here -- run the scripts directly with QVAL/env vars set if needed.)")
