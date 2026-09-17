# ============================================================
# STAGE 02 -- HISTORICAL TRANSMISSION FITTING (orchestrator only)
# ============================================================
#
# THE ONLY current-pipeline stage that can refit Stan. Refits are
# expensive and are NEVER triggered automatically just because this
# runner is sourced: each of the three required fits is checked for
# existence first, and a fit is only (re)run if it is missing AND the
# environment variable CEARA_ALLOW_REFIT is "TRUE". This is deliberate --
# see Section 15/16 of the refactor task: "a full refit must always be an
# explicit user decision."
#
# To force a refit of something that already exists, delete (or rename)
# the corresponding .rds file first, set:
#   Sys.setenv(CEARA_ALLOW_REFIT = "TRUE")
# then re-source this runner.
#
# No scientific change: sources the exact frozen fitting scripts/functions
# unchanged (17_v4_0, 28_v4_9_hierarchical_seasonality) and the moved (but
# logically identical) 01_fit_ceara_climate_v49.R.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}

PIPELINE_DIR <- project_path("02_Script/02_ceara_pipeline")
STAGE_DIR <- file.path(PIPELINE_DIR, "02_transmission_fitting")
ALLOW_REFIT <- identical(Sys.getenv("CEARA_ALLOW_REFIT", "FALSE"), "TRUE")

# ---- Stage 01 prerequisite ----
climate_covariates_path <- project_path("03_Output/02_ceara_pipeline/tables/climate_forced_v4_9/climate_anomaly_covariates.csv")
if (!file.exists(climate_covariates_path)) {
  stop("Missing prerequisite:\n  ", climate_covariates_path,
       "\n\nRun Stage 01 first:\n  01_input_preparation (00_run_stage_01.R)", call. = FALSE)
}
message("[Stage 02] Prerequisite OK: Stage 01 climate covariates found.")

# ---- STEP 02A: td14 sampler-initialisation posterior ----
message("\nSTEP 02A: checking td14 initialisation posterior")
td14_path <- project_path("03_Output/02_ceara_pipeline/model_fits/initialization/renewal_ceara_v4_0_minimal_no_vaccine_fit_td14.rds")
if (file.exists(td14_path)) {
  message("  Found: ", td14_path, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  message("  NOT found -- CEARA_ALLOW_REFIT=TRUE, generating with the frozen v4.0 fitter's existing td14 settings.")
  source(project_path("02_Script/00_shared/legacy_model_functions/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R"), local = .GlobalEnv) # defines run_v4_0_minimal_no_vaccine(); does NOT auto-run (sys.nframe() != 0 under source())
  Sys.setenv(RENEWAL_V4_RUN_TAG = "td14", RENEWAL_V4_MAX_TREEDEPTH = "14")
  run_v4_0_minimal_no_vaccine()
  Sys.unsetenv(c("RENEWAL_V4_RUN_TAG", "RENEWAL_V4_MAX_TREEDEPTH"))
  if (!file.exists(td14_path)) stop("td14 fit did not produce the expected file: ", td14_path, call. = FALSE)
} else {
  stop("Missing prerequisite:\n  ", td14_path,
       "\n\nThis is expensive to regenerate and is not done automatically.\n",
       "To allow it: Sys.setenv(CEARA_ALLOW_REFIT = \"TRUE\") then re-source this stage.", call. = FALSE)
}

# ---- STEP 02B: baseline (non-climate) v4.9 q=0.05 ----
message("\nSTEP 02B: fitting baseline v4.9 q=0.05")
v49_baseline_path <- project_path("03_Output/02_ceara_pipeline/model_fits/baseline/v4_9_hierarchical_seasonality/outputs/q0.05/renewal_ceara_v4_9_fit_q0.05.rds")
if (file.exists(v49_baseline_path)) {
  message("  Found: ", v49_baseline_path, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  message("  NOT found -- CEARA_ALLOW_REFIT=TRUE, fitting now (expensive).")
  source(project_path("02_Script/00_shared/legacy_model_functions/28_v4_9_hierarchical_seasonality/scripts/01_fit_v4_9.R"), local = .GlobalEnv) # defines run_fit_v4_9(); does NOT auto-run
  run_fit_v4_9(qval = 0.05)
  if (!file.exists(v49_baseline_path)) stop("v4.9 baseline fit did not produce the expected file: ", v49_baseline_path, call. = FALSE)
} else {
  stop("Missing prerequisite:\n  ", v49_baseline_path,
       "\n\nThis is expensive to regenerate and is not done automatically.\n",
       "To allow it: Sys.setenv(CEARA_ALLOW_REFIT = \"TRUE\") then re-source this stage.", call. = FALSE)
}

# ---- STEP 02C: climate-forced v4.9 q=0.05 (the accepted canary) ----
message("\nSTEP 02C: fitting climate-forced v4.9 q=0.05")
climate_fit_path <- project_path("03_Output/02_ceara_pipeline/model_fits/climate_forced/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds")
if (file.exists(climate_fit_path)) {
  message("  Found: ", climate_fit_path, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  message("  NOT found -- CEARA_ALLOW_REFIT=TRUE, fitting now (expensive).")
  source(file.path(STAGE_DIR, "01_fit_ceara_climate_v49.R"), local = .GlobalEnv) # defines run_full_canary_fit(); does NOT auto-run under source()
  run_full_canary_fit()
  if (!file.exists(climate_fit_path)) stop("Climate-forced v4.9 fit did not produce the expected file: ", climate_fit_path, call. = FALSE)
} else {
  stop("Missing prerequisite:\n  ", climate_fit_path,
       "\n\nThis is expensive to regenerate and is not done automatically.\n",
       "To allow it: Sys.setenv(CEARA_ALLOW_REFIT = \"TRUE\") then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] Complete. All three required fits are present (existing or freshly generated).")
