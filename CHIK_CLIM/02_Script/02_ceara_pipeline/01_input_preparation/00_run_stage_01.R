# ============================================================
# STAGE 01 -- INPUT / CLIMATE PREPARATION (orchestrator only)
# ============================================================
#
# Sources the three climate-preparation scripts in order. No new analysis
# logic here -- each script's scientific content is unchanged from its
# original location under 36_climate_forced_v4_9/scripts/.
#
# Note (structural, not changed by this refactor): script 02 reads the
# baseline v4.9 q=0.05 fit (Stage 02B's output) to obtain the exact
# harmonic/week alignment used at fitting time. That fit already existed
# in this project before the climate work began, so Stage 01 has always
# implicitly depended on it. See REFACTOR_NOTES_SCIENTIFIC_ISSUES_NOT_CHANGED.md.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}

STAGE_DIR <- project_path("02_Script/02_ceara_pipeline/01_input_preparation")

national_climate_path <- project_path("03_Output/07_national_pipeline/tables/national_wave_analysis/brazil_chik_uf_weekly_climate.csv")
if (!file.exists(national_climate_path)) {
  stop("Missing prerequisite:\n  ", national_climate_path,
       "\n\nThis is built by the national pipeline's provisional climate-susceptibility checkpoint (02_Script/07_national_pipeline/07_provisional_climate_susceptibility_checkpoint/), not by this stage.", call. = FALSE)
}
v49_baseline_path <- project_path("03_Output/02_ceara_pipeline/model_fits/baseline/v4_9_hierarchical_seasonality/outputs/q0.05/renewal_ceara_v4_9_fit_q0.05.rds")
if (!file.exists(v49_baseline_path)) {
  stop("Missing prerequisite (needed by script 02 of this stage):\n  ", v49_baseline_path,
       "\n\nRun Stage 02 (STEP 02B) first:\n  02_transmission_fitting (00_run_stage_02.R)", call. = FALSE)
}
message("[Stage 01] Prerequisites OK.")

message("\n[Stage 01] 01/03: build_climate_anomaly_covariates ...")
source(file.path(STAGE_DIR, "01_build_climate_anomaly_covariates.R"), local = .GlobalEnv)

message("\n[Stage 01] 02/03: check_climate_orthogonality_and_priors ...")
source(file.path(STAGE_DIR, "02_check_climate_orthogonality_and_priors.R"), local = .GlobalEnv)

message("\n[Stage 01] 03/03: check_climate_multiplier_prior ...")
source(file.path(STAGE_DIR, "03_check_climate_multiplier_prior.R"), local = .GlobalEnv)

final_product <- project_path("03_Output/02_ceara_pipeline/tables/climate_forced_v4_9/climate_anomaly_covariates.csv")
if (!file.exists(final_product)) stop("Stage 01 ran but did not produce the expected product: ", final_product, call. = FALSE)
message("\n[Stage 01] Complete. Required product present: ", final_product)
