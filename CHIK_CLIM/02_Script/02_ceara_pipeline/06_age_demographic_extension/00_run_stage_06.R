# ============================================================
# STAGE 06 -- AGE / DEMOGRAPHIC EXTENSION (orchestrator only)
# ============================================================
#
# Merges what were previously two separate stages in 38_ceara_current_pipeline
# (04_age_reconstruction, 05_closed_loop_validation) into one, since both are
# read-only post-processing of the single accepted climate-forced fit and
# neither refits Stan. Execution order: 01, 03, 04, 05 (age reconstruction),
# then 06 (closed-loop replay validation, which itself needs 01's age-shares
# product). Script 02 (02_age_cohort_simulator.R) is a HELPER (defines
# simulate_age_cohort()) sourced internally by 03/04/06 where needed -- it is
# never run as its own standalone analysis step.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}

STAGE_DIR <- project_path("02_Script/02_ceara_pipeline/06_age_demographic_extension")

climate_fit_path <- project_path("03_Output/model_fits/ceara/climate_forced_v4_9/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds")
age_pop_source_path <- project_path("01_Data/ibge_pop_uf_single_age_expanded_2015_2024.rds")
missing <- c(climate_fit_path, age_pop_source_path)[!file.exists(c(climate_fit_path, age_pop_source_path))]
if (length(missing) > 0) {
  stop("Missing prerequisite(s):\n  ", paste(missing, collapse = "\n  "),
       "\n\nIf the climate fit is missing, run Stage 02 first:\n  02_transmission_fitting (00_run_stage_02.R)", call. = FALSE)
}
message("[Stage 06] Prerequisites OK: accepted climate fit and age population source both found.")

message("\n[Stage 06] 01/05: build_age_population_shares ...")
source(file.path(STAGE_DIR, "01_build_age_population_shares.R"), local = .GlobalEnv)

message("\n[Stage 06] 02/05: validate_age_aggregation_identity ...")
source(file.path(STAGE_DIR, "03_validate_age_aggregation_identity.R"), local = .GlobalEnv) # sources 02_age_cohort_simulator.R internally as a helper

message("\n[Stage 06] 03/05: reconstruct_age_specific_susceptibility ...")
source(file.path(STAGE_DIR, "04_reconstruct_age_specific_susceptibility.R"), local = .GlobalEnv) # sources 02_age_cohort_simulator.R internally as a helper

message("\n[Stage 06] 04/05: plot_age_specific_susceptibility ...")
source(file.path(STAGE_DIR, "05_plot_age_specific_susceptibility.R"), local = .GlobalEnv)

message("\n[Stage 06] 05/05: validate_closed_loop_replay -- does the closed-loop age simulator reproduce the accepted fit with vaccination OFF? ...")
source(file.path(STAGE_DIR, "06_validate_closed_loop_replay.R"), local = .GlobalEnv) # needs 01's CE_age_population_shares.csv, already produced above

message("\n[Stage 06] Complete.")
