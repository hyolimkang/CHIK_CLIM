# ============================================================
# STAGE 08 -- IMPACT + NNV SUMMARIES (orchestrator only)
# ============================================================
#
# Order: 01 summarise vaccine impacts -> 02 calculate NNV -> 03 validate
# NNV. Reads the already-saved, VALIDATED three-arm result only -- does
# not rerun transmission or vaccine simulation, and does not make
# publication figures (that is Stage 09).

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}

STAGE_DIR <- project_path("02_Script/02_ceara_pipeline/08_impact_summaries")

three_arm_results_path <- project_path("03_Output/02_ceara_pipeline/results/ceara_historical_age12_vaccine_counterfactual/three_arm_results.rds")
if (!file.exists(three_arm_results_path)) {
  stop("Missing prerequisite:\n  ", three_arm_results_path,
       "\n\nRun Stage 07 first:\n  07_counterfactuals (00_run_stage_07.R)", call. = FALSE)
}
qa_summary_path <- project_path("03_Output/02_ceara_pipeline/diagnostics/ceara_historical_age12_vaccine_counterfactual/QA_overall_summary.csv")
if (!file.exists(qa_summary_path)) {
  stop("Missing prerequisite (three-arm result has not been validated):\n  ", qa_summary_path,
       "\n\nRun Stage 07 first (it performs the counterfactual simulation AND its validation):\n  07_counterfactuals (00_run_stage_07.R)", call. = FALSE)
}
qa_summary <- read.csv(qa_summary_path)
# NA rows are informational (not pass/fail) items, e.g. qa9's first-divergence WEEK NUMBER is stored under
# "pass" as NA because it is not a logical value -- na.rm=TRUE correctly ignores those, not "pass by default".
if (!all(as.logical(qa_summary$pass), na.rm = TRUE)) {
  stop("Stage 07 QA did not fully pass (see ", qa_summary_path, ") -- re-run and fix Stage 07 before computing impact/NNV summaries.", call. = FALSE)
}
message("[Stage 08] Prerequisites OK: three-arm result exists and Stage 07 QA passed.")

message("\n[Stage 08] 01/03: summarise_vaccine_impacts ...")
source(file.path(STAGE_DIR, "01_summarise_vaccine_impacts.R"), local = .GlobalEnv)

message("\n[Stage 08] 02/03: calculate_nnv ...")
source(file.path(STAGE_DIR, "02_calculate_nnv.R"), local = .GlobalEnv)

message("\n[Stage 08] 03/03: validate_nnv (fails loudly on any QA failure) ...")
source(file.path(STAGE_DIR, "03_validate_nnv.R"), local = .GlobalEnv)

message("\n[Stage 08] Complete.")
