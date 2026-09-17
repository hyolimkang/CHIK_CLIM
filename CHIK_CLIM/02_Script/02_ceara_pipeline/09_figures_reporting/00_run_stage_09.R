# ============================================================
# STAGE 09 -- PUBLICATION FIGURES (orchestrator only)
# ============================================================
#
# Order: 01 main vaccine figures -> 02 supplementary vaccine figures ->
# 03 NNV figures. Reads saved Stage 06 (three-arm results) and Stage 08
# (impact + NNV summary) results only. No simulation, no fitting.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}

STAGE_DIR <- project_path("02_Script/02_ceara_pipeline/09_figures_reporting")
RESULTS_DIR_HISTVAX <- project_path("03_Output/02_ceara_pipeline/results/ceara_historical_age12_vaccine_counterfactual")
TABLES_DIR_HISTVAX <- project_path("03_Output/02_ceara_pipeline/tables/ceara_historical_age12_vaccine_counterfactual")

required_files <- c(
  file.path(RESULTS_DIR_HISTVAX, "three_arm_results.rds"),
  file.path(TABLES_DIR_HISTVAX, "TABLE3_effect_decomposition.csv"),
  file.path(TABLES_DIR_HISTVAX, "TABLE3_decomposition_by_draw.csv"),
  file.path(RESULTS_DIR_HISTVAX, "nnv/nnv_posterior_draws.rds"),
  file.path(RESULTS_DIR_HISTVAX, "nnv/nnv_year_end_cumulative_by_draw.rds")
)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Missing prerequisite(s):\n  ", paste(missing, collapse = "\n  "),
       "\n\nRun Stage 08 first:\n  08_impact_summaries (00_run_stage_08.R)", call. = FALSE)
}
message("[Stage 09] Prerequisites OK: Stage 06/08 result and summary files found.")

message("\n[Stage 09] 01/03: make_main_vaccine_figures ...")
source(file.path(STAGE_DIR, "01_make_main_vaccine_figures.R"), local = .GlobalEnv)

message("\n[Stage 09] 02/03: make_supplementary_vaccine_figures ...")
source(file.path(STAGE_DIR, "02_make_supplementary_vaccine_figures.R"), local = .GlobalEnv)

message("\n[Stage 09] 03/03: make_nnv_figures ...")
source(file.path(STAGE_DIR, "03_make_nnv_figures.R"), local = .GlobalEnv)

message("\n[Stage 09] Complete. All publication figures regenerated.")
