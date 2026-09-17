# ============================================================
# CEARA CURRENT PRODUCTION PIPELINE -- TOP-LEVEL ORCHESTRATOR
# ============================================================
#
# ORCHESTRATOR ONLY -- contains no scientific equations. Sources the
# 00_run_stage_0X.R runner for each stage in [RUN_FROM_STAGE, RUN_TO_STAGE],
# in order. Each stage runner checks its own prerequisites and stops with
# a clear message if something upstream is missing.
#
# Safe default: RUN_FROM_STAGE = 3, RUN_TO_STAGE = 9 -- this NEVER touches
# Stage 02 (Stan fitting), so opening the project and sourcing this file
# with no configuration will not accidentally trigger an expensive refit,
# and will not even attempt Stage 01 (which stage 02 also does not need
# unless it is itself being (re)run). Stage 02 is the ONLY stage that can
# fit Stan, and even when explicitly included in the range, an individual
# fit only actually (re)runs if its output is missing (see
# 02_historical_transmission_fitting/00_run_stage_02.R).
#
# ---- Example usage ----
#
#   # full rebuild including Stan (only refits what is actually missing)
#   Sys.setenv(RUN_FROM_STAGE = "1", RUN_TO_STAGE = "9")
#   source("02_Script/40_renewal_model/38_ceara_current_pipeline/99_run_ceara_pipeline.R")
#
#   # reuse accepted fit; rerun downstream analysis only
#   Sys.setenv(RUN_FROM_STAGE = "4", RUN_TO_STAGE = "9")
#   source("02_Script/40_renewal_model/38_ceara_current_pipeline/99_run_ceara_pipeline.R")
#
#   # figures only
#   Sys.setenv(RUN_FROM_STAGE = "9", RUN_TO_STAGE = "9")
#   source("02_Script/40_renewal_model/38_ceara_current_pipeline/99_run_ceara_pipeline.R")
#
#   # default (no env vars set): Stage 03 through Stage 09
#   source("02_Script/40_renewal_model/38_ceara_current_pipeline/99_run_ceara_pipeline.R")

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this orchestrator.")
  source(setup_path, local = .GlobalEnv)
}

RUN_FROM_STAGE <- as.integer(Sys.getenv("RUN_FROM_STAGE", "3"))
RUN_TO_STAGE   <- as.integer(Sys.getenv("RUN_TO_STAGE", "9"))
if (is.na(RUN_FROM_STAGE) || is.na(RUN_TO_STAGE) || RUN_FROM_STAGE < 1L || RUN_TO_STAGE > 9L || RUN_FROM_STAGE > RUN_TO_STAGE) {
  stop("RUN_FROM_STAGE/RUN_TO_STAGE must be integers with 1 <= RUN_FROM_STAGE <= RUN_TO_STAGE <= 9.", call. = FALSE)
}

STAGE_FOLDERS <- c(
  "01_input_climate_preparation", "02_historical_transmission_fitting", "03_model_diagnostics",
  "04_age_reconstruction", "05_closed_loop_validation", "06_vaccine_counterfactual",
  "07_vaccine_validation", "08_impact_nnv_summaries", "09_publication_figures"
)
PIPELINE_DIR <- project_path("02_Script/40_renewal_model/38_ceara_current_pipeline")

# Including Stage 2 in an EXPLICIT run range is the "explicit user decision"
# required before any Stan refit is even considered; Stage 02's own runner
# still only refits whatever is actually missing (never routinely).
if (RUN_FROM_STAGE <= 2L && RUN_TO_STAGE >= 2L) {
  Sys.setenv(CEARA_ALLOW_REFIT = "TRUE")
  message("[pipeline] Stage 02 is in range -- Stan refitting is ALLOWED for any missing fit (existing fits are still reused).")
} else {
  Sys.setenv(CEARA_ALLOW_REFIT = "FALSE")
}

message(sprintf("[pipeline] Running Ceara current production pipeline, Stage %02d -> Stage %02d\n", RUN_FROM_STAGE, RUN_TO_STAGE))

for (stage_num in RUN_FROM_STAGE:RUN_TO_STAGE) {
  stage_folder <- STAGE_FOLDERS[stage_num]
  runner_path <- file.path(PIPELINE_DIR, stage_folder, sprintf("00_run_stage_%02d.R", stage_num))
  message(strrep("=", 60))
  message(sprintf("STAGE %02d: %s", stage_num, stage_folder))
  message(strrep("=", 60))
  source(runner_path, local = .GlobalEnv)
}

message("\n[pipeline] Done. Stages ", RUN_FROM_STAGE, " through ", RUN_TO_STAGE, " completed successfully.")
