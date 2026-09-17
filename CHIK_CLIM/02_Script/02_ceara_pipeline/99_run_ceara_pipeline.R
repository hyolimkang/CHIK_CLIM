# ============================================================
# CEARA PIPELINE -- TOP-LEVEL ORCHESTRATOR
# ============================================================
#
# ORCHESTRATOR ONLY -- contains no scientific equations. Sources the
# 00_run_stage_0X.R runner for each EXISTING stage in [RUN_FROM_STAGE,
# RUN_TO_STAGE], in order. Stages 04 (sensitivity_identification) and 05
# (spatial_validation) do not exist for Ceara (no such analysis was done
# here) and are silently skipped with a message -- the stage NUMBERS stay
# fixed across all state pipelines (02_ceara_pipeline, 03_bahia_pipeline,
# ...) even where a given state has no content for a slot, so "stage 6 =
# age/demographic extension" means the same thing everywhere.
#
# Safe default: RUN_FROM_STAGE = 3, RUN_TO_STAGE = 9 -- this NEVER touches
# Stage 02 (Stan fitting), so opening the project and sourcing this file
# with no configuration will not accidentally trigger an expensive refit.
# Stage 02 is the ONLY stage that can fit Stan, and even when explicitly
# included in the range, an individual fit only actually (re)runs if its
# output is missing (see 02_transmission_fitting/00_run_stage_02.R).
#
# ---- Example usage ----
#
#   # full rebuild including Stan (only refits what is actually missing)
#   Sys.setenv(RUN_FROM_STAGE = "1", RUN_TO_STAGE = "9")
#   source("02_Script/02_ceara_pipeline/99_run_ceara_pipeline.R")
#
#   # reuse accepted fit; rerun downstream analysis only (the safe default)
#   source("02_Script/02_ceara_pipeline/99_run_ceara_pipeline.R")
#
#   # figures only
#   Sys.setenv(RUN_FROM_STAGE = "9", RUN_TO_STAGE = "9")
#   source("02_Script/02_ceara_pipeline/99_run_ceara_pipeline.R")

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this orchestrator.")
  source(setup_path, local = .GlobalEnv)
}

RUN_FROM_STAGE <- as.integer(Sys.getenv("RUN_FROM_STAGE", "3"))
RUN_TO_STAGE   <- as.integer(Sys.getenv("RUN_TO_STAGE", "9"))
if (is.na(RUN_FROM_STAGE) || is.na(RUN_TO_STAGE) || RUN_FROM_STAGE < 1L || RUN_TO_STAGE > 9L || RUN_FROM_STAGE > RUN_TO_STAGE) {
  stop("RUN_FROM_STAGE/RUN_TO_STAGE must be integers with 1 <= RUN_FROM_STAGE <= RUN_TO_STAGE <= 9.", call. = FALSE)
}

STAGE_FOLDERS <- c(
  "01" = "01_input_preparation", "02" = "02_transmission_fitting", "03" = "03_model_diagnostics",
  "04" = NA, "05" = NA, # not applicable to Ceara
  "06" = "06_age_demographic_extension", "07" = "07_counterfactuals",
  "08" = "08_impact_summaries", "09" = "09_figures_reporting"
)
PIPELINE_DIR <- project_path("02_Script/02_ceara_pipeline")

# Including Stage 2 in an EXPLICIT run range is the "explicit user decision"
# required before any Stan refit is even considered; Stage 02's own runner
# still only refits whatever is actually missing (never routinely).
if (RUN_FROM_STAGE <= 2L && RUN_TO_STAGE >= 2L) {
  Sys.setenv(CEARA_ALLOW_REFIT = "TRUE")
  message("[pipeline] Stage 02 is in range -- Stan refitting is ALLOWED for any missing fit (existing fits are still reused).")
} else {
  Sys.setenv(CEARA_ALLOW_REFIT = "FALSE")
}

message(sprintf("[pipeline] Running Ceara pipeline, Stage %02d -> Stage %02d\n", RUN_FROM_STAGE, RUN_TO_STAGE))

for (stage_num in RUN_FROM_STAGE:RUN_TO_STAGE) {
  key <- sprintf("%02d", stage_num)
  stage_folder <- STAGE_FOLDERS[[key]]
  if (is.na(stage_folder)) {
    message(sprintf("STAGE %s: (not applicable to Ceara -- skipped)", key))
    next
  }
  runner_path <- file.path(PIPELINE_DIR, stage_folder, sprintf("00_run_stage_%s.R", key))
  message(strrep("=", 60))
  message(sprintf("STAGE %s: %s", key, stage_folder))
  message(strrep("=", 60))
  source(runner_path, local = .GlobalEnv)
}

message("\n[pipeline] Done. Stages ", RUN_FROM_STAGE, " through ", RUN_TO_STAGE, " completed successfully.")
