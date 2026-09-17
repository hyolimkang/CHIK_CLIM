# ============================================================
# MATO GROSSO PIPELINE -- TOP-LEVEL ORCHESTRATOR
# ============================================================
# ORCHESTRATOR ONLY. Safe default: RUN_FROM_STAGE=3, RUN_TO_STAGE=3 (never
# touches Stage 02 Stan fitting unless explicitly included in range, in
# which case MATO_GROSSO_ALLOW_REFIT is set automatically -- see Stage 02's
# own runner for per-fit reuse-if-exists logic). Stages 04-09 do not exist
# for Mato Grosso (no sensitivity/spatial/age/vaccine/counterfactual/impact
# work was done here).

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this orchestrator.")
  source(setup_path, local = .GlobalEnv)
}

RUN_FROM_STAGE <- as.integer(Sys.getenv("RUN_FROM_STAGE", "3"))
RUN_TO_STAGE   <- as.integer(Sys.getenv("RUN_TO_STAGE", "3"))
if (is.na(RUN_FROM_STAGE) || is.na(RUN_TO_STAGE) || RUN_FROM_STAGE < 1L || RUN_TO_STAGE > 9L || RUN_FROM_STAGE > RUN_TO_STAGE) {
  stop("RUN_FROM_STAGE/RUN_TO_STAGE must be integers with 1 <= RUN_FROM_STAGE <= RUN_TO_STAGE <= 9.", call. = FALSE)
}

STAGE_FOLDERS <- c(
  "01" = "01_input_preparation", "02" = "02_transmission_fitting", "03" = "03_model_diagnostics",
  "04" = NA, "05" = NA, "06" = NA, "07" = NA, "08" = NA, "09" = NA # not applicable to Mato Grosso
)
PIPELINE_DIR <- project_path("02_Script/06_mato_grosso_pipeline")

if (RUN_FROM_STAGE <= 2L && RUN_TO_STAGE >= 2L) {
  Sys.setenv(MATO_GROSSO_ALLOW_REFIT = "TRUE")
  message("[pipeline] Stage 02 is in range -- Stan refitting is ALLOWED for any missing fit (existing fits are still reused).")
} else {
  Sys.setenv(MATO_GROSSO_ALLOW_REFIT = "FALSE")
}

message(sprintf("[pipeline] Running Mato Grosso pipeline, Stage %02d -> Stage %02d\n", RUN_FROM_STAGE, RUN_TO_STAGE))
for (stage_num in RUN_FROM_STAGE:RUN_TO_STAGE) {
  key <- sprintf("%02d", stage_num)
  stage_folder <- STAGE_FOLDERS[[key]]
  if (is.na(stage_folder)) { message(sprintf("STAGE %s: (not applicable to Mato Grosso -- skipped)", key)); next }
  runner_path <- file.path(PIPELINE_DIR, stage_folder, sprintf("00_run_stage_%s.R", key))
  message(strrep("=", 60)); message(sprintf("STAGE %s: %s", key, stage_folder)); message(strrep("=", 60))
  source(runner_path, local = .GlobalEnv)
}
message("\n[pipeline] Done. Stages ", RUN_FROM_STAGE, " through ", RUN_TO_STAGE, " completed successfully.")
