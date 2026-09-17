# ============================================================
# PERNAMBUCO PIPELINE -- TOP-LEVEL ORCHESTRATOR
# ============================================================
# ORCHESTRATOR ONLY. Safe default: RUN_FROM_STAGE=3, RUN_TO_STAGE=5. This
# NEVER auto-enables refitting, even though Stage 04 (modelC/D ridge/curved
# reparam) and Stage 05 (5-strata pilot) also contain fit steps -- unlike
# Stage 02, those fits are all either already on disk (reused as-is) or, in
# the one case that is missing (the 5-strata pilot), a documented HALT
# (chain-1 warmup pathology; see 05_spatial_validation/00_run_stage_05.R)
# that must never be silently re-attempted just because Stage 05 is in
# range. PERNAMBUCO_ALLOW_REFIT is auto-enabled ONLY when Stage 02 itself is
# in range; a refit inside Stage 04/05 requires the user to set
# Sys.setenv(PERNAMBUCO_ALLOW_REFIT = "TRUE") themselves beforehand -- an
# explicit, separate decision, per the project-wide no-silent-refit rule.
# Stages 06-09 do not exist for Pernambuco (no age/vaccine/counterfactual/
# impact work was done here).

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this orchestrator.")
  source(setup_path, local = .GlobalEnv)
}

RUN_FROM_STAGE <- as.integer(Sys.getenv("RUN_FROM_STAGE", "3"))
RUN_TO_STAGE   <- as.integer(Sys.getenv("RUN_TO_STAGE", "5"))
if (is.na(RUN_FROM_STAGE) || is.na(RUN_TO_STAGE) || RUN_FROM_STAGE < 1L || RUN_TO_STAGE > 9L || RUN_FROM_STAGE > RUN_TO_STAGE) {
  stop("RUN_FROM_STAGE/RUN_TO_STAGE must be integers with 1 <= RUN_FROM_STAGE <= RUN_TO_STAGE <= 9.", call. = FALSE)
}

STAGE_FOLDERS <- c(
  "01" = "01_input_preparation", "02" = "02_transmission_fitting", "03" = "03_model_diagnostics",
  "04" = "04_sensitivity_identification", "05" = "05_spatial_validation",
  "06" = NA, "07" = NA, "08" = NA, "09" = NA # not applicable to Pernambuco
)
PIPELINE_DIR <- project_path("02_Script/04_pernambuco_pipeline")

if (RUN_FROM_STAGE <= 2L && RUN_TO_STAGE >= 2L) {
  Sys.setenv(PERNAMBUCO_ALLOW_REFIT = "TRUE")
  message("[pipeline] Stage 02 is in range -- Stan refitting is ALLOWED for any missing fit (existing fits are still reused).")
} else {
  Sys.setenv(PERNAMBUCO_ALLOW_REFIT = "FALSE")
  if (RUN_FROM_STAGE <= 5L && RUN_TO_STAGE >= 4L) {
    message("[pipeline] Stage 04/05 fit steps will reuse existing fits only (PERNAMBUCO_ALLOW_REFIT=FALSE). The one fit currently missing on disk (the 5-strata pilot, a documented halt) will be skipped, not refit -- set Sys.setenv(PERNAMBUCO_ALLOW_REFIT = \"TRUE\") yourself first if you deliberately want to re-attempt it.")
  }
}

message(sprintf("[pipeline] Running Pernambuco pipeline, Stage %02d -> Stage %02d\n", RUN_FROM_STAGE, RUN_TO_STAGE))
for (stage_num in RUN_FROM_STAGE:RUN_TO_STAGE) {
  key <- sprintf("%02d", stage_num)
  stage_folder <- STAGE_FOLDERS[[key]]
  if (is.na(stage_folder)) { message(sprintf("STAGE %s: (not applicable to Pernambuco -- skipped)", key)); next }
  runner_path <- file.path(PIPELINE_DIR, stage_folder, sprintf("00_run_stage_%s.R", key))
  message(strrep("=", 60)); message(sprintf("STAGE %s: %s", key, stage_folder)); message(strrep("=", 60))
  source(runner_path, local = .GlobalEnv)
}
message("\n[pipeline] Done. Stages ", RUN_FROM_STAGE, " through ", RUN_TO_STAGE, " completed successfully.")
