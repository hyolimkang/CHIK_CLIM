# ============================================================
# STAGE 02 -- TRANSMISSION FITTING (Mato Grosso, orchestrator only)
# ============================================================
# Can refit Stan (global-q case-only + ad098 rescue, fixed-q profile
# q=0.020-0.070). Reuses each fit if it already exists; only (re)fits
# when missing AND MATO_GROSSO_ALLOW_REFIT="TRUE" -- a full refit is
# never triggered just because this file is sourced.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/06_mato_grosso_pipeline/02_transmission_fitting")
ALLOW_REFIT <- identical(Sys.getenv("MATO_GROSSO_ALLOW_REFIT", "FALSE"), "TRUE")
BASE_OUT <- project_path("03_Output/06_mato_grosso_pipeline/model_fits/v4_9_replication")

message("\n[Stage 02] 01/03: fit_mt_global_q_case_only ...")
fp <- file.path(BASE_OUT, "outputs/caseonly/mt_global_q_case_only.rds")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(STAGE_DIR, "01_fit_mt_global_q_case_only.R"), local = .GlobalEnv) # defines run_fit_mt_global_q_case_only(); does not auto-run under source()
  run_fit_mt_global_q_case_only()
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(MATO_GROSSO_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] 02/03: fit_mt_global_q_case_only_ad098 (rescue) ...")
fp <- file.path(BASE_OUT, "outputs/caseonly_ad098/mt_global_q_case_only_ad098.rds")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(STAGE_DIR, "02_fit_mt_global_q_case_only_ad098.R"), local = .GlobalEnv) # defines run_rescue(); does not auto-run under source()
  run_rescue()
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(MATO_GROSSO_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] 03/03: fit_mt_fixed_q_profile (q = 0.020-0.070) ...")
q_grid <- c(0.020, 0.030, 0.040, 0.050, 0.070)
q_tags <- sprintf("q%.3f", q_grid)
fixedq_paths <- file.path(BASE_OUT, "outputs/fixedq", paste0("mt_fixedq_", q_tags, ".rds"))
missing_fixedq <- fixedq_paths[!file.exists(fixedq_paths)]
if (length(missing_fixedq) == 0) {
  message("  Found all q = ", paste(q_grid, collapse = ", "), " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(STAGE_DIR, "03_fit_mt_fixed_q_profile.R"), local = .GlobalEnv) # defines fit_one_q(); does not auto-run under source()
  for (qv in q_grid) {
    tag <- sprintf("q%.3f", qv)
    fp <- file.path(BASE_OUT, "outputs/fixedq", paste0("mt_fixedq_", tag, ".rds"))
    if (!file.exists(fp)) fit_one_q(qv)
  }
} else {
  stop("Missing prerequisite(s):\n  ", paste(missing_fixedq, collapse = "\n  "),
       "\n\nSet Sys.setenv(MATO_GROSSO_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] Complete.")
