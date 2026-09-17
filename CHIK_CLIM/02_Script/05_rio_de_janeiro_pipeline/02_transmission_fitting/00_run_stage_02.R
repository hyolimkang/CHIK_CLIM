# ============================================================
# STAGE 02 -- TRANSMISSION FITTING (Rio de Janeiro, orchestrator only)
# ============================================================
# Can refit Stan (state case-only + ad098 rescue, targeted fixed-q sweep,
# city case-only + U10 + ad099 rescue). Reuses each fit if it already
# exists; only (re)fits when missing AND RIO_DE_JANEIRO_ALLOW_REFIT="TRUE"
# -- a full refit is never triggered just because this file is sourced.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/05_rio_de_janeiro_pipeline/02_transmission_fitting")
ALLOW_REFIT <- identical(Sys.getenv("RIO_DE_JANEIRO_ALLOW_REFIT", "FALSE"), "TRUE")
BASE_OUT <- project_path("03_Output/05_rio_de_janeiro_pipeline/model_fits/v4_9_replication")

message("\n[Stage 02] 01/06: fit_rj_global_q_case_only (state) ...")
fp <- file.path(BASE_OUT, "outputs/caseonly/rj_global_q_case_only.rds")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(STAGE_DIR, "01_fit_rj_global_q_case_only.R"), local = .GlobalEnv) # defines run_fit_rj_global_q_case_only(); does not auto-run under source()
  run_fit_rj_global_q_case_only()
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(RIO_DE_JANEIRO_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] 02/06: fit_rj_global_q_case_only_ad098 (state rescue) ...")
fp <- file.path(BASE_OUT, "outputs/caseonly_ad098/rj_global_q_case_only_ad098.rds")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(STAGE_DIR, "02_fit_rj_global_q_case_only_ad098.R"), local = .GlobalEnv) # defines run_rescue(); does not auto-run under source()
  run_rescue()
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(RIO_DE_JANEIRO_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] 03/06: fit_rj_targeted_fixed_q (state, q = 0.010-0.020) ...")
q_grid <- c(0.0100, 0.0125, 0.0150, 0.0175, 0.0200)
q_tags <- sprintf("q%.4f", q_grid)
fixedq_paths <- file.path(BASE_OUT, "outputs/fixedq", paste0("rj_fixedq_", q_tags, ".rds"))
missing_fixedq <- fixedq_paths[!file.exists(fixedq_paths)]
if (length(missing_fixedq) == 0) {
  message("  Found all q = ", paste(q_grid, collapse = ", "), " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(STAGE_DIR, "03_fit_rj_targeted_fixed_q.R"), local = .GlobalEnv) # defines fit_one_q(); does not auto-run under source()
  for (qv in q_grid) {
    tag <- sprintf("q%.4f", qv)
    fp <- file.path(BASE_OUT, "outputs/fixedq", paste0("rj_fixedq_", tag, ".rds"))
    if (!file.exists(fp)) fit_one_q(qv)
  }
} else {
  stop("Missing prerequisite(s):\n  ", paste(missing_fixedq, collapse = "\n  "),
       "\n\nSet Sys.setenv(RIO_DE_JANEIRO_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] 04/06: fit_rj_city_global_q_case_only (city) ...")
fp <- file.path(BASE_OUT, "outputs/city_caseonly/rj_city_global_q_case_only.rds")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(STAGE_DIR, "04_fit_rj_city_global_q_case_only.R"), local = .GlobalEnv) # defines run_fit_rj_city_case_only(); does not auto-run under source()
  run_fit_rj_city_case_only()
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(RIO_DE_JANEIRO_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] 05/06: fit_rj_city_u10 (city + U10 serology) ...")
fp <- file.path(BASE_OUT, "outputs/city_u10/rj_city_u10.rds")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(STAGE_DIR, "05_fit_rj_city_u10.R"), local = .GlobalEnv) # defines run_fit_rj_city_u10(); does not auto-run under source()
  run_fit_rj_city_u10()
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(RIO_DE_JANEIRO_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] 06/06: fit_rj_city_u10_ad099 (city + U10 rescue) ...")
fp <- file.path(BASE_OUT, "outputs/city_u10_ad099/rj_city_u10_ad099.rds")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(STAGE_DIR, "06_fit_rj_city_u10_ad099.R"), local = .GlobalEnv) # defines run_fit_rj_city_u10_ad099(); does not auto-run under source()
  run_fit_rj_city_u10_ad099()
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(RIO_DE_JANEIRO_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] Complete.")
