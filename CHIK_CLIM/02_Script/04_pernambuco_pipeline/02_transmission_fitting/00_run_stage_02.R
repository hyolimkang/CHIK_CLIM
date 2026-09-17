# ============================================================
# STAGE 02 -- TRANSMISSION FITTING (Pernambuco, orchestrator only)
# ============================================================
#
# Can refit Stan (global-q model A/B/B-full variants, fixed-q sweep).
# Reuses each fit if it already exists; only (re)fits when missing AND
# PERNAMBUCO_ALLOW_REFIT="TRUE" -- a full refit is never triggered just
# because this file is sourced. modelA_caseonly/modelB_U14/modelB_U14_full
# and the full q0.05-0.30 fixed-q sweep are all already fit and present on
# disk as of this refactor; this runner's fit branches are a documented
# fallback, not something exercised in normal validation.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this stage runner.")
  source(setup_path, local = .GlobalEnv)
}
STAGE_DIR <- project_path("02_Script/04_pernambuco_pipeline/02_transmission_fitting")
ALLOW_REFIT <- identical(Sys.getenv("PERNAMBUCO_ALLOW_REFIT", "FALSE"), "TRUE")
BASE_OUT <- project_path("03_Output/04_pernambuco_pipeline/model_fits/v4_9_replication")

message("\n[Stage 02] 01/02: fit_pe_global_q (modelA_caseonly / modelB_U14 / modelB_U14_full) ...")
required_tags <- c("modelA_caseonly", "modelB_U14", "modelB_U14_full")
fit_paths <- file.path(BASE_OUT, "outputs", required_tags, paste0("renewal_pe_global_q_fit_", required_tags, ".rds"))
missing_fits <- fit_paths[!file.exists(fit_paths)]
if (length(missing_fits) == 0) {
  message("  Found all of: ", paste(basename(fit_paths), collapse = ", "), " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(STAGE_DIR, "01_build_pe_serology_window.R"), local = .GlobalEnv) # helper, no auto-run under source()
  source(file.path(STAGE_DIR, "02_fit_pe_global_q.R"), local = .GlobalEnv) # defines run_fit_pe_global_q(); does not auto-run under source()
  stop("Missing fit(s):\n  ", paste(missing_fits, collapse = "\n  "),
       "\n\nCall run_fit_pe_global_q(tag = <one of modelA_caseonly/modelB_U14/modelB_U14_full>) with the matching",
       " sero_on/serology-full settings for that tag (see the script's own header for each model's definition),",
       " then re-source this stage.", call. = FALSE)
} else {
  stop("Missing prerequisite(s):\n  ", paste(missing_fits, collapse = "\n  "),
       "\n\nSet Sys.setenv(PERNAMBUCO_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] 02/02: fit_pe_fixed_q_sweep (q = 0.05-0.30) ...")
q_grid <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
q_tags <- sprintf("q%.2f", q_grid)
sweep_paths <- file.path(BASE_OUT, "outputs/fixed_q_sweep", q_tags, paste0("renewal_pe_fixedq_fit_", q_tags, ".rds"))
missing_sweep <- sweep_paths[!file.exists(sweep_paths)]
if (length(missing_sweep) == 0) {
  message("  Found all q = ", paste(q_grid, collapse = ", "), " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  if (!exists("build_pe_sero_window_index")) source(file.path(STAGE_DIR, "01_build_pe_serology_window.R"), local = .GlobalEnv)
  source(file.path(STAGE_DIR, "03_fit_pe_fixed_q_sweep.R"), local = .GlobalEnv) # defines run_fit_pe_fixedq(); does not auto-run under source()
  for (qv in q_grid) {
    tag <- sprintf("q%.2f", qv)
    fp <- file.path(BASE_OUT, "outputs/fixed_q_sweep", tag, paste0("renewal_pe_fixedq_fit_", tag, ".rds"))
    if (!file.exists(fp)) run_fit_pe_fixedq(qval = qv)
  }
} else {
  stop("Missing prerequisite(s):\n  ", paste(missing_sweep, collapse = "\n  "),
       "\n\nSet Sys.setenv(PERNAMBUCO_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this stage.", call. = FALSE)
}

message("\n[Stage 02] Complete.")
