# ============================================================
# 06_national_wave_analysis (National pipeline, orchestrator only)
# ============================================================
# Nationwide early-Re (per major wave) + annual-SHAPE retrospective
# susceptibility: fit (01, 02 -- Stan HMC, existence-gated), deterministic
# weekly reconstruction + master-table build (03, 04 -- always re-run,
# fast/non-Stan), and diagnostics/plots (05-08 -- always re-run).

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this runner.")
  source(setup_path, local = .GlobalEnv)
}
FOLDER_DIR <- project_path("02_Script/07_national_pipeline/06_national_wave_analysis")
ALLOW_REFIT <- identical(Sys.getenv("NATIONAL_ALLOW_REFIT", "FALSE"), "TRUE")
TABLE_DIR <- project_path("03_Output/07_national_pipeline/tables/national_wave_analysis")

message("\n[06_national_wave_analysis] 01/08: fit_brazil_major_wave_early_re ...")
fp <- file.path(TABLE_DIR, "brazil_chik_major_wave_early_re.csv")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(FOLDER_DIR, "01_fit_brazil_major_wave_early_re.R"), local = .GlobalEnv) # defines run_national_early_re(); does not auto-run under source()
  run_national_early_re()
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(NATIONAL_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this runner.", call. = FALSE)
}

message("\n[06_national_wave_analysis] 02/08: fit_brazil_annual_shape_susceptibility ...")
fp <- file.path(TABLE_DIR, "brazil_chik_annual_shape_hmc_gate.csv")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not refit).")
} else if (ALLOW_REFIT) {
  source(file.path(FOLDER_DIR, "02_fit_brazil_annual_shape_susceptibility.R"), local = .GlobalEnv) # defines run_national_annual_shape(); does not auto-run under source()
  run_national_annual_shape()
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(NATIONAL_ALLOW_REFIT = \"TRUE\") to allow fitting, then re-source this runner.", call. = FALSE)
}

message("\n[06_national_wave_analysis] 03/08: reconstruct_brazil_weekly_susceptibility ...")
source(file.path(FOLDER_DIR, "03_reconstruct_brazil_weekly_susceptibility.R"), local = .GlobalEnv) # defines run_national_weekly_susceptibility(); does not auto-run under source()
run_national_weekly_susceptibility()

message("\n[06_national_wave_analysis] 04/08: build_brazil_wave_analysis_master ...")
source(file.path(FOLDER_DIR, "04_build_brazil_wave_analysis_master.R"), local = .GlobalEnv) # defines run_national_wave_master(); does not auto-run under source()
run_national_wave_master()

message("\n[06_national_wave_analysis] 05/08: plot_brazil_early_re_diagnostics ...")
source(file.path(FOLDER_DIR, "05_plot_brazil_early_re_diagnostics.R"), local = .GlobalEnv) # defines run_early_re_diagnostics(); does not auto-run under source()
run_early_re_diagnostics()

message("\n[06_national_wave_analysis] 06/08: plot_brazil_annual_early_re_maps ...")
source(file.path(FOLDER_DIR, "06_plot_brazil_annual_early_re_maps.R"), local = .GlobalEnv) # defines run_annual_early_re_maps(); does not auto-run under source()
run_annual_early_re_maps()

message("\n[06_national_wave_analysis] 07/08: plot_brazil_susceptibility_diagnostics ...")
source(file.path(FOLDER_DIR, "07_plot_brazil_susceptibility_diagnostics.R"), local = .GlobalEnv) # defines run_susceptibility_diagnostics(); does not auto-run under source()
run_susceptibility_diagnostics()

message("\n[06_national_wave_analysis] 08/08: diagnose_brazil_annual_q_implied ...")
source(file.path(FOLDER_DIR, "08_diagnose_brazil_annual_q_implied.R"), local = .GlobalEnv) # defines run_annual_q_implied_diagnostic(); does not auto-run under source()
run_annual_q_implied_diagnostic()

message("\n[06_national_wave_analysis] Complete.")
