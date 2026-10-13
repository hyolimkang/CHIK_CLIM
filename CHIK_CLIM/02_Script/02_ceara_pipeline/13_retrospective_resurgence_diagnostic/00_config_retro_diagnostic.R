# Retrospective resurgence diagnostic: central configuration.
#
# Diagnoses WHY the H1/H2 retrospective holdout validation (built in
# 12_future_year_effect_validation) under- or over-predicts the 2022
# Ceara resurgence, decomposing susceptibility / transmission /
# importation / observation-process contributions. NO NEW STAN FITTING
# in this module -- reuses the training-only tau_year posteriors
# (ce_tau_year_training_H1.rds / H2.rds) already fitted and accepted in
# that module. All outputs live under this module's own
# *_retrospective_resurgence_diagnostic/ trees; nothing accepted is
# modified. NO SSP/ISIMIP forcing, NO vaccination.

if (!exists("ROOT")) stop("ROOT not found -- source CHIK_CLIM/00_project_setup.R first (opening CHIK_CLIM.Rproj does this automatically).")

# ---- Read-only accepted/prior-module sources (never written to) ----
DIR_YV      <- file.path(ROOT, "02_Script/02_ceara_pipeline/12_future_year_effect_validation")
DIR_FITS_YV <- file.path(ROOT, "03_Output/02_ceara_pipeline/model_fits/future_year_effect_validation")
DIR_RESULTS_YV_PRIOR <- file.path(ROOT, "03_Output/02_ceara_pipeline/results/future_year_effect_validation")

for (f in c(file.path(DIR_FITS_YV, "ce_tau_year_training_H1.rds"), file.path(DIR_FITS_YV, "ce_tau_year_training_H2.rds"))) {
  if (!file.exists(f)) stop("Missing required upstream training-only fit, cannot proceed: ", f)
}

# Reuse Step-2 future-renewal primitives AND the year-effect-validation
# config (HOLDOUTS, IMPORTATION_SENSITIVITY_GRID, ENSEMBLE_RETRO, QVAL,
# ACCEPTED_ANOMALY_CSV, build_training_only_tau_data(), WAVE_DETECTION_YV),
# all UNMODIFIED.
source(file.path(DIR_YV, "00_config_year_effect_validation.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R")) # build_state_weekly()
source(file.path(DIR_YV, "03_build_training_only_data.R"))

# ---- This module's own output locations (new, additive only) ----
DIR_SCRIPTS_RD <- file.path(ROOT, "02_Script/02_ceara_pipeline/13_retrospective_resurgence_diagnostic")
DIR_REPORTS_RD <- file.path(ROOT, "03_Output/02_ceara_pipeline/reports/retrospective_resurgence_diagnostic")
DIR_TABLES_RD  <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/retrospective_resurgence_diagnostic")
DIR_FIGURES_RD <- file.path(ROOT, "03_Output/02_ceara_pipeline/figures/retrospective_resurgence_diagnostic")
DIR_RESULTS_RD <- file.path(ROOT, "03_Output/02_ceara_pipeline/results/retrospective_resurgence_diagnostic")
for (d in c(DIR_REPORTS_RD, DIR_TABLES_RD, DIR_FIGURES_RD, DIR_RESULTS_RD)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Diagnostic-simulation settings ----
# Same reference importation level, sensitivity grid and holdout
# definitions as 12_future_year_effect_validation (NOT re-chosen here).
# Diagnostic re-simulation uses IDENTICAL seeds to the saved
# PHASE3_4_holdout_simulation_results.rds cells so that S/X/R0/R_eff
# reproduce bit-for-bit; only ADDITIONAL columns (realised imports M,
# U, lambda_t, per-week year_effect/A_year) are newly captured.
DIAG_ENSEMBLE <- ENSEMBLE_RETRO$full   # n_draws=100, n_replicates=20 -- matches the reference cells being decomposed
DIAG_SEED_OFFSET <- BASE_SEED_YV + 777L # identical to 06_phase3_4_retrospective_validation.R's reference-cell call

# Experiment C (initial-susceptibility decomposition) needs its own,
# clearly-separate seed stream since it perturbs the ensemble
# construction itself (posterior-mean-fixed S/U at cutoff) rather than
# reproducing an existing cell.
EXPERIMENT_C_SEED_OFFSET <- BASE_SEED_YV + 8888L
EXPERIMENT_C_ENSEMBLE <- ENSEMBLE_RETRO$smoke # n_draws=10, n_replicates=5 -- a controlled, single-mechanism contrast; smoke scale, expanded only if MC noise is too large to distinguish arms

# ---- Phase 3: outbreak-threshold audit ----
# The existing WAVE_DETECTION_YV threshold (>=100 absolute cases AND
# >=5/100k cumulative incidence, at CE-state population scale ~9M) is
# evaluated against an alternative, empirically derived from PRE-HOLDOUT
# (2015-2019 training-period) historical epidemic sizes only -- fixed
# BEFORE any holdout comparison is computed, never tuned to validation outcomes.
THRESHOLD_SENSITIVITY <- list(
  current  = list(abs_min = 100,  incidence_min_per_100k = 5),
  stricter = list(abs_min = 1000, incidence_min_per_100k = 20)  # ~4x current incidence floor; justified empirically in 01_phase0_design_audit.R / Phase 3 against 2015-2019 wave sizes
)

message(sprintf("[retro_diagnostic config] module ready. Reusing training fits: %s | holdouts: %s",
                 paste(c("H1","H2"), collapse=", "), paste(names(HOLDOUTS), collapse=", ")))
