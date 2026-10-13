# STEP 6 retrospective failure diagnosis: isolated, read-only configuration.
#
# This module performs no fitting and no calibration. It reuses the accepted
# H1/H2 training-only posterior bundles and the bit-for-bit diagnostic
# re-simulations produced by 13_retrospective_resurgence_diagnostic. All new
# outputs are written under step6_retrospective_failure_diagnosis/.

if (!exists("ROOT")) stop("ROOT not found -- source CHIK_CLIM/00_project_setup.R first.")

DIR_STEP6 <- file.path(ROOT, "02_Script/02_ceara_pipeline/15_step6_retrospective_failure_diagnosis")
DIR_YV <- file.path(ROOT, "02_Script/02_ceara_pipeline/12_future_year_effect_validation")
DIR_RD <- file.path(ROOT, "02_Script/02_ceara_pipeline/13_retrospective_resurgence_diagnostic")

DIR_FITS_YV <- file.path(ROOT, "03_Output/02_ceara_pipeline/model_fits/future_year_effect_validation")
DIR_RESULTS_YV <- file.path(ROOT, "03_Output/02_ceara_pipeline/results/future_year_effect_validation")
DIR_RESULTS_RD <- file.path(ROOT, "03_Output/02_ceara_pipeline/results/retrospective_resurgence_diagnostic")
DIR_TABLES_RD <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/retrospective_resurgence_diagnostic")
DIR_TABLES_HMC <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/hmc_reliability_audit")

DIR_REPORTS_STEP6 <- file.path(ROOT, "03_Output/02_ceara_pipeline/reports/step6_retrospective_failure_diagnosis")
DIR_TABLES_STEP6 <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/step6_retrospective_failure_diagnosis")
DIR_FIGURES_STEP6 <- file.path(ROOT, "03_Output/02_ceara_pipeline/figures/step6_retrospective_failure_diagnosis")
DIR_RESULTS_STEP6 <- file.path(ROOT, "03_Output/02_ceara_pipeline/results/step6_retrospective_failure_diagnosis")
for (d in c(DIR_REPORTS_STEP6, DIR_TABLES_STEP6, DIR_FIGURES_STEP6, DIR_RESULTS_STEP6)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

DIAG_RDS <- file.path(DIR_RESULTS_RD, "PHASE1_diagnostic_simulation_results.rds")
ORIG_RDS <- file.path(DIR_RESULTS_YV, "PHASE3_4_holdout_simulation_results.rds")
OBS_RDS <- file.path(DIR_RESULTS_YV, "PHASE3_4_observed_holdout_summaries.rds")
FIT_PATHS <- c(
  H1 = file.path(DIR_FITS_YV, "ce_tau_year_training_H1.rds"),
  H2 = file.path(DIR_FITS_YV, "ce_tau_year_training_H2.rds")
)
required_inputs <- c(DIAG_RDS, ORIG_RDS, OBS_RDS, FIT_PATHS,
                     file.path(DIR_TABLES_RD, "PHASE0_design_audit.csv"),
                     file.path(DIR_TABLES_RD, "PHASE0_climatology_leakage_magnitude.csv"),
                     file.path(DIR_TABLES_RD, "PHASE1_integrity_check_vs_original.csv"),
                     file.path(DIR_TABLES_HMC, "HMC_AUDIT_model_level.csv"))
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Missing required Step 6 input(s): ", paste(missing_inputs, collapse = ", "))

source(file.path(DIR_YV, "00_config_year_effect_validation.R"))

# Fixed before inspecting outcome classifications. This is a mechanistic,
# validation-independent establishment rule, not a calibrated epidemic rule:
# local-transmission hazard must exceed introduction-related hazard for at least
# four consecutive weeks, with positive latent infections, following at least
# one week with Reff>1 and M_t>0 in the same calendar year.
ESTABLISHMENT_RUN_WEEKS <- 4L

# Reuse the project's already accepted major-epidemic thresholds. They were not
# selected for this Step 6 analysis or tuned to the 2022 observed total.
MAJOR_ABS_CASE_MINIMUM <- WAVE_DETECTION_YV$major_abs_case_minimum
MAJOR_INCIDENCE_MINIMUM_PER_100K <- WAVE_DETECTION_YV$major_incidence_minimum_per_100k

PROVISIONAL_NOTE <- paste(
  "PROVISIONAL: H1/H2 have low but nonzero divergences and incomplete mixing;",
  "tau_year upper-tail reliability is not fully established."
)

