# STEP 7 -- importation / persistence sensitivity: isolated, read-only
# configuration. Reuses the H1/H2 training-only posteriors (NO Stan
# refitting) and the pure-R future-renewal simulator UNMODIFIED. All new
# outputs live under step7_importation_persistence/. PROVISIONAL per the
# HMC reliability audit.

if (!exists("ROOT")) stop("ROOT not found -- source CHIK_CLIM/00_project_setup.R first.")

DIR_STEP6 <- file.path(ROOT, "02_Script/02_ceara_pipeline/15_step6_retrospective_failure_diagnosis")
DIR_YV    <- file.path(ROOT, "02_Script/02_ceara_pipeline/12_future_year_effect_validation")

DIR_FITS_YV    <- file.path(ROOT, "03_Output/02_ceara_pipeline/model_fits/future_year_effect_validation")
DIR_TABLES_STEP6  <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/step6_retrospective_failure_diagnosis")
DIR_RESULTS_STEP6 <- file.path(ROOT, "03_Output/02_ceara_pipeline/results/step6_retrospective_failure_diagnosis")

DIR_SCRIPTS_STEP7 <- file.path(ROOT, "02_Script/02_ceara_pipeline/16_step7_importation_persistence")
DIR_REPORTS_STEP7 <- file.path(ROOT, "03_Output/02_ceara_pipeline/reports/step7_importation_persistence")
DIR_TABLES_STEP7  <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/step7_importation_persistence")
DIR_FIGURES_STEP7 <- file.path(ROOT, "03_Output/02_ceara_pipeline/figures/step7_importation_persistence")
DIR_RESULTS_STEP7 <- file.path(ROOT, "03_Output/02_ceara_pipeline/results/step7_importation_persistence")
for (d in c(DIR_REPORTS_STEP7, DIR_TABLES_STEP7, DIR_FIGURES_STEP7, DIR_RESULTS_STEP7)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

FIT_PATHS <- c(H1 = file.path(DIR_FITS_YV, "ce_tau_year_training_H1.rds"),
               H2 = file.path(DIR_FITS_YV, "ce_tau_year_training_H2.rds"))
required_inputs <- c(FIT_PATHS, file.path(DIR_TABLES_STEP6, "STEP6_2022_headline_mechanistic_metrics.csv"),
                      file.path(DIR_TABLES_STEP6, "STEP6_entry2022_infection_history.csv"),
                      file.path(DIR_TABLES_STEP6, "STEP6_2022_exceedance_probability_mc.csv"))
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Missing required Step 7 input(s): ", paste(missing_inputs, collapse = ", "))

source(file.path(DIR_YV, "00_config_year_effect_validation.R")) # HOLDOUTS, QVAL, ACCEPTED_ANOMALY_CSV, generation_weights_future(), etc.
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R")) # build_state_weekly()

# ---- Step 6 reference metrics (for the Phase 0 reproduction check) ----
STEP6_REFERENCE <- list(
  H1 = list(entry_weighted_infectiousness_median = 4.77, S_prop_median = 0.741, peak_Reff_median = 1.41,
            prob_major_epidemic = 0.1055, prob_exceed_2022 = 0.003),
  H2 = list(entry_weighted_infectiousness_median = 886.4, S_prop_median = 0.744, peak_Reff_median = 1.40,
            prob_exceed_2022 = 0.119)
)
REPRODUCTION_TOLERANCE <- 0.03 # relative tolerance (3%) for the Phase 0 numerical-agreement check; the comparison is diagnostic, not a gate that silently replaces the baseline on failure.

# ---- Prespecified importation grid (fixed here, BEFORE any arm is run; never tuned to 2022 outcomes) ----
LAMBDA_GRID <- c(I0 = 0, I1 = 0.25, I2 = 1, I3 = 4)
LAMBDA_BASELINE <- 1 # matches the accepted fit's own fixed imports_per_week and Step 6's reference arm

# ---- Ensemble sizes (Phase 6: smoke test first, then >=100 draws x 20 replicates per main arm) ----
ENSEMBLE_STEP7 <- list(smoke = list(n_draws = 10L, n_replicates = 5L), full = list(n_draws = 100L, n_replicates = 20L))
BASE_SEED_STEP7 <- 20261013L + 7000L # distinct stream from prior modules' BASE_SEED_YV (+777) and the HMC-audit module's offsets

PROVISIONAL_NOTE <- paste0(
  "PROVISIONAL: H1/H2 training-only posteriors have low but nonzero HMC divergences\n",
  "and incomplete within-chain mixing (tau_year upper tail); all results below are conditional on this."
)
KNOWN_LIMITATIONS <- c(
  "Future year effects are re-centred (sum to zero) within each finite holdout horizon, shrinking the marginal SD from tau_year to tau_year*sqrt((K-1)/K) -- quantified in Step 6, not corrected here except in the explicit, isolated Phase 5 arm.",
  "Climate anomaly standardisation (z_T_anom/z_P_anom) uses scaling constants estimated over the FULL 2015-2025 period, not training-only -- a small, disclosed, uncorrected preprocessing leakage common to every arm in this module."
)

message(sprintf("[step7 config] ready. Lambda grid: %s | ensembles: smoke=%dx%d, full=%dx%d",
                 paste(names(LAMBDA_GRID), LAMBDA_GRID, sep = "=", collapse = ", "),
                 ENSEMBLE_STEP7$smoke$n_draws, ENSEMBLE_STEP7$smoke$n_replicates,
                 ENSEMBLE_STEP7$full$n_draws, ENSEMBLE_STEP7$full$n_replicates))
