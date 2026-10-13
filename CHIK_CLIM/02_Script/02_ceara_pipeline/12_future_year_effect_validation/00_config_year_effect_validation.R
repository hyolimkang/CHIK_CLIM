# Future year-effect validation: central configuration.
#
# Historical inference + retrospective (held-out) validation only. NO
# SSP/ISIMIP forcing, NO vaccination. Does not modify any accepted
# model, fit, or previous diagnostic/robustness output -- all new fits
# and outputs live under this module's own *_year_effect_validation/
# trees.

if (!exists("ROOT")) stop("ROOT not found -- source CHIK_CLIM/00_project_setup.R first (opening CHIK_CLIM.Rproj does this automatically).")

# ---- Read-only accepted upstream sources ----
ACCEPTED_CLIMATE_FIT_M0 <- file.path(ROOT, "03_Output/02_ceara_pipeline/model_fits/climate_forced/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds")
ACCEPTED_ANOMALY_CSV    <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/climate_forced_v4_9/climate_anomaly_covariates.csv")
TD14_FIT_PATH           <- file.path(ROOT, "03_Output/02_ceara_pipeline/model_fits/initialization/renewal_ceara_v4_0_minimal_no_vaccine_fit_td14.rds")
STAN_SOURCE_TAU_YEAR    <- file.path(ROOT, "02_Script/00_shared/stan/current/ceara/renewal_ceara_v4_9_climate_forced_tau_year.stan")
STAN_SOURCE_M0          <- file.path(ROOT, "02_Script/00_shared/stan/current/ceara/renewal_ceara_v4_9_climate_forced_canary.stan")
DIR_PILOT               <- file.path(ROOT, "02_Script/02_ceara_pipeline/10_future_projection_pilot")

for (f in c(ACCEPTED_CLIMATE_FIT_M0, ACCEPTED_ANOMALY_CSV, TD14_FIT_PATH, STAN_SOURCE_TAU_YEAR, STAN_SOURCE_M0)) {
  if (!file.exists(f)) stop("Missing required upstream input, cannot proceed: ", f)
}

# Reuse the Step-2 future-renewal simulation primitives, UNMODIFIED
source(file.path(DIR_PILOT, "00_config_future_pilot.R"))
source(file.path(DIR_PILOT, "03_future_renewal_functions.R"))

# ---- This module's own output locations ----
DIR_SCRIPTS_YV <- file.path(ROOT, "02_Script/02_ceara_pipeline/12_future_year_effect_validation")
DIR_REPORTS_YV <- file.path(ROOT, "03_Output/02_ceara_pipeline/reports/future_year_effect_validation")
DIR_TABLES_YV  <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/future_year_effect_validation")
DIR_FIGURES_YV <- file.path(ROOT, "03_Output/02_ceara_pipeline/figures/future_year_effect_validation")
DIR_RESULTS_YV <- file.path(ROOT, "03_Output/02_ceara_pipeline/results/future_year_effect_validation")
DIR_FITS_YV    <- file.path(ROOT, "03_Output/02_ceara_pipeline/model_fits/future_year_effect_validation")
for (d in c(DIR_REPORTS_YV, DIR_TABLES_YV, DIR_FIGURES_YV, DIR_RESULTS_YV, DIR_FITS_YV)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Fixed scientific settings, byte-identical to the accepted fit ----
QVAL <- 0.05
SEED_START_DATE_FULL <- as.Date("2022-01-02")
SEED_LENGTH_FULL <- 8L
FIXED_FULL <- list(date_start = as.Date("2015-01-04"), date_end = as.Date("2025-12-21"),
                    imports_per_week = 1, sero_pos = 103L, sero_n = 404L, kappa_sero = 50)
HMC_SETTINGS <- list(warmup = 800L, sampling = 400L, chains = 4L, adapt_delta = 0.95, max_treedepth = 14L, metric = "dense_e")

# ---- Phase 1: tau_year prior sensitivity grid (full 2015-2025 period, retains 2022 seed) ----
TAU_PRIORS <- list(P1 = 0.20, P2 = 0.40, P3 = 0.80)

# ---- Phase 3: retrospective holdout experiments ----
# Training windows deliberately EXCLUDE the 2022 seed date entirely (H1
# ends 2019, H2 ends 2021) -- no seed-forcing is needed or used in EITHER
# training fit; a dedicated seed-free training-data builder is used
# throughout (04_build_training_only_data.R), never the accepted
# seed-requiring builder.
HOLDOUTS <- list(
  H1 = list(train_start = as.Date("2015-01-04"), train_end = as.Date("2019-12-29"),
            holdout_start = as.Date("2020-01-05"), holdout_end = as.Date("2025-12-21"), label = "H1: train 2015-2019, holdout 2020-2025"),
  H2 = list(train_start = as.Date("2015-01-04"), train_end = as.Date("2021-12-26"),
            holdout_start = as.Date("2022-01-02"), holdout_end = as.Date("2025-12-21"), label = "H2: train 2015-2021, holdout 2022-2025")
)
# Serology (2018-06-01 to 2018-12-31 survey window) precedes BOTH training cutoffs -- included in both, as permitted.

IMPORTATION_SENSITIVITY_GRID <- c(0.25, 1, 4) # weekly Poisson rate, prespecified, NOT tuned to holdout outcomes

ENSEMBLE_RETRO <- list(smoke = list(n_draws = 10L, n_replicates = 5L), full = list(n_draws = 100L, n_replicates = 20L))
BASE_SEED_YV <- 20261013L

WAVE_DETECTION_YV <- WAVE_DETECTION # reuse the Step-2 adapted peak-trough definition unchanged, for consistency with earlier modules

message(sprintf("[year_effect_validation config] tau priors: %s | holdouts: %s",
                 paste(names(TAU_PRIORS), unlist(TAU_PRIORS), sep = "=", collapse = ", "), paste(names(HOLDOUTS), collapse = ", ")))
