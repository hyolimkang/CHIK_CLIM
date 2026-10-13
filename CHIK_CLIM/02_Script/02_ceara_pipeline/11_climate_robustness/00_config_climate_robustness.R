# Climate robustness assessment: central configuration.
#
# Historical inference only (2015-2025). No future SSP scenarios, no
# vaccination, no new epidemic generator. Preserves all accepted
# historical models, fitted objects and output files -- new Stan fits
# (M1-M4) are saved under a separate model_fits/climate_robustness/ tree
# and M0 is READ (never re-fit) from the existing accepted location.

if (!exists("ROOT")) stop("ROOT not found -- source CHIK_CLIM/00_project_setup.R first (opening CHIK_CLIM.Rproj does this automatically).")

# ---- Read-only accepted upstream sources ----
ACCEPTED_CLIMATE_FIT_M0 <- file.path(ROOT, "03_Output/02_ceara_pipeline/model_fits/climate_forced/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds")
ACCEPTED_ANOMALY_CSV    <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/climate_forced_v4_9/climate_anomaly_covariates.csv")
ACCEPTED_SCALING_CSV    <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/climate_forced_v4_9/climate_anomaly_scaling_metadata.csv")
ACCEPTED_LAGGED_RAW_CSV <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/climate_forced_v4_9/climate_lagged_raw.csv")
NATIONAL_CLIMATE_CSV    <- file.path(ROOT, "03_Output/07_national_pipeline/tables/national_wave_analysis/brazil_chik_uf_weekly_climate.csv")
STAN_SOURCE_CLIMATE_M0  <- file.path(ROOT, "02_Script/00_shared/stan/current/ceara/renewal_ceara_v4_9_climate_forced_canary.stan")
TD14_FIT_PATH           <- file.path(ROOT, "03_Output/02_ceara_pipeline/model_fits/initialization/renewal_ceara_v4_0_minimal_no_vaccine_fit_td14.rds")

for (f in c(ACCEPTED_CLIMATE_FIT_M0, ACCEPTED_ANOMALY_CSV, ACCEPTED_SCALING_CSV, ACCEPTED_LAGGED_RAW_CSV, NATIONAL_CLIMATE_CSV, STAN_SOURCE_CLIMATE_M0, TD14_FIT_PATH)) {
  if (!file.exists(f)) stop("Missing required upstream input, cannot proceed: ", f)
}

# ---- This module's own output locations (new, additive only) ----
DIR_SCRIPTS_CR <- file.path(ROOT, "02_Script/02_ceara_pipeline/11_climate_robustness")
DIR_REPORTS_CR <- file.path(ROOT, "03_Output/02_ceara_pipeline/reports/climate_robustness")
DIR_TABLES_CR  <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/climate_robustness")
DIR_FIGURES_CR <- file.path(ROOT, "03_Output/02_ceara_pipeline/figures/climate_robustness")
DIR_FITS_CR    <- file.path(ROOT, "03_Output/02_ceara_pipeline/model_fits/climate_robustness")
for (d in c(DIR_REPORTS_CR, DIR_TABLES_CR, DIR_FIGURES_CR, DIR_FITS_CR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Fixed scientific settings, byte-identical to M0's accepted fit ----
QVAL <- 0.05
SEED_START_DATE <- as.Date("2022-01-02")
SEED_LENGTH <- 8L
FIXED <- list(date_start = as.Date("2015-01-04"), date_end = as.Date("2025-12-21"),
              imports_per_week = 1, year_effect_prior_sd = 0.40,
              sero_pos = 103L, sero_n = 404L, kappa_sero = 50)
HMC_SETTINGS <- list(warmup = 800L, sampling = 400L, chains = 4L, adapt_delta = 0.95, max_treedepth = 14L, metric = "dense_e")
PRIMARY_STATES <- c("CE", "BA", "RJ", "MT")

# ---- The 5 model specifications (M0 reused, M1-M4 newly fit) ----
MODEL_SPECS <- list(
  M0 = list(lag_spec = "primary",       year_effect_prior_sd = 0.40, label = "M0: primary lag, prior SD=0.40 (ACCEPTED, reused)"),
  M1 = list(lag_spec = "sensitivity_A", year_effect_prior_sd = 0.40, label = "M1: sensitivity_A lag, prior SD=0.40"),
  M2 = list(lag_spec = "sensitivity_B", year_effect_prior_sd = 0.40, label = "M2: sensitivity_B lag, prior SD=0.40"),
  M3 = list(lag_spec = "primary",       year_effect_prior_sd = 0.20, label = "M3: primary lag, prior SD=0.20"),
  M4 = list(lag_spec = "primary",       year_effect_prior_sd = 0.60, label = "M4: primary lag, prior SD=0.60")
)

message(sprintf("[climate_robustness config] %d model specifications (M0 reused + %d new fits)", length(MODEL_SPECS), length(MODEL_SPECS) - 1))
