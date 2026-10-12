# Future projection pilot (Steps 1-2 only): central configuration.
#
# Scope reminder (see README.md in this folder): Step 1 re-analyses the
# ALREADY-FITTED climate-forced v4.9 canary posterior (no refit). Step 2 is a
# pure-R stochastic simulation engine for 2026-2035 under stationary/
# resampled historical climate (NOT SSP/CMIP6) and NO vaccination. Nothing
# under 02_ceara_pipeline/{01..09} is read-write here except as a read-only
# accepted input; every output of this module lives under the
# *_pilot/ subfolders created below.

if (!exists("ROOT")) stop("ROOT not found -- source CHIK_CLIM/00_project_setup.R first (opening CHIK_CLIM.Rproj does this automatically).")

# ---- Read-only accepted upstream sources (never written to) ----
ACCEPTED_CLIMATE_FIT   <- file.path(ROOT, "03_Output/02_ceara_pipeline/model_fits/climate_forced/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds")
ACCEPTED_BASELINE_FIT  <- file.path(ROOT, "03_Output/02_ceara_pipeline/model_fits/baseline/v4_9_hierarchical_seasonality/outputs/q0.05/renewal_ceara_v4_9_fit_q0.05.rds")
ACCEPTED_ANOMALY_CSV   <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/climate_forced_v4_9/climate_anomaly_covariates.csv")
ACCEPTED_SCALING_CSV   <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/climate_forced_v4_9/climate_anomaly_scaling_metadata.csv")
ACCEPTED_IBGE_UF_PROJ  <- file.path(ROOT, "01_Data/ibge_population_projection_uf_2024revision.rds")
STAN_SOURCE_CLIMATE    <- file.path(ROOT, "02_Script/00_shared/stan/current/ceara/renewal_ceara_v4_9_climate_forced_canary.stan")

for (f in c(ACCEPTED_CLIMATE_FIT, ACCEPTED_BASELINE_FIT, ACCEPTED_ANOMALY_CSV, ACCEPTED_SCALING_CSV, ACCEPTED_IBGE_UF_PROJ)) {
  if (!file.exists(f)) stop("Missing required upstream input, cannot proceed: ", f)
}

# ---- This module's own output locations (new, additive only) ----
DIR_SCRIPTS    <- file.path(ROOT, "02_Script/02_ceara_pipeline/10_future_projection_pilot")
DIR_REPORTS    <- file.path(ROOT, "03_Output/02_ceara_pipeline/reports/future_projection_pilot")
DIR_TABLES     <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/future_projection_pilot")
DIR_FIGURES    <- file.path(ROOT, "03_Output/02_ceara_pipeline/figures/future_projection_pilot")
DIR_RESULTS    <- file.path(ROOT, "03_Output/02_ceara_pipeline/results/future_projection_pilot")
for (d in c(DIR_REPORTS, DIR_TABLES, DIR_FIGURES, DIR_RESULTS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Ceara identifiers ----
UF_CODE_CEARA <- 23L
STATE_ABBR_CEARA <- "CE"

# ---- Historical window (must match the accepted fit exactly) ----
HIST_DATE_START <- as.Date("2015-01-04")
HIST_DATE_END   <- as.Date("2025-12-21")

# ---- Future test horizon ----
FUTURE_YEAR_START <- 2026L
FUTURE_YEAR_END   <- 2035L
# First future week continues the exact historical weekly sequence (+7 days);
# trimmed below to the first/last ISO week whose calendar year falls in
# [FUTURE_YEAR_START, FUTURE_YEAR_END].
build_future_dates <- function() {
  candidate <- seq(HIST_DATE_START, by = "week", length.out = 1200)
  candidate <- candidate[candidate > HIST_DATE_END]
  yrs <- as.integer(format(candidate, "%Y"))
  candidate[yrs >= FUTURE_YEAR_START & yrs <= FUTURE_YEAR_END]
}

# ---- Generation interval / seasonal period constants (byte-identical to the accepted model -- DO NOT change) ----
WEEKS_PER_YEAR <- 52.1775
generation_weights_future <- function() {
  w <- diff(pgamma(0:8, shape = 4, rate = 2))
  w / sum(w)
}

# ---- Structural arms (Section 2.2-2.4 of the task) ----
YEAR_EFFECT_ARMS <- c("Y0_reference", "Y1_resampled_distributional")
CLIMATE_ARMS     <- c("C0_stationary_neutral", "C1_resampled_weather_variability", "C2_resampled_weather_zero_beta")
DEMOGRAPHY_ARMS  <- c("D1_ibge_projection", "D0_stationary_2025_levels")

# Weekly external introduction-pressure grid (Poisson rate, M_t ~ Poisson(lambda_import)).
# lambda_import = 1 matches the accepted fit's fixed imports_per_week. This
# grid is a prespecified STRUCTURAL SENSITIVITY range, not an estimated
# parameter -- lambda_import has NOT been fit to any data.
IMPORTATION_GRID <- c(0, 0.25, 1, 4)
IMPORTATION_REFERENCE <- 1

# ---- Simulation ensemble sizes (Section 2.8) ----
ENSEMBLE_SETTINGS <- list(
  smoke = list(n_draws = 10L,  n_replicates = 5L),
  pilot = list(n_draws = 100L, n_replicates = 20L)
)

# ---- Reproducibility ----
BASE_SEED <- 20261012L

# ---- Outbreak/wave detection settings, adapted from the accepted national
# wave-census "primary" definition (07_national_pipeline/05_national_wave_census_v2);
# applied to SIMULATED REPORTED CASES only, never to latent infections directly. ----
WAVE_DETECTION <- list(
  smooth_weeks = 3L, trough_fraction = 0.30, trough_persistence = 3L,
  onset_fraction = 0.05, peak_neighborhood_weeks = 3L,
  major_abs_case_minimum = 100, major_incidence_minimum_per_100k = 5
)

message(sprintf(
  "[future_pilot config] future horizon: %d weeks (%s to %s) | %d year-effect arms x %d climate arms x %d demography arms x %d importation levels",
  length(build_future_dates()), format(min(build_future_dates())), format(max(build_future_dates())),
  length(YEAR_EFFECT_ARMS), length(CLIMATE_ARMS), length(DEMOGRAPHY_ARMS), length(IMPORTATION_GRID)
))
