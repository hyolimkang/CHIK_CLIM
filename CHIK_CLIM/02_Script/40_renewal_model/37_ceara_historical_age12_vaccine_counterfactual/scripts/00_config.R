# Historical age-12 routine vaccination counterfactual (Ceara, 2015-2021).
# Central configuration -- ALL new outputs stay under this analysis
# directory. Does NOT modify any frozen v4.9 / climate-canary / age-
# reconstruction file or output.

ROOT <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
ANALYSIS_DIR <- file.path(ROOT, "02_Script/40_renewal_model/37_ceara_historical_age12_vaccine_counterfactual")
DIR_SCRIPTS <- file.path(ANALYSIS_DIR, "scripts")
DIR_RESULTS <- file.path(ANALYSIS_DIR, "results")
DIR_TABLES <- file.path(ANALYSIS_DIR, "tables")
DIR_DIAGNOSTICS <- file.path(ANALYSIS_DIR, "diagnostics")
DIR_REPORTS <- file.path(ANALYSIS_DIR, "reports")
# Figures follow the project-wide convention (03_Output/figures/<analysis>/),
# not the script-local tree -- moved here so all analyses are discoverable
# in one place; consolidated from 37_.../figures/ on 2026-09-17.
DIR_FIG_ROOT <- file.path(ROOT, "03_Output/figures/ceara_historical_age12_vaccine_counterfactual")
DIR_FIG_MAIN <- file.path(DIR_FIG_ROOT, "main")
DIR_FIG_SUPP <- file.path(DIR_FIG_ROOT, "supplementary")
DIR_FIG_NNV <- file.path(DIR_FIG_ROOT, "nnv")
for (d in c(DIR_RESULTS, DIR_TABLES, DIR_FIG_MAIN, DIR_FIG_SUPP, DIR_FIG_NNV, DIR_DIAGNOSTICS, DIR_REPORTS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Upstream (READ-ONLY) accepted sources -- never written to by this analysis.
ACCEPTED_CLIMATE_FIT <- file.path(ROOT, "02_Script/40_renewal_model/36_climate_forced_v4_9/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds")
ACCEPTED_AGE_SHARES <- file.path(ROOT, "03_Output/tables/climate_forced_v4_9/CE_age_population_shares.csv")
ACCEPTED_SIMULATOR <- file.path(ROOT, "02_Script/40_renewal_model/36_climate_forced_v4_9/scripts/07_age_cohort_simulator.R")

# Historical analysis window: 2015 through 2021 INCLUSIVE, stopping before
# the externally conditioned 2022 seed period (per instruction).
WINDOW_START <- as.Date("2015-01-04")
WINDOW_END <- as.Date("2021-12-26") # last week fully within calendar year 2021

MAX_AGE <- 100L
TARGET_AGE <- 12L

# Scenario 0 / Scenario 1 (mechanism stress test -- NOT a policy value)
SCENARIOS <- list(
  novaccine = list(vaccination_on = FALSE),
  stresstest = list(vaccination_on = TRUE, coverage = 1.0, ve_infection = 1.0)
)

N_DRAWS <- 200L # "all usable draws unless computationally prohibitive" -- see reports/ for the justification
RANDOM_SEED <- 20260920L

REPORTING_AGE_BREAKS <- c(0, 12, 13, 18, 65, MAX_AGE + 1)
REPORTING_AGE_LABELS <- c("0-11", "12", "13-17", "18-64", "65+")

message(sprintf("[config] Analysis dir: %s", ANALYSIS_DIR))
message(sprintf("[config] Window: %s to %s | N_DRAWS=%d | seed=%d", WINDOW_START, WINDOW_END, N_DRAWS, RANDOM_SEED))
