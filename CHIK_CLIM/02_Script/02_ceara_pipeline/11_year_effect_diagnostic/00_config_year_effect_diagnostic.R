# Year-effect validation diagnostic: central configuration.
#
# ONE targeted diagnostic step following up on the future_projection_pilot
# (Steps 1-2): why does Y1 (resampled year effects) produce ~3.26M vs Y0's
# ~8,100 median 10-year cumulative infections? This module does NOT modify
# or re-run Steps 1-2, does NOT refit Stan, and does NOT proceed to SSP
# projection or vaccination. It only READS the existing accepted fit and
# SOURCES (unmodified) the existing Step 2 simulation functions.
#
# Writes to its own *_diagnostic/ output subfolders only.

if (!exists("ROOT")) stop("ROOT not found -- source CHIK_CLIM/00_project_setup.R first (opening CHIK_CLIM.Rproj does this automatically).")

DIR_PILOT <- file.path(ROOT, "02_Script/02_ceara_pipeline/10_future_projection_pilot")
source(file.path(DIR_PILOT, "00_config_future_pilot.R"))          # read-only: paths, ACCEPTED_CLIMATE_FIT, generation_weights_future(), WEEKS_PER_YEAR, build_future_dates()
source(file.path(DIR_PILOT, "03_future_renewal_functions.R"))     # read-only: simulate_renewal_core(), compute_R0_future(), advance_terminal_state(), build_future_demography(), detect_simulated_waves(), simulate_reported_cases()

DIR_SCRIPTS_DIAG <- file.path(ROOT, "02_Script/02_ceara_pipeline/11_year_effect_diagnostic")
DIR_REPORTS_DIAG  <- file.path(ROOT, "03_Output/02_ceara_pipeline/reports/year_effect_diagnostic")
DIR_TABLES_DIAG   <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/year_effect_diagnostic")
DIR_FIGURES_DIAG  <- file.path(ROOT, "03_Output/02_ceara_pipeline/figures/year_effect_diagnostic")
DIR_RESULTS_DIAG  <- file.path(ROOT, "03_Output/02_ceara_pipeline/results/year_effect_diagnostic")
for (d in c(DIR_REPORTS_DIAG, DIR_TABLES_DIAG, DIR_FIGURES_DIAG, DIR_RESULTS_DIAG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Four structural arms for the paired decomposition (Section 3) ----
# A: neutral: year_effect=0, A_year=1
# B: year-level variation only: year_effect=Y1-sampled, A_year=1
# C: seasonal-amplitude variation only: year_effect=0, A_year=Y1-sampled
# D: both (== the existing Y1 arm)
DECOMP_ARMS <- c("A_neutral", "B_year_only", "C_amplitude_only", "D_both")
DECOMP_ARM_LABELS <- c(
  A_neutral = "A: neutral (year_effect=0, A_year=1)",
  B_year_only = "B: year-level only (year_effect=Y1, A_year=1)",
  C_amplitude_only = "C: seasonal-amplitude only (year_effect=0, A_year=Y1)",
  D_both = "D: both (= Y1)"
)

DIAG_LAMBDA_IMPORT <- 1          # fixed, per instruction, to isolate year-effect mechanisms
DIAG_CLIMATE_ARM <- "C0_stationary_neutral" # fixed per instruction

DIAG_ENSEMBLE <- list(
  smoke = list(n_draws = 10L, n_replicates = 5L),
  full  = list(n_draws = 100L, n_replicates = 20L)
)
DIAG_BASE_SEED <- 20261012L + 7L # distinct stream from the Step-2 BASE_SEED, still reproducible

message(sprintf("[year_effect_diagnostic config] 4 arms x lambda_import=%s x climate=%s | future horizon %d weeks",
                 DIAG_LAMBDA_IMPORT, DIAG_CLIMATE_ARM, length(build_future_dates())))
