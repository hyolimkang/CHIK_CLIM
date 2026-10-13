# Year-effect diagnostic -- standalone runner.
#
# ONE targeted diagnostic following up on the future_projection_pilot
# (Steps 1-2): why does Y1 produce ~8,100 (Y0) vs ~3.26M (Y1) median
# 10-year cumulative infections? Does NOT modify Steps 1-2, does NOT
# refit Stan, does NOT proceed to SSP projection or vaccination.
#
# Usage:
#   source("02_Script/02_ceara_pipeline/11_year_effect_diagnostic/00_run_year_effect_diagnostic.R")
#
# Outputs under 03_Output/02_ceara_pipeline/{reports,tables,figures,results}/year_effect_diagnostic/

if (!exists("ROOT")) stop("ROOT not found -- open CHIK_CLIM.Rproj or source CHIK_CLIM/00_project_setup.R first.")
DIR_DIAG <- file.path(ROOT, "02_Script/02_ceara_pipeline/11_year_effect_diagnostic")

message("=== Year-effect diagnostic Part 2: historical vs. future variability ===")
source(file.path(DIR_DIAG, "01_historical_vs_future_variability.R"))

message("\n=== Year-effect diagnostic Part 3/4: paired 4-arm simulation ===")
source(file.path(DIR_DIAG, "02_paired_arm_simulation.R"))

message("\n=== Year-effect diagnostic Part 5: figures ===")
source(file.path(DIR_DIAG, "03_make_diagnostic_figures.R"))

message("\n=== Year-effect diagnostic complete. See 03_Output/02_ceara_pipeline/reports/year_effect_diagnostic/FUTURE_YEAR_EFFECT_VALIDATION.md ===")
