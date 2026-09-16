# Bahia data audit -- run BEFORE any fitting, per the external-replication
# design (BAHIA_V4_9_FROZEN_MODEL_ASSESSMENT.md Section 3). Checks the same
# identities/continuity properties that make_v4_3_data() checks for Ceara.

required_packages <- c("here", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
source(file.path(root, "02_Script/40_renewal_model/00_shared/01_build_ceara_state_weekly.R"))

DATE_START <- as.Date("2015-01-04")
DATE_END <- as.Date("2025-12-21") # matches v4.9's Ceara window exactly

panel <- readRDS(file.path(root, "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
weekly_cases <- build_state_weekly(panel, "29") |> filter(week_start >= DATE_START, week_start <= DATE_END)
weekly_demography <- readRDS(file.path(root, "01_Data/bahia_weekly_demography_2015_2025.rds")) |>
  filter(week_start >= DATE_START, week_start <= DATE_END)

report <- list()
report$n_weeks_cases <- nrow(weekly_cases)
report$n_weeks_demography <- nrow(weekly_demography)
report$date_range_cases <- paste(range(weekly_cases$week_start), collapse = " to ")
report$date_range_demography <- paste(range(weekly_demography$week_start), collapse = " to ")

expected_dates <- seq(DATE_START, DATE_END, by = "week")
report$n_expected_weeks <- length(expected_dates)
report$cases_dates_continuous <- identical(as.Date(weekly_cases$week_start), expected_dates)
report$demography_dates_continuous <- identical(as.Date(weekly_demography$week_start), expected_dates)
report$cases_duplicate_weeks <- sum(duplicated(weekly_cases$week_start))
report$demography_duplicate_weeks <- sum(duplicated(weekly_demography$week_start))

weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
report$n_weeks_joined <- nrow(weekly)
report$join_matches_both_inputs <- nrow(weekly) == nrow(weekly_cases) && nrow(weekly) == nrow(weekly_demography)

accounting_error <- weekly$N_end - weekly$N_start - weekly$births + weekly$all_cause_deaths - weekly$net_population_reconciliation
report$max_abs_demographic_accounting_error <- max(abs(accounting_error))

report$total_cases <- sum(weekly_cases$cases)
report$population_range <- paste(round(range(weekly$N_start)), collapse = " to ")
report$mean_weekly_cases <- mean(weekly_cases$cases)
report$zero_case_weeks <- sum(weekly_cases$cases == 0)
report$max_weekly_cases <- max(weekly_cases$cases)
report$max_weekly_cases_week <- as.character(weekly_cases$week_start[which.max(weekly_cases$cases)])

cat("=== Bahia data audit ===\n")
for (nm in names(report)) cat(sprintf("%-38s: %s\n", nm, report[[nm]]))

audit_df <- tibble::tibble(check = names(report), value = as.character(unlist(report)))
out_dir <- file.path(root, "03_Output/tables/renewal_bahia_v4_9_replication")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(audit_df, file.path(out_dir, "bahia_data_audit.csv"))
message("\n[saved] ", file.path(out_dir, "bahia_data_audit.csv"))

stopifnot(
  "cases dates not continuous" = report$cases_dates_continuous,
  "demography dates not continuous" = report$demography_dates_continuous,
  "duplicate weeks in cases" = report$cases_duplicate_weeks == 0,
  "duplicate weeks in demography" = report$demography_duplicate_weeks == 0,
  "join dropped rows" = report$join_matches_both_inputs,
  "demographic accounting identity violated" = report$max_abs_demographic_accounting_error < 1e-6
)
message("[audit] ALL CHECKS PASSED")
