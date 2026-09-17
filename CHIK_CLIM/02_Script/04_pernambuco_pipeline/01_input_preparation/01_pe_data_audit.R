# Pernambuco external replication -- Section 2 & 3 of the design spec.
# Build and audit the PE state-weekly case series (existing validated
# municipality/week panel, UF prefix "26"), and determine whether PE
# qualifies for a conditioned reintroduction seed under the SAME
# pre-specified long-gap rule already used for Ceara/Bahia.

required_packages <- c("here", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
setwd(root)
source(file.path(root, "02_Script/00_shared/functions/01_build_ceara_state_weekly.R"))

table_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/tables/pernambuco_v4_9_replication")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

DATE_START <- as.Date("2015-01-04")
DATE_END <- as.Date("2025-12-21") # same window convention as Ceara/Bahia

panel <- readRDS(file.path(root, "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
weekly_cases <- build_state_weekly(panel, "26") |> dplyr::filter(week_start >= DATE_START, week_start <= DATE_END)

demography_path <- file.path(root, "01_Data/pernambuco_weekly_demography_2015_2025.rds")
if (!file.exists(demography_path)) stop("Pernambuco weekly demography not yet built -- run 09c/10c in 02_Script/00_data_prep first.")
weekly_demography <- readRDS(demography_path) |> dplyr::filter(week_start >= DATE_START, week_start <= DATE_END)

weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch -- PE case series and demography weeks do not align.")

# ---- Audit -----------------------------------------------------------------
weeks_expected <- seq(DATE_START, DATE_END, by = "week")
accounting_error <- max(abs(weekly$N_end - weekly$N_start - weekly$births + weekly$all_cause_deaths - weekly$net_population_reconciliation))

audit <- tibble(
  state = "PE", n_weeks = nrow(weekly), n_weeks_expected = length(weeks_expected),
  week_continuity_ok = identical(weekly$week_start, weeks_expected),
  n_zero_case_weeks = sum(weekly$cases == 0), n_positive_case_weeks = sum(weekly$cases > 0),
  total_cases = sum(weekly$cases), n_na_any = sum(!complete.cases(weekly)),
  N_start_2015 = weekly$N_start[1], N_end_2025 = weekly$N_end[nrow(weekly)],
  year_min = min(weekly$year), year_max = max(weekly$year),
  demographic_accounting_max_abs_error = accounting_error,
  data_source_confirmed = "chik_dlnm_panel_muni_week_2015_2025.rds (UF prefix 26) + pernambuco_weekly_demography_2015_2025.rds -- NOT inherited Bahia objects/paths"
)

message("=== Pernambuco v4.9 input audit ===")
print(as.data.frame(audit))
write_csv(audit, file.path(table_dir, "pernambuco_v49_input_audit.csv"))
message("[saved] ", file.path(table_dir, "pernambuco_v49_input_audit.csv"))

if (!audit$week_continuity_ok || audit$n_na_any > 0 || accounting_error > 1e-6) {
  stop("STOP: Pernambuco input audit failed -- do not proceed to fitting.")
}
saveRDS(weekly, file.path(table_dir, "pernambuco_weekly_input.rds"))
message("\n[saved] ", file.path(table_dir, "pernambuco_weekly_input.rds"), " (audited weekly input for Stan)")

# ---- Section 3: seed-qualification check, SAME rule as Ceara/Bahia --------
wave_census <- read_csv(file.path(root, "03_Output/07_national_pipeline/tables/national_wave_census_v2/brazil_chik_wave_census_v2.csv"), show_col_types = FALSE)
pe_major_waves <- wave_census |> dplyr::filter(state == "PE", major_epidemic_primary == TRUE) |>
  transmute(wave_id, wave_order,
            onset_week = as.Date(onset_week, format = "%m/%d/%Y"),
            end_week = as.Date(end_week, format = "%m/%d/%Y")) |>
  arrange(onset_week)

gaps <- pe_major_waves |>
  mutate(previous_end = lag(end_week), weeks_since_previous_wave = as.numeric(onset_week - previous_end) / 7) |>
  dplyr::filter(!is.na(weeks_since_previous_wave))

message("\n=== Pernambuco major-wave gaps (Section 3 seed-qualification check) ===")
print(as.data.frame(gaps |> select(wave_id, previous_end, onset_week, weeks_since_previous_wave)))

CEARA_2022_GAP_WEEKS <- 399 # the established long-gap standard (Ceara's own only seed-qualifying transition)
qualifying <- gaps |> dplyr::filter(weeks_since_previous_wave >= CEARA_2022_GAP_WEEKS)

message(sprintf("\nMax PE inter-wave gap = %.0f weeks (vs the %.0f-week Ceara-2022 standard). Qualifying transitions: %d.",
                 max(gaps$weeks_since_previous_wave), CEARA_2022_GAP_WEEKS, nrow(qualifying)))
if (nrow(qualifying) == 0) {
  message("CONCLUSION: No PE inter-wave gap qualifies under the established long-gap rule. NO SEED MECHANISM will be used for Pernambuco (is_seed all-zero), matching the Bahia outcome under the SAME objective rule -- not copied, independently re-derived.")
} else {
  message("CONCLUSION: The following PE transition(s) QUALIFY for a conditioned seed under the established rule -- reporting before any fitting, per instruction:")
  print(as.data.frame(qualifying))
}
write_csv(gaps, file.path(table_dir, "pernambuco_wave_gap_seed_check.csv"))
message("[saved] ", file.path(table_dir, "pernambuco_wave_gap_seed_check.csv"))
