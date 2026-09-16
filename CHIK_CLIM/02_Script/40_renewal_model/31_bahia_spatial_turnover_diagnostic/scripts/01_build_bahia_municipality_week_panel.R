# Bahia spatial-turnover diagnostic -- Section 1.
#
# Descriptive, model-assumption diagnostic. Does NOT fit or modify any Stan
# model, does NOT estimate q, does NOT use serology. Purpose: locate and
# audit the existing, already-validated municipality-by-week chikungunya
# panel for Bahia, and confirm that summing it up reproduces the exact
# Bahia state-level weekly case series already used by the frozen renewal
# model (build_state_weekly(panel, "29")).

required_packages <- c("dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
setwd(root)
source(file.path(root, "02_Script/40_renewal_model/00_shared/01_build_ceara_state_weekly.R"))

table_dir <- file.path(root, "03_Output/tables/bahia_spatial_turnover_diagnostic")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

DATE_START <- as.Date("2015-01-04")
DATE_END <- as.Date("2025-12-21") # matches the frozen Bahia renewal model's exact window (573 weeks)

# ---- 1. Load the EXISTING validated muni-week panel (no SINAN rebuild) ----
panel <- readRDS(file.path(root, "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
name_lookup <- readRDS(file.path(root, "01_Data/ibge_muni_name_lookup.rds")) |>
  filter(uf == "BA") |>
  select(muni6, name_muni, uf)

bahia_panel <- panel |>
  filter(substr(muni6, 1, 2) == "29", week_start >= DATE_START, week_start <= DATE_END) |>
  left_join(name_lookup, by = "muni6") |>
  select(muni6, name_muni, week_start, year, week_of_year, t, cases_confirmed, population) |>
  arrange(muni6, week_start)

if (any(is.na(bahia_panel$name_muni))) {
  stop("Some Bahia muni6 codes in the panel have no name-lookup match -- investigate before proceeding.")
}

saveRDS(bahia_panel, file.path(table_dir, "bahia_municipality_week_panel.rds"))

# ---- 2. Audit ----
n_munis <- length(unique(bahia_panel$muni6))
weeks_expected <- seq(DATE_START, DATE_END, by = "week")
n_weeks_expected <- length(weeks_expected)

per_muni_week_count <- bahia_panel |> count(muni6, name = "n_weeks")
incomplete_munis <- per_muni_week_count |> filter(n_weeks != n_weeks_expected)

zero_vs_missing <- bahia_panel |>
  summarise(
    n_muni_week_rows = n(),
    n_expected_rows = n_munis * n_weeks_expected,
    n_zero_case_weeks = sum(cases_confirmed == 0L, na.rm = TRUE),
    n_positive_case_weeks = sum(cases_confirmed > 0L, na.rm = TRUE),
    n_na_cases = sum(is.na(cases_confirmed)),
    n_na_population = sum(is.na(population))
  )

total_weekly_cases <- bahia_panel |>
  group_by(week_start) |>
  summarise(cases = sum(cases_confirmed, na.rm = TRUE), .groups = "drop")

# ---- 3. Confirm aggregation reproduces the EXACT state series used by the ----
#         frozen Bahia renewal model (build_state_weekly, same panel/window).
official_state_weekly <- build_state_weekly(panel, "29") |>
  filter(week_start >= DATE_START, week_start <= DATE_END)

reconciled <- total_weekly_cases |>
  inner_join(official_state_weekly |> select(week_start, cases_official = cases), by = "week_start")

max_abs_diff <- max(abs(reconciled$cases - reconciled$cases_official))
n_week_mismatch <- sum(reconciled$cases != reconciled$cases_official)
weeks_covered_match <- nrow(reconciled) == n_weeks_expected && nrow(reconciled) == nrow(official_state_weekly)

audit <- tibble(
  n_municipalities = n_munis,
  n_weeks_expected = n_weeks_expected,
  n_municipalities_with_incomplete_weeks = nrow(incomplete_munis),
  n_muni_week_rows = zero_vs_missing$n_muni_week_rows,
  n_expected_rows = zero_vs_missing$n_expected_rows,
  n_zero_case_weeks = zero_vs_missing$n_zero_case_weeks,
  n_positive_case_weeks = zero_vs_missing$n_positive_case_weeks,
  n_na_cases = zero_vs_missing$n_na_cases,
  n_na_population = zero_vs_missing$n_na_population,
  total_bahia_cases_2015_2025 = sum(bahia_panel$cases_confirmed, na.rm = TRUE),
  population_coverage_min = min(bahia_panel$population, na.rm = TRUE),
  population_coverage_max = max(bahia_panel$population, na.rm = TRUE),
  weeks_covered_match_expected_and_official = weeks_covered_match,
  max_abs_diff_vs_official_state_series = max_abs_diff,
  n_weeks_mismatched_vs_official = n_week_mismatch
)

message("=== Bahia municipality-week panel audit ===")
print(as.data.frame(audit))
if (nrow(incomplete_munis) > 0L) {
  message("\nMunicipalities with incomplete week coverage:")
  print(as.data.frame(incomplete_munis))
}

write_csv(audit, file.path(table_dir, "bahia_municipality_week_panel_audit.csv"))
write_csv(incomplete_munis, file.path(table_dir, "bahia_municipality_incomplete_weeks.csv"))
message("\n[saved] ", file.path(table_dir, "bahia_municipality_week_panel.rds"))
message("[saved] ", file.path(table_dir, "bahia_municipality_week_panel_audit.csv"))

if (max_abs_diff > 0) {
  stop("STOP: municipality-level aggregation does NOT exactly reproduce the state-level series used by the frozen renewal model. Investigate before proceeding to wave burden analysis.")
} else {
  message("\nCONFIRMED: summing the municipality-week panel reproduces the exact state-level case series used by the frozen Bahia renewal model (max abs diff = 0 across ", nrow(reconciled), " weeks).")
}
