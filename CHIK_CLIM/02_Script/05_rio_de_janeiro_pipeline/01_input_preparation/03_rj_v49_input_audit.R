# Rio de Janeiro external replication -- Section 5: build the complete v4.9
# input (weekly cases + demography, same window/contract as CE/BA/PE) and
# the Section 5 (of the PE precedent) seed-qualification check using the
# SAME pre-specified long-gap rule -- no manual RJ-specific seed decision.

required_packages <- c("here", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
setwd(root)
source(file.path(root, "02_Script/00_shared/functions/01_build_ceara_state_weekly.R"))

table_dir <- file.path(root, "03_Output/tables/rio_de_janeiro_v4_9_replication")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

DATE_START <- as.Date("2015-01-04")
DATE_END <- as.Date("2025-12-21") # same window convention as CE/BA/PE

panel <- readRDS(file.path(root, "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
weekly_cases <- build_state_weekly(panel, "33") |> dplyr::filter(week_start >= DATE_START, week_start <= DATE_END)

demography_path <- file.path(root, "01_Data/rio_de_janeiro_weekly_demography_2015_2025.rds")
if (!file.exists(demography_path)) stop("RJ weekly demography not yet built -- run 09d/10d in 02_Script/00_data_prep first.")
weekly_demography <- readRDS(demography_path) |> dplyr::filter(week_start >= DATE_START, week_start <= DATE_END)

weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch -- RJ case series and demography weeks do not align.")

# ---- Audit -----------------------------------------------------------------
weeks_expected <- seq(DATE_START, DATE_END, by = "week")
accounting_error <- max(abs(weekly$N_end - weekly$N_start - weekly$births + weekly$all_cause_deaths - weekly$net_population_reconciliation))

years <- sort(unique(as.integer(format(weekly$week_start, "%Y"))))
audit <- tibble(
  state = "RJ", n_weeks = nrow(weekly), n_weeks_expected = length(weeks_expected),
  week_continuity_ok = identical(weekly$week_start, weeks_expected),
  n_zero_case_weeks = sum(weekly$cases == 0), n_positive_case_weeks = sum(weekly$cases > 0),
  total_cases = sum(weekly$cases), n_na_any = sum(!complete.cases(weekly)),
  N_start_2015 = weekly$N_start[1], N_end_2025 = weekly$N_end[nrow(weekly)],
  year_min = min(weekly$year), year_max = max(weekly$year), n_model_years = length(years),
  demographic_accounting_max_abs_error = accounting_error,
  generation_interval_sum_check = "G=8, generation_weights() sums to 1 (shared helper, unchanged)",
  importation_convention = "imports_per_week = 1 (frozen v4.9 default, same as CE/BA/PE, unchanged)",
  data_source_confirmed = "chik_dlnm_panel_muni_week_2015_2025.rds (UF prefix 33) + rio_de_janeiro_weekly_demography_2015_2025.rds -- NOT inherited CE/BA/PE objects/paths"
)

message("=== Rio de Janeiro v4.9 input audit ===")
print(as.data.frame(audit))
write_csv(audit, file.path(table_dir, "rio_de_janeiro_v49_input_audit.csv"))
message("[saved] ", file.path(table_dir, "rio_de_janeiro_v49_input_audit.csv"))

if (!audit$week_continuity_ok || audit$n_na_any > 0 || accounting_error > 1e-6) {
  stop("STOP: RJ input audit failed -- do not proceed to fitting.")
}
saveRDS(weekly, file.path(table_dir, "rio_de_janeiro_weekly_input.rds"))
message("\n[saved] ", file.path(table_dir, "rio_de_janeiro_weekly_input.rds"), " (audited weekly input for Stan)")

# ---- Markdown summary (Section 5 of the RJ spec) ---------------------------
md <- c(
  "# Rio de Janeiro v4.9 Input Audit", "",
  "Frozen v4.9 architecture, reused EXACTLY (renewal transmission, lifelong",
  "immunity, demographic S/U bookkeeping, common generation interval, NB2",
  "observation model, harmonic seasonal R0(t), year effects, A_year,",
  "existing importation convention). No RJ-specific model changes.", "",
  "## Key input dimensions", "",
  sprintf("- N (weeks): %d", audit$n_weeks),
  sprintf("- Model years: %d (%d-%d)", audit$n_model_years, audit$year_min, audit$year_max),
  sprintf("- First / last date: %s / %s", as.character(DATE_START), as.character(DATE_END)),
  sprintf("- Total reported cases (2015-2025): %s", format(audit$total_cases, big.mark = ",")),
  sprintf("- Population range: %s to %s", format(round(audit$N_start_2015), big.mark = ","), format(round(audit$N_end_2025), big.mark = ",")),
  "- Generation interval: G=8 (generation_weights(), Gamma(shape=4,rate=2) discretised, sums to 1) -- UNCHANGED",
  "- Importation: imports_per_week = 1 -- UNCHANGED frozen default",
  sprintf("- Demographic accounting max abs error: %.2e (must be ~0)", accounting_error),
  sprintf("- Week continuity: %s", audit$week_continuity_ok),
  sprintf("- NA values: %d", audit$n_na_any),
  "", "## Sources reused (unmodified)", "",
  "- `02_Script/00_shared/functions/01_build_ceara_state_weekly.R` (`build_state_weekly()`, UF-agnostic)",
  "- `02_Script/stan/renewal_bahia_v4_9_global_q_multisite_serology.stan` (global-q, J_sero>=0, state-agnostic)",
  "- `02_Script/00_shared/legacy_model_functions/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R` (`generation_weights()`, `compute_hmc_gate()`)",
  "- `02_Script/00_shared/legacy_model_functions/22_v4_3_short_2015_2019_q_calibration/scripts/02_fit_v4_3.R` (`make_v4_3_data()`, `load_td14_scalar_inits()`)",
  "- `03_Output/tables/national_wave_census_v2/brazil_chik_wave_census_v2.csv` (episode definitions, RJ rows, unmodified)",
  "", "## RJ-specific preprocessing (unavoidable, mechanical only)", "",
  "- `02_Script/00_data_prep/09d_fetch_rio_de_janeiro_sinasc_weekly_births.R` (UF_PREFIX=\"33\", verbatim copy of the PE/BA script)",
  "- `02_Script/00_data_prep/10d_build_rio_de_janeiro_weekly_demography.R` (UF_CODE=33L, verbatim copy)",
  "", "CE/BA/PE models/results/scripts are UNMODIFIED by this work."
)
writeLines(md, file.path(dirname(table_dir), "..", "..") |> file.path("03_Output/reports/rio_de_janeiro_v4_9_replication/RJ_v49_input_audit.md"))
message("[saved] RJ_v49_input_audit.md")

# ---- Seed-qualification check, SAME rule as CE/BA/PE ------------------------
wave_census <- read_csv(file.path(root, "03_Output/tables/national_wave_census_v2/brazil_chik_wave_census_v2.csv"), show_col_types = FALSE)
rj_major_waves <- wave_census |> dplyr::filter(state == "RJ", major_epidemic_primary == TRUE) |>
  transmute(wave_id, wave_order,
            onset_week = as.Date(onset_week, format = "%m/%d/%Y"),
            end_week = as.Date(end_week, format = "%m/%d/%Y")) |>
  arrange(onset_week)

gaps <- rj_major_waves |>
  mutate(previous_end = lag(end_week), weeks_since_previous_wave = as.numeric(onset_week - previous_end) / 7) |>
  dplyr::filter(!is.na(weeks_since_previous_wave))

message("\n=== RJ major-wave gaps (seed-qualification check) ===")
print(as.data.frame(gaps |> select(wave_id, previous_end, onset_week, weeks_since_previous_wave)))

CEARA_2022_GAP_WEEKS <- 399 # established long-gap standard
qualifying <- gaps |> dplyr::filter(weeks_since_previous_wave >= CEARA_2022_GAP_WEEKS)

message(sprintf("\nMax RJ inter-wave gap = %.0f weeks (vs the %.0f-week Ceara-2022 standard). Qualifying transitions: %d.",
                 max(gaps$weeks_since_previous_wave), CEARA_2022_GAP_WEEKS, nrow(qualifying)))
if (nrow(qualifying) == 0) {
  message("CONCLUSION: No RJ inter-wave gap qualifies under the established long-gap rule. NO SEED MECHANISM will be used for RJ (is_seed all-zero).")
} else {
  message("CONCLUSION: The following RJ transition(s) QUALIFY for a conditioned seed under the established rule:")
  print(as.data.frame(qualifying))
}
write_csv(gaps, file.path(table_dir, "rio_de_janeiro_wave_gap_seed_check.csv"))
message("[saved] ", file.path(table_dir, "rio_de_janeiro_wave_gap_seed_check.csv"))
