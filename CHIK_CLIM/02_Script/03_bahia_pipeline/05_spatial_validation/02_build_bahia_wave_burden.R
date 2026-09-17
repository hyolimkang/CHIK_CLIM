# Bahia spatial-turnover diagnostic -- Sections 2 & 3.
#
# Section 2: reuse the EXISTING, objectively-defined Bahia major-wave
# census (brazil_chik_wave_census_v2.csv, major_epidemic_primary == TRUE)
# -- the SAME 9 waves already used for this project's wave-level PPC. No
# hand-picked windows, no threshold changes.
#
# Section 3: per-municipality, per-wave burden (reported cases, incidence
# per 100k, share of state cases) plus cumulative PRIOR burden (all weeks
# strictly before the wave's onset week).

required_packages <- c("dplyr", "readr", "tibble", "lubridate")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tibble); library(lubridate) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
table_dir <- file.path(root, "03_Output/03_bahia_pipeline/tables/bahia_spatial_turnover_diagnostic")

bahia_panel <- readRDS(file.path(table_dir, "bahia_municipality_week_panel.rds"))

# ---- Section 2: existing objective Bahia major-wave definitions ----------
wave_census <- read_csv(file.path(root, "03_Output/07_national_pipeline/tables/national_wave_census_v2/brazil_chik_wave_census_v2.csv"),
                         show_col_types = FALSE)

bahia_waves <- wave_census |>
  dplyr::filter(state == "BA", major_epidemic_primary == TRUE) |>
  transmute(
    wave_id, wave_order,
    start_week = as.Date(onset_week, format = "%m/%d/%Y"),
    peak_week = as.Date(peak_week, format = "%m/%d/%Y"),
    end_week = as.Date(end_week, format = "%m/%d/%Y"),
    total_cases_census = total_cases
  ) |>
  arrange(start_week)

message("=== Bahia major waves (existing objective census, major_epidemic_primary == TRUE) ===")
print(as.data.frame(bahia_waves))

# Cross-check: does summing the muni-week panel over [start_week, end_week]
# reproduce the census's own total_cases for each wave?
wave_case_check <- bahia_waves |>
  rowwise() |>
  mutate(
    total_cases_recomputed = sum(bahia_panel$cases_confirmed[
      bahia_panel$week_start >= start_week & bahia_panel$week_start <= end_week
    ], na.rm = TRUE)
  ) |>
  ungroup() |>
  mutate(diff = total_cases_recomputed - total_cases_census)

message("\n=== Cross-check: recomputed vs census total_cases per wave ===")
print(as.data.frame(wave_case_check |> select(wave_id, total_cases_census, total_cases_recomputed, diff)))

if (any(abs(wave_case_check$diff) > 1)) {
  warning("Recomputed wave totals differ from the census by >1 case for at least one wave -- likely a window/inclusivity convention difference (e.g. end_week boundary); using RECOMPUTED totals (from the exact same panel used everywhere else in this diagnostic) as total_cases going forward.")
}
bahia_waves <- wave_case_check |> mutate(total_cases = total_cases_recomputed) |> select(-total_cases_recomputed, -diff)

write_csv(bahia_waves, file.path(table_dir, "bahia_wave_definitions.csv"))
message("\n[saved] ", file.path(table_dir, "bahia_wave_definitions.csv"))

# ---- Section 3: municipality-level wave burden ----------------------------
compute_wave_burden <- function(wave_row, panel) {
  in_wave <- panel |> dplyr::filter(week_start >= wave_row$start_week, week_start <= wave_row$end_week)
  before_wave <- panel |> dplyr::filter(week_start < wave_row$start_week)

  burden <- in_wave |>
    group_by(muni6, name_muni) |>
    summarise(
      reported_cases_iw = sum(cases_confirmed, na.rm = TRUE),
      population_mean_iw = mean(population, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      wave_id = wave_row$wave_id,
      incidence_iw = reported_cases_iw / population_mean_iw * 100000,
      total_state_cases_w = sum(reported_cases_iw)
    ) |>
    mutate(share_of_state_cases_iw = if_else(total_state_cases_w > 0, reported_cases_iw / total_state_cases_w, NA_real_))

  prior <- before_wave |>
    group_by(muni6) |>
    summarise(
      cumulative_previous_cases_iw = sum(cases_confirmed, na.rm = TRUE),
      population_mean_prior = mean(population, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(cumulative_previous_incidence_iw = cumulative_previous_cases_iw / population_mean_prior * 100000)

  burden |>
    left_join(prior |> select(muni6, cumulative_previous_cases_iw, cumulative_previous_incidence_iw), by = "muni6") |>
    select(wave_id, muni6, name_muni, reported_cases_iw, population_mean_iw, incidence_iw,
           share_of_state_cases_iw, cumulative_previous_cases_iw, cumulative_previous_incidence_iw)
}

muni_wave_burden <- bind_rows(lapply(seq_len(nrow(bahia_waves)), function(i) {
  compute_wave_burden(bahia_waves[i, ], bahia_panel)
}))

# Municipalities present in the panel but with zero cases in a given wave
# will already appear via the group_by (only if reported_cases_iw==0 rows
# survive summarise, which they do since group_by uses all muni6 present in
# in_wave -- but muni6 with ALL-zero weeks in a wave still appear because
# they are present in `panel` for every week). Confirm full 414-muni coverage
# per wave.
coverage_check <- muni_wave_burden |> count(wave_id, name = "n_munis_covered")
message("\n=== Municipality coverage per wave (expect 414 for all) ===")
print(as.data.frame(coverage_check))

write_csv(muni_wave_burden, file.path(table_dir, "bahia_municipality_wave_burden.csv"))
message("\n[saved] ", file.path(table_dir, "bahia_municipality_wave_burden.csv"))
