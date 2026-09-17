# Pernambuco spatial-turnover diagnostic -- Sections 2 & 3.
#
# Section 2: reuse the EXISTING, objectively-defined PE major-wave census
# (brazil_chik_wave_census_v2.csv, major_epidemic_primary == TRUE) -- the
# SAME 9 waves already used for this project's wave-level PPC (PE_wave_03
# through PE_wave_11). No hand-picked windows, no threshold changes.
#
# Section 3: per-municipality, per-wave burden (reported cases, incidence
# per 100k, share of state cases) plus cumulative PRIOR burden (all weeks
# strictly before the wave's onset week).

required_packages <- c("dplyr", "readr", "tibble", "lubridate")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tibble); library(lubridate) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_turnover_diagnostic")

pe_panel <- readRDS(file.path(table_dir, "pe_municipality_week_panel.rds"))

# ---- Section 2: existing objective PE major-wave definitions --------------
wave_census <- read_csv(file.path(root, "03_Output/tables/national_wave_census_v2/brazil_chik_wave_census_v2.csv"),
                         show_col_types = FALSE)

pe_waves <- wave_census |>
  dplyr::filter(state == "PE", major_epidemic_primary == TRUE) |>
  transmute(
    wave_id, wave_order,
    start_week = as.Date(onset_week, format = "%m/%d/%Y"),
    peak_week = as.Date(peak_week, format = "%m/%d/%Y"),
    end_week = as.Date(end_week, format = "%m/%d/%Y"),
    total_cases_census = total_cases
  ) |>
  arrange(start_week)

message("=== PE major waves (existing objective census, major_epidemic_primary == TRUE) ===")
print(as.data.frame(pe_waves))

# Cross-check: does summing the muni-week panel over [start_week, end_week]
# reproduce the census's own total_cases for each wave?
wave_case_check <- pe_waves |>
  rowwise() |>
  mutate(
    total_cases_recomputed = sum(pe_panel$cases_confirmed[
      pe_panel$week_start >= start_week & pe_panel$week_start <= end_week
    ], na.rm = TRUE)
  ) |>
  ungroup() |>
  mutate(diff = total_cases_recomputed - total_cases_census)

message("\n=== Cross-check: recomputed vs census total_cases per wave ===")
print(as.data.frame(wave_case_check |> select(wave_id, total_cases_census, total_cases_recomputed, diff)))

if (any(abs(wave_case_check$diff) > 1)) {
  warning("Recomputed wave totals differ from the census by >1 case for at least one wave -- likely a window/inclusivity convention difference (e.g. end_week boundary); using RECOMPUTED totals (from the exact same panel used everywhere else in this diagnostic) as total_cases going forward.")
}
pe_waves <- wave_case_check |> mutate(total_cases = total_cases_recomputed) |> select(-total_cases_recomputed, -diff)

write_csv(pe_waves, file.path(table_dir, "pe_wave_definitions.csv"))
message("\n[saved] ", file.path(table_dir, "pe_wave_definitions.csv"))

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

muni_wave_burden <- bind_rows(lapply(seq_len(nrow(pe_waves)), function(i) {
  compute_wave_burden(pe_waves[i, ], pe_panel)
}))

# Municipalities present in the panel but with zero cases in a given wave
# still appear via group_by (present in `panel` for every week). Confirm
# full 184-muni coverage per wave (PE has 184 municipalities in the case
# panel, vs Bahia's 414).
coverage_check <- muni_wave_burden |> count(wave_id, name = "n_munis_covered")
message("\n=== Municipality coverage per wave (expect 184 for all) ===")
print(as.data.frame(coverage_check))

write_csv(muni_wave_burden, file.path(table_dir, "pe_municipality_wave_burden.csv"))
message("\n[saved] ", file.path(table_dir, "pe_municipality_wave_burden.csv"))
