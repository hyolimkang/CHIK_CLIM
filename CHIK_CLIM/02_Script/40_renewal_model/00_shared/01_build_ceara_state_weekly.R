# ===========================================================================
# 01_build_ceara_state_weekly.R
#
# Purpose
# -------
# Aggregate the muni-week DLNM panel up to state (UF) level for the renewal
# model track. Cases are summed and temperature/precip are population-
# weighted means, same logic as 90_exploratory/11_plot_case_climate_timeseries.R.
# The STATE TOTAL population, however, is not the sum of the muni-level
# annual estimates used elsewhere in this repo <U+2014> those mix vintages (SIDRA
# annual estimates, the 2022 Census, and 2023 interpolation) and produce a
# spurious ~450k (~5%) dip-and-recover in Ceara's population across 2021-23
# that has nothing to do with real demographic change (see
# 00_data_prep/07_fetch_ibge_population_projection_uf.R for the comparison).
# Since the renewal model drives weekly susceptible-pool changes directly
# off population deltas, that artifact would read as a mass migration
# event. Instead, state population is joined in from IBGE's population
# projection revision 2024 <U+2014> a single internally consistent UF-level series
# for the whole 2015-2025 span.
#
# Output
# ------
# A data.frame with one row per week: week_start, week_of_year, t (time
# index), year, cases, population, Tmean, PRCP <U+2014> not saved to disk, this is
# meant to be sourced by the fit scripts in this folder.
# ===========================================================================

for (p in c("here", "dplyr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
})

build_state_weekly <- function(panel, uf_code) {
  target_uf_code <- as.integer(uf_code)
  uf_population <- readRDS(here::here("01_Data/ibge_population_projection_uf_2024revision.rds")) |>
    dplyr::filter(.data$uf_code == target_uf_code) |>
    dplyr::transmute(year = .data$year, population = .data$pop_total)

  weekly <- panel |>
    dplyr::filter(substr(muni6, 1, 2) == uf_code) |>
    dplyr::group_by(week_start, week_of_year, t, year) |>
    dplyr::summarise(
      cases = sum(cases_confirmed, na.rm = TRUE),
      # Climate is still weighted by each municipality's own annual
      # population share within the state -- that relative weighting is
      # far less sensitive to which population vintage is used than the
      # state TOTAL below, which is what the renewal model actually reads.
      Tmean = weighted.mean(Tmean, w = population, na.rm = TRUE),
      PRCP = weighted.mean(PRCP, w = population, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(uf_population, by = "year") |>
    dplyr::arrange(week_start)

  if (anyNA(weekly$population)) {
    stop(
      "State population projection (revision 2024) has no coverage for some ",
      "years present in the panel -- check 01_Data/ibge_population_projection_uf_2024revision.rds"
    )
  }
  weekly
}

if (sys.nframe() == 0) {
  panel <- readRDS(here::here("01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
  ceara_weekly <- build_state_weekly(panel, "23")
  cat("Ceara state-weekly series:", nrow(ceara_weekly), "weeks\n")
  print(head(ceara_weekly))
}
