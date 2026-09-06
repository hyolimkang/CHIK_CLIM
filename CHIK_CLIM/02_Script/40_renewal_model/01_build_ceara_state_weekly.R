# ===========================================================================
# 01_build_ceara_state_weekly.R
#
# Purpose
# -------
# Aggregate the muni-week DLNM panel up to state (UF) level for the renewal
# model track. Same aggregation logic as
# 90_exploratory/11_plot_case_climate_timeseries.R (cases summed,
# temperature/precip population-weighted mean), factored out here so both
# the exploratory plots and the renewal model read from one function.
#
# Output
# ------
# A data.frame with one row per week: week_start, week_of_year, t (time
# index), cases, population, Tmean, PRCP — not saved to disk, this is meant
# to be sourced by the fit scripts in this folder.
# ===========================================================================

for (p in c("here", "dplyr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({ library(here); library(dplyr) })

build_state_weekly <- function(panel, uf_code) {
  panel |>
    dplyr::filter(substr(muni6, 1, 2) == uf_code) |>
    dplyr::group_by(week_start, week_of_year, t) |>
    dplyr::summarise(
      cases      = sum(cases_confirmed, na.rm = TRUE),
      # weighted.mean() must run before `population` is overwritten by its
      # own sum below — dplyr::summarise() evaluates expressions in order
      # and a later expression sees the already-summarised scalar, not the
      # original per-muni vector.
      Tmean      = weighted.mean(Tmean, w = population, na.rm = TRUE),
      PRCP       = weighted.mean(PRCP, w = population, na.rm = TRUE),
      population = sum(population, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(week_start)
}

if (sys.nframe() == 0) {
  panel <- readRDS(here::here("01_Data/chik_dlnm_panel_muni_week_2015_2024.rds"))
  ceara_weekly <- build_state_weekly(panel, "23")
  cat("Ceara state-weekly series:", nrow(ceara_weekly), "weeks\n")
  print(head(ceara_weekly))
}
