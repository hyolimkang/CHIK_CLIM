# ===========================================================================
# 04_build_dlnm_panel.R
#
# Purpose
# -------
# Join the case panel (muni x week) + population (muni x year) + national
# climate (muni x week, from 01/02_fetch_era5land_*.R +
# 03_build_national_climate_weekly.R) into a single DLNM-ready analysis panel.
#
# Handling
# --------
# - muni6 with no climate at all (~30 SINAN "unknown residence" placeholder
#   codes not present in the IBGE centroid/polygon data) are excluded.
# - muni6 with NA climate values (e.g. Fernando de Noronha — tiny island
#   municipalities that barely overlap the ERA5-Land grid since it only
#   covers land — see script 03's comments) are excluded the same way.
# - Population only exists for 2015-2024 (IBGE hasn't published a 2014
#   estimate). The panel's first week (2014-12-28, epi_year=2014) has its
#   population-join year substituted with 2015 (carry-back) — the same
#   logic build_uf_weekly_panel.R uses to carry the last known year
#   forward, just applied in the opposite direction.
# - t: a shared national calendar-time index (the same week_start gets the
#   same t across every municipality) — used for the DLNM model's s(t)
#   long-term trend term.
#
# Output
# ------
#   01_Data/chik_dlnm_panel_muni_week_2015_2024.rds
#     columns: muni6, week_start, year, week_of_year, t,
#              cases_confirmed, population, PRCP, Tmin, Tmax, Tmean
# ===========================================================================

for (p in c("here", "dplyr", "lubridate", "readr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(dplyr); library(lubridate); library(readr)
})

build_dlnm_panel <- function(case_panel, climate, population) {
  n_case_muni <- dplyr::n_distinct(case_panel$muni6)

  no_climate_at_all <- setdiff(unique(case_panel$muni6), unique(climate$muni6))

  na_muni <- climate |>
    dplyr::filter(if_any(c(PRCP, Tmin, Tmax, Tmean), is.na)) |>
    dplyr::distinct(muni6) |>
    dplyr::pull(muni6)

  dropped <- union(no_climate_at_all, na_muni)
  if (length(dropped) > 0) {
    message(sprintf(
      "[drop] %d/%d case-panel muni6 excluded (no climate data or NA climate values): %s%s",
      length(dropped), n_case_muni,
      paste(head(dropped, 10), collapse = ", "),
      if (length(dropped) > 10) ", ..." else ""
    ))
  }

  climate <- climate |> dplyr::filter(!muni6 %in% dropped)

  panel <- case_panel |>
    dplyr::filter(muni6 %in% unique(climate$muni6)) |>
    dplyr::mutate(
      week_start = as.Date(week_start),
      year         = lubridate::year(week_start),
      week_of_year = lubridate::isoweek(week_start),
      year_for_pop = pmin(pmax(year, min(population$year)), max(population$year))
    )

  panel <- panel |>
    dplyr::inner_join(
      climate |> dplyr::mutate(week_start = as.Date(week_start)),
      by = c("muni6", "week_start")
    )

  panel <- panel |>
    dplyr::left_join(
      population |> dplyr::select(muni6, year, population),
      by = c("muni6", "year_for_pop" = "year")
    )

  week_lookup <- tibble::tibble(week_start = sort(unique(panel$week_start))) |>
    dplyr::mutate(t = dplyr::row_number())
  panel <- panel |> dplyr::left_join(week_lookup, by = "week_start")

  panel <- panel |>
    dplyr::select(
      muni6, week_start, year, week_of_year, t,
      cases_confirmed, population, PRCP, Tmin, Tmax, Tmean
    ) |>
    dplyr::arrange(muni6, week_start)

  # ---- fail-fast checks -----------------------------------------------
  expected_rows <- case_panel |>
    dplyr::filter(muni6 %in% unique(climate$muni6)) |>
    nrow()
  if (nrow(panel) != expected_rows) {
    stop(sprintf(
      "Row count mismatch after join: got %d, expected %d (case panel restricted to munis with climate).",
      nrow(panel), expected_rows
    ))
  }

  na_counts <- colSums(is.na(panel))
  if (any(na_counts[c("cases_confirmed", "population", "PRCP", "Tmin", "Tmax", "Tmean")] > 0)) {
    stop("Missing values found after join:\n", paste(capture.output(print(na_counts)), collapse = "\n"))
  }

  panel
}

if (sys.nframe() == 0) {
  case_panel <- readRDS(here::here("01_Data/chik_brazil_muni_week_2015_2024.rds"))
  climate    <- readRDS(here::here("01_Data/chik_muni_week_climate_2015_2024.rds"))
  population <- readRDS(here::here("01_Data/ibge_pop_muni_year_2015_2024.rds"))

  dlnm_panel <- build_dlnm_panel(case_panel, climate, population)

  out_path <- here::here("01_Data/chik_dlnm_panel_muni_week_2015_2024.rds")
  saveRDS(dlnm_panel, out_path)
  readr::write_csv(dlnm_panel, sub("\\.rds$", ".csv", out_path))

  message(sprintf(
    "[save] %s\n  %s rows | %d municipalities | %d weeks",
    out_path,
    format(nrow(dlnm_panel), big.mark = ","),
    dplyr::n_distinct(dlnm_panel$muni6),
    dplyr::n_distinct(dlnm_panel$week_start)
  ))
}
