# =============================================================================
# 08_plot_ceara_muni_incidence_map.R
#
# Diagnostic map to decide the v3.0 municipality subset: raw case totals
# (used for the >250-case threshold so far) conflate population size with
# real transmission signal -- a big city can clear a raw-case threshold on
# population alone while a small municipality with a real, sharp epidemic
# stays under it. This maps per-capita incidence (per 100,000) directly,
# with population shown alongside so the two can be told apart at a glance.
#
# Output
# ------
#   03_Output/figures/renewal_v3_0/ceara_muni_incidence_per100k_map.png
#   03_Output/tables/renewal_v3_0/ceara_muni_incidence_per100k.csv
# =============================================================================

required_packages <- c("here", "dplyr", "sf", "ggplot2", "patchwork", "scales", "lubridate")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(sf)
  library(ggplot2)
  library(patchwork)
})

if (sys.nframe() == 0) {
  UF_CODE <- "23"
  DATE_START <- as.Date("2015-01-04")
  DATE_END <- as.Date("2019-12-29")

  panel <- readRDS(here::here("01_Data/chik_dlnm_panel_muni_week_2015_2025.rds")) |>
    dplyr::filter(
      substr(muni6, 1, 2) == UF_CODE,
      week_start >= DATE_START, week_start <= DATE_END
    )
  total_cases <- panel |>
    dplyr::group_by(muni6) |>
    dplyr::summarise(total_cases = sum(cases_confirmed), .groups = "drop")

  # Average annual population over the same window as the case totals, not a
  # single reference year -- avoids an arbitrary single-year denominator.
  population <- readRDS(here::here("01_Data/ibge_pop_muni_year_2015_2025.rds")) |>
    dplyr::filter(
      muni6 %in% total_cases$muni6,
      year >= lubridate::year(DATE_START), year <= lubridate::year(DATE_END)
    ) |>
    dplyr::group_by(muni6) |>
    dplyr::summarise(mean_population = mean(population), .groups = "drop")

  muni_summary <- total_cases |>
    dplyr::inner_join(population, by = "muni6") |>
    dplyr::mutate(
      incidence_per_100k = total_cases / mean_population * 1e5,
      has_epidemic_signal_1000 = incidence_per_100k >= 1000,
      has_epidemic_signal_500 = incidence_per_100k >= 500,
      has_epidemic_signal_100 = incidence_per_100k >= 100
    ) |>
    dplyr::arrange(dplyr::desc(incidence_per_100k))

  if (nrow(muni_summary) != nrow(total_cases)) {
    warning(sprintf(
      "%d of %d municipalities had no matching population and are excluded",
      nrow(total_cases) - nrow(muni_summary), nrow(total_cases)
    ))
  }

  message(sprintf(
    "[signal] municipalities with incidence >=100/100k: %d | >=500/100k: %d | >=1000/100k: %d (of %d)",
    sum(muni_summary$has_epidemic_signal_100), sum(muni_summary$has_epidemic_signal_500),
    sum(muni_summary$has_epidemic_signal_1000), nrow(muni_summary)
  ))

  polygons <- readRDS(here::here("01_Data/ibge_muni_polygons.rds"))
  ceara_outline <- polygons |> dplyr::filter(substr(muni6, 1, 2) == UF_CODE)
  map_data <- ceara_outline |> dplyr::inner_join(muni_summary, by = "muni6")

  p_incidence <- ggplot() +
    geom_sf(data = ceara_outline, fill = "#F3F4F6", colour = "#D1D5DB", linewidth = 0.1) +
    geom_sf(data = map_data, aes(fill = incidence_per_100k), colour = "white", linewidth = 0.1) +
    scale_fill_viridis_c(
      option = "inferno", trans = "sqrt",
      breaks = c(0, 100, 500, 1000, 2500, 5000, 10000),
      labels = scales::label_number(big.mark = ","),
      name = "Cases per\n100,000\n(2015-2019)"
    ) +
    labs(title = "Per-capita incidence") +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5), legend.position = "right")

  p_raw_cases <- ggplot() +
    geom_sf(data = ceara_outline, fill = "#F3F4F6", colour = "#D1D5DB", linewidth = 0.1) +
    geom_sf(data = map_data, aes(fill = total_cases), colour = "white", linewidth = 0.1) +
    scale_fill_viridis_c(
      option = "viridis", trans = "sqrt",
      breaks = c(0, 100, 500, 1000, 2500, 5000, 10000, 80000),
      labels = scales::label_number(big.mark = ","),
      name = "Raw cases\n(2015-2019)"
    ) +
    labs(title = "Raw case count") +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5), legend.position = "right")

  p_population <- ggplot() +
    geom_sf(data = ceara_outline, fill = "#F3F4F6", colour = "#D1D5DB", linewidth = 0.1) +
    geom_sf(data = map_data, aes(fill = mean_population), colour = "white", linewidth = 0.1) +
    scale_fill_viridis_c(
      option = "mako", direction = -1, trans = "log10",
      labels = scales::label_number(big.mark = ",", scale_cut = scales::cut_short_scale()),
      name = "Population\n(2015-2019\nmean)"
    ) +
    labs(title = "Raw population size") +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5), legend.position = "right")

  figure <- (p_incidence | p_raw_cases | p_population) +
    patchwork::plot_annotation(
      title = "Cear<U+00E1> municipalities: case burden vs. population size, 2015-2019",
      subtitle = sprintf(
        "%d/%d municipalities reach >=100 cases per 100,000; %d reach >=500; %d reach >=1,000",
        sum(muni_summary$has_epidemic_signal_100), nrow(muni_summary),
        sum(muni_summary$has_epidemic_signal_500), sum(muni_summary$has_epidemic_signal_1000)
      ),
      theme = theme(
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9.5, colour = "#6B7280")
      )
    )

  figure_dir <- here::here("03_Output/figures/renewal_v3_0")
  table_dir <- here::here("03_Output/tables/renewal_v3_0")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  figure_path <- file.path(figure_dir, "ceara_muni_incidence_per100k_map.png")
  table_path <- file.path(table_dir, "ceara_muni_incidence_per100k.csv")

  ggsave(figure_path, figure, width = 340, height = 130, units = "mm", dpi = 300, bg = "white")
  utils::write.csv(muni_summary, table_path, row.names = FALSE)

  message("[save] ", figure_path)
  message("[save] ", table_path)
}
