# ===========================================================================
# 11_plot_case_climate_timeseries.R
#
# Purpose
# -------
# Descriptive time-series plot: weekly confirmed cases, temperature, and
# precipitation for one state (UF), aggregated up from the muni-week DLNM
# panel. Meant as an eyeball check for synchrony/lag between climate and
# cases before formal DLNM modelling — not a substitute for it.
#
# Why 3 stacked panels instead of one dual/triple-axis plot
# -----------------------------------------------------------
# Cases, temperature (°C), and precipitation (mm) have unrelated scales.
# Overlaying them on a shared or dual y-axis lets an arbitrary axis-scaling
# choice make the series look more or less correlated than they really are
# (the classic "dual-axis chart" pitfall). Instead each variable gets its
# own panel on its own y-axis, sharing only the x-axis (week_start) — you
# can still compare timing/lag visually without the scale distortion.
#
# Aggregation from muni to state (UF)
# ------------------------------------
# - cases_confirmed: summed across municipalities in the UF (a real count).
# - Tmean / PRCP: population-weighted mean across municipalities (mm of
#   rain or degrees C don't sum meaningfully across places the way case
#   counts do; weighting by population approximates "what the average
#   resident of this state experienced that week").
# ===========================================================================

for (p in c("here", "dplyr", "ggplot2", "patchwork", "scales")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(dplyr); library(ggplot2); library(patchwork); library(scales)
})

# Palette slots 1-3 (blue/orange/aqua) — validated colorblind-safe, all-pairs.
COL_CASES <- "#2a78d6"
COL_TEMP  <- "#eb6834"
COL_PRCP  <- "#1baf7a"
COL_GRID  <- "#e1e0d9"
COL_AXIS  <- "#898781"
COL_INK   <- "#0b0b0b"

theme_stack <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.4),
      axis.line.x = element_line(color = COL_AXIS, linewidth = 0.4),
      axis.ticks = element_blank(),
      axis.text = element_text(color = COL_AXIS),
      axis.title.y = element_text(color = COL_INK, size = 10),
      plot.title = element_text(color = COL_INK, size = 11, face = "bold"),
      plot.margin = margin(4, 8, 2, 8)
    )
}

plot_case_climate_timeseries <- function(panel, uf_code, uf_name = uf_code) {
  state_weekly <- panel |>
    dplyr::filter(substr(muni6, 1, 2) == uf_code) |>
    dplyr::group_by(week_start) |>
    dplyr::summarise(
      cases = sum(cases_confirmed, na.rm = TRUE),
      Tmean = weighted.mean(Tmean, w = population, na.rm = TRUE),
      PRCP  = weighted.mean(PRCP, w = population, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(week_start)

  p_cases <- ggplot(state_weekly, aes(week_start, cases)) +
    geom_col(fill = COL_CASES, width = 6) +
    labs(title = sprintf("Chikungunya — %s", uf_name), x = NULL, y = "Confirmed cases / week") +
    theme_stack()

  p_temp <- ggplot(state_weekly, aes(week_start, Tmean)) +
    geom_line(color = COL_TEMP, linewidth = 0.7) +
    labs(x = NULL, y = "Mean temperature (°C)") +
    theme_stack()

  p_prcp <- ggplot(state_weekly, aes(week_start, PRCP)) +
    geom_line(color = COL_PRCP, linewidth = 0.7) +
    labs(x = "Week", y = "Precipitation (mm/week)") +
    theme_stack()

  (p_cases / p_temp / p_prcp) +
    patchwork::plot_layout(heights = c(1.3, 1, 1))
}

if (sys.nframe() == 0) {
  # ============================================================
  UF_CODE <- "23"    # Ceara
  UF_NAME <- "Ceara (CE)"
  # ============================================================

  panel <- readRDS(here::here("01_Data/chik_dlnm_panel_muni_week_2015_2024.rds"))
  fig <- plot_case_climate_timeseries(panel, UF_CODE, UF_NAME)

  out_path <- here::here(sprintf("03_Output/figures/case_climate_timeseries_%s.png", tolower(UF_CODE)))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_path, fig, width = 10, height = 7, dpi = 150)
  message("[save] ", out_path)
}
