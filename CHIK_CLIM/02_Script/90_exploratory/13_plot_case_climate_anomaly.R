# ===========================================================================
# 13_plot_case_climate_anomaly.R
#
# Purpose
# -------
# Builds a week-of-year climatology (the 2015-2024 average Tmean/PRCP for
# each of week-of-year 1-52/53) and expresses each week's observed value as
# an anomaly relative to that typical value:
#   anomaly = observed - climatology[week_of_year]
# This isolates "was this particular year's climate unusual for this time
# of year" from the deterministic seasonal cycle that repeats every year
# regardless of whether an outbreak occurs — the confound flagged when
# looking at 12_plot_case_climate_ccf.R's raw-value CCF (temperature's
# lag-0 correlation there is a plausible seasonality artifact, since both
# series have their own annual cycle independent of any causal link).
#
# Produces two figures:
#   1. Stacked panel: cases + Tmean anomaly + PRCP anomaly over time
#      (diverging fill: warmer/wetter-than-normal vs cooler/drier-than-normal)
#   2. The same CCF diagnostic as script 12, but on anomalies instead of
#      raw values — if a lag effect survives after removing the seasonal
#      cycle, that's a much stronger signal than the raw-value version.
# ===========================================================================

for (p in c("here", "dplyr", "ggplot2", "patchwork", "scales")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(dplyr); library(ggplot2); library(patchwork)
})

COL_CASES <- "#2a78d6"
COL_WARM  <- "#e34948"   # diverging pair: warmer/wetter than normal
COL_COOL  <- "#2a78d6"   # diverging pair: cooler/drier than normal
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

# Aggregates the muni-week DLNM panel up to one state, then attaches a
# week-of-year climatology and each week's anomaly relative to it.
build_state_anomaly_weekly <- function(panel, uf_code) {
  state_weekly <- panel |>
    dplyr::filter(substr(muni6, 1, 2) == uf_code) |>
    dplyr::group_by(week_start, week_of_year) |>
    dplyr::summarise(
      cases = sum(cases_confirmed, na.rm = TRUE),
      Tmean = weighted.mean(Tmean, w = population, na.rm = TRUE),
      PRCP  = weighted.mean(PRCP, w = population, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(week_start)

  # Week-of-year climatology: the multi-year average for each ISO
  # week-of-year (1-52/53), pooling all years in the panel.
  climatology <- state_weekly |>
    dplyr::group_by(week_of_year) |>
    dplyr::summarise(
      clim_Tmean = mean(Tmean, na.rm = TRUE),
      clim_PRCP  = mean(PRCP, na.rm = TRUE),
      .groups = "drop"
    )

  state_weekly |>
    dplyr::left_join(climatology, by = "week_of_year") |>
    dplyr::mutate(
      Tmean_anomaly = Tmean - clim_Tmean,
      PRCP_anomaly  = PRCP - clim_PRCP,
      log_cases     = log1p(cases)
    )
}

plot_anomaly_timeseries <- function(state_weekly, uf_name) {
  p_cases <- ggplot(state_weekly, aes(week_start, cases)) +
    geom_col(fill = COL_CASES, width = 6) +
    labs(title = sprintf("Chikungunya — %s (climate anomaly vs week-of-year norm)", uf_name),
         x = NULL, y = "Confirmed cases / week") +
    theme_stack()

  p_temp <- ggplot(state_weekly, aes(week_start, Tmean_anomaly, fill = Tmean_anomaly > 0)) +
    geom_col(width = 6, show.legend = FALSE) +
    geom_hline(yintercept = 0, color = COL_AXIS, linewidth = 0.4) +
    scale_fill_manual(values = c(`TRUE` = COL_WARM, `FALSE` = COL_COOL)) +
    labs(x = NULL, y = "Temp anomaly (°C)") +
    theme_stack()

  p_prcp <- ggplot(state_weekly, aes(week_start, PRCP_anomaly, fill = PRCP_anomaly > 0)) +
    geom_col(width = 6, show.legend = FALSE) +
    geom_hline(yintercept = 0, color = COL_AXIS, linewidth = 0.4) +
    scale_fill_manual(values = c(`TRUE` = COL_WARM, `FALSE` = COL_COOL)) +
    labs(x = "Week", y = "Precip anomaly (mm)") +
    theme_stack()

  (p_cases / p_temp / p_prcp) + patchwork::plot_layout(heights = c(1.3, 1, 1))
}

compute_ccf_df <- function(x, y, max_lag) {
  cc <- ccf(x, y, lag.max = max_lag, plot = FALSE, na.action = na.pass)
  n <- length(x)
  ci <- qnorm(0.975) / sqrt(n)
  tibble::tibble(
    lag = as.integer(cc$lag),
    corr = as.numeric(cc$acf),
    significant = abs(corr) > ci
  ) |>
    dplyr::mutate(ci = ci)
}

plot_anomaly_ccf <- function(state_weekly, uf_name, max_lag = 16) {
  ccf_temp <- compute_ccf_df(state_weekly$Tmean_anomaly, state_weekly$log_cases, max_lag) |>
    dplyr::filter(lag >= 0)
  ccf_prcp <- compute_ccf_df(state_weekly$PRCP_anomaly, state_weekly$log_cases, max_lag) |>
    dplyr::filter(lag >= 0)

  make_ccf_panel <- function(df, col, title) {
    ggplot(df, aes(lag, corr)) +
      geom_hline(yintercept = 0, color = COL_AXIS, linewidth = 0.4) +
      geom_hline(aes(yintercept = ci), linetype = "dashed", color = COL_AXIS, linewidth = 0.4) +
      geom_hline(aes(yintercept = -ci), linetype = "dashed", color = COL_AXIS, linewidth = 0.4) +
      geom_col(aes(fill = significant), width = 0.6, show.legend = FALSE) +
      scale_fill_manual(values = c(`TRUE` = col, `FALSE` = COL_GRID)) +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      labs(title = title, x = "Lag (weeks, climate anomaly leads cases)", y = "Cross-correlation") +
      theme_stack()
  }

  p_temp <- make_ccf_panel(ccf_temp, "#eb6834", sprintf("%s — Temperature ANOMALY vs log(cases)", uf_name))
  p_prcp <- make_ccf_panel(ccf_prcp, "#1baf7a", sprintf("%s — Precipitation ANOMALY vs log(cases)", uf_name))

  p_temp / p_prcp
}

if (sys.nframe() == 0) {
  UF_CODE <- "23"
  UF_NAME <- "Ceara (CE)"

  panel <- readRDS(here::here("01_Data/chik_dlnm_panel_muni_week_2015_2024.rds"))
  state_weekly <- build_state_anomaly_weekly(panel, UF_CODE)

  fig_ts  <- plot_anomaly_timeseries(state_weekly, UF_NAME)
  fig_ccf <- plot_anomaly_ccf(state_weekly, UF_NAME)

  out_ts  <- here::here(sprintf("03_Output/figures/case_climate_anomaly_timeseries_%s.png", tolower(UF_CODE)))
  out_ccf <- here::here(sprintf("03_Output/figures/case_climate_anomaly_ccf_%s.png", tolower(UF_CODE)))
  dir.create(dirname(out_ts), recursive = TRUE, showWarnings = FALSE)

  ggsave(out_ts, fig_ts, width = 10, height = 7, dpi = 150)
  ggsave(out_ccf, fig_ccf, width = 9, height = 7, dpi = 150)

  message("[save] ", out_ts)
  message("[save] ", out_ccf)
}
