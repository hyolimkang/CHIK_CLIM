# ===========================================================================
# 12_plot_case_climate_ccf.R
#
# Purpose
# -------
# Cross-correlation function (CCF) between weekly climate (Tmean, PRCP) and
# chikungunya cases for one state — the standard exploratory step before
# DLNM: it asks "at what lag is climate most correlated with cases?" rather
# than relying on eyeballing overlaid time series (see
# 11_plot_case_climate_timeseries.R, which motivated this — raw overlays
# can't reveal a lagged relationship, especially when outbreaks are episodic
# (2016-17, 2022) while the climate seasonal cycle repeats every year).
#
# Positive lag = climate leads cases by that many weeks (climate in week
# t-lag correlated with cases in week t) — the direction that matters
# biologically (mosquito life-cycle response to climate, then transmission).
#
# log1p(cases) is used instead of raw counts so the huge 2016-17/2022 spikes
# don't completely dominate the correlation at the expense of the smaller,
# more frequent fluctuations.
# ===========================================================================

for (p in c("here", "dplyr", "ggplot2", "patchwork")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(dplyr); library(ggplot2); library(patchwork)
})

COL_TEMP <- "#eb6834"
COL_PRCP <- "#1baf7a"
COL_SIG  <- "#d03b3b"
COL_GRID <- "#e1e0d9"
COL_AXIS <- "#898781"
COL_INK  <- "#0b0b0b"

theme_ccf <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.4),
      axis.line.x = element_line(color = COL_AXIS, linewidth = 0.4),
      axis.ticks = element_blank(),
      axis.text = element_text(color = COL_AXIS),
      axis.title = element_text(color = COL_INK, size = 10),
      plot.title = element_text(color = COL_INK, size = 11, face = "bold"),
      plot.margin = margin(4, 8, 4, 8)
    )
}

compute_ccf_df <- function(x, y, max_lag) {
  cc <- ccf(x, y, lag.max = max_lag, plot = FALSE, na.action = na.pass)
  n <- length(x)
  ci <- qnorm(0.975) / sqrt(n)  # standard 95% CI for a white-noise CCF
  tibble::tibble(
    lag = as.integer(cc$lag),
    corr = as.numeric(cc$acf),
    significant = abs(corr) > ci
  ) |>
    dplyr::mutate(ci = ci)
}

plot_case_climate_ccf <- function(panel, uf_code, uf_name = uf_code, max_lag = 16) {
  state_weekly <- panel |>
    dplyr::filter(substr(muni6, 1, 2) == uf_code) |>
    dplyr::group_by(week_start) |>
    dplyr::summarise(
      cases = sum(cases_confirmed, na.rm = TRUE),
      Tmean = weighted.mean(Tmean, w = population, na.rm = TRUE),
      PRCP  = weighted.mean(PRCP, w = population, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(week_start) |>
    dplyr::mutate(log_cases = log1p(cases))

  # ccf(x, y): positive lag = x leads y. We want "climate leads cases", so
  # x = climate, y = log_cases, and we keep only lag >= 0 in the plot
  # (climate leading or simultaneous) — a negative lag here would mean
  # cases predicting future climate, which isn't the mechanism of interest.
  ccf_temp <- compute_ccf_df(state_weekly$Tmean, state_weekly$log_cases, max_lag) |>
    dplyr::filter(lag >= 0)
  ccf_prcp <- compute_ccf_df(state_weekly$PRCP, state_weekly$log_cases, max_lag) |>
    dplyr::filter(lag >= 0)

  make_ccf_panel <- function(df, col, title) {
    ggplot(df, aes(lag, corr)) +
      geom_hline(yintercept = 0, color = COL_AXIS, linewidth = 0.4) +
      geom_hline(aes(yintercept = ci), linetype = "dashed", color = COL_AXIS, linewidth = 0.4) +
      geom_hline(aes(yintercept = -ci), linetype = "dashed", color = COL_AXIS, linewidth = 0.4) +
      geom_col(aes(fill = significant), width = 0.6, show.legend = FALSE) +
      scale_fill_manual(values = c(`TRUE` = col, `FALSE` = COL_GRID)) +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      labs(title = title, x = "Lag (weeks, climate leads cases)", y = "Cross-correlation") +
      theme_ccf()
  }

  p_temp <- make_ccf_panel(ccf_temp, COL_TEMP, sprintf("%s — Temperature vs log(cases)", uf_name))
  p_prcp <- make_ccf_panel(ccf_prcp, COL_PRCP, sprintf("%s — Precipitation vs log(cases)", uf_name))

  p_temp / p_prcp
}

if (sys.nframe() == 0) {
  UF_CODE <- "23"
  UF_NAME <- "Ceara (CE)"
  MAX_LAG <- 16

  panel <- readRDS(here::here("01_Data/chik_dlnm_panel_muni_week_2015_2024.rds"))
  fig <- plot_case_climate_ccf(panel, UF_CODE, UF_NAME, MAX_LAG)

  out_path <- here::here(sprintf("03_Output/figures/case_climate_ccf_%s.png", tolower(UF_CODE)))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_path, fig, width = 9, height = 7, dpi = 150)
  message("[save] ", out_path)
}
