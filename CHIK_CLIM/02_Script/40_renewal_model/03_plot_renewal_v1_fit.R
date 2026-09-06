# ===========================================================================
# 03_plot_renewal_v1_fit.R
#
# Purpose
# -------
# Visual check for the Stage 1 renewal fit: observed weekly cases vs the
# model's posterior median I_t (+ 95% credible band), and the posterior for
# the single constant R. Confirms in a figure what the diagnostics already
# suggested — a constant R can't reproduce the epidemic's rise/peak/decline
# shape on its own, so the large process-noise term is doing that work
# instead (motivating Stage 2's time-varying R).
# ===========================================================================

for (p in c("here", "rstan", "ggplot2", "dplyr", "patchwork")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(rstan); library(ggplot2); library(dplyr); library(patchwork)
})

source(here::here("02_Script/40_renewal_model/01_build_ceara_state_weekly.R"))

COL_OBS <- "#2a78d6"
COL_FIT <- "#eb6834"
COL_GRID <- "#e1e0d9"
COL_AXIS <- "#898781"
COL_INK <- "#0b0b0b"

theme_fit <- function() {
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
      legend.position = "top",
      legend.title = element_blank()
    )
}

if (sys.nframe() == 0) {
  DATE_START <- as.Date("2015-06-01")
  DATE_END   <- as.Date("2018-06-30")

  panel <- readRDS(here::here("01_Data/chik_dlnm_panel_muni_week_2015_2024.rds"))
  ceara_weekly <- build_state_weekly(panel, "23") |>
    dplyr::filter(week_start >= DATE_START, week_start <= DATE_END)

  fit <- readRDS(here::here("02_Script/stan/renewal_ceara_v1.rds"))
  I_draws <- rstan::extract(fit, pars = "I")$I  # [draws, N]

  fit_summary <- tibble::tibble(
    week_start = ceara_weekly$week_start,
    observed = ceara_weekly$cases,
    I_median = apply(I_draws, 2, median),
    I_lo = apply(I_draws, 2, quantile, 0.025),
    I_hi = apply(I_draws, 2, quantile, 0.975)
  )

  p_fit <- ggplot(fit_summary, aes(week_start)) +
    geom_ribbon(aes(ymin = I_lo, ymax = I_hi), fill = COL_FIT, alpha = 0.2) +
    geom_line(aes(y = I_median, color = "Model I_t (median, 95% CrI)"), linewidth = 0.8) +
    geom_point(aes(y = observed, color = "Observed cases"), size = 1) +
    scale_color_manual(values = c("Observed cases" = COL_OBS, "Model I_t (median, 95% CrI)" = COL_FIT)) +
    labs(title = "Stage 1 renewal fit — Ceara 2016-17 wave (constant R, no depletion)",
         x = NULL, y = "Weekly cases") +
    theme_fit()

  R_draws <- rstan::extract(fit, pars = "R")$R
  p_R <- ggplot(tibble::tibble(R = R_draws), aes(R)) +
    geom_histogram(bins = 40, fill = COL_FIT, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = COL_AXIS) +
    labs(title = "Posterior: constant R", x = "R", y = "Draws") +
    theme_fit()

  fig <- p_fit / p_R + patchwork::plot_layout(heights = c(1.5, 1))

  out_path <- here::here("03_Output/figures/renewal_v1_fit_ceara.png")
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_path, fig, width = 9, height = 8, dpi = 150)
  message("[save] ", out_path)
}
