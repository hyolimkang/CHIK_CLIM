# Figures 1-6 for the provisional Re x climate x susceptibility checkpoint
# (task Section 15). Visualisation-only -- reads already-produced tables and
# does not recompute or alter any posterior/model quantity.

required_fig_packages <- c("here", "dplyr", "readr", "tidyr", "tibble", "ggplot2", "mgcv", "scales")
missing_fig_packages <- required_fig_packages[
  !vapply(required_fig_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_fig_packages)) stop("Missing package(s): ", paste(missing_fig_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(tibble); library(ggplot2); library(mgcv); library(scales) })

fig_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}
source(file.path(fig_root(), "02_Script", "07_national_pipeline", "06_national_wave_analysis", "00_national_wave_analysis_helpers.R"))

theme_episode <- function() {
  theme_minimal(base_size = 12) + theme(strip.text = element_text(face = "bold", size = 8), legend.position = "bottom")
}

run_plot_provisional_episode_figures <- function() {
  paths <- ensure_national_output_dirs()
  episode <- read_csv(file.path(paths$table, "brazil_chik_Re_S_climate_episode_analysis_scale1_0.csv"), show_col_types = FALSE) |>
    mutate(wave_type = if_else(is_recurrence, "Recurrent wave", "First wave"), log1p_precip = log1p(precip_pre6_2))

  m_climate <- gam(log(Re_early_6_median) ~ s(temp_pre6_2, k = 4) + s(log1p_precip, k = 4) + s(state, bs = "re"),
                    data = episode |> mutate(state = factor(state)), method = "REML")
  m_climate_s <- gam(log(Re_early_6_median) ~ s(temp_pre6_2, k = 4) + s(log1p_precip, k = 4) + s(S_prop, k = 3) + s(state, bs = "re"),
                      data = episode |> mutate(state = factor(state)), method = "REML")
  m_gate <- gam(log(Re_early_6_median) ~ offset(log(S_prop)) + s(temp_pre6_2, k = 4) + s(log1p_precip, k = 4) + s(state, bs = "re"),
                data = episode |> mutate(state = factor(state)), method = "REML")

  # ---- Figure 1: Re vs S_prop ----
  f1 <- ggplot(episode, aes(x = S_prop, y = Re_early_6_median, colour = wave_type)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = .7) +
    geom_smooth(aes(group = 1), method = "loess", colour = "black", se = TRUE, linewidth = .6) +
    scale_y_log10() +
    scale_colour_manual(values = c("First wave" = "#1B9E77", "Recurrent wave" = "#D95F02")) +
    labs(title = "Figure 1: Early Re vs. wave-onset susceptible fraction",
         subtitle = "Each point = one major epidemic wave (n=117); log10 y-axis; dashed line = Re=1",
         x = "S_prop at onset (M0, weekly)", y = "Re_early_6 (median)", colour = NULL) +
    theme_episode()
  ggsave(file.path(paths$figure, "brazil_chik_episode_figure1_Re_vs_S.png"), f1, width = 190, height = 150, units = "mm", dpi = 300, bg = "white")

  # ---- Figure 2: Re vs temperature and precipitation ----
  f2a <- ggplot(episode, aes(x = temp_pre6_2, y = Re_early_6_median, colour = wave_type)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = .7) + geom_smooth(aes(group = 1), method = "loess", colour = "black", se = TRUE, linewidth = .6) +
    scale_y_log10() + scale_colour_manual(values = c("First wave" = "#1B9E77", "Recurrent wave" = "#D95F02")) +
    labs(x = "Pre-onset mean temperature (weeks -6:-2, C)", y = "Re_early_6 (median)", colour = NULL) + theme_episode()
  f2b <- ggplot(episode, aes(x = precip_pre6_2, y = Re_early_6_median, colour = wave_type)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = .7) + geom_smooth(aes(group = 1), method = "loess", colour = "black", se = TRUE, linewidth = .6) +
    scale_y_log10() + scale_colour_manual(values = c("First wave" = "#1B9E77", "Recurrent wave" = "#D95F02")) +
    labs(x = "Pre-onset summed precipitation (weeks -6:-2, mm)", y = NULL, colour = NULL) + theme_episode()
  f2 <- patchwork_or_grid(f2a, f2b, title = "Figure 2: Early Re vs. pre-onset climate")
  ggsave(file.path(paths$figure, "brazil_chik_episode_figure2_Re_vs_climate.png"), f2, width = 260, height = 150, units = "mm", dpi = 300, bg = "white")

  # ---- Figure 3: partial effects from M_climate_S ----
  # plot.gam() always draws when called; redirect to a null device so this
  # batch script doesn't leave a stray Rplots.pdf behind. The function
  # returns the per-smooth plot data (term, fitted partial effect, se)
  # regardless of the drawing target.
  pdf(NULL)
  plot_data <- plot(m_climate_s, n = 200)
  dev.off()
  # s(state, bs="re") has a character/factor $x (state levels), not a
  # continuous smooth -- exclude it INSIDE the loop (before bind_rows) since
  # combining its character $x with the other terms' numeric $x errors.
  keep_terms <- c("temp_pre6_2", "log1p_precip", "S_prop")
  partials <- bind_rows(lapply(plot_data, function(pd) {
    if (is.null(pd$xlab) || !(pd$xlab %in% keep_terms)) return(NULL)
    tibble(term = pd$xlab, x = as.numeric(pd$x), fit = pd$fit, se = pd$se)
  }))
  f3 <- ggplot(partials, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), alpha = .2) +
    geom_line() +
    facet_wrap(~term, scales = "free_x", nrow = 1,
               labeller = as_labeller(c(temp_pre6_2 = "Temperature (pre-onset)", log1p_precip = "log1p(Precipitation)", S_prop = "S_prop"))) +
    labs(title = "Figure 3: Partial effects from M_climate_S", subtitle = "95% CI (grey band); y = partial effect on log(Re_early_6)",
         x = NULL, y = "Partial effect") +
    theme_episode()
  ggsave(file.path(paths$figure, "brazil_chik_episode_figure3_partial_effects.png"), f3, width = 260, height = 120, units = "mm", dpi = 300, bg = "white")

  # ---- Figure 4: observed vs predicted, three models ----
  obs_pred <- bind_rows(
    tibble(model = "M_climate", observed = episode$Re_early_6_median, predicted = exp(fitted(m_climate))),
    tibble(model = "M_climate_S", observed = episode$Re_early_6_median, predicted = exp(fitted(m_climate_s))),
    tibble(model = "M_mechanistic_gate", observed = episode$Re_early_6_median, predicted = exp(fitted(m_gate)))
  )
  f4 <- ggplot(obs_pred, aes(x = observed, y = predicted)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = .6) +
    scale_x_log10() + scale_y_log10() +
    facet_wrap(~model, nrow = 1) +
    labs(title = "Figure 4: Observed vs. fitted Re_early_6 (in-sample)", x = "Observed", y = "Fitted") +
    theme_episode()
  ggsave(file.path(paths$figure, "brazil_chik_episode_figure4_obs_vs_pred.png"), f4, width = 260, height = 110, units = "mm", dpi = 300, bg = "white")

  # ---- Figure 5: R0_climate and predicted Re = S_prop * R0_climate ----
  episode_gate <- episode |> mutate(
    R0_climate = exp(predict(m_gate, newdata = episode |> mutate(state = factor(state)), exclude = grep("^s\\(state\\)", vapply(m_gate$smooth, function(s) s$label, character(1)), value = TRUE)) - log(S_prop)),
    predicted_Re = S_prop * R0_climate
  )
  f5a <- ggplot(episode_gate, aes(x = temp_pre6_2, y = R0_climate, colour = precip_pre6_2)) +
    geom_point(alpha = .8) + scale_colour_viridis_c(option = "D", name = "Precip\n(mm)") +
    labs(x = "Pre-onset temperature (C)", y = expression(Inferred~R[0]^climate)) + theme_episode()
  f5b <- ggplot(episode_gate, aes(x = predicted_Re, y = Re_early_6_median)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(aes(colour = S_prop), alpha = .8) + scale_colour_viridis_c(option = "C", name = "S_prop") +
    scale_x_log10() + scale_y_log10() +
    labs(x = expression(Predicted~Re == S[prop]%*%R[0]^climate), y = "Observed Re_early_6") + theme_episode()
  f5 <- patchwork_or_grid(f5a, f5b, title = "Figure 5: Mechanistic gate -- R0_climate and S-modified prediction")
  ggsave(file.path(paths$figure, "brazil_chik_episode_figure5_mechanistic_gate.png"), f5, width = 260, height = 150, units = "mm", dpi = 300, bg = "white")

  # ---- Figure 6: model comparison across FOI scale scenarios ----
  comparison <- read_csv(file.path(paths$table, "brazil_chik_provisional_gam_comparison.csv"), show_col_types = FALSE) |>
    dplyr::filter(climate_window == "pre6_2_PRIMARY") |>
    mutate(FOI_scale_scenario = factor(FOI_scale_scenario, levels = c("0_5", "1_0", "2_0"),
                                        labels = c("0.5x (low)", "1.0x (baseline)", "2.0x (high)")))
  f6 <- ggplot(comparison, aes(x = FOI_scale_scenario, y = leave_one_state_out_RMSE, colour = model_name, group = model_name)) +
    geom_line() + geom_point(size = 2) +
    labs(title = "Figure 6: Model comparison across long-term FOI-scale scenarios",
         subtitle = "Primary climate window (-6:-2 weeks); leave-one-state-out RMSE",
         x = "Long-term FOI scale scenario", y = "Leave-one-state-out RMSE (log Re scale)", colour = NULL) +
    theme_episode()
  ggsave(file.path(paths$figure, "brazil_chik_episode_figure6_foi_scale_comparison.png"), f6, width = 190, height = 150, units = "mm", dpi = 300, bg = "white")

  manifest <- tibble(
    figure = c("1_Re_vs_S", "2_Re_vs_climate", "3_partial_effects", "4_obs_vs_pred", "5_mechanistic_gate", "6_foi_scale_comparison"),
    path = file.path(paths$figure, c(
      "brazil_chik_episode_figure1_Re_vs_S.png", "brazil_chik_episode_figure2_Re_vs_climate.png",
      "brazil_chik_episode_figure3_partial_effects.png", "brazil_chik_episode_figure4_obs_vs_pred.png",
      "brazil_chik_episode_figure5_mechanistic_gate.png", "brazil_chik_episode_figure6_foi_scale_comparison.png"
    ))
  )
  write_csv(manifest, file.path(paths$figure, "brazil_chik_episode_figure_manifest.csv"))
  message("[episode-figures] done, manifest: ", file.path(paths$figure, "brazil_chik_episode_figure_manifest.csv"))
  invisible(manifest)
}

# Minimal 2-panel combiner avoiding a hard dependency on the patchwork
# package (not confirmed installed elsewhere in this repo).
patchwork_or_grid <- function(p1, p2, title) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    (p1 + p2) + patchwork::plot_annotation(title = title)
  } else {
    if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
    gridExtra::grid.arrange(p1, p2, ncol = 2, top = title)
  }
}

if (sys.nframe() == 0L) run_plot_provisional_episode_figures()
