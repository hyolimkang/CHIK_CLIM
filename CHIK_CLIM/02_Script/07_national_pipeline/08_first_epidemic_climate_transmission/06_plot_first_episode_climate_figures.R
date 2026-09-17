# First-epidemic climate-transmission figures.
#
# These are calibration figures for the first-epidemic GAM. The plotted
# climate response is a population-level transmission-potential prediction,
# not a newly observed state-specific R0.

required_packages <- c("here", "dplyr", "readr", "ggplot2", "mgcv", "scales", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(mgcv)
  library(scales)
})

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
source(file.path(root, "02_Script/07_national_pipeline/06_national_wave_analysis/00_national_wave_analysis_helpers.R"))
paths <- ensure_national_output_dirs()

theme_nm <- theme_classic(base_size = 9.5) +
  theme(
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_text(colour = "grey15"),
    axis.text = element_text(colour = "grey25"),
    plot.title = element_text(face = "bold", size = 10, margin = margin(b = 3)),
    plot.subtitle = element_text(colour = "grey35", size = 8, margin = margin(b = 6)),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", colour = "grey20")
  )

training <- read_csv(file.path(paths$table, "first_epidemic_climate_training_data.csv"), show_col_types = FALSE) |>
  mutate(
    UF = factor(UF),
    window_weeks = factor(window_weeks, levels = c(4, 6, 8), labels = c("4", "6", "8"))
  )
gam_fit <- readRDS(file.path(paths$table, "first_epidemic_climate_gam.rds"))
grid <- read_csv(file.path(paths$table, "climate_R0_prediction_grid.csv"), show_col_types = FALSE)

# ================================================================
# A. Input data: use window length, not overcrowded date axes, as the
# within-UF comparison. This is a calibration-data audit, not a main result.
# ================================================================
p_a <- ggplot(training, aes(window_weeks, Re_median)) +
  geom_hline(yintercept = 1, linetype = "22", colour = "grey55", linewidth = 0.35) +
  geom_linerange(aes(ymin = Re_lower, ymax = Re_upper), colour = "#486B87", linewidth = 0.45) +
  geom_point(colour = "#183B56", size = 1.45) +
  facet_wrap(~ UF, ncol = 6) +
  scale_y_continuous(breaks = c(1, 3, 5, 7), limits = c(0.7, 7.1), expand = c(0, 0)) +
  labs(
    title = "Early-growth estimates used for first-epidemic climate calibration",
    subtitle = "Each state contributes 4-, 6-, and 8-week estimates; vertical bars are 95% credible intervals",
    x = "Early-growth window (weeks)", y = expression(Early~R[e])
  ) +
  theme_nm + theme(strip.text = element_text(size = 7.5))
ggsave(file.path(paths$figure, "FIRST_EPISODE_RE_BY_UF.png"), p_a,
       width = 190, height = 150, units = "mm", dpi = 400, bg = "white")
message("[saved] FIRST_EPISODE_RE_BY_UF.png")

# ================================================================
# B-C. Marginal climate responses. Both curves use the same y scale so that
# the visibly wider precipitation uncertainty is directly comparable.
# ================================================================
precip_med <- median(training$precipitation)
temp_med <- median(training$temperature)

temp_seq <- seq(min(training$temperature), max(training$temperature), length.out = 150)
temp_prediction <- predict(
  gam_fit,
  newdata = data.frame(temperature = temp_seq, precipitation = precip_med),
  se.fit = TRUE, exclude = "s(UF)", newdata.guaranteed = TRUE
)
temp_curve <- tibble(
  temperature = temp_seq,
  fit = exp(as.numeric(temp_prediction$fit)),
  lo = exp(as.numeric(temp_prediction$fit) - 1.96 * as.numeric(temp_prediction$se.fit)),
  hi = exp(as.numeric(temp_prediction$fit) + 1.96 * as.numeric(temp_prediction$se.fit))
)

precip_seq <- seq(min(training$precipitation), max(training$precipitation), length.out = 150)
precip_prediction <- predict(
  gam_fit,
  newdata = data.frame(temperature = temp_med, precipitation = precip_seq),
  se.fit = TRUE, exclude = "s(UF)", newdata.guaranteed = TRUE
)
precip_curve <- tibble(
  precipitation = precip_seq,
  fit = exp(as.numeric(precip_prediction$fit)),
  lo = exp(as.numeric(precip_prediction$fit) - 1.96 * as.numeric(precip_prediction$se.fit)),
  hi = exp(as.numeric(precip_prediction$fit) + 1.96 * as.numeric(precip_prediction$se.fit))
)

response_limits <- range(c(temp_curve$lo, temp_curve$hi, precip_curve$lo, precip_curve$hi))
response_limits <- response_limits + c(-0.04, 0.04) * diff(response_limits)

p_b <- ggplot(temp_curve, aes(temperature, fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#C76B3B", alpha = 0.20) +
  geom_line(colour = "#A9471C", linewidth = 0.8) +
  geom_rug(data = training, aes(x = temperature), inherit.aes = FALSE,
           sides = "b", colour = "grey30", alpha = 0.45, length = unit(0.025, "npc")) +
  scale_y_continuous(limits = response_limits, expand = c(0, 0)) +
  labs(
    title = "Temperature response",
    subtitle = sprintf("Precipitation held at median (%.0f mm); shaded band: 95%% CI", precip_med),
    x = "Mean temperature over weeks t-2 to t (degrees C)",
    y = expression(widehat(R)[climate])
  ) +
  theme_nm
ggsave(file.path(paths$figure, "CLIMATE_R0_TEMPERATURE_RESPONSE.png"), p_b,
       width = 130, height = 100, units = "mm", dpi = 400, bg = "white")
message("[saved] CLIMATE_R0_TEMPERATURE_RESPONSE.png")

p_c <- ggplot(precip_curve, aes(precipitation, fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#4B8CB7", alpha = 0.20) +
  geom_line(colour = "#1F6F9E", linewidth = 0.8) +
  geom_rug(data = training, aes(x = precipitation), inherit.aes = FALSE,
           sides = "b", colour = "grey30", alpha = 0.45, length = unit(0.025, "npc")) +
  scale_y_continuous(limits = response_limits, expand = c(0, 0)) +
  labs(
    title = "Precipitation response",
    subtitle = sprintf("Temperature held at median (%.1f degrees C); shaded band: 95%% CI", temp_med),
    x = "Cumulative precipitation over weeks t-6 to t-2 (mm)",
    y = expression(widehat(R)[climate])
  ) +
  theme_nm
ggsave(file.path(paths$figure, "CLIMATE_R0_PRECIPITATION_RESPONSE.png"), p_c,
       width = 130, height = 100, units = "mm", dpi = 400, bg = "white")
message("[saved] CLIMATE_R0_PRECIPITATION_RESPONSE.png")

# ================================================================
# D. Two-dimensional transmission-potential surface.
# ================================================================
p_d <- ggplot(grid, aes(temperature, precipitation, fill = R0_climate_hat)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(data = grid, aes(temperature, precipitation, z = R0_climate_hat), inherit.aes = FALSE,
               colour = "white", alpha = 0.55, linewidth = 0.28) +
  geom_point(data = training, aes(temperature, precipitation), inherit.aes = FALSE,
             colour = "grey10", fill = "white", shape = 21, size = 1.65, stroke = 0.55) +
  scale_fill_gradientn(colours = c("#482878", "#2D708E", "#22A884", "#FDE725"),
                       name = expression(widehat(R)[climate])) +
  guides(fill = guide_colourbar(barheight = unit(42, "mm"), barwidth = unit(4, "mm"))) +
  labs(
    title = "First-epidemic climate transmission-potential surface",
    subtitle = "Open circles mark observed temperature-precipitation combinations used for calibration",
    x = "Mean temperature over weeks t-2 to t (degrees C)",
    y = "Cumulative precipitation over weeks t-6 to t-2 (mm)"
  ) +
  theme_nm + theme(legend.position = "right")
ggsave(file.path(paths$figure, "CLIMATE_TRANSMISSION_SURFACE.png"), p_d,
       width = 155, height = 122, units = "mm", dpi = 400, bg = "white")
message("[saved] CLIMATE_TRANSMISSION_SURFACE.png")

# ================================================================
# E. In-sample calibration. Window length is the only grouping retained;
# state labels are omitted because they obscure the pooled calibration cloud.
# ================================================================
pred_train <- predict(gam_fit, newdata = training, exclude = "s(UF)", newdata.guaranteed = TRUE)
training$R_climate_hat <- exp(as.numeric(pred_train))

p_e <- ggplot(training, aes(R_climate_hat, Re_median, shape = window_weeks)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey45", linewidth = 0.5) +
  geom_linerange(aes(ymin = Re_lower, ymax = Re_upper), colour = "#6B7D88", alpha = 0.48, linewidth = 0.35) +
  geom_point(colour = "#183B56", fill = "white", size = 2.25, stroke = 0.6) +
  scale_shape_manual(values = c("4" = 21, "6" = 22, "8" = 24), name = "Window (weeks)") +
  scale_x_continuous(limits = c(min(training$R_climate_hat) - 0.03, max(training$R_climate_hat) + 0.03)) +
  scale_y_continuous(limits = c(0.7, max(training$Re_upper) + 0.35), expand = c(0, 0)) +
  labs(
    title = "In-sample climate calibration",
    subtitle = "Population-level GAM prediction; vertical bars are 95% credible intervals; dashed line denotes y = x",
    x = expression(Climate-predicted~transmission~potential~~widehat(R)[climate]),
    y = expression(Observed~early~R[e])
  ) +
  theme_nm
ggsave(file.path(paths$figure, "FIRST_EPISODE_RE_PREDICTED_VS_OBSERVED.png"), p_e,
       width = 145, height = 122, units = "mm", dpi = 400, bg = "white")
message("[saved] FIRST_EPISODE_RE_PREDICTED_VS_OBSERVED.png")
