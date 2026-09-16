# Recurrent-epidemic climate x susceptibility figures -- Sections 6-10.
#
# This script only visualises existing episode-level estimates. It does not
# refit either the first-epidemic climate GAM or any long-term reconstruction.
# Primary figures are restricted to BA, RJ, and MT. CE remains a conditional
# fixed-q sensitivity analysis and is deliberately shown in its own figure.

required_packages <- c("here", "dplyr", "readr", "ggplot2", "tibble", "scales", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(scales)
  library(patchwork)
})

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
source(file.path(root, "02_Script/40_renewal_model/15_national_wave_analysis/00_national_wave_analysis_helpers.R"))
paths <- ensure_national_output_dirs()

# A quiet, journal-style theme: legends encode states; individual UF-year
# labels are intentionally omitted to keep the scientific pattern visible.
theme_nm <- theme_classic(base_size = 9.5) +
  theme(
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_text(colour = "grey15"),
    axis.text = element_text(colour = "grey25"),
    plot.title = element_text(face = "bold", size = 10, margin = margin(b = 3)),
    plot.subtitle = element_text(colour = "grey35", size = 8, margin = margin(b = 6)),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", colour = "grey20")
  )

state_labels <- c(BA = "Bahia", RJ = "Rio de Janeiro", MT = "Mato Grosso")
state_colours <- c("Bahia" = "#2F5D8A", "Rio de Janeiro" = "#B45F2B", "Mato Grosso" = "#2D7F6F")
state_shapes <- c("Bahia" = 21, "Rio de Janeiro" = 22, "Mato Grosso" = 24)

combined <- read_csv(file.path(paths$table, "RECURRENT_CLIMATE_SUSCEPTIBILITY_EPISODES.csv"), show_col_types = FALSE) |>
  mutate(
    episode_start = as.Date(episode_start),
    year = as.integer(format(episode_start, "%Y")),
    state_name = recode(UF, !!!state_labels),
    state_name = factor(state_name, levels = unname(state_labels))
  )

# Use reported cumulative incidence, rather than raw case totals, for the
# visual definition of epidemic size. Raw totals are not comparable across
# UFs with very different population sizes.
episode_size <- read_csv(file.path(paths$table, "brazil_chik_wave_analysis_master.csv"), show_col_types = FALSE) |>
  select(UF = state, episode_id = wave_id, reported_incidence_per_100k = cumulative_incidence_per_100k)

combined <- combined |>
  left_join(episode_size, by = c("UF", "episode_id"))
if (anyNA(combined$reported_incidence_per_100k)) {
  stop("Missing reported cumulative incidence for one or more recurrent episodes.")
}

# CE is conditional on a fixed reporting fraction. It is retained in the
# sensitivity panel below, but not mixed with the identified/partially
# identified primary-state comparison or its validation metrics.
primary_points <- combined |>
  filter(UF %in% names(state_labels)) |>
  mutate(
    Re_obs_median = as.numeric(Re_obs_median),
    reported_incidence_per_100k = as.numeric(reported_incidence_per_100k)
  )

threshold_curve <- tibble(
  R_climate_hat = seq(min(primary_points$R_climate_hat) * 0.94,
                      max(primary_points$R_climate_hat) * 1.04,
                      length.out = 300)
) |>
  mutate(S_threshold = 1 / R_climate_hat)

# Iso-Re_pred curves make the phase-plane interpretation explicit:
# Re_pred = R_climate_hat x S_pre/N.  The original threshold at one remains
# the shaded growth boundary, while higher curves show increasingly favourable
# predicted-growth regimes without changing any fitted quantity.
re_pred_contour_levels <- c(1, 1.5, 2, 2.5)
re_pred_contour_labels <- c("Predicted R[e] = 1", "Predicted R[e] = 1.5",
                            "Predicted R[e] = 2", "Predicted R[e] = 2.5")
re_pred_contours <- bind_rows(lapply(re_pred_contour_levels, function(level) {
  tibble(
    R_climate_hat = threshold_curve$R_climate_hat,
    S_pre_prop = level / threshold_curve$R_climate_hat,
    contour = factor(
      re_pred_contour_labels[match(level, re_pred_contour_levels)],
      levels = re_pred_contour_labels
    )
  )
})) |>
  filter(S_pre_prop >= 0.35, S_pre_prop <= 1.02)

# ================================================================
# Primary phase plane
# ================================================================
p_phase <- ggplot(primary_points, aes(R_climate_hat, S_pre_median)) +
  geom_ribbon(
    data = threshold_curve,
    aes(x = R_climate_hat, ymin = 0.35, ymax = S_threshold),
    inherit.aes = FALSE, fill = "#E7D8CE", alpha = 0.65
  ) +
  geom_ribbon(
    data = threshold_curve,
    aes(x = R_climate_hat, ymin = S_threshold, ymax = 1.02),
    inherit.aes = FALSE, fill = "#E7F0EC", alpha = 0.65
  ) +
  geom_line(
    data = re_pred_contours,
    aes(R_climate_hat, S_pre_prop, linetype = contour),
    inherit.aes = FALSE, colour = "grey25", linewidth = 0.55
  ) +
  geom_point(
    aes(size = reported_incidence_per_100k, fill = Re_obs_median, shape = state_name),
    colour = "grey15", stroke = 0.55, alpha = 0.95
  ) +
  scale_shape_manual(values = state_shapes, name = NULL) +
  scale_linetype_manual(
    values = c("Predicted R[e] = 1" = "22", "Predicted R[e] = 1.5" = "dashed",
               "Predicted R[e] = 2" = "solid", "Predicted R[e] = 2.5" = "dotdash"),
    name = expression(Predicted~R[e]), labels = c("1", "1.5", "2", "2.5")
  ) +
  scale_fill_gradient(low = "#E5EEF1", high = "#A63C1F", name = expression(Observed~early~R[e])) +
  scale_size_continuous(
    range = c(2.8, 10), name = "Reported cumulative\nincidence per 100,000",
    breaks = c(25, 100, 500, 1000), labels = label_number(accuracy = 1)
  ) +
  scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0.35, 1.02), expand = c(0, 0)) +
  labs(
    title = "Climate-susceptibility phase plane",
    subtitle = "Primary recurrent episodes only; curves denote climate x susceptibility predicted R[e]",
    x = expression(Climate-predicted~transmission~potential~~widehat(R)[climate]),
    y = expression(Pre-epidemic~susceptibility~~S[pre]/N)
  ) +
  theme_nm + theme(legend.box = "vertical")
ggsave(file.path(paths$figure, "RECURRENT_CLIMATE_SUSCEPTIBILITY_PHASE_PLANE.png"),
       p_phase, width = 174, height = 132, units = "mm", dpi = 400, bg = "white")
message("[saved] RECURRENT_CLIMATE_SUSCEPTIBILITY_PHASE_PLANE.png")

# ================================================================
# Direct outcome views: Re vs S/N and epidemic size vs S/N
# ================================================================
p_re_s <- ggplot(primary_points, aes(S_pre_median, Re_obs_median, colour = state_name, shape = state_name)) +
  geom_hline(yintercept = 1, colour = "grey55", linetype = "22", linewidth = 0.45) +
  geom_smooth(data = primary_points, aes(S_pre_median, Re_obs_median, group = 1), inherit.aes = FALSE,
              method = "lm", se = TRUE, colour = "grey30", fill = "grey75",
              alpha = 0.32, linewidth = 0.5, show.legend = FALSE) +
  geom_errorbar(aes(ymin = Re_obs_lower, ymax = Re_obs_upper), width = 0,
                alpha = 0.42, linewidth = 0.38, show.legend = FALSE) +
  geom_errorbarh(aes(xmin = S_pre_lower, xmax = S_pre_upper), height = 0,
                 alpha = 0.42, linewidth = 0.38, show.legend = FALSE) +
  geom_point(size = 2.9, stroke = 0.65, fill = "white") +
  scale_colour_manual(values = state_colours, name = NULL) +
  scale_shape_manual(values = state_shapes, name = NULL) +
  scale_x_continuous(labels = label_percent(accuracy = 1), limits = c(0.50, 1.01), expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_continuous(breaks = 1:4, limits = c(0.75, 4.65), expand = c(0, 0)) +
  labs(
    title = expression(Observed~early~R[e]~" versus "~S[pre]/N),
    subtitle = "Grey line: unadjusted linear summary; bars: 95% intervals",
    x = expression(Pre-epidemic~susceptibility~~S[pre]/N),
    y = expression(Observed~early~R[e])
  ) +
  theme_nm

p_size_s <- ggplot(primary_points, aes(S_pre_median, reported_incidence_per_100k, colour = state_name, shape = state_name)) +
  geom_smooth(data = primary_points, aes(S_pre_median, reported_incidence_per_100k, group = 1), inherit.aes = FALSE,
              method = "lm", se = TRUE, colour = "grey30", fill = "grey75",
              alpha = 0.32, linewidth = 0.5, show.legend = FALSE) +
  geom_errorbarh(aes(xmin = S_pre_lower, xmax = S_pre_upper), height = 0,
                 alpha = 0.42, linewidth = 0.38, show.legend = FALSE) +
  geom_point(size = 2.9, stroke = 0.65, fill = "white") +
  scale_colour_manual(values = state_colours, name = NULL) +
  scale_shape_manual(values = state_shapes, name = NULL) +
  scale_x_continuous(labels = label_percent(accuracy = 1), limits = c(0.50, 1.01), expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_log10(labels = label_number(accuracy = 1), breaks = c(10, 30, 100, 300, 1000)) +
  labs(
    title = expression(Reported~epidemic~size~" versus "~S[pre]/N),
    subtitle = "Cumulative reported incidence; logarithmic scale",
    x = expression(Pre-epidemic~susceptibility~~S[pre]/N),
    y = "Reported cumulative incidence per 100,000"
  ) +
  theme_nm

p_outcomes <- p_re_s + p_size_s +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Recurrent epidemic expression and pre-epidemic susceptibility",
    subtitle = "Primary states only (Bahia, Rio de Janeiro, Mato Grosso)"
  ) &
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 11))

ggsave(file.path(paths$figure, "RECURRENT_CLIMATE_SUSCEPTIBILITY_PHASE_OUTCOMES.png"),
       p_outcomes, width = 220, height = 110, units = "mm", dpi = 400, bg = "white")
message("[saved] RECURRENT_CLIMATE_SUSCEPTIBILITY_PHASE_OUTCOMES.png")

# ================================================================
# CE fixed-q sensitivity: separate conditional display, faceted by episode
# ================================================================
ce_points <- combined |>
  filter(UF == "CE") |>
  mutate(
    q_scenario = factor(sprintf("q = %.2f", as.numeric(q_scenario)), levels = sprintf("q = %.2f", c(0.05, 0.10, 0.15, 0.20))),
    episode_year = format(episode_start, "%Y")
  )

p_ce_sensitivity <- ggplot(ce_points, aes(R_climate_hat, S_pre_median, colour = q_scenario, group = episode_id)) +
  geom_line(data = threshold_curve, aes(R_climate_hat, S_threshold), inherit.aes = FALSE,
            colour = "grey30", linetype = "22", linewidth = 0.45) +
  geom_line(colour = "grey70", linewidth = 0.35, show.legend = FALSE) +
  geom_point(size = 2.15) +
  facet_wrap(~ episode_year, ncol = 3) +
  scale_colour_manual(values = c("q = 0.05" = "#3B528B", "q = 0.10" = "#21918C", "q = 0.15" = "#5EC962", "q = 0.20" = "#FDE725"), name = "Fixed reporting fraction") +
  scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0.35, 1.02), expand = c(0, 0)) +
  labs(
    title = "Ceara conditional sensitivity across the fixed-q grid",
    subtitle = "Each panel is one episode; grey paths connect the four conditional reconstructions",
    x = expression(widehat(R)[climate]), y = expression(S[pre]/N)
  ) +
  theme_nm + theme(legend.position = "bottom")
ggsave(file.path(paths$figure, "RECURRENT_CE_FIXED_Q_SENSITIVITY_PANEL.png"),
       p_ce_sensitivity, width = 174, height = 125, units = "mm", dpi = 400, bg = "white")
message("[saved] RECURRENT_CE_FIXED_Q_SENSITIVITY_PANEL.png")

# ================================================================
# Direct mechanistic validation (primary states only)
# ================================================================
val_data <- primary_points
p_validation <- ggplot(val_data, aes(Re_pred_median, Re_obs_median, colour = state_name, shape = state_name)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey35", linewidth = 0.55) +
  geom_point(size = 3.1, stroke = 0.65, fill = "white") +
  scale_colour_manual(values = state_colours, name = NULL) +
  scale_shape_manual(values = state_shapes, name = NULL) +
  coord_equal(xlim = c(0.8, 4.5), ylim = c(0.8, 4.5), expand = FALSE) +
  labs(
    title = "Mechanistic prediction versus observed early growth",
    subtitle = expression("Points are posterior medians; dashed line denotes "~y == x),
    x = expression(Predicted~R[e]~~widehat(R)[climate]~times~S[pre]/N),
    y = expression(Observed~early~R[e])
  ) +
  theme_nm
ggsave(file.path(paths$figure, "RECURRENT_RE_PREDICTED_VS_OBSERVED.png"),
       p_validation, width = 128, height = 128, units = "mm", dpi = 400, bg = "white")
message("[saved] RECURRENT_RE_PREDICTED_VS_OBSERVED.png")

# Metrics are descriptive and now match the declared primary analysis. CE is
# intentionally not folded into these values because it is conditional on q.
cor_val <- cor(val_data$Re_pred_median, val_data$Re_obs_median)
calib_fit <- lm(Re_obs_median ~ Re_pred_median, data = val_data)
mae <- mean(abs(val_data$Re_pred_median - val_data$Re_obs_median))
rmse <- sqrt(mean((val_data$Re_pred_median - val_data$Re_obs_median)^2))
classification_agree <- mean((val_data$Re_pred_median > 1) == (val_data$Re_obs_median > 1))

mae_A <- mean(abs(val_data$R_climate_hat - val_data$Re_obs_median))
rmse_A <- sqrt(mean((val_data$R_climate_hat - val_data$Re_obs_median)^2))
cor_A <- cor(val_data$R_climate_hat, val_data$Re_obs_median)
classification_agree_A <- mean((val_data$R_climate_hat > 1) == (val_data$Re_obs_median > 1))

validation_summary <- tibble(
  model = c("A: climate-only (R_climate_hat)", "B: climate x susceptibility (Re_pred)"),
  n_episodes = nrow(val_data),
  correlation_with_Re_obs = c(cor_A, cor_val),
  MAE = c(mae_A, mae),
  RMSE = c(rmse_A, rmse),
  classification_agreement_Re_gt1 = c(classification_agree_A, classification_agree)
)
print(as.data.frame(validation_summary), digits = 3)
message(sprintf("Calibration slope (Re_obs ~ Re_pred): %.3f (intercept %.3f)", coef(calib_fit)[2], coef(calib_fit)[1]))
write_csv(validation_summary, file.path(paths$table, "RECURRENT_RE_VALIDATION_METRICS.csv"))
message("[saved] RECURRENT_RE_VALIDATION_METRICS.csv")
