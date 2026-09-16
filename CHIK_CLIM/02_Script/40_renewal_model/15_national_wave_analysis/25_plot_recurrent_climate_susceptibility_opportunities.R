# All-opportunities recurrent climate x susceptibility figures.
#
# This script consumes the additive weekly opportunity panel built in 24.
# It does not refit the first-epidemic climate GAM or any state-level
# susceptibility reconstruction.  The outcome is deliberately a frozen
# major-wave onset, not an inferred threshold chosen after examining these
# figures.  See the accompanying report for the exact non-outbreak definition
# and limitations of retrospective S/N.

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

root <- here::here()
if (!file.exists(file.path(root, "CHIK_CLIM.Rproj"))) root <- file.path(root, "CHIK_CLIM")
if (!file.exists(file.path(root, "CHIK_CLIM.Rproj"))) stop("Could not identify the inner CHIK_CLIM R-project root.")
source(file.path(root, "02_Script/40_renewal_model/15_national_wave_analysis/00_national_wave_analysis_helpers.R"))
paths <- ensure_national_output_dirs()

PRIMARY_STATES <- c("BA", "RJ", "MT")
state_labels <- c(BA = "Bahia", RJ = "Rio de Janeiro", MT = "Mato Grosso")
state_shapes <- c("Bahia" = 21, "Rio de Janeiro" = 22, "Mato Grosso" = 24)
state_colours <- c("Bahia" = "#2F5D8A", "Rio de Janeiro" = "#B45F2B", "Mato Grosso" = "#2D7F6F")

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

opportunity_path <- file.path(paths$table, "RECURRENT_CLIMATE_SUSCEPTIBILITY_OPPORTUNITIES_WEEKLY.csv")
if (!file.exists(opportunity_path)) {
  stop("Weekly opportunity panel is missing. Run 24_build_recurrent_weekly_opportunity_panel.R first.")
}

opportunities <- read_csv(opportunity_path, show_col_types = FALSE) |>
  mutate(
    anchor_week = as.Date(anchor_week),
    S_week = as.Date(S_week),
    outcome_followup_start = as.Date(outcome_followup_start),
    outcome_followup_end_6w = as.Date(outcome_followup_end_6w),
    outcome_followup_end_8w = as.Date(outcome_followup_end_8w),
    outcome_episode_onset = as.Date(outcome_episode_onset),
    outcome_episode_end = as.Date(outcome_episode_end),
    state_name = factor(recode(state, !!!state_labels), levels = unname(state_labels)),
    outcome_status_8w = factor(
      outcome_status_8w,
      levels = c("No major wave within 8 weeks", "Major wave within 8 weeks", "Right-censored follow-up")
    )
  )

if (!all(PRIMARY_STATES %in% unique(opportunities$state))) {
  stop("The opportunity panel does not contain all three primary states (BA, RJ, MT).")
}

# All phase-plane points have a complete 8-week ascertainment window.  Thus
# no right-censored row is implicitly treated as a non-outbreak control.
phase_points <- opportunities |>
  filter(
    state %in% PRIMARY_STATES,
    outcome_followup_complete_8w,
    !is.na(R_climate_hat), !is.na(S_pre_median), !is.na(R_pred_median)
  ) |>
  arrange(major_wave_onset_next_8w, state, anchor_week)
if (!nrow(phase_points)) stop("No complete 8-week opportunity windows are available for the phase plane.")

# Sliding windows linked to one future wave are deliberately not plotted as
# repeated large symbols.  All non-outbreak opportunities remain in the
# two-dimensional count surface; each outcome is represented once, using the
# latest eligible anchor before that wave's onset.  This preserves the full
# negative-opportunity denominator while avoiding pseudo-replication in the
# visual encoding of a single observed epidemic.
non_outbreak_points <- phase_points |>
  filter(!major_wave_onset_next_8w)
future_wave_points <- phase_points |>
  filter(major_wave_onset_next_8w) |>
  group_by(outcome_episode_id) |>
  arrange(desc(anchor_week), .by_group = TRUE) |>
  slice_head(n = 1L) |>
  ungroup()
if (anyNA(future_wave_points$outcome_Re_early_6_median) ||
    anyNA(future_wave_points$outcome_reported_cumulative_incidence_per_100k)) {
  stop("A representative future major-wave point lacks early-Re or incidence information.")
}

# Iso-Re_pred curves: S/N = c / R_climate.  Values are analytical contours,
# not estimates fitted to the outcome-coded panel.
x_seq <- seq(min(phase_points$R_climate_hat) * 0.985,
             max(phase_points$R_climate_hat) * 1.015, length.out = 400)
y_min <- max(0, min(0.35, floor(min(phase_points$S_pre_median) * 20) / 20 - 0.03))
y_max <- min(1.02, max(1.02, ceiling(max(phase_points$S_pre_median) * 20) / 20 + 0.03))
contour_levels <- c(1, 1.5, 2, 2.5)
contour_labels <- c("Predicted R[e] = 1", "Predicted R[e] = 1.5",
                    "Predicted R[e] = 2", "Predicted R[e] = 2.5")
contours <- bind_rows(lapply(contour_levels, function(level) {
  tibble(
    R_climate_hat = x_seq,
    S_pre_median = level / x_seq,
    contour = factor(
      contour_labels[match(level, contour_levels)],
      levels = contour_labels
    )
  )
})) |>
  filter(S_pre_median >= y_min, S_pre_median <= y_max)
threshold_curve <- tibble(R_climate_hat = x_seq, S_threshold = 1 / x_seq)

phase_background <- list(
  geom_ribbon(
    data = threshold_curve,
    aes(x = R_climate_hat, ymin = y_min, ymax = S_threshold),
    inherit.aes = FALSE, fill = "#E7D8CE", alpha = 0.65
  ),
  geom_ribbon(
    data = threshold_curve,
    aes(x = R_climate_hat, ymin = S_threshold, ymax = y_max),
    inherit.aes = FALSE, fill = "#E7F0EC", alpha = 0.65
  )
)

phase_scales <- list(
  scale_shape_manual(values = state_shapes, name = NULL),
  scale_linetype_manual(
    values = c("Predicted R[e] = 1" = "22", "Predicted R[e] = 1.5" = "dashed",
               "Predicted R[e] = 2" = "solid", "Predicted R[e] = 2.5" = "dotdash"),
    name = expression(Predicted~R[e]), labels = c("1", "1.5", "2", "2.5")
  ),
  scale_size_continuous(
    trans = "sqrt", range = c(3, 10),
    name = "Subsequent reported cumulative\nincidence per 100,000",
    breaks = c(25, 100, 500, 1000), labels = label_number(accuracy = 1)
  ),
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))),
  scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(y_min, y_max), expand = c(0, 0))
)

p_outbreak <- ggplot() +
  phase_background +
  stat_bin2d(
    data = non_outbreak_points,
    aes(R_climate_hat, S_pre_median, fill = after_stat(count)),
    bins = 16, colour = "white", linewidth = 0.12, alpha = 0.74
  ) +
  geom_line(
    data = contours,
    aes(R_climate_hat, S_pre_median, linetype = contour),
    colour = "grey25", linewidth = 0.5
  ) +
  geom_point(
    data = future_wave_points,
    aes(
      R_climate_hat, S_pre_median,
      size = outcome_reported_cumulative_incidence_per_100k,
      colour = outcome_Re_early_6_median,
      shape = state_name
    ),
    fill = "white", stroke = 1.05, alpha = 0.96
  ) +
  scale_fill_gradient(
    low = "#F7FBFF", high = "#6A8FA7",
    name = "Non-outbreak\nopportunity count"
  ) +
  scale_colour_gradient(
    low = "#F4A582", high = "#B2182B",
    name = expression(Observed~early~R[e])
  ) +
  phase_scales +
  labs(
    title = "All recurrent climate-susceptibility opportunities",
    subtitle = sprintf(
      "Blue bins: %d complete non-outbreak windows; outlined points: %d distinct future major waves (one latest eligible window each)",
      nrow(non_outbreak_points), nrow(future_wave_points)
    ),
    x = expression(Climate-predicted~transmission~potential~~widehat(R)[climate]),
    y = expression(Pre-window~susceptibility~~S[pre]/N)
  ) +
  theme_nm

p_phase <- p_outbreak +
  labs(
    title = "Recurrent epidemic opportunities: climate potential, susceptibility, and subsequent expression",
    subtitle = sprintf(
      "Blue bins: %d complete non-outbreak windows; outlines: %d future waves; contours show predicted R[e] = R_climate x S_pre/N",
      nrow(non_outbreak_points), nrow(future_wave_points)
    ),
    caption = paste(
      "How to read: a wave marker above an iso-R[e] curve has predicted R[e] greater than that label.",
      "Blue density in the same region shows that climate x susceptibility alone is not sufficient for a major wave.",
      sep = "\n"
    )
  ) +
  guides(
    colour = guide_colourbar(order = 1, barheight = unit(20, "mm"), barwidth = unit(3.5, "mm")),
    size = guide_legend(order = 2),
    shape = guide_legend(order = 3),
    fill = guide_colourbar(order = 4, barheight = unit(20, "mm"), barwidth = unit(3.5, "mm")),
    linetype = guide_legend(order = 5)
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.caption = element_text(size = 7.8, colour = "grey30", hjust = 0, margin = margin(t = 7)),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(3, "pt")
  )

ggsave(
  file.path(paths$figure, "RECURRENT_ALL_OPPORTUNITIES_PHASE_PLANE_FINAL.png"),
  p_phase, width = 235, height = 160, units = "mm", dpi = 400, bg = "white"
)
message("[saved] RECURRENT_ALL_OPPORTUNITIES_PHASE_PLANE_FINAL.png")

# Time trajectories retain all supported recurrent windows, including the
# right-censored tail, because they are descriptive trajectories rather than
# outcome-classified controls.  Major-wave bands use the same frozen census as
# the phase-plane outcome; they are not re-detected from the case series here.
trajectory_points <- opportunities |>
  filter(state %in% PRIMARY_STATES, !is.na(R_climate_hat), !is.na(S_pre_median), !is.na(R_pred_median)) |>
  select(state, state_name, anchor_week, R_climate_hat, S_pre_median, R_pred_median)

trajectory_long <- bind_rows(
  trajectory_points |>
    transmute(state, state_name, anchor_week, metric = "Climate transmission potential", value = R_climate_hat),
  trajectory_points |>
    transmute(state, state_name, anchor_week, metric = "Susceptible fraction", value = S_pre_median),
  trajectory_points |>
    transmute(state, state_name, anchor_week, metric = "Predicted effective reproduction number", value = R_pred_median)
) |>
  mutate(
    metric = factor(
      metric,
      levels = c("Climate transmission potential", "Susceptible fraction", "Predicted effective reproduction number")
    )
  )

master <- read_csv(file.path(paths$table, "brazil_chik_wave_analysis_master.csv"), show_col_types = FALSE) |>
  mutate(across(c(onset_week, end_week), as.Date))
trajectory_end <- max(trajectory_points$anchor_week)
major_wave_bands <- master |>
  filter(state %in% PRIMARY_STATES) |>
  transmute(
    state,
    state_name = factor(recode(state, !!!state_labels), levels = unname(state_labels)),
    onset_week,
    end_week = coalesce(end_week, trajectory_end)
  )
metric_levels <- levels(trajectory_long$metric)
major_wave_bands <- bind_rows(lapply(metric_levels, function(metric) {
  major_wave_bands |> mutate(metric = factor(metric, levels = metric_levels))
}))

p_trajectory <- ggplot(trajectory_long, aes(anchor_week, value, colour = state_name)) +
  geom_rect(
    data = major_wave_bands,
    aes(xmin = onset_week, xmax = end_week, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "grey70", alpha = 0.30
  ) +
  geom_hline(
    data = tibble(metric = factor("Predicted effective reproduction number", levels = metric_levels), value = 1),
    aes(yintercept = value), colour = "grey28", linetype = "22", linewidth = 0.45
  ) +
  geom_line(linewidth = 0.38, show.legend = FALSE) +
  facet_grid(metric ~ state_name, scales = "free_y", switch = "y") +
  scale_colour_manual(values = state_colours, guide = "none") +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "Climate potential, susceptibility, and predicted epidemic expression through time",
    subtitle = "Grey bands: frozen major-wave intervals; S_pre is before each six-week window; lower dashed line: predicted R[e] = 1.",
    x = NULL, y = NULL
  ) +
  theme_nm +
  theme(
    legend.position = "none",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 8),
    panel.spacing = unit(4, "pt")
  )

ggsave(
  file.path(paths$figure, "RECURRENT_CLIMATE_SUSCEPTIBILITY_TRAJECTORIES.png"),
  p_trajectory, width = 220, height = 175, units = "mm", dpi = 400, bg = "white"
)
message("[saved] RECURRENT_CLIMATE_SUSCEPTIBILITY_TRAJECTORIES.png")

outcome_summary <- phase_points |>
  group_by(state, state_name, outcome_status_8w) |>
  summarise(
    n_opportunities = n(),
    median_R_climate = median(R_climate_hat),
    median_S_pre_prop = median(S_pre_median),
    median_R_pred = median(R_pred_median),
    .groups = "drop"
  ) |>
  arrange(state, outcome_status_8w)
write_csv(outcome_summary, file.path(paths$table, "RECURRENT_ALL_OPPORTUNITIES_OUTCOME_SUMMARY.csv"))
message("[saved] RECURRENT_ALL_OPPORTUNITIES_OUTCOME_SUMMARY.csv")
