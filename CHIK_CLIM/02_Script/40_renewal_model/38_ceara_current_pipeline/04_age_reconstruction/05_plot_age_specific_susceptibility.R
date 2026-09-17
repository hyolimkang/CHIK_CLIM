# Age-structured vaccine extension -- Figures A-D + summary table for the
# age-specific susceptibility/infection diagnostic (descriptive only; no
# vaccination). Uses the arrays saved by
# 09_age_specific_susceptibility_diagnostic.R.

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(tidyr); library(ggplot2); library(patchwork) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
base_dir <- file.path(root, "02_Script/40_renewal_model/36_climate_forced_v4_9")
table_dir <- file.path(root, "03_Output/tables/climate_forced_v4_9")
figure_dir <- file.path(root, "03_Output/figures/climate_forced_v4_9")
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

arr <- readRDS(file.path(base_dir, "outputs/ce_canary/age_diagnostic_arrays.rds"))
dates <- arr$dates; ages <- arr$ages; age_group_labels <- arr$age_group_labels
checkpoint_dates <- arr$checkpoint_dates; checkpoint_labels <- arr$checkpoint_labels

CE_2017_ONSET <- as.Date("2017-02-12"); CE_2017_END <- as.Date("2018-12-16")
CE_2022_SEED_START <- as.Date("2022-01-02"); CE_2022_WAVE_END <- as.Date("2023-12-17")

# ================================================================
# Figure A: age x time heatmaps (posterior median)
# ================================================================
S_prop_median <- apply(arr$S_prop_arr, c(2, 3), median)
immune_prop_median <- apply(arr$immune_prop_arr, c(2, 3), median)
heat_S <- expand.grid(t = seq_along(dates), age = ages) |>
  mutate(week_start = dates[t], S_prop = S_prop_median[cbind(t, age + 1)])
heat_U <- expand.grid(t = seq_along(dates), age = ages) |>
  mutate(week_start = dates[t], immune_prop = immune_prop_median[cbind(t, age + 1)])

epidemic_lines <- tibble(week_start = c(CE_2017_ONSET, CE_2022_SEED_START), label = c("2017 epidemic onset", "2022 seed period start"))

p_heat_S <- ggplot(heat_S, aes(week_start, age, fill = S_prop)) +
  geom_raster() +
  geom_vline(data = epidemic_lines, aes(xintercept = week_start), colour = "white", linetype = "dashed", linewidth = 0.4) +
  scale_fill_viridis_c(name = "Susceptible\nfraction", labels = scales::label_percent(accuracy = 1), option = "viridis") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Figure A1: age x time susceptible-fraction heatmap (posterior median)",
       subtitle = "Dashed lines: 2017 epidemic onset, 2022 seed period start (actual model dates)",
       x = NULL, y = "Age (years)") + theme_v4
ggsave(file.path(figure_dir, "AGE_DIAG_FigA1_susceptibility_heatmap.png"), p_heat_S, width = 220, height = 140, units = "mm", dpi = 300, bg = "white")

p_heat_U <- ggplot(heat_U, aes(week_start, age, fill = immune_prop)) +
  geom_raster() +
  geom_vline(data = epidemic_lines, aes(xintercept = week_start), colour = "white", linetype = "dashed", linewidth = 0.4) +
  scale_fill_viridis_c(name = "Cumulative\nimmune fraction", labels = scales::label_percent(accuracy = 1), option = "magma") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Figure A2: age x time cumulative immune-fraction heatmap (posterior median)",
       subtitle = "Dashed lines: 2017 epidemic onset, 2022 seed period start (actual model dates)",
       x = NULL, y = "Age (years)") + theme_v4
ggsave(file.path(figure_dir, "AGE_DIAG_FigA2_immunity_heatmap.png"), p_heat_U, width = 220, height = 140, units = "mm", dpi = 300, bg = "white")
message("[saved] Figure A1, A2 (heatmaps)")

# ================================================================
# Figure B: age-group susceptibility trajectories (median + 95% CrI)
# ================================================================
S_prop_group_arr <- arr$S_group_arr / arr$N_group_arr # [draw, t, group]
group_traj <- bind_rows(lapply(seq_along(age_group_labels), function(g) {
  tibble(week_start = dates, age_group = age_group_labels[g],
         median = apply(S_prop_group_arr[, , g], 2, median),
         lo95 = apply(S_prop_group_arr[, , g], 2, quantile, .025),
         hi95 = apply(S_prop_group_arr[, , g], 2, quantile, .975))
})) |> mutate(age_group = factor(age_group, levels = age_group_labels))

p_group_traj <- ggplot(group_traj, aes(week_start, median, colour = age_group, fill = age_group)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.6) +
  scale_colour_viridis_d(name = "Age group") + scale_fill_viridis_d(name = "Age group") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Figure B: age-group susceptible-fraction trajectories", subtitle = "Population-weighted; median + 95% posterior interval",
       x = NULL, y = "Susceptible fraction (S_prop)") + theme_v4
ggsave(file.path(figure_dir, "AGE_DIAG_FigB_group_susceptibility_trajectories.png"), p_group_traj, width = 220, height = 130, units = "mm", dpi = 300, bg = "white")
message("[saved] Figure B")

# ================================================================
# Figure C: age-specific infection curves (weekly incidence/100k + cumulative)
# ================================================================
incidence_group <- bind_rows(lapply(seq_along(age_group_labels), function(g) {
  inc <- arr$Inew_group_arr[, , g] / arr$N_group_arr[, , g] * 1e5
  tibble(week_start = dates, age_group = age_group_labels[g],
         median = apply(inc, 2, median), lo95 = apply(inc, 2, quantile, .025), hi95 = apply(inc, 2, quantile, .975))
})) |> mutate(age_group = factor(age_group, levels = age_group_labels))

p_incidence <- ggplot(incidence_group, aes(week_start, median, colour = age_group)) +
  geom_line(linewidth = 0.5) +
  scale_colour_viridis_d(name = "Age group") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Figure C1: age-group weekly infection incidence", subtitle = "Per 100,000 population in that age group; median across posterior draws",
       x = NULL, y = "Weekly infections per 100,000") + theme_v4

# NOTE (bug fixed): cumulative infection proportion for an age-BAND must
# use the CURRENT immune stock Uinf_group(t)/N_group(t) -- i.e. the same
# quantity as immune_prop_group, consistent with the summary table and
# with this project's established "immune_prop = cumulative proportion
# infected" convention. An earlier version of this figure used
# cumsum(Inew_group(t))/N_group(t) instead, which is NOT a valid cumulative
# proportion for an age band: it only counts infections that occurred
# WHILE a person was already in that band, so anyone ageing in after being
# infected at a younger age is missing from the numerator while still
# counted in the (growing) denominator -- producing an artificial decline
# between epidemics that has nothing to do with waning immunity (there is
# none in this lifelong-immunity model; Uinf/N is confirmed non-decreasing
# for every age band, to floating-point precision).
cumulative_group <- bind_rows(lapply(seq_along(age_group_labels), function(g) {
  immune_prop_g <- arr$Uinf_group_arr[, , g] / arr$N_group_arr[, , g]
  tibble(week_start = dates, age_group = age_group_labels[g],
         median = apply(immune_prop_g, 2, median), lo95 = apply(immune_prop_g, 2, quantile, .025), hi95 = apply(immune_prop_g, 2, quantile, .975))
})) |> mutate(age_group = factor(age_group, levels = age_group_labels))

p_cumulative <- ggplot(cumulative_group, aes(week_start, median, colour = age_group, fill = age_group)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.1, colour = NA) +
  geom_line(linewidth = 0.6) +
  scale_colour_viridis_d(name = "Age group") + scale_fill_viridis_d(name = "Age group") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Figure C2: age-group cumulative infection proportion", subtitle = "= Uinf_group(t)/N_group(t): share of the CURRENTLY LIVING group ever infected (not an individual lifetime probability).\nCan decline in fast-turnover young bands via birth dilution -- see AGE_SPECIFIC_SUSCEPTIBILITY_DIAGNOSTIC.md",
       x = NULL, y = "Cumulative proportion infected") + theme_v4

figure_C <- p_incidence / p_cumulative
ggsave(file.path(figure_dir, "AGE_DIAG_FigC_infection_curves.png"), figure_C, width = 220, height = 200, units = "mm", dpi = 300, bg = "white")
message("[saved] Figure C")

# ================================================================
# Figure D: cross-sections at key dates (single-year age resolution)
# ================================================================
cross_section_S <- bind_rows(lapply(seq_along(checkpoint_dates), function(k) {
  tibble(checkpoint = checkpoint_labels[k], checkpoint_date = checkpoint_dates[k], age = ages,
         median = apply(arr$S_prop_checkpoint_arr[, k, ], 2, median),
         lo95 = apply(arr$S_prop_checkpoint_arr[, k, ], 2, quantile, .025),
         hi95 = apply(arr$S_prop_checkpoint_arr[, k, ], 2, quantile, .975))
})) |> mutate(checkpoint = factor(checkpoint, levels = checkpoint_labels))

p_cross_S <- ggplot(cross_section_S, aes(age, median, colour = checkpoint)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = checkpoint), alpha = 0.08, colour = NA) +
  geom_line(linewidth = 0.6) +
  scale_colour_viridis_d(name = NULL) + scale_fill_viridis_d(name = NULL) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Figure D: age profile of susceptible fraction at key dates",
       subtitle = paste(sprintf("%s (%s)", checkpoint_labels, format(checkpoint_dates, "%Y-%m-%d")), collapse = " | "),
       x = "Age (years)", y = "Susceptible fraction (S_prop_age)") + theme_v4 +
  theme(plot.subtitle = element_text(size = 6.2))
ggsave(file.path(figure_dir, "AGE_DIAG_FigD_cross_sections.png"), p_cross_S, width = 200, height = 130, units = "mm", dpi = 300, bg = "white")
message("[saved] Figure D")

# ================================================================
# Table: age-group summary at checkpoints (population-weighted, using
# S_group_arr/N_group_arr at the checkpoint weeks)
# ================================================================
checkpoint_idx <- arr$checkpoint_idx
group_checkpoint_tbl <- bind_rows(lapply(seq_along(checkpoint_dates), function(k) {
  t_idx <- checkpoint_idx[k]
  bind_rows(lapply(seq_along(age_group_labels), function(g) {
    S_prop_draws <- arr$S_group_arr[, t_idx, g] / arr$N_group_arr[, t_idx, g]
    immune_prop_draws <- arr$Uinf_group_arr[, t_idx, g] / arr$N_group_arr[, t_idx, g]
    tibble(checkpoint = checkpoint_labels[k], checkpoint_date = checkpoint_dates[k], age_group = age_group_labels[g],
           population = round(median(arr$N_group_arr[, t_idx, g])),
           S_prop_median = median(S_prop_draws), S_prop_lo95 = quantile(S_prop_draws, .025), S_prop_hi95 = quantile(S_prop_draws, .975),
           immune_prop_median = median(immune_prop_draws), cumulative_infection_prop_median = median(immune_prop_draws))
  }))
})) |> mutate(age_group = factor(age_group, levels = age_group_labels))

range_by_checkpoint <- group_checkpoint_tbl |> group_by(checkpoint, checkpoint_date) |>
  summarise(S_prop_range_min = min(S_prop_median), S_prop_range_max = max(S_prop_median),
            S_prop_range_span = max(S_prop_median) - min(S_prop_median), .groups = "drop")

message("\n=== Table: age-group summary at checkpoints ===")
print(as.data.frame(group_checkpoint_tbl), digits = 3)
message("\n=== Range of S_prop across age groups at each checkpoint ===")
print(as.data.frame(range_by_checkpoint), digits = 3)

write_csv(group_checkpoint_tbl, file.path(table_dir, "AGE_DIAGNOSTIC_group_summary_table.csv"))
write_csv(range_by_checkpoint, file.path(table_dir, "AGE_DIAGNOSTIC_group_range_by_checkpoint.csv"))
message("\n[saved] AGE_DIAGNOSTIC_group_summary_table.csv, AGE_DIAGNOSTIC_group_range_by_checkpoint.csv")
