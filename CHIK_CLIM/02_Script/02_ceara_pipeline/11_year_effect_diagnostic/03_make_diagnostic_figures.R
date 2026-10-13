# Year-effect diagnostic -- FIG10-12.
#
# Reads the tables/results from 01_historical_vs_future_variability.R and
# 02_paired_arm_simulation.R (run both first). Colour convention matches
# the rest of this repo's figures: the Okabe-Ito colourblind-safe
# qualitative palette.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr", "ggplot2", "patchwork", "scales", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr); library(ggplot2); library(patchwork) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/11_year_effect_diagnostic/00_config_year_effect_diagnostic.R"))

OKABE <- c(blue = "#0072B2", orange = "#D55E00", yellow = "#E69F00", green = "#009E73", purple = "#CC79A7", grey = "#56B4E9")
theme_pub <- theme_classic(base_size = 9) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey30"), legend.position = "top",
        plot.caption = element_text(size = 7, colour = "grey40", face = "italic"))
DIAG_NOTE <- "Diagnostic simulation (paired design, common random numbers across arms); not a scenario projection."

save_fig <- function(plot, name, width_mm = 180, height_mm = 120) {
  ggsave(file.path(DIR_FIGURES_DIAG, paste0(name, ".pdf")), plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(file.path(DIR_FIGURES_DIAG, paste0(name, ".png")), plot, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
  ggsave(file.path(DIR_FIGURES_DIAG, paste0(name, ".svg")), plot, width = width_mm, height = height_mm, units = "mm")
  message("[saved] ", name, " (.pdf/.png@600dpi/.svg)")
}

message("[Diag 3] Loading accepted fit + diagnostic results...")
climate <- readRDS(ACCEPTED_CLIMATE_FIT)
draws <- rstan::extract(climate$fit, pars = c("year_effect", "A_year", "sigma_season_year"), permuted = TRUE)
YEAR_EFFECT_PRIOR_SD <- 0.40

# ================================================================
# FIG10 -- Historical posterior vs. future Y1 pooled distributions
# ================================================================
set.seed(DIAG_BASE_SEED + 2L)
n_draws_hist <- nrow(draws$year_effect)
N_REAL_FIG <- 2000L
draw_for_fig <- sample.int(n_draws_hist, N_REAL_FIG, replace = TRUE)
future_individual <- bind_rows(lapply(seq_len(N_REAL_FIG), function(i) {
  d <- draw_for_fig[i]
  z_year_f <- rnorm(10); z_season_f <- rnorm(10)
  tibble(year_effect = YEAR_EFFECT_PRIOR_SD * (z_year_f - mean(z_year_f)), A_year = exp(draws$sigma_season_year[d] * z_season_f))
}))

hist_individual <- tibble(year_effect = as.vector(draws$year_effect), A_year = as.vector(draws$A_year))

fig10_df <- bind_rows(
  hist_individual |> mutate(source = "Historical (2015-2025 posterior)"),
  future_individual |> mutate(source = "Future Y1 (2026-2035 resampled)")
) |> mutate(source = factor(source, levels = c("Historical (2015-2025 posterior)", "Future Y1 (2026-2035 resampled)")))

p_year_effect <- ggplot(fig10_df, aes(year_effect, fill = source, colour = source)) +
  geom_density(alpha = 0.45, linewidth = 0.5) + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(name = NULL, values = c("Historical (2015-2025 posterior)" = OKABE[["blue"]], "Future Y1 (2026-2035 resampled)" = OKABE[["orange"]])) +
  scale_colour_manual(name = NULL, values = c("Historical (2015-2025 posterior)" = OKABE[["blue"]], "Future Y1 (2026-2035 resampled)" = OKABE[["orange"]])) +
  labs(subtitle = "A. year_effect (pooled across years x draws/realisations)", x = "year_effect (log R0 scale)", y = "Density") + theme_pub
p_A_year <- ggplot(fig10_df, aes(A_year, fill = source, colour = source)) +
  geom_density(alpha = 0.45, linewidth = 0.5) + geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(name = NULL, values = c("Historical (2015-2025 posterior)" = OKABE[["blue"]], "Future Y1 (2026-2035 resampled)" = OKABE[["orange"]])) +
  scale_colour_manual(name = NULL, values = c("Historical (2015-2025 posterior)" = OKABE[["blue"]], "Future Y1 (2026-2035 resampled)" = OKABE[["orange"]])) +
  scale_x_continuous(limits = c(0, 6)) +
  labs(subtitle = "B. A_year (seasonal-amplitude multiplier)", x = "A_year", y = "Density") + theme_pub

fig10 <- p_year_effect / p_A_year
fig10 <- fig10 + plot_annotation(
  title = "FIG10: Historical posterior vs. future Y1 pooled distributions",
  subtitle = "Historical pools across-year AND across-draw variation (posterior uncertainty + interannual\nvariability mixed). Future pools 2,000 independent 10-year realisations. See\nDIAG_future_vs_historical_spread_ratios.csv for the cleaner within-draw-only comparison.",
  caption = DIAG_NOTE, theme = theme(plot.subtitle = element_text(size = 8, colour = "grey30"), plot.caption = element_text(size = 7, colour = "grey40", face = "italic")))
save_fig(fig10, "FIG10_historical_vs_future_distributions", width_mm = 170, height_mm = 150)

# ================================================================
# FIG11 -- Annual R0 peak distributions across 4 arms + historical reference
# ================================================================
annual_R0_future <- read_csv(file.path(DIR_RESULTS_DIAG, "DIAG_full_annual_R0_by_arm.csv"), show_col_types = FALSE)
hist_annual_sample <- read_csv(file.path(DIR_RESULTS_DIAG, "DIAG_historical_annual_R0_sample.csv"), show_col_types = FALSE)

SHORT_ARM_LABELS <- c(A_neutral = "A: neutral", B_year_only = "B: year-level\nonly", C_amplitude_only = "C: amplitude\nonly", D_both = "D: both\n(= Y1)")
fig11_df <- bind_rows(
  hist_annual_sample |> transmute(group = "Historical\n(2015-2025)", peak_R0 = peak_R0),
  annual_R0_future |> transmute(group = recode(arm, !!!SHORT_ARM_LABELS), peak_R0 = peak_R0)
) |> mutate(group = factor(group, levels = c("Historical\n(2015-2025)", SHORT_ARM_LABELS[DECOMP_ARMS])))

ARM_COLOURS <- c(OKABE[["grey"]], OKABE[["blue"]], OKABE[["green"]], OKABE[["purple"]], OKABE[["orange"]])
names(ARM_COLOURS) <- levels(fig11_df$group)

fig11 <- ggplot(fig11_df, aes(group, peak_R0, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.7, linewidth = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = ARM_COLOURS, guide = "none") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(title = "FIG11: Annual peak R0(t) -- four structural arms vs. historical reference",
       subtitle = "A: year_effect=0,A_year=1 | B: year_effect=Y1,A_year=1 | C: year_effect=0,A_year=Y1 | D: both (=Y1).\nHistorical: posterior x year pooled (sampled). Arms: 100 draws x 10 future years (C0 climate-neutral, lambda_import=1).",
       x = NULL, y = "Annual peak R0(t)", caption = DIAG_NOTE) + theme_pub +
  theme(axis.text.x = element_text(size = 7.5))
save_fig(fig11, "FIG11_annual_peak_R0_by_arm", width_mm = 180, height_mm = 120)

# ================================================================
# FIG12 -- Paired ten-year cumulative infection outcomes
# ================================================================
per_traj <- read_csv(file.path(DIR_RESULTS_DIAG, "DIAG_full_per_trajectory_summary.csv"), show_col_types = FALSE) |>
  mutate(arm_label = factor(recode(arm, !!!SHORT_ARM_LABELS), levels = SHORT_ARM_LABELS[DECOMP_ARMS]))
paired_summary <- read_csv(file.path(DIR_TABLES_DIAG, "DIAG_paired_contrast_summary.csv"), show_col_types = FALSE) |>
  mutate(contrast = factor(contrast, levels = rev(contrast)))

ARM_COLOURS2 <- c(OKABE[["blue"]], OKABE[["green"]], OKABE[["purple"]], OKABE[["orange"]]); names(ARM_COLOURS2) <- SHORT_ARM_LABELS[DECOMP_ARMS]

p_dist <- ggplot(per_traj, aes(arm_label, cum_10y + 1, fill = arm_label)) +
  geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.7, linewidth = 0.3) +
  scale_y_log10(labels = scales::label_number(big.mark = ",")) +
  scale_fill_manual(values = ARM_COLOURS2, guide = "none") +
  labs(subtitle = "A. 10-year cumulative latent infections by arm (log scale; n=2,000 draw x replicate pairs each).\nA: neutral | B: year-level only | C: amplitude only | D: both (=Y1)",
       x = NULL, y = "Cumulative infections (log scale)") +
  theme_pub + theme(axis.text.x = element_text(size = 8))

p_paired <- ggplot(paired_summary, aes(contrast, median_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_pointrange(aes(ymin = lo95, ymax = hi95), colour = OKABE[["orange"]], linewidth = 0.6, fatten = 3) +
  coord_flip() +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(subtitle = "B. Paired differences in 10-year cumulative infections (same posterior draw x replicate across arms)",
       x = NULL, y = "Difference in cumulative infections (95% interval across paired draws)", caption = DIAG_NOTE) + theme_pub

fig12 <- p_dist / p_paired + plot_layout(heights = c(1.3, 1))
fig12 <- fig12 + plot_annotation(title = "FIG12: Paired 10-year cumulative infection outcomes across four arms")
save_fig(fig12, "FIG12_paired_cumulative_infections", width_mm = 180, height_mm = 170)

message("\n[Diag 3] FIG10-12 complete -> ", DIR_FIGURES_DIAG)
