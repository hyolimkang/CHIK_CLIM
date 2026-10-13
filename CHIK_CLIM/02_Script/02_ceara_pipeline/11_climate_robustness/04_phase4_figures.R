# PHASE 4 -- Publication-quality figures FIG13-16.
#
# Reads the tables written by 03_phase3_scientific_comparison.R (run that
# first). Okabe-Ito colourblind-safe palette, matching this repo's other
# figures (dataviz skill consulted for prior figures in this project;
# convention reused here for consistency).

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr", "ggplot2", "patchwork", "scales", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(tidyr); library(ggplot2); library(patchwork) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/11_climate_robustness/00_config_climate_robustness.R"))

OKABE <- c(blue = "#0072B2", orange = "#D55E00", yellow = "#E69F00", green = "#009E73", purple = "#CC79A7", grey = "#56B4E9")
MODEL_COLOURS <- c(M0 = OKABE[["blue"]], M1 = OKABE[["orange"]], M2 = OKABE[["green"]], M3 = OKABE[["purple"]], M4 = OKABE[["yellow"]])
MODEL_LABELS <- c(M0 = "M0: primary, SD=0.40\n(accepted)", M1 = "M1: sens_A, SD=0.40", M2 = "M2: sens_B, SD=0.40", M3 = "M3: primary, SD=0.20", M4 = "M4: primary, SD=0.60")
theme_pub <- theme_classic(base_size = 9) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey30"), legend.position = "top",
        plot.caption = element_text(size = 7, colour = "grey40", face = "italic"))
HIST_NOTE <- "Historical inference only (2015-2025); not a future climate-change projection."

save_fig <- function(plot, name, width_mm = 180, height_mm = 120) {
  ggsave(file.path(DIR_FIGURES_CR, paste0(name, ".pdf")), plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(file.path(DIR_FIGURES_CR, paste0(name, ".png")), plot, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
  ggsave(file.path(DIR_FIGURES_CR, paste0(name, ".svg")), plot, width = width_mm, height = height_mm, units = "mm")
  message("[saved] ", name, " (.pdf/.png@600dpi/.svg)")
}

message("[Phase4] Loading Phase 3 tables...")
beta_summary <- read_csv(file.path(DIR_TABLES_CR, "PHASE3_beta_posterior_summary_M0_M4.csv"), show_col_types = FALSE) |>
  mutate(model_label = factor(recode(model, !!!MODEL_LABELS), levels = rev(MODEL_LABELS)))

# ================================================================
# FIG13 -- Climate coefficient robustness across M0-M4
# ================================================================
fig13 <- ggplot(beta_summary, aes(model_label, posterior_median, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_pointrange(aes(ymin = lo95, ymax = hi95), fatten = 3, linewidth = 0.7) +
  coord_flip() +
  facet_wrap(~ parameter, labeller = as_labeller(c(beta_T = "beta[T]~(temperature)", beta_P = "beta[P]~(precipitation)"), label_parsed)) +
  scale_colour_manual(values = MODEL_COLOURS, guide = "none") +
  labs(title = "FIG13: Climate coefficient robustness across 5 historical specifications",
       subtitle = "Posterior median and 95% credible interval. Dashed line: zero effect.\nClimate-coefficient prior N(0,0.15) unchanged across all 5 specifications.",
       x = NULL, y = "Coefficient value (own standardised z-scale; see PHASE3_beta_physical_scale_comparison.csv for physical units)",
       caption = HIST_NOTE) + theme_pub + theme(axis.title.x = element_text(size = 7.5))
save_fig(fig13, "FIG13_climate_coefficient_robustness", width_mm = 185, height_mm = 115)

# ================================================================
# FIG14 -- Climate forcing multiplier under M0-M4 (small multiples --
# different lag transformations make a single overlaid panel misleading)
# ================================================================
multiplier_weekly <- read_csv(file.path(DIR_TABLES_CR, "PHASE3_climate_multiplier_weekly_M0_M4.csv"), show_col_types = FALSE) |>
  mutate(model_label = factor(recode(model, !!!MODEL_LABELS), levels = MODEL_LABELS))

fig14 <- ggplot(multiplier_weekly, aes(week_start, colour = model, fill = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.25, colour = NA) +
  geom_line(aes(y = median), linewidth = 0.4) +
  facet_wrap(~ model_label, ncol = 1) +
  scale_colour_manual(values = MODEL_COLOURS, guide = "none") + scale_fill_manual(values = MODEL_COLOURS, guide = "none") +
  labs(title = "FIG14: Historical weekly climate multiplier exp(beta_T*z_T + beta_P*z_P) under M0-M4",
       subtitle = "Each panel uses that model's OWN lag transformation and standardisation -- multiplier values are\ncomparable as a forcing on R0(t), NOT as directly comparable raw covariate scales (see FIG13 caption).",
       x = NULL, y = "Multiplier on R0(t)", caption = HIST_NOTE) + theme_pub + theme(strip.text = element_text(size = 7.5))
save_fig(fig14, "FIG14_climate_forcing_multiplier_by_model", width_mm = 180, height_mm = 210)

# ================================================================
# FIG15 -- Climate-year-effect confounding, especially 2022
# ================================================================
year_cor <- read_csv(file.path(DIR_TABLES_CR, "PHASE3_beta_vs_year_effect_correlations_M0_M4.csv"), show_col_types = FALSE) |>
  pivot_longer(c(cor_beta_T, cor_beta_P), names_to = "parameter", values_to = "correlation") |>
  mutate(parameter = recode(parameter, cor_beta_T = "beta[T]", cor_beta_P = "beta[P]"),
         model = factor(model, levels = c("M0", "M1", "M2", "M3", "M4")))

p_heatmap <- ggplot(year_cor, aes(factor(year), model, fill = correlation)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(data = year_cor |> dplyr::filter(year == 2022), aes(label = sprintf("%.2f", correlation)), size = 2.6, colour = "black") +
  facet_wrap(~ parameter, labeller = label_parsed) +
  scale_fill_gradient2(name = "Posterior\ncorrelation", low = OKABE[["blue"]], mid = "white", high = OKABE[["orange"]], midpoint = 0, limits = c(-1, 1)) +
  labs(subtitle = "A. cor(climate coefficient, year_effect[y]) by calendar year and model (2022 value annotated)", x = "Calendar year", y = NULL) +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

cor_2022 <- read_csv(file.path(DIR_TABLES_CR, "PHASE3_beta_vs_2022_year_effect_correlation.csv"), show_col_types = FALSE) |>
  pivot_longer(c(cor_beta_T, cor_beta_P), names_to = "parameter", values_to = "correlation") |>
  mutate(parameter = recode(parameter, cor_beta_T = "beta[T]", cor_beta_P = "beta[P]"), model = factor(model, levels = c("M0", "M1", "M2", "M3", "M4")))
p_2022 <- ggplot(cor_2022, aes(model, correlation, fill = model)) +
  geom_col(width = 0.6) + geom_hline(yintercept = 0, colour = "grey40") +
  facet_wrap(~ parameter, labeller = label_parsed) +
  scale_fill_manual(values = MODEL_COLOURS, guide = "none") + scale_y_continuous(limits = c(-1, 1)) +
  labs(subtitle = "B. cor(climate coefficient, year_effect[2022]) by model", x = NULL, y = "Posterior correlation") + theme_pub

fig15 <- p_heatmap / p_2022 + plot_layout(heights = c(1.4, 1))
fig15 <- fig15 + plot_annotation(title = "FIG15: Climate-year-effect confounding across model specifications",
                                  caption = HIST_NOTE, theme = theme(plot.caption = element_text(size = 7, colour = "grey40", face = "italic")))
save_fig(fig15, "FIG15_climate_year_effect_confounding", width_mm = 180, height_mm = 190)

# ================================================================
# FIG16 -- Historical reconstruction sensitivity: R0(t), S_prop(t), PPC
# ================================================================
recon_weekly <- read_csv(file.path(DIR_TABLES_CR, "PHASE3_R0_Sprop_weekly_M0_M4.csv"), show_col_types = FALSE) |>
  mutate(model = factor(model, levels = c("M0", "M1", "M2", "M3", "M4")))
annual_ppc <- read_csv(file.path(DIR_TABLES_CR, "PHASE3_annual_PPC_M0_M4.csv"), show_col_types = FALSE) |>
  mutate(model = factor(model, levels = c("M0", "M1", "M2", "M3", "M4")))

p_r0 <- ggplot(recon_weekly, aes(week_start, R0_median, colour = model)) + geom_line(linewidth = 0.3, alpha = 0.85) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_colour_manual(name = NULL, values = MODEL_COLOURS, labels = MODEL_LABELS) +
  labs(subtitle = "A. Posterior median R0(t)", x = NULL, y = "R0(t)") + theme_pub + theme(legend.text = element_text(size = 6.5))

p_s <- ggplot(recon_weekly, aes(week_start, S_prop_median, colour = model)) + geom_line(linewidth = 0.4, alpha = 0.85) +
  scale_colour_manual(values = MODEL_COLOURS, guide = "none") + scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(subtitle = "B. Susceptible proportion S(t)/N(t)", x = NULL, y = "S / N") + theme_pub

p_ppc <- ggplot(annual_ppc, aes(factor(year), ratio, fill = model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(name = NULL, values = MODEL_COLOURS, labels = MODEL_LABELS) +
  labs(subtitle = "C. Annual posterior-predictive median / observed reported cases", x = NULL, y = "Predicted / observed ratio") +
  theme_pub + theme(legend.text = element_text(size = 6.5), axis.text.x = element_text(angle = 45, hjust = 1))

fig16 <- p_r0 / p_s / p_ppc
fig16 <- fig16 + plot_annotation(title = "FIG16: Historical reconstruction sensitivity across M0-M4", caption = HIST_NOTE,
                                  theme = theme(plot.caption = element_text(size = 7, colour = "grey40", face = "italic")))
save_fig(fig16, "FIG16_historical_reconstruction_sensitivity", width_mm = 180, height_mm = 220)

message("\n[Phase4] FIG13-16 complete -> ", DIR_FIGURES_CR)
