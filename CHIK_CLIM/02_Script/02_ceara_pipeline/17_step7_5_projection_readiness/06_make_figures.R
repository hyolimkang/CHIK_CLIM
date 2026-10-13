# STEP 7.5 -- Figures 1-3. Reads only already-produced Phase 0-3
# tables/results (no new simulation).

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(tidyr); library(ggplot2); library(patchwork) })
filter <- dplyr::filter

source(file.path(ROOT, "02_Script/02_ceara_pipeline/17_step7_5_projection_readiness/00_config_step7_5.R"))

OKABE <- c(blue = "#0072B2", orange = "#D55E00", yellow = "#E69F00", green = "#009E73", purple = "#CC79A7", grey = "#56B4E9")
theme_pub <- theme_classic(base_size = 9) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey30"), legend.position = "top",
        plot.caption = element_text(size = 7, colour = "grey40", face = "italic"))
save_fig <- function(plot, name, width_mm = 190, height_mm = 150) {
  ggsave(file.path(DIR_FIGURES_STEP75, paste0(name, ".pdf")), plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(file.path(DIR_FIGURES_STEP75, paste0(name, ".png")), plot, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
  message("[saved] ", name)
}
PROV_CAPTION <- paste0(PROVISIONAL_NOTE, "\nDiagnostic structural arms (p=0.10/0.25 prespecified, NOT fit); no lambda or process is chosen as a new reference.")

# ================================================================
# FIG1 -- H1 2020-2021 infection persistence and observed cases
# ================================================================
weekly <- read_csv(file.path(DIR_RESULTS_STEP75, "PHASE1_weekly_2020_2021_predictive_vs_observed.csv"), show_col_types = FALSE)

p1a <- ggplot(weekly, aes(week_start)) +
  geom_ribbon(aes(ymin = pmax(X_latent_lo95, 0.1), ymax = X_latent_hi95), fill = OKABE[["blue"]], alpha = 0.15) +
  geom_ribbon(aes(ymin = pmax(X_latent_lo50, 0.1), ymax = X_latent_hi50), fill = OKABE[["blue"]], alpha = 0.3) +
  geom_line(aes(y = pmax(X_latent_median, 0.1)), colour = OKABE[["blue"]], linewidth = 0.4) +
  scale_y_log10(labels = scales::comma) +
  labs(x = NULL, y = "Simulated latent\ninfections X(t) (log)",
       subtitle = "A: simulated weekly latent infections (H1/G2, lambda=1) -- median, 50%/95% intervals.") + theme_pub

p1b <- ggplot(weekly, aes(week_start)) +
  geom_ribbon(aes(ymin = pmax(Y_reported_lo95, 0.1), ymax = Y_reported_hi95), fill = OKABE[["orange"]], alpha = 0.15) +
  geom_ribbon(aes(ymin = pmax(Y_reported_lo50, 0.1), ymax = Y_reported_hi50), fill = OKABE[["orange"]], alpha = 0.3) +
  geom_line(aes(y = pmax(Y_reported_median, 0.1)), colour = OKABE[["orange"]], linewidth = 0.4) +
  geom_point(aes(y = pmax(observed_cases, 0.1)), colour = "black", size = 0.8) +
  geom_line(aes(y = pmax(observed_cases, 0.1)), colour = "black", linewidth = 0.3, linetype = "dotted") +
  scale_y_log10(labels = scales::comma) +
  labs(x = NULL, y = "Reported cases (log)\nblack = observed",
       subtitle = "B: simulated reported cases vs. OBSERVED weekly surveillance (used only for evaluation, never to calibrate).") + theme_pub

p1c <- ggplot(weekly, aes(week_start)) +
  geom_line(aes(y = R_eff_median), colour = OKABE[["purple"]], linewidth = 0.4) +
  geom_ribbon(aes(ymin = R_eff_lo50, ymax = R_eff_hi50), fill = OKABE[["purple"]], alpha = 0.25) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  labs(x = NULL, y = "Reff(t)", subtitle = "C: Reff(t) -- crosses 1 repeatedly; the bottleneck is not insufficient transmission potential.") + theme_pub

p1d <- ggplot(weekly, aes(week_start)) +
  geom_col(aes(y = M_import_median), fill = OKABE[["grey"]], width = 6) +
  labs(x = "Week", y = "Median weekly\nintroduction pressure M_t",
       subtitle = "D: median introduction pressure M_t (lambda=1 throughout) -- small and constant.", caption = PROV_CAPTION) + theme_pub

fig1 <- p1a / p1b / p1c / p1d
save_fig(fig1, "FIG_S75_1_H1_2020_2021_persistence", height_mm = 240)

# ================================================================
# FIG2 -- constant vs. episodic introduction, identical expected pressure
# ================================================================
arms <- readRDS(file.path(DIR_RESULTS_STEP75, "PHASE2_structural_arms.rds"))
ARM_LABS <- c(ModelA_constant = "Model A: constant\nPoisson(lambda=1)", ModelB_p10 = "Model B: episodic\np=0.10 (mean=1)", ModelB_p25 = "Model B: episodic\np=0.25 (mean=1)")
ARM_COLOURS <- c("Model A: constant\nPoisson(lambda=1)" = unname(OKABE["grey"]), "Model B: episodic\np=0.10 (mean=1)" = unname(OKABE["blue"]), "Model B: episodic\np=0.25 (mean=1)" = unname(OKABE["orange"]))

entry_by_arm <- bind_rows(lapply(names(arms), function(nm) {
  arms[[nm]]$trajectories |> filter(week_start < as.Date("2022-01-02")) |> arrange(posterior_draw_id, stochastic_replicate_id, week_start) |>
    group_by(posterior_draw_id, stochastic_replicate_id) |> summarise(entry_weighted_infectiousness = sum(W_GEN * rev(tail(X_latent, G_GEN))), .groups = "drop") |>
    mutate(arm_lab = ARM_LABS[[nm]])
}))
p2a <- ggplot(entry_by_arm, aes(arm_lab, pmax(entry_weighted_infectiousness, 0.01), fill = arm_lab)) +
  geom_violin(trim = TRUE, alpha = 0.75, linewidth = 0.2) + geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  scale_y_log10() + scale_fill_manual(values = ARM_COLOURS, guide = "none") +
  labs(x = NULL, y = "Entry-2022 weighted\ninfectiousness (log)",
       subtitle = "A: entry-2022 infectiousness is nearly identical across structures at equal mean pressure.") + theme_pub

mtraj <- bind_rows(lapply(names(arms), function(nm) arms[[nm]]$trajectories |> filter(posterior_draw_id == arms[[nm]]$trajectories$posterior_draw_id[1], stochastic_replicate_id == 1,
                                                                                         week_start >= as.Date("2021-06-01"), week_start <= as.Date("2022-02-01")) |> mutate(arm_lab = ARM_LABS[[nm]])))
p2b <- ggplot(mtraj, aes(week_start, M_import, colour = arm_lab)) + geom_point(size = 1) + geom_line(linewidth = 0.3, alpha = 0.5) +
  facet_wrap(~arm_lab, ncol = 1, scales = "free_y") + scale_colour_manual(values = ARM_COLOURS, guide = "none") +
  labs(x = "Week", y = "Realised M_t\n(one example draw/replicate)",
       subtitle = "B: realised weekly introduction pressure -- Model A steady and small; Model B zero most weeks, occasional larger bursts.",
       caption = PROV_CAPTION) + theme_pub + theme(strip.text = element_text(size = 7))

fig2 <- p2a / p2b + plot_layout(heights = c(1, 1.3)) +
  plot_annotation(title = "FIG_S75_2. Constant vs. episodic introduction, identical expected pressure (H1/G2)", theme = theme_pub)
save_fig(fig2, "FIG_S75_2_constant_vs_episodic", height_mm = 220)

# ================================================================
# FIG3 -- local amplification, entry-2022 infectiousness, epidemic outcomes
# ================================================================
structural_summary <- read_csv(file.path(DIR_TABLES_STEP75, "PHASE3_structural_comparison_summary.csv"), show_col_types = FALSE) |>
  mutate(arm_lab = factor(ARM_LABS[arm], levels = unname(ARM_LABS)))

p3a <- ggplot(structural_summary, aes(arm_lab, prob_major_epidemic_2022, fill = arm_lab)) +
  geom_col(width = 0.6) + scale_y_continuous(labels = scales::percent) + scale_fill_manual(values = ARM_COLOURS, guide = "none") +
  labs(x = NULL, y = "P(major epidemic, 2022)", subtitle = "A: P(major epidemic)") + theme_pub + theme(axis.text.x = element_text(size = 6.5))

p3b <- ggplot(structural_summary, aes(arm_lab, prob_reaches_observed_2022, fill = arm_lab)) +
  geom_col(width = 0.6) + scale_y_continuous(labels = scales::percent) + scale_fill_manual(values = ARM_COLOURS, guide = "none") +
  labs(x = NULL, y = "P(reaches observed\n2022 burden)", subtitle = "B: P(>=50,986 reported cases)") + theme_pub + theme(axis.text.x = element_text(size = 6.5))

p3c <- ggplot(structural_summary, aes(arm_lab, median_cumulative_reported_full_period, fill = arm_lab)) +
  geom_col(width = 0.6) + geom_errorbar(aes(ymin = lo95_cumulative, ymax = hi95_cumulative), width = 0.15) +
  scale_y_log10(labels = scales::comma) + scale_fill_manual(values = ARM_COLOURS, guide = "none") +
  labs(x = NULL, y = "Cumulative reported\ncases, full period (log)", subtitle = "C: full-period cumulative reported cases, 95% interval") + theme_pub + theme(axis.text.x = element_text(size = 6.5))

p3d <- ggplot(structural_summary, aes(arm_lab, WIS_cumulative_full_period, fill = arm_lab)) +
  geom_col(width = 0.6) + scale_fill_manual(values = ARM_COLOURS, guide = "none") +
  labs(x = NULL, y = "Cumulative WIS\n(lower = better)", subtitle = "D: full-period WIS (near-identical across structures)", caption = PROV_CAPTION) +
  theme_pub + theme(axis.text.x = element_text(size = 6.5))

fig3 <- (p3a | p3b) / (p3c | p3d) +
  plot_annotation(title = "FIG_S75_3. Local amplification, entry-2022 infectiousness, and epidemic outcomes",
                   subtitle = "At equal mean introduction pressure, structure (constant vs. episodic burst) changes outcomes only marginally --\nmagnitude, not clustering, is the binding constraint.",
                   theme = theme_pub)
save_fig(fig3, "FIG_S75_3_amplification_and_outcomes", height_mm = 200)

message("\n[Step7.5] Figures 1-3 complete -> ", DIR_FIGURES_STEP75)
