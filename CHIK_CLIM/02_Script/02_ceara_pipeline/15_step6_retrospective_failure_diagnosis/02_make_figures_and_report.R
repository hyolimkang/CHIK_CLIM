# STEP 6 -- Figures FIG21-23 and the final decision report. Reads only
# the tables/results already produced by 01_build_diagnostics.R (no new
# simulation, no new fitting). PROVISIONAL throughout per the HMC audit.

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(tidyr); library(ggplot2); library(patchwork) })
filter <- dplyr::filter

source(file.path(ROOT, "02_Script/02_ceara_pipeline/15_step6_retrospective_failure_diagnosis/00_config_step6.R"))

OKABE <- c(blue = "#0072B2", orange = "#D55E00", yellow = "#E69F00", green = "#009E73", purple = "#CC79A7", grey = "#56B4E9")
theme_pub <- theme_classic(base_size = 9) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey30"), legend.position = "top",
        plot.caption = element_text(size = 7, colour = "grey40", face = "italic"))
save_fig <- function(plot, name, width_mm = 180, height_mm = 150) {
  ggsave(file.path(DIR_FIGURES_STEP6, paste0(name, ".pdf")), plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(file.path(DIR_FIGURES_STEP6, paste0(name, ".png")), plot, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
  message("[saved] ", name, " (.pdf/.png@600dpi)")
}
PROV_CAPTION <- paste0(PROVISIONAL_NOTE, "\nTraining-derived quantities (posterior draws) are distinguished from stochastic holdout simulations\n(replicate-level trajectories) and from observed validation data (never used to fit/calibrate the holdout).")

weekly_q <- read_csv(file.path(DIR_RESULTS_STEP6, "STEP6_weekly_trajectory_quantiles.csv"), show_col_types = FALSE)
epidemic_period <- read_csv(file.path(DIR_TABLES_STEP6, "STEP6_observed_2022_epidemic_period.csv"), show_col_types = FALSE)
entry_draws <- read_csv(file.path(DIR_RESULTS_STEP6, "STEP6_entry2022_infection_history_draws.csv"), show_col_types = FALSE)
mech_class <- read_csv(file.path(DIR_TABLES_STEP6, "STEP6_2022_mechanism_classification.csv"), show_col_types = FALSE)
establishment <- read_csv(file.path(DIR_TABLES_STEP6, "STEP6_2022_establishment_summary.csv"), show_col_types = FALSE)
outcome_contrasts <- read_csv(file.path(DIR_TABLES_STEP6, "STEP6_G2_2022_outcome_contrasts.csv"), show_col_types = FALSE)
annual_traj <- readRDS(file.path(DIR_RESULTS_STEP6, "STEP6_trajectory_annual_metrics.rds"))
HOLDOUT_LABS <- c(H1 = "H1: train 2015-2019, holdout 2020-2025", H2 = "H2: train 2015-2021, holdout 2022-2025")
GEN_LABS <- c(G1_fixed_SD_040 = "G1: fixed SD=0.40", G2_estimated_tau = "G2: estimated tau_year")

# ================================================================
# FIG21 -- Susceptibility and transmission diagnostics (G2, both holdouts)
# ================================================================
wq_g2 <- weekly_q |> filter(generator == "G2_estimated_tau") |> mutate(holdout_lab = HOLDOUT_LABS[holdout])
rect_2022 <- tibble(xmin = epidemic_period$start, xmax = epidemic_period$end)

panel_s <- ggplot(wq_g2, aes(week_start)) +
  geom_rect(data = rect_2022, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey85", alpha = 0.6) +
  geom_ribbon(aes(ymin = S_prop_lo95, ymax = S_prop_hi95), fill = OKABE[["blue"]], alpha = 0.15) +
  geom_ribbon(aes(ymin = S_prop_lo50, ymax = S_prop_hi50), fill = OKABE[["blue"]], alpha = 0.3) +
  geom_line(aes(y = S_prop_median), colour = OKABE[["blue"]], linewidth = 0.4) +
  facet_wrap(~holdout_lab, scales = "free_x") +
  labs(y = "S(t)/N(t)", x = NULL, title = "FIG21. Susceptibility and transmission diagnostics (G2, estimated tau_year)",
       subtitle = "Grey band = observed central-90% 2022 epidemic period (2022-02-27 to 2022-08-07).\nRibbons: 50%/95% across posterior draws x stochastic replicates.") +
  theme_pub

panel_r0 <- ggplot(wq_g2, aes(week_start)) +
  geom_rect(data = rect_2022, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey85", alpha = 0.6) +
  geom_ribbon(aes(ymin = R0_lo95, ymax = R0_hi95), fill = OKABE[["purple"]], alpha = 0.15) +
  geom_ribbon(aes(ymin = R0_lo50, ymax = R0_hi50), fill = OKABE[["purple"]], alpha = 0.3) +
  geom_line(aes(y = R0_median), colour = OKABE[["purple"]], linewidth = 0.4) +
  facet_wrap(~holdout_lab, scales = "free_x") + coord_cartesian(ylim = c(0, 4)) +
  labs(y = "R0(t)", x = NULL) + theme_pub + theme(strip.text = element_blank())

panel_reff <- ggplot(wq_g2, aes(week_start)) +
  geom_rect(data = rect_2022, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey85", alpha = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_ribbon(aes(ymin = R_eff_lo95, ymax = R_eff_hi95), fill = OKABE[["orange"]], alpha = 0.15) +
  geom_ribbon(aes(ymin = R_eff_lo50, ymax = R_eff_hi50), fill = OKABE[["orange"]], alpha = 0.3) +
  geom_line(aes(y = R_eff_median), colour = OKABE[["orange"]], linewidth = 0.4) +
  facet_wrap(~holdout_lab, scales = "free_x") + coord_cartesian(ylim = c(0, 3)) +
  labs(y = "Reff(t)", x = "Week", caption = PROV_CAPTION) + theme_pub + theme(strip.text = element_blank())

fig21 <- panel_s / panel_r0 / panel_reff
save_fig(fig21, "FIG21_susceptibility_transmission_diagnostics", height_mm = 200)

# ================================================================
# FIG22 -- Introduction and establishment diagnostics
# ================================================================
entry_plot_df <- entry_draws |> mutate(holdout_lab = HOLDOUT_LABS[holdout],
                                        x_plot = pmax(entry2022_weighted_infectiousness, 0.01)) # floor for log scale display only
panel_entry <- ggplot(entry_plot_df, aes(holdout_lab, x_plot, fill = holdout)) +
  geom_violin(trim = TRUE, alpha = 0.7, linewidth = 0.2) +
  geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white") +
  scale_y_log10(labels = scales::comma) + scale_fill_manual(values = c(H1 = unname(OKABE["blue"]), H2 = unname(OKABE["green"])), guide = "none") +
  labs(title = "FIG22. Introduction and establishment diagnostics (2022)", x = NULL,
       y = "Entry-2022 weighted infectiousness\n(gen.-weighted X, 8wk before 2022, log)",
       subtitle = "A: H1's entry infectiousness comes from its OWN 2-year stochastic holdout simulation (2020-2021, no real data);\nH2's comes directly from the training-only posterior fitted to REAL observed 2015-2021 surveillance.") +
  theme_pub

mech_plot_df <- mech_class |> mutate(holdout_lab = HOLDOUT_LABS[holdout], generator_lab = GEN_LABS[generator])
panel_mech <- ggplot(mech_plot_df, aes(generator_lab, probability, fill = mechanism_class)) +
  geom_col(position = "stack", width = 0.7) +
  facet_wrap(~holdout_lab, scales = "free_x") +
  scale_fill_manual(name = "2022 outcome class", values = c(
    "1 No positive M during Reff>1" = unname(OKABE["grey"]), "2 Pressure present; no sustained local dominance" = unname(OKABE["yellow"]),
    "3 Established; no major annual epidemic" = unname(OKABE["blue"]), "4 Major annual epidemic" = unname(OKABE["orange"]))) +
  scale_y_continuous(labels = scales::percent) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  labs(x = NULL, y = "Share of trajectories (draws x replicates)",
       subtitle = "B: mechanism classification -- establishment (class 3+4) is common in H1/G2 (90.9%), but crossing the major-epidemic\nthreshold (class 4) is rare (10.6%); in H2/G2 the two are nearly identical (94.0% vs 99.6%).",
       caption = PROV_CAPTION) +
  theme_pub + theme(legend.position = "bottom", legend.text = element_text(size = 6.5))

fig22 <- panel_entry / panel_mech + plot_layout(heights = c(1, 1.1))
save_fig(fig22, "FIG22_introduction_establishment_diagnostics", height_mm = 220)

# ================================================================
# FIG23 -- Mechanistic explanation of the 2022 prediction failure
# ================================================================
oc_vars <- c("cumulative_latent_infections", "peak_Reff", "year_effect", "A_year")
oc_df <- outcome_contrasts |> filter(holdout == "H1", quantity %in% oc_vars) |>
  mutate(quantity_lab = recode(quantity, cumulative_latent_infections = "Cumulative latent\ninfections (log10)",
                                 peak_Reff = "Peak Reff", year_effect = "Realised year_effect[2022]", A_year = "Realised A_year[2022]"))
panel_outcome <- ggplot(oc_df, aes(outcome, median, colour = outcome)) +
  geom_pointrange(aes(ymin = lo80, ymax = hi80), linewidth = 0.8, fatten = 2.5) +
  geom_pointrange(aes(ymin = lo95, ymax = hi95), linewidth = 0.35, fatten = 1) +
  facet_wrap(~quantity_lab, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c("No major annual epidemic" = unname(OKABE["grey"]), "Major annual epidemic" = unname(OKABE["orange"])), guide = "none") +
  labs(title = "FIG23. Mechanistic explanation of the 2022 prediction failure (H1/G2)", x = NULL, y = NULL,
       subtitle = "A: H1/G2 trajectories split by 2022 outcome (descriptive; conditioning on simulated outcome, not a causal contrast).\nPoint=median, thick bar=80% range, thin bar=95% range across draws x replicates.") +
  theme_pub + theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6.5))

annual_2022 <- annual_traj |> filter(year == 2022) |> mutate(holdout_lab = HOLDOUT_LABS[holdout], generator_lab = GEN_LABS[generator])
obs_2022 <- tibble(holdout_lab = HOLDOUT_LABS[c("H1", "H2")], observed = epidemic_period$observed_2022_cases)
panel_predictive <- ggplot(annual_2022, aes(generator_lab, pmax(annual_reported_cases, 1), fill = generator_lab)) +
  geom_violin(trim = TRUE, alpha = 0.7, linewidth = 0.2) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  geom_hline(data = obs_2022, aes(yintercept = observed), colour = "black", linetype = "dotted", linewidth = 0.6) +
  facet_wrap(~holdout_lab, scales = "free_x") +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = c("G1: fixed SD=0.40" = unname(OKABE["purple"]), "G2: estimated tau_year" = unname(OKABE["blue"])), guide = "none") +
  labs(x = NULL, y = "Simulated 2022 annual reported cases (log scale)",
       subtitle = "B: simulated 2022 predictive distributions vs. the observed total (dotted line = 50,986 cases, never used to fit/calibrate).",
       caption = PROV_CAPTION) +
  theme_pub + theme(axis.text.x = element_text(angle = 15, hjust = 1))

fig23 <- panel_outcome / panel_predictive + plot_layout(heights = c(1, 1.2))
save_fig(fig23, "FIG23_mechanistic_explanation_2022_failure", height_mm = 210)

message("\n[Step6] Figures FIG21-23 complete -> ", DIR_FIGURES_STEP6)
