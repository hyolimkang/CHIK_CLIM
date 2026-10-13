# STEP 7 -- Figures FIG24-27. Reads only already-produced Phase 0-6
# tables/results (no new simulation).

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(tidyr); library(ggplot2); library(patchwork) })
filter <- dplyr::filter

source(file.path(ROOT, "02_Script/02_ceara_pipeline/16_step7_importation_persistence/00_config_step7.R"))

OKABE <- c(blue = "#0072B2", orange = "#D55E00", yellow = "#E69F00", green = "#009E73", purple = "#CC79A7", grey = "#56B4E9")
theme_pub <- theme_classic(base_size = 9) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey30"), legend.position = "top",
        plot.caption = element_text(size = 7, colour = "grey40", face = "italic"))
save_fig <- function(plot, name, width_mm = 190, height_mm = 150) {
  ggsave(file.path(DIR_FIGURES_STEP7, paste0(name, ".pdf")), plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(file.path(DIR_FIGURES_STEP7, paste0(name, ".png")), plot, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
  message("[saved] ", name)
}
PROV_CAPTION <- paste0(PROVISIONAL_NOTE, "\nSensitivity ARMS (diagnostic, prespecified lambda grid), not calibrated model selection -- no lambda is chosen as a new reference.")
LAMBDA_LABS <- c(I0 = "I0: lambda=0", I1 = "I1: lambda=0.25", I2 = "I2: lambda=1 (baseline)", I3 = "I3: lambda=4")

# ================================================================
# FIG24 -- importation intensity vs entry-2022 weighted infectiousness
# ================================================================
entry_draws <- read_csv(file.path(DIR_RESULTS_STEP7, "PHASE1_entry_infectiousness_draws.csv"), show_col_types = FALSE)
h1_entry <- entry_draws |> filter(holdout == "H1") |> mutate(lambda_lab = factor(LAMBDA_LABS[lambda_id], levels = LAMBDA_LABS))
h2_ref <- median(entry_draws$entry_weighted_infectiousness[entry_draws$holdout == "H2"]) # invariant across lambda by construction (training-posterior quantity)

fig24 <- ggplot(h1_entry, aes(lambda_lab, pmax(entry_weighted_infectiousness, 0.01), fill = lambda_lab)) +
  geom_violin(trim = TRUE, alpha = 0.75, linewidth = 0.2) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  geom_hline(aes(yintercept = h2_ref, linetype = "H2 entry infectiousness (training-posterior,\ninvariant to holdout lambda -- diagnostic reference)"), colour = "black", linewidth = 0.6) +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = unname(OKABE[c("grey", "yellow", "blue", "orange")]), guide = "none") +
  scale_linetype_manual(name = NULL, values = "dotted") +
  labs(title = "FIG24. Importation intensity vs. entry-2022 weighted infectiousness (H1/G2)",
       subtitle = "Even lambda=4 (4x baseline) leaves H1 ~50x below H2's real-data-informed entry level.",
       x = NULL, y = "Entry-2022 weighted infectiousness (log scale)", caption = PROV_CAPTION) +
  theme_pub + theme(legend.position = "bottom")
save_fig(fig24, "FIG24_importation_vs_entry_infectiousness")

# ================================================================
# FIG25 -- pre-2022 vs 2022-onward importation experiments
# ================================================================
phase2 <- read_csv(file.path(DIR_TABLES_STEP7, "PHASE2_pre_post_2022_split_summary.csv"), show_col_types = FALSE) |>
  mutate(varied_lambda = ifelse(varied_period == "pre-2022 (2020-2021)", lambda_pre, lambda_post),
         experiment_lab = ifelse(grepl("ExperimentA", experiment), "Exp A: vary PRE-2022 (post fixed=1)", "Exp B: vary 2022-ONWARD (pre fixed=1)"))
EXP_COLOURS <- c("Exp A: vary PRE-2022 (post fixed=1)" = unname(OKABE["blue"]), "Exp B: vary 2022-ONWARD (pre fixed=1)" = unname(OKABE["orange"]))

p25a <- ggplot(phase2, aes(varied_lambda, pmax(median_entry_weighted_infectiousness, 0.01), colour = experiment_lab)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) + scale_y_log10() +
  scale_colour_manual(name = NULL, values = EXP_COLOURS) +
  labs(x = "Varied lambda_import", y = "Median entry-2022\ninfectiousness (log)") + theme_pub
p25b <- ggplot(phase2, aes(varied_lambda, prob_established, colour = experiment_lab)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_colour_manual(values = EXP_COLOURS, guide = "none") +
  labs(x = "Varied lambda_import", y = "P(establishment)") + theme_pub + theme(legend.position = "none")
p25c <- ggplot(phase2, aes(varied_lambda, prob_major_epidemic, colour = experiment_lab)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) + scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = EXP_COLOURS, guide = "none") +
  labs(x = "Varied lambda_import", y = "P(major epidemic)") + theme_pub + theme(legend.position = "none")

fig25 <- (p25a | p25b | p25c) + plot_layout(guides = "collect") & theme(legend.position = "top")
fig25 <- fig25 +
  plot_annotation(title = "FIG25. Pre-2022 vs. 2022-onward importation: which matters more for the 2022 outcome? (H1/G2)",
                   subtitle = "Pre-2022 pressure (blue) moves P(major epidemic) ~2.7x more than equal-sized changes to 2022-onward pressure (orange).",
                   caption = PROV_CAPTION, theme = theme_pub)
save_fig(fig25, "FIG25_pre_vs_post_2022_importation", width_mm = 230, height_mm = 120)

# ================================================================
# FIG26 -- introduction timing, local transmission, and fadeout
# ================================================================
phase0 <- readRDS(file.path(DIR_RESULTS_STEP7, "PHASE0_baseline_lambda1_simulations.rds"))
T1 <- readRDS(file.path(DIR_RESULTS_STEP7, "PHASE3_T1_permuted_timing_trajectories.rds"))
example_ids <- phase0$H1$trajectories |> distinct(posterior_draw_id, stochastic_replicate_id) |> slice_sample(n = 2, replace = FALSE)

make_example_df <- function(traj, label, ids) {
  traj |> inner_join(ids, by = c("posterior_draw_id", "stochastic_replicate_id")) |>
    mutate(pair_id = paste0("draw ", posterior_draw_id, " / rep ", stochastic_replicate_id), arm = label,
           year = as.integer(format(week_start, "%Y"))) |> filter(year %in% 2021:2022)
}
example_T0 <- make_example_df(phase0$H1$trajectories, "T0: original timing", example_ids)
example_T1 <- make_example_df(T1$trajectories, "T1: permuted timing (same yearly total)", example_ids)
example_df <- bind_rows(example_T0, example_T1)

p26a <- ggplot(example_df, aes(week_start, X_latent, colour = arm)) +
  geom_line(linewidth = 0.4) + geom_point(aes(y = M_import * 5), shape = 1, size = 1, alpha = 0.5) +
  facet_wrap(~pair_id, scales = "free_y", ncol = 2) +
  scale_colour_manual(name = NULL, values = c("T0: original timing" = unname(OKABE["blue"]), "T1: permuted timing (same yearly total)" = unname(OKABE["orange"]))) +
  labs(x = NULL, y = "Latent infections X(t)\n(open circles: M_t x5, same scale)",
       subtitle = "A: two representative paired (draw, replicate) trajectories -- local infections (lines) vs. introduction pressure (circles).") +
  theme_pub

timing_summary <- read_csv(file.path(DIR_TABLES_STEP7, "PHASE3_timing_permutation_summary.csv"), show_col_types = FALSE)
paired_delta <- read_csv(file.path(DIR_TABLES_STEP7, "PHASE3_paired_timing_delta_2022_cases.csv"), show_col_types = FALSE)
p26b <- ggplot(timing_summary, aes(arm, median_annual_reported_cases_2022, fill = arm)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = c("T0_original_timing" = unname(OKABE["blue"]), "T1_permuted_timing" = unname(OKABE["orange"])), guide = "none") +
  labs(x = NULL, y = "Median 2022\nreported cases",
       subtitle = sprintf("B: ensemble summary -- near-identical (paired delta median=%.0f [%.0f, %.0f]).", paired_delta$median_delta, paired_delta$lo95_delta, paired_delta$hi95_delta),
       caption = PROV_CAPTION) +
  theme_pub

fig26 <- p26a / p26b + plot_layout(heights = c(1.4, 1)) + plot_annotation(title = "FIG26. Introduction timing, local transmission, and fadeout (H1/G2, lambda=1)", theme = theme_pub)
save_fig(fig26, "FIG26_timing_local_transmission_fadeout", height_mm = 210)

# ================================================================
# FIG27 -- predictive-performance sensitivity across the importation grid
# ================================================================
annual_perf <- read_csv(file.path(DIR_TABLES_STEP7, "PHASE6_annual_predictive_performance_grid.csv"), show_col_types = FALSE) |>
  mutate(lambda_lab = factor(LAMBDA_LABS[lambda_id], levels = LAMBDA_LABS), holdout_lab = ifelse(holdout == "H1", "H1 (2020-2025)", "H2 (2022-2025)"))
exceed_mc <- read_csv(file.path(DIR_TABLES_STEP7, "PHASE6_exceedance_probability_mc_grid.csv"), show_col_types = FALSE) |>
  mutate(lambda_lab = factor(LAMBDA_LABS[lambda_id], levels = LAMBDA_LABS), holdout_lab = ifelse(holdout == "H1", "H1 (2020-2025)", "H2 (2022-2025)"))
cum_perf <- read_csv(file.path(DIR_TABLES_STEP7, "PHASE6_cumulative_predictive_performance_grid.csv"), show_col_types = FALSE) |>
  mutate(lambda_lab = factor(LAMBDA_LABS[lambda_id], levels = LAMBDA_LABS), holdout_lab = ifelse(holdout == "H1", "H1 (2020-2025)", "H2 (2022-2025)"))

p27a <- ggplot(annual_perf, aes(factor(year), pmax(pred_median, 0.5), colour = lambda_lab, group = lambda_lab)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.5) +
  geom_point(aes(y = pmax(observed_cases, 0.5)), colour = "black", shape = 4, size = 2) +
  facet_wrap(~holdout_lab, scales = "free_x") + scale_y_log10(labels = scales::comma) +
  scale_colour_manual(name = NULL, values = unname(OKABE[c("grey", "yellow", "blue", "orange")])) +
  labs(x = "Year", y = "Annual reported cases\n(median, log scale; X = observed)",
       subtitle = "A: annual predictive medians across the importation grid vs. observed (black X).") + theme_pub

p27b <- ggplot(exceed_mc, aes(lambda_lab, point_estimate, colour = holdout_lab, group = holdout_lab)) +
  geom_line(linewidth = 0.5) + geom_pointrange(aes(ymin = cluster_bootstrap_lo95, ymax = cluster_bootstrap_hi95), fatten = 2) +
  scale_y_continuous(labels = scales::percent) + scale_colour_manual(name = NULL, values = c("H1 (2020-2025)" = unname(OKABE["blue"]), "H2 (2022-2025)" = unname(OKABE["orange"]))) +
  labs(x = NULL, y = "P(2022 cases >= observed)", subtitle = "B: 2022 exceedance probability\n(cluster-bootstrap 95% CI).") +
  theme_pub + theme(axis.text.x = element_text(angle = 20, hjust = 1))

p27c <- ggplot(cum_perf, aes(lambda_lab, WIS, colour = holdout_lab, group = holdout_lab)) +
  geom_line(linewidth = 0.5) + geom_point(size = 2) +
  scale_colour_manual(values = c("H1 (2020-2025)" = unname(OKABE["blue"]), "H2 (2022-2025)" = unname(OKABE["orange"])), guide = "none") +
  labs(x = NULL, y = "Multi-year cumulative WIS\n(lower = better)", subtitle = "C: multi-year cumulative WIS\n(proper scoring rule).") +
  theme_pub + theme(axis.text.x = element_text(angle = 20, hjust = 1))

fig27 <- p27a / (p27b | p27c) + plot_layout(heights = c(1.1, 1)) +
  plot_annotation(title = "FIG27. Predictive-performance sensitivity across the importation grid",
                   subtitle = "H1's WIS improves with higher lambda (panel C) but its 2022 exceedance probability stays near zero except at lambda=4 (panel B);\nH2 is comparatively insensitive to holdout-period lambda throughout.",
                   caption = PROV_CAPTION, theme = theme_pub)
save_fig(fig27, "FIG27_predictive_performance_sensitivity", height_mm = 230)

message("\n[Step7] Figures FIG24-27 complete -> ", DIR_FIGURES_STEP7)
