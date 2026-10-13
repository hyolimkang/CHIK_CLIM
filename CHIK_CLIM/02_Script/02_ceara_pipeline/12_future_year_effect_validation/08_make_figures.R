# PHASE 5 -- Publication-quality figures FIG17-20.

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr", "ggplot2", "patchwork", "scales", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(tidyr); library(ggplot2); library(patchwork) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/12_future_year_effect_validation/00_config_year_effect_validation.R"))

OKABE <- c(blue = "#0072B2", orange = "#D55E00", yellow = "#E69F00", green = "#009E73", purple = "#CC79A7", grey = "#56B4E9")
theme_pub <- theme_classic(base_size = 9) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey30"), legend.position = "top",
        plot.caption = element_text(size = 7, colour = "grey40", face = "italic"))
RETRO_NOTE <- "Retrospective, observed-climate-conditional simulation from a TRAINING-ONLY posterior; no holdout cases or 2022 seed used."

save_fig <- function(plot, name, width_mm = 180, height_mm = 120) {
  ggsave(file.path(DIR_FIGURES_YV, paste0(name, ".pdf")), plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(file.path(DIR_FIGURES_YV, paste0(name, ".png")), plot, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
  ggsave(file.path(DIR_FIGURES_YV, paste0(name, ".svg")), plot, width = width_mm, height = height_mm, units = "mm")
  message("[saved] ", name, " (.pdf/.png@600dpi/.svg)")
}

PRIOR_COLOURS <- c(P1 = OKABE[["blue"]], P2 = OKABE[["orange"]], P3 = OKABE[["green"]])
PRIOR_LABELS <- c(P1 = "P1: HalfNormal(0,0.20)", P2 = "P2: HalfNormal(0,0.40)", P3 = "P3: HalfNormal(0,0.80)")

# ================================================================
# FIG17 -- Prior and posterior distributions of tau_year under P1-P3
# ================================================================
tau_draws <- readRDS(file.path(DIR_RESULTS_YV, "PHASE1_tau_year_draws.rds"))
tau_df <- bind_rows(lapply(names(tau_draws), function(p) tibble(prior_id = p, tau_year = tau_draws[[p]])))
prior_x <- seq(0, 1.6, length.out = 400)
prior_curves <- bind_rows(lapply(names(TAU_PRIORS), function(p) tibble(prior_id = p, x = prior_x, density = 2 * dnorm(prior_x, 0, TAU_PRIORS[[p]])))) # HalfNormal density = 2*Normal density for x>=0

fig17 <- ggplot() +
  geom_line(data = prior_curves, aes(x, density, colour = prior_id), linetype = "dashed", linewidth = 0.5) +
  geom_density(data = tau_df, aes(tau_year, fill = prior_id, colour = prior_id), alpha = 0.35, linewidth = 0.6) +
  geom_vline(xintercept = 0.40, linetype = "dotted", colour = "black", linewidth = 0.6) +
  annotate("text", x = 0.40, y = Inf, label = "historical fixed SD=0.40", vjust = 1.3, hjust = -0.05, size = 2.8, angle = 90) +
  facet_wrap(~ prior_id, labeller = as_labeller(PRIOR_LABELS)) +
  scale_colour_manual(values = PRIOR_COLOURS, guide = "none") + scale_fill_manual(values = PRIOR_COLOURS, guide = "none") +
  coord_cartesian(xlim = c(0, 1.4)) +
  labs(title = "FIG17: tau_year prior (dashed) vs. posterior (filled) under P1-P3",
       subtitle = "Dotted vertical line: the historical fixed value (year_effect_prior_sd=0.40) this study replaces with an estimated parameter.",
       x = "tau_year", y = "Density") + theme_pub
save_fig(fig17, "FIG17_tau_year_prior_posterior", width_mm = 180, height_mm = 100)

# ================================================================
# FIG18 -- Historical year effects vs. posterior-predictive new-year effects (G1 vs G2)
# ================================================================
hist_ye <- read_csv(file.path(DIR_TABLES_YV, "PHASE1_implied_year_effect_by_prior.csv"), show_col_types = FALSE) |> dplyr::filter(prior_id == "P2")
gen_traj <- readRDS(file.path(DIR_RESULTS_YV, "PHASE2_generator_trajectories.rds"))
# Reconstruct G1/G2 future year-effect realisations is not directly stored in trajectories; use the paired-differences/annual_R0 proxy instead via a lightweight re-derivation is avoided here -- read from a dedicated per-generator year-effect sample saved by Phase 2 if present, else skip gracefully.
ye_future_path <- file.path(DIR_RESULTS_YV, "PHASE2_future_year_effect_samples.csv")
if (file.exists(ye_future_path)) {
  future_ye <- read_csv(ye_future_path, show_col_types = FALSE)
  fig18 <- ggplot() +
    geom_density(data = future_ye |> dplyr::filter(generator != "G0_neutral"), aes(year_effect, fill = generator, colour = generator), alpha = 0.4) +
    geom_point(data = hist_ye, aes(x = median, y = 0), shape = 17, size = 2, colour = "black") +
    geom_errorbarh(data = hist_ye, aes(xmin = lo95, xmax = hi95, y = 0), height = 0.3, colour = "black", alpha = 0.5) +
    scale_fill_manual(name = NULL, values = c(G1_fixed_SD_040 = OKABE[["blue"]], G2_estimated_tau = OKABE[["orange"]]), labels = c(G1_fixed_SD_040 = "G1: fixed SD=0.40", G2_estimated_tau = "G2: estimated tau_year")) +
    scale_colour_manual(name = NULL, values = c(G1_fixed_SD_040 = OKABE[["blue"]], G2_estimated_tau = OKABE[["orange"]]), labels = c(G1_fixed_SD_040 = "G1: fixed SD=0.40", G2_estimated_tau = "G2: estimated tau_year")) +
    labs(title = "FIG18: Historical vs. posterior-predictive future year effects",
         subtitle = "Black triangles+bars: historical year_effect[y] posterior medians/95% CrI (2015-2025, P2 fit). Density: future-year-effect generators.",
         x = "year_effect (log R0 scale)", y = "Density", caption = RETRO_NOTE) + theme_pub
  save_fig(fig18, "FIG18_historical_vs_future_year_effects", width_mm = 180, height_mm = 110)
} else {
  message("[Phase5] PHASE2_future_year_effect_samples.csv not found -- FIG18 skipped (see Phase2 script note).")
}

# ================================================================
# FIG19 -- Retrospective epidemic trajectories, H1 & H2
# ================================================================
annual_pi <- read_csv(file.path(DIR_TABLES_YV, "PHASE4_annual_PI_vs_observed.csv"), show_col_types = FALSE)
obs_cells <- annual_pi |> distinct(cell) |> pull(cell)
# Primary cells: observed-climate, reference importation, both generators, both holdouts
primary_cells <- obs_cells[grepl(sprintf("_observed_%s$", IMPORTATION_REFERENCE), obs_cells)]
fig19_df <- annual_pi |> dplyr::filter(cell %in% primary_cells) |>
  mutate(holdout = sub("_.*", "", cell), generator = ifelse(grepl("G1_fixed_SD_040", cell), "G1: fixed SD=0.40", "G2: estimated tau_year"))

fig19 <- ggplot(fig19_df, aes(year)) +
  geom_ribbon(aes(ymin = lo025, ymax = hi975, fill = generator), alpha = 0.18) +
  geom_ribbon(aes(ymin = lo10, ymax = hi90, fill = generator), alpha = 0.28) +
  geom_ribbon(aes(ymin = lo25, ymax = hi75, fill = generator), alpha = 0.40) +
  geom_line(aes(y = median, colour = generator), linewidth = 0.6) +
  geom_point(aes(y = cases), colour = "black", size = 1.8) +
  geom_line(aes(y = cases), colour = "black", linewidth = 0.4, linetype = "dotted") +
  facet_grid(holdout ~ generator, scales = "free_y") +
  scale_fill_manual(values = c("G1: fixed SD=0.40" = OKABE[["blue"]], "G2: estimated tau_year" = OKABE[["orange"]]), guide = "none") +
  scale_colour_manual(values = c("G1: fixed SD=0.40" = OKABE[["blue"]], "G2: estimated tau_year" = OKABE[["orange"]]), guide = "none") +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 100), labels = scales::label_number(big.mark = ",")) +
  labs(title = "FIG19: Retrospective holdout validation -- annual reported cases",
       subtitle = "Black dotted line + points: OBSERVED (never used to fit/calibrate the holdout simulation).\nColoured bands: 50%/80%/95% predictive intervals; line: predictive median.",
       x = NULL, y = "Annual reported cases (pseudo-log scale)", caption = RETRO_NOTE) + theme_pub
save_fig(fig19, "FIG19_retrospective_holdout_trajectories", width_mm = 190, height_mm = 150)

# ================================================================
# FIG20 -- Fixed-SD vs. estimated-tau generators: cum infections, outbreak probability, annual peak R0
# ================================================================
gen_summary <- read_csv(file.path(DIR_TABLES_YV, "PHASE2_generator_summary.csv"), show_col_types = FALSE)
peak_extreme <- read_csv(file.path(DIR_TABLES_YV, "PHASE2_peak_R0_extremes_by_generator.csv"), show_col_types = FALSE)
SHORT_GEN <- c(G0_neutral = "G0", G1_fixed_SD_040 = "G1", G2_estimated_tau = "G2")
GEN_COLOURS <- c(G0 = OKABE[["grey"]], G1 = OKABE[["blue"]], G2 = OKABE[["orange"]])
per_traj <- read_csv(file.path(DIR_RESULTS_YV, "PHASE2_per_trajectory_summary.csv"), show_col_types = FALSE) |>
  mutate(gen_label = recode(generator, !!!SHORT_GEN))

p_cum <- ggplot(per_traj, aes(gen_label, cum_10y + 1, fill = gen_label)) +
  geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3) + geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.7, linewidth = 0.3) +
  scale_y_log10(labels = scales::label_number(big.mark = ",")) +
  scale_fill_manual(values = GEN_COLOURS, guide = "none") +
  labs(subtitle = "A. 10-year cumulative\nlatent infections (log scale)", x = NULL, y = "Cumulative infections") + theme_pub

p_major <- ggplot(gen_summary |> mutate(gen_label = recode(generator, !!!SHORT_GEN)), aes(gen_label, proportion_major_epidemic, fill = gen_label)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = GEN_COLOURS, guide = "none") +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(subtitle = "B. P(>=1 major epidemic,\n2026-2035)", x = NULL, y = "Probability") + theme_pub

p_peak <- ggplot(peak_extreme |> mutate(gen_label = recode(generator, !!!SHORT_GEN)), aes(gen_label, median_peak_R0, fill = gen_label)) +
  geom_col(width = 0.6) + geom_errorbar(aes(ymin = median_peak_R0, ymax = p99_peak_R0), width = 0.15) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = GEN_COLOURS, guide = "none") +
  labs(subtitle = "C. Annual peak R0(t): median\n(bar) to 99th pct (whisker)", x = NULL, y = "Peak R0(t)") + theme_pub

fig20 <- (p_cum | p_major | p_peak) + plot_annotation(
  title = "FIG20: Fixed-SD vs. estimated-tau future generators",
  subtitle = "G0: neutral | G1: fixed SD=0.40 (existing Y1) | G2: estimated tau_year",
  caption = "Illustrative scenario simulation; no SSP forcing or vaccination applied.",
  theme = theme(plot.subtitle = element_text(size = 8, colour = "grey30"), plot.caption = element_text(size = 7, colour = "grey40", face = "italic")))
save_fig(fig20, "FIG20_generator_comparison", width_mm = 190, height_mm = 115)

message("\n[Phase5] Figures complete -> ", DIR_FIGURES_YV)
