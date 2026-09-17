# Figures for the Bahia global-q identification pilot (Sections
# "IDENTIFICATION DIAGNOSTICS" / "CRITICAL COMPARISONS" / "SEROLOGY
# INFORMATION FLOW" of the design spec). Diagnostic only -- does not alter
# the model. Compares:
#   A = Bahia v4.9 case-only, fixed q=0.05
#   B = Bahia v4.9 multisite serology, fixed q=0.05
#   C = Bahia v4.9 multisite serology, GLOBAL q estimated (this pilot)

required_packages <- c("rstan", "dplyr", "tibble", "ggplot2", "tidyr", "gridExtra")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(ggplot2); library(tidyr); library(gridExtra) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
base_dir <- file.path(root, "03_Output/03_bahia_pipeline/model_fits/v4_9_replication")
figure_dir <- file.path(root, "03_Output/03_bahia_pipeline/figures/renewal_bahia_v4_9_global_q_multisite_serology/comparison")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

A <- readRDS(file.path(base_dir, "outputs/q0.05/renewal_bahia_v4_9_fit_q0.05.rds"))
Bm <- readRDS(file.path(base_dir, "outputs_multisite_serology/q0.05/renewal_bahia_multisite_fit_q0.05.rds"))
Cg <- readRDS(file.path(base_dir, "outputs_global_q_multisite_serology/on_full2/renewal_bahia_global_q_fit_on_full2.rds"))

# ---- Figure 1: q prior vs posterior ----------------------------------
q_post <- as.vector(rstan::extract(Cg$fit, "q")$q)
q_prior <- as.vector(rstan::extract(Cg$fit, "q_prior_draw")$q_prior_draw)
q_df <- bind_rows(
  tibble(q = q_prior, dist = "Prior"),
  tibble(q = q_post, dist = "Posterior")
)
fixed_q_grid <- c(0.05, 0.10, 0.15, 0.20)
p1 <- ggplot(q_df, aes(q, fill = dist, colour = dist)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  geom_vline(xintercept = fixed_q_grid, linetype = 3, colour = "grey50") +
  annotate("text", x = fixed_q_grid, y = Inf, label = sprintf("q=%.2f", fixed_q_grid),
           angle = 90, vjust = 1.2, hjust = 1.1, size = 2.5, colour = "grey40") +
  scale_x_continuous(limits = c(0, 0.6)) +
  scale_fill_manual(values = c(Prior = "grey60", Posterior = "#76558F")) +
  scale_colour_manual(values = c(Prior = "grey40", Posterior = "#4B2E6B")) +
  labs(title = "Bahia global-q (stabilized, HMC PASS): prior vs posterior of q",
       subtitle = "Dotted lines = previously explored fixed-q sensitivity grid",
       x = "q", y = "density") +
  theme_v4 + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(figure_dir, "fig1_q_prior_vs_posterior.png"), p1, width = 160, height = 110, units = "mm", dpi = 300)

# ---- Figure 2: immune fraction trajectory, A vs B vs C ----------------
immune_traj <- function(bundle, label) {
  dates <- as.Date(bundle$weekly_data$week_start)
  ip <- rstan::extract(bundle$fit, "immune_prop")$immune_prop
  tibble(week_start = dates, model = label,
         median = apply(ip, 2, median), lo = apply(ip, 2, quantile, .025), hi = apply(ip, 2, quantile, .975))
}
traj_df <- bind_rows(
  immune_traj(A, "A: case-only, fixed q=0.05"),
  immune_traj(Bm, "B: multisite serology, fixed q=0.05"),
  immune_traj(Cg, "C: multisite serology, global q (stabilized)")
)
p2 <- ggplot(traj_df, aes(week_start, median, colour = model, fill = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.5) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_colour_manual(values = c("A: case-only, fixed q=0.05" = "grey40",
                                  "B: multisite serology, fixed q=0.05" = "#009E73",
                                  "C: multisite serology, global q (stabilized)" = "#76558F")) +
  scale_fill_manual(values = c("A: case-only, fixed q=0.05" = "grey40",
                                "B: multisite serology, fixed q=0.05" = "#009E73",
                                "C: multisite serology, global q (stabilized)" = "#76558F")) +
  labs(title = "Bahia cumulative immune fraction (state-level): A vs B vs C",
       subtitle = "C's much lower median and far wider 95% CrI reflect q floating instead of fixed at 0.05",
       x = NULL, y = "immune fraction (U/N)") +
  theme_v4 + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(figure_dir, "fig2_immune_fraction_ABC_comparison.png"), p2, width = 200, height = 120, units = "mm", dpi = 300)

# ---- Figure 3: eta_geo comparison, global-q vs fixed q=0.05 -----------
eta_C <- rstan::extract(Cg$fit, "eta_geo")$eta_geo
eta_B <- rstan::extract(Bm$fit, "eta_geo")$eta_geo
eta_df <- bind_rows(
  tibble(sero_id = Cg$sero_audit$sero_id, model = "C: global q (stabilized)",
         median = apply(eta_C, 2, median), lo = apply(eta_C, 2, quantile, .025), hi = apply(eta_C, 2, quantile, .975)),
  tibble(sero_id = Bm$sero_audit$sero_id, model = "B: fixed q=0.05",
         median = apply(eta_B, 2, median), lo = apply(eta_B, 2, quantile, .025), hi = apply(eta_B, 2, quantile, .975))
)
p3 <- ggplot(eta_df, aes(sero_id, median, colour = model)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.4), fatten = 1.8) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  scale_colour_manual(values = c("B: fixed q=0.05" = "#009E73", "C: global q (stabilized)" = "#76558F")) +
  labs(title = "Geographic offset (eta_geo) by survey: fixed q vs global q",
       subtitle = "If serology pinned down state-scale via q, eta_geo under C should shrink toward 0 -- it does not",
       x = "survey", y = "eta_geo (log-odds)") +
  theme_v4 + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(figure_dir, "fig3_eta_geo_comparison_global_vs_fixed.png"), p3, width = 160, height = 110, units = "mm", dpi = 300)

# ---- Figure 4: posterior scatter -- q vs alpha_R and q vs immune_2025 --
alpha_R <- as.vector(rstan::extract(Cg$fit, "alpha_R")$alpha_R)
dates_C <- as.Date(Cg$weekly_data$week_start)
idx_2025 <- length(dates_C)
immune_2025 <- rstan::extract(Cg$fit, "immune_prop")$immune_prop[, idx_2025]
scatter_df <- tibble(q = q_post, alpha_R = alpha_R, immune_2025 = immune_2025)

p4a <- ggplot(scatter_df, aes(q, alpha_R)) +
  geom_point(alpha = 0.25, size = 0.7, colour = "#76558F") +
  labs(title = sprintf("q vs alpha_R (cor = %.2f)", cor(q_post, alpha_R)), x = "q", y = "alpha_R") +
  theme_v4
p4b <- ggplot(scatter_df, aes(q, immune_2025)) +
  geom_point(alpha = 0.25, size = 0.7, colour = "#009E73") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(title = sprintf("q vs immune fraction, 2025 (cor = %.2f)", cor(q_post, immune_2025)),
       x = "q", y = "immune fraction (2025)") +
  theme_v4
p4 <- gridExtra::grid.arrange(p4a, p4b, ncol = 2,
                               top = grid::textGrob("Posterior ridge diagnostics: q vs key transmission/immunity quantities", gp = grid::gpar(fontsize = 10)))
ggsave(file.path(figure_dir, "fig4_q_ridge_scatter.png"), p4, width = 200, height = 100, units = "mm", dpi = 300)

message("[saved] figures under: ", figure_dir)
