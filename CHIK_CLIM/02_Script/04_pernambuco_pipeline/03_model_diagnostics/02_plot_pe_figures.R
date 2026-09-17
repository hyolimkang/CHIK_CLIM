# Pernambuco external replication -- required figures (Section 15).

required_packages <- c("rstan", "dplyr", "tibble", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(ggplot2); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/model_fits/v4_9_replication")
# NOTE: kept at the TOP level (not a subfolder) -- these are the exact
# named deliverables from the design spec's Section 15 output list.
figure_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/figures/pernambuco_v4_9_replication")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

A <- readRDS(file.path(base_dir, "outputs/modelA_caseonly/renewal_pe_global_q_fit_modelA_caseonly.rds"))
B <- readRDS(file.path(base_dir, "outputs/modelB_U14/renewal_pe_global_q_fit_modelB_U14.rds"))

plot_case_fit <- function(bundle, title) {
  w <- bundle$weekly_data; dates <- as.Date(w$week_start)
  C_pred <- rstan::extract(bundle$fit, "C_pred")$C_pred
  med <- apply(C_pred, 2, median); lo <- apply(C_pred, 2, quantile, .025); hi <- apply(C_pred, 2, quantile, .975)
  df <- tibble(week_start = dates, observed = w$cases, median = med, lo = lo, hi = hi)
  ggplot(df, aes(week_start)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#56B4E9", alpha = 0.2) +
    geom_line(aes(y = median), colour = "#315A7D") +
    geom_point(aes(y = observed), size = 0.4, colour = "black", alpha = 0.6) +
    labs(title = title, x = NULL, y = "weekly cases") + theme_v4
}
ggsave(file.path(figure_dir, "pernambuco_global_q_case_only.png"), plot_case_fit(A, "Pernambuco Model A: case-only, global q"), width = 200, height = 100, units = "mm", dpi = 300)
ggsave(file.path(figure_dir, "pernambuco_global_q_U14_serology.png"), plot_case_fit(B, "Pernambuco Model B: global q + U14 Recife serology"), width = 200, height = 100, units = "mm", dpi = 300)
message("[saved] pernambuco_global_q_case_only.png, pernambuco_global_q_U14_serology.png")

# ---- q prior vs Model A vs Model B -----------------------------------------
q_A <- as.vector(rstan::extract(A$fit, "q")$q)
q_B <- as.vector(rstan::extract(B$fit, "q")$q)
q_prior <- as.vector(rstan::extract(B$fit, "q_prior_draw")$q_prior_draw)
q_df <- bind_rows(tibble(q = q_prior, dist = "Prior"), tibble(q = q_A, dist = "A: case-only"), tibble(q = q_B, dist = "B: U14 serology"))
p_q <- ggplot(q_df, aes(q, fill = dist, colour = dist)) +
  geom_density(alpha = 0.3) +
  scale_x_continuous(limits = c(0, 0.6)) +
  scale_fill_manual(values = c(Prior = "grey60", "A: case-only" = "#D55E00", "B: U14 serology" = "#76558F")) +
  scale_colour_manual(values = c(Prior = "grey40", "A: case-only" = "#A34700", "B: U14 serology" = "#4B2E6B")) +
  labs(title = "Pernambuco: q prior vs case-only vs U14-serology posterior", x = "q", y = "density") +
  theme_v4 + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(figure_dir, "pernambuco_q_prior_vs_case_vs_serology.png"), p_q, width = 190, height = 120, units = "mm", dpi = 300)
message("[saved] pernambuco_q_prior_vs_case_vs_serology.png")

# ---- Susceptibility trajectory: A vs B -------------------------------------
immune_traj <- function(bundle, label) {
  dates <- as.Date(bundle$weekly_data$week_start)
  ip <- rstan::extract(bundle$fit, "immune_prop")$immune_prop
  tibble(week_start = dates, model = label, median = apply(ip, 2, median), lo = apply(ip, 2, quantile, .025), hi = apply(ip, 2, quantile, .975))
}
traj_df <- bind_rows(immune_traj(A, "A: case-only"), immune_traj(B, "B: U14 serology"))
p_susc <- ggplot(traj_df, aes(week_start, median, colour = model, fill = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_colour_manual(values = c("A: case-only" = "#D55E00", "B: U14 serology" = "#76558F")) +
  scale_fill_manual(values = c("A: case-only" = "#D55E00", "B: U14 serology" = "#76558F")) +
  labs(title = "Pernambuco cumulative immune fraction: case-only vs U14 serology", x = NULL, y = "immune fraction") +
  theme_v4 + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(figure_dir, "pernambuco_susceptibility_case_vs_serology.png"), p_susc, width = 200, height = 120, units = "mm", dpi = 300)
message("[saved] pernambuco_susceptibility_case_vs_serology.png")

# ---- U14 serology fit -------------------------------------------------------
p_state_window <- rstan::extract(B$fit, "p_state_window")$p_state_window[, 1]
p_site_window <- rstan::extract(B$fit, "p_site_window")$p_site_window[, 1]
sero_df <- tibble(
  quantity = c("Observed (Recife)", "PE state-level (no geo offset)", "Geo-adjusted Recife"),
  median = c(770/2070, median(p_state_window), median(p_site_window)),
  lo = c(0.340, quantile(p_state_window, .025), quantile(p_site_window, .025)),
  hi = c(0.404, quantile(p_state_window, .975), quantile(p_site_window, .975))
)
p_sero <- ggplot(sero_df, aes(quantity, median, colour = quantity)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), fatten = 3) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_colour_manual(values = c("Observed (Recife)" = "black", "PE state-level (no geo offset)" = "#009E73", "Geo-adjusted Recife" = "#76558F")) +
  labs(title = "Pernambuco U14 (Recife) serology fit", x = NULL, y = "seroprevalence") +
  theme_v4 + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(figure_dir, "pernambuco_U14_serology_fit.png"), p_sero, width = 160, height = 130, units = "mm", dpi = 300)
message("[saved] pernambuco_U14_serology_fit.png")
