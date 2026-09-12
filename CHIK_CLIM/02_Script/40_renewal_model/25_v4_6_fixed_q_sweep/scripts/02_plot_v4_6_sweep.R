# v4.6 fixed-q sweep: summary figures across the pre-registered q grid
# (0.05, 0.10, 0.15, 0.20, 0.25, 0.30). Each q value is an independently
# HMC-passing fit (q is data, not estimated); this script compares them.

required_packages <- c("rstan", "dplyr", "tibble", "readr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(readr); library(ggplot2); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "02_Script/40_renewal_model/25_v4_6_fixed_q_sweep")
Q_GRID <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

load_q_bundle <- function(q) {
  tag <- paste0("q", sprintf("%.2f", q))
  readRDS(file.path(base_dir, "outputs", tag, paste0("renewal_ceara_v4_6_fit_", tag, ".rds")))
}

bundles <- lapply(Q_GRID, load_q_bundle)
names(bundles) <- sprintf("%.2f", Q_GRID)

# --- Per-q summary table -----------------------------------------------------
summary_rows <- lapply(seq_along(Q_GRID), function(i) {
  q <- Q_GRID[i]; b <- bundles[[i]]; fit <- b$fit
  sero_idx <- b$sero_week$index
  immune_prop <- rstan::extract(fit, "immune_prop")$immune_prop
  sero_pred <- rstan::extract(fit, "sero_pred")$sero_pred
  C_pred <- rstan::extract(fit, "C_pred")$C_pred
  observed <- b$weekly_data$cases
  lo <- apply(C_pred, 2, quantile, .025); hi <- apply(C_pred, 2, quantile, .975)
  tibble(
    q = q, hmc_pass = b$hmc$hmc_pass, divergences = b$hmc$divergences,
    max_rhat = b$hmc$maximum_rhat, min_bulk_ess = b$hmc$minimum_bulk_ess,
    immune_2018_median = median(immune_prop[, sero_idx]),
    immune_2018_lo = quantile(immune_prop[, sero_idx], .025), immune_2018_hi = quantile(immune_prop[, sero_idx], .975),
    sero_pred_median = median(sero_pred), sero_pred_lo = quantile(sero_pred, .025), sero_pred_hi = quantile(sero_pred, .975),
    sero_observed_covered = 103 >= quantile(sero_pred, .025) && 103 <= quantile(sero_pred, .975),
    total_case_observed = sum(observed), total_case_pred_median = median(rowSums(C_pred)),
    total_case_pred_lo = quantile(rowSums(C_pred), .025), total_case_pred_hi = quantile(rowSums(C_pred), .975),
    weekly_95pct_coverage = mean(observed >= lo & observed <= hi)
  )
})
summary_df <- bind_rows(summary_rows)
write_csv(summary_df, file.path(base_dir, "fixed_q_sweep_comparison.csv"))
message("[v4.6] summary table saved: ", file.path(base_dir, "fixed_q_sweep_comparison.csv"))
print(as.data.frame(summary_df))

theme_v4 <- theme_classic(base_size = 8.5) + theme(panel.grid.major.y = element_line(colour = "grey90"))

# --- Panel A: weekly case PPC time series, faceted by q ----------------------
weekly_ppc <- bind_rows(lapply(seq_along(Q_GRID), function(i) {
  q <- Q_GRID[i]; b <- bundles[[i]]; fit <- b$fit
  C_pred <- rstan::extract(fit, "C_pred")$C_pred
  tibble(q = q, week_start = as.Date(b$weekly_data$week_start), observed = b$weekly_data$cases,
         pred_median = apply(C_pred, 2, median), pred_lo = apply(C_pred, 2, quantile, .025), pred_hi = apply(C_pred, 2, quantile, .975))
}))
weekly_ppc$q_label <- factor(sprintf("q = %.2f", weekly_ppc$q), levels = sprintf("q = %.2f", Q_GRID))

pA <- ggplot(weekly_ppc, aes(week_start)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), fill = "steelblue", alpha = .25) +
  geom_line(aes(y = pred_median), colour = "steelblue", linewidth = .3) +
  geom_line(aes(y = observed), colour = "black", linewidth = .3) +
  facet_wrap(~q_label, ncol = 2) +
  labs(title = "A. Weekly case PPC by fixed q (black = observed, blue = predicted 95% CI)", x = NULL, y = "weekly cases") +
  theme_v4 + theme(strip.background = element_blank())

# --- Panel B: model-implied 2018 immune fraction vs q -------------------------
pB <- ggplot(summary_df, aes(q)) +
  geom_ribbon(aes(ymin = immune_2018_lo, ymax = immune_2018_hi), fill = "darkorange", alpha = .25) +
  geom_line(aes(y = immune_2018_median), colour = "darkorange") +
  geom_point(aes(y = immune_2018_median), colour = "darkorange") +
  geom_hline(yintercept = 103 / 404, linetype = 2, colour = "black") +
  annotate("text", x = max(Q_GRID), y = 103 / 404, label = "observed Juazeiro (25.5%)", hjust = 1, vjust = -0.5, size = 2.6) +
  labs(title = "B. Model-implied 2018 immune fraction vs fixed q", x = "q", y = "immune fraction") +
  theme_v4

# --- Panel C: total case-count PPC vs q ---------------------------------------
pC <- ggplot(summary_df, aes(q)) +
  geom_ribbon(aes(ymin = total_case_pred_lo, ymax = total_case_pred_hi), fill = "seagreen", alpha = .25) +
  geom_line(aes(y = total_case_pred_median), colour = "seagreen") +
  geom_point(aes(y = total_case_pred_median), colour = "seagreen") +
  geom_hline(aes(yintercept = total_case_observed), linetype = 2, colour = "black") +
  annotate("text", x = min(Q_GRID), y = summary_df$total_case_observed[1], label = "observed total cases", hjust = 0, vjust = -0.5, size = 2.6) +
  labs(title = "C. Total 2015-2019 case-count PPC vs fixed q", x = "q", y = "total cases (2015-2019)") +
  theme_v4

# --- Panel D: HMC gate summary vs q -------------------------------------------
pD <- ggplot(summary_df, aes(q, max_rhat, colour = hmc_pass)) +
  geom_hline(yintercept = 1.01, linetype = 2, colour = "grey50") +
  geom_point(size = 2) + geom_line(aes(group = 1), colour = "grey60") +
  labs(title = "D. HMC gate (max Rhat) across the q grid", x = "q", y = "max Rhat", colour = "HMC pass") +
  theme_v4 + theme(legend.position = "bottom")

figure <- (pA) / (pB | pC | pD) + plot_layout(heights = c(2, 1))
figure_dir <- file.path(root, "03_Output/figures/renewal_v4_6_fixed_q_sweep")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(figure_dir, "v4_6_fixed_q_sweep_summary.png")
ggsave(out_path, figure, width = 220, height = 260, units = "mm", dpi = 300)
message("[v4.6] figure saved: ", out_path)
