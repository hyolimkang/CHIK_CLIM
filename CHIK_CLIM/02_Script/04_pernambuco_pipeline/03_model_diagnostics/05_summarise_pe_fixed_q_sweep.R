# Pernambuco fixed-q sweep -- consolidate HMC, immune fraction, and U14
# serology fit across q = 0.05-0.30 (U14 serology ON throughout).

required_packages <- c("rstan", "dplyr", "readr", "tibble", "ggplot2", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr); library(tibble); library(ggplot2); library(tidyr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/model_fits/v4_9_replication")
table_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/tables/pernambuco_v4_9_replication/fixed_q_sweep")
figure_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/figures/pernambuco_v4_9_replication/fixed_q_sweep")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

q_grid <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

summarise_one <- function(qval) {
  tag <- paste0("q", sprintf("%.2f", qval))
  b <- readRDS(file.path(base_dir, "outputs/fixed_q_sweep", tag, paste0("renewal_pe_fixedq_fit_", tag, ".rds")))
  fit <- b$fit; w <- b$weekly_data; dates <- as.Date(w$week_start)

  immune_prop <- rstan::extract(fit, "immune_prop")$immune_prop
  idx_2018 <- max(which(format(dates, "%Y") == "2018")); idx_2025 <- nrow(w)
  p_state_window <- rstan::extract(fit, "p_state_window")$p_state_window[, 1]
  p_site_window <- rstan::extract(fit, "p_site_window")$p_site_window[, 1]
  eta_geo <- rstan::extract(fit, "eta_geo")$eta_geo[, 1]

  C_pred <- rstan::extract(fit, "C_pred")$C_pred
  idx_2021 <- which(format(dates, "%Y") == "2021") # PE's largest wave (wave_07)

  tibble(
    q = qval, divergences = b$hmc$divergences, max_treedepth_hits = b$hmc$max_treedepth_hits,
    max_rhat = b$hmc$maximum_rhat, min_bulk_ess = b$hmc$minimum_bulk_ess, bfmi_min = min(b$hmc$bfmi), hmc_pass = b$hmc$hmc_pass,
    immune_2018_median = median(immune_prop[, idx_2018]), immune_2018_lo95 = quantile(immune_prop[, idx_2018], .025), immune_2018_hi95 = quantile(immune_prop[, idx_2018], .975),
    immune_2025_median = median(immune_prop[, idx_2025]), immune_2025_lo95 = quantile(immune_prop[, idx_2025], .025), immune_2025_hi95 = quantile(immune_prop[, idx_2025], .975),
    pe_state_window_median = median(p_state_window), geo_adjusted_recife_median = median(p_site_window),
    eta_geo_median = median(eta_geo), eta_geo_lo95 = quantile(eta_geo, .025), eta_geo_hi95 = quantile(eta_geo, .975),
    obs_2021_total = sum(w$cases[idx_2021]), pred_2021_total_median = median(rowSums(C_pred[, idx_2021]))
  )
}

sweep_summary <- bind_rows(lapply(q_grid, summarise_one))
message("=== Pernambuco fixed-q sweep summary (U14 serology ON throughout) ===")
print(as.data.frame(sweep_summary), digits = 4)
write_csv(sweep_summary, file.path(table_dir, "pernambuco_fixed_q_sweep_summary.csv"))
message("\n[saved] ", file.path(table_dir, "pernambuco_fixed_q_sweep_summary.csv"))

# ---- Figure: immune fraction 2025 + U14 fit + HMC status across q --------
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

p1 <- ggplot(sweep_summary, aes(factor(q), immune_2025_median, colour = hmc_pass)) +
  geom_pointrange(aes(ymin = immune_2025_lo95, ymax = immune_2025_hi95)) +
  scale_colour_manual(values = c(`TRUE` = "#009E73", `FALSE` = "#D55E00"), labels = c(`TRUE` = "HMC PASS", `FALSE` = "HMC FAIL"), name = NULL) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(title = "Pernambuco fixed-q sweep: 2025 cumulative immune fraction by q", x = "q (fixed)", y = "Immune fraction (2025)") +
  theme_v4 + theme(legend.position = "bottom")

p2 <- ggplot(sweep_summary |> pivot_longer(c(pe_state_window_median, geo_adjusted_recife_median), names_to = "quantity", values_to = "value"),
             aes(factor(q), value, colour = quantity)) +
  geom_point(size = 2.5) + geom_line(aes(group = quantity)) +
  geom_hline(yintercept = 0.372, linetype = 2, colour = "black") +
  annotate("text", x = 1, y = 0.372, label = "Observed U14 = 37.2%", vjust = -0.5, hjust = 0, size = 2.8) +
  scale_colour_manual(values = c(pe_state_window_median = "#009E73", geo_adjusted_recife_median = "#76558F"),
                       labels = c(pe_state_window_median = "PE state-level (no offset)", geo_adjusted_recife_median = "Geo-adjusted Recife"), name = NULL) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(title = "U14 fit across the fixed-q sweep", x = "q (fixed)", y = "seroprevalence") +
  theme_v4 + theme(legend.position = "bottom")

library(patchwork)
combined <- p1 / p2
ggsave(file.path(figure_dir, "pernambuco_fixed_q_sweep_comparison.png"), combined, width = 190, height = 220, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "pernambuco_fixed_q_sweep_comparison.png"))
