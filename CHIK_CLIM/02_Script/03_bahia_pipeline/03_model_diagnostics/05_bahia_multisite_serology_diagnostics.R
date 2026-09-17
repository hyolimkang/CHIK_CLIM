# Required diagnostics for the Bahia multi-site geo-adjusted serology
# model: per-survey p_state_window/p_site_window/eta_geo posterior
# summaries, HMC diagnostics, and headline case/wave PPC, for the q=0.05
# and q=0.10 pilots (Section 10 of the design spec).

required_packages <- c("rstan", "dplyr", "tibble", "readr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(readr); library(ggplot2) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
base_dir <- file.path(root, "03_Output/model_fits/bahia/v4_9_replication")

load_multisite <- function(q, stage = "pilot") {
  tag <- paste0("q", sprintf("%.2f", q))
  readRDS(file.path(base_dir, "outputs_multisite_serology", tag, paste0("renewal_bahia_multisite_fit_", tag, ".rds")))
}

summarise_one_q <- function(q) {
  b <- load_multisite(q)
  fit <- b$fit
  audit <- b$sero_audit

  p_state_window <- rstan::extract(fit, "p_state_window")$p_state_window
  p_site_window <- rstan::extract(fit, "p_site_window")$p_site_window
  eta_geo <- rstan::extract(fit, "eta_geo")$eta_geo

  sero_table <- bind_rows(lapply(seq_len(nrow(audit)), function(j) {
    tibble(
      q = q, sero_id = audit$sero_id[j], location = audit$location[j],
      n_positive = audit$n_positive[j], n_tested = audit$n_tested[j],
      observed_prevalence = audit$observed_prevalence[j],
      p_state_window_median = median(p_state_window[, j]),
      p_state_window_lo = quantile(p_state_window[, j], .025), p_state_window_hi = quantile(p_state_window[, j], .975),
      p_site_window_median = median(p_site_window[, j]),
      p_site_window_lo = quantile(p_site_window[, j], .025), p_site_window_hi = quantile(p_site_window[, j], .975),
      eta_geo_median = median(eta_geo[, j]),
      eta_geo_lo = quantile(eta_geo[, j], .025), eta_geo_hi = quantile(eta_geo[, j], .975)
    )
  }))

  w <- b$weekly_data; dates <- as.Date(w$week_start)
  C_pred <- rstan::extract(fit, "C_pred")$C_pred
  idx_2017 <- which(format(dates, "%Y") == "2017")
  case_2017 <- list(obs_total = sum(w$cases[idx_2017]), pred_total_median = median(rowSums(C_pred[, idx_2017])),
                     obs_peak = max(w$cases[idx_2017]), pred_peak_median = median(apply(C_pred[, idx_2017], 1, max)))

  S_prop <- rstan::extract(fit, "S_prop")$S_prop
  immune_prop <- rstan::extract(fit, "immune_prop")$immune_prop
  idx_2018_end <- max(which(format(dates, "%Y") == "2018"))
  idx_2025_end <- nrow(w)
  immune_summary <- tibble(
    q = q,
    immune_2018_median = median(immune_prop[, idx_2018_end]), immune_2018_lo = quantile(immune_prop[, idx_2018_end], .025), immune_2018_hi = quantile(immune_prop[, idx_2018_end], .975),
    immune_2025_median = median(immune_prop[, idx_2025_end]), immune_2025_lo = quantile(immune_prop[, idx_2025_end], .025), immune_2025_hi = quantile(immune_prop[, idx_2025_end], .975)
  )

  list(sero_table = sero_table, case_2017 = case_2017, immune_summary = immune_summary, hmc = b$hmc)
}

Q_PILOT <- c(0.05, 0.10)
results <- lapply(Q_PILOT, summarise_one_q)

sero_all <- bind_rows(lapply(results, `[[`, "sero_table"))
immune_all <- bind_rows(lapply(results, `[[`, "immune_summary"))
hmc_all <- bind_rows(lapply(seq_along(Q_PILOT), function(i) {
  h <- results[[i]]$hmc
  tibble(q = Q_PILOT[i], divergences = h$divergences, max_treedepth_hits = h$max_treedepth_hits,
         max_rhat = h$maximum_rhat, min_bulk_ess = h$minimum_bulk_ess, bfmi_min = min(h$bfmi), hmc_pass = h$hmc_pass)
}))
case_all <- bind_rows(lapply(seq_along(Q_PILOT), function(i) {
  c2017 <- results[[i]]$case_2017
  tibble(q = Q_PILOT[i], obs_2017_total = c2017$obs_total, pred_2017_total_median = c2017$pred_total_median,
         obs_2017_peak = c2017$obs_peak, pred_2017_peak_median = c2017$pred_peak_median)
}))

message("=== HMC diagnostics (pilots) ==="); print(as.data.frame(hmc_all))
message("\n=== 2017 case total/peak ==="); print(as.data.frame(case_all))
message("\n=== 2018/2025 immune fraction ==="); print(as.data.frame(immune_all))
message("\n=== 6-survey serology comparison ==="); print(as.data.frame(sero_table <- sero_all))

table_dir <- file.path(root, "03_Output/tables/renewal_bahia_v4_9_multisite_serology")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(sero_all, file.path(table_dir, "bahia_multisite_serology_ppc_by_q.csv"))
write_csv(hmc_all, file.path(table_dir, "bahia_multisite_hmc_diagnostics_by_q.csv"))
write_csv(case_all, file.path(table_dir, "bahia_multisite_2017_case_ppc_by_q.csv"))
write_csv(immune_all, file.path(table_dir, "bahia_multisite_immune_fraction_by_q.csv"))
message("\n[saved] tables under: ", table_dir)

theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
sero_plot_df <- bind_rows(
  sero_all |> transmute(q, sero_id, quantity = "State-level (no geo offset)", median = p_state_window_median, lo = p_state_window_lo, hi = p_state_window_hi),
  sero_all |> transmute(q, sero_id, quantity = "Site-adjusted (with geo offset)", median = p_site_window_median, lo = p_site_window_lo, hi = p_site_window_hi),
  sero_all |> distinct(sero_id, observed_prevalence) |> transmute(q = NA_real_, sero_id, quantity = "Observed", median = observed_prevalence, lo = NA_real_, hi = NA_real_)
)
p_compare <- ggplot(sero_plot_df |> dplyr::filter(!is.na(q)), aes(factor(q), median, colour = quantity)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = .4), fatten = 1.8) +
  geom_hline(data = sero_plot_df |> dplyr::filter(is.na(q)), aes(yintercept = median), linetype = 2, colour = "black") +
  facet_wrap(~sero_id, nrow = 2) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_colour_manual(values = c("State-level (no geo offset)" = "#009E73", "Site-adjusted (with geo offset)" = "#76558F")) +
  labs(title = "Bahia 6-site serology: state vs geo-adjusted prediction by q (dashed = observed)", x = "q", y = "seroprevalence") +
  theme_v4 + theme(legend.position = "bottom", strip.background = element_blank())

figure_dir <- file.path(root, "03_Output/figures/renewal_bahia_v4_9_multisite_serology")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(figure_dir, "bahia_multisite_serology_comparison_by_q.png"), p_compare, width = 220, height = 160, units = "mm", dpi = 300)
message("[saved] ", file.path(figure_dir, "bahia_multisite_serology_comparison_by_q.png"))
