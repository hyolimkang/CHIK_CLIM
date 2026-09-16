# Mato Grosso -- fixed-q identification profile analysis. Loads the 5
# fixed-q fits from 08_fit_mt_fixed_q_profile.R and computes: case
# log-likelihood profile (total + episode-specific), case PPC, latent
# susceptibility consequences, and the q-alpha_R/R0/immune trade-off. This
# is an IDENTIFICATION DIAGNOSTIC only -- it does not select a new primary
# model.

required_packages <- c("rstan", "dplyr", "readr", "tibble", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr); library(tibble); library(ggplot2); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "02_Script/40_renewal_model/35_mato_grosso_v4_9_replication")
out_dir <- file.path(base_dir, "outputs", "fixedq")
table_dir <- file.path(root, "03_Output/tables/mato_grosso_v4_9_replication")
figure_dir <- file.path(root, "03_Output/figures/mato_grosso_v4_9_replication")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

Q_GRID <- c(0.020, 0.030, 0.040, 0.050, 0.070)
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

bundles <- lapply(Q_GRID, function(q) readRDS(file.path(out_dir, sprintf("mt_fixedq_q%.3f.rds", q))))
names(bundles) <- sprintf("q%.3f", Q_GRID)

wave_census <- read_csv(file.path(root, "03_Output/tables/national_wave_census_v2/brazil_chik_wave_census_v2.csv"), show_col_types = FALSE)
mt_waves <- wave_census |> filter(state == "MT", major_epidemic_primary == TRUE) |>
  transmute(wave_id, start_week = as.Date(onset_week, format = "%m/%d/%Y"), end_week = as.Date(end_week, format = "%m/%d/%Y")) |>
  mutate(era = if_else(start_week < as.Date("2021-01-01"), "early", "late"))
early_bounds <- range(c(mt_waves$start_week[mt_waves$era == "early"], mt_waves$end_week[mt_waves$era == "early"]))
late_bounds <- range(c(mt_waves$start_week[mt_waves$era == "late"], mt_waves$end_week[mt_waves$era == "late"]))
message(sprintf("Early episode window: %s to %s | Late episode window: %s to %s",
                 early_bounds[1], early_bounds[2], late_bounds[1], late_bounds[2]))

checkpoint_dates <- as.Date(c("2018-12-31", "2020-12-31", "2022-12-31", "2025-12-21"))

# ================================================================
# Per-q extraction: log-lik (total/early/late), PPC, latent summaries
# ================================================================
per_q <- lapply(names(bundles), function(nm) {
  b <- bundles[[nm]]
  fit <- b$fit
  weekly <- b$weekly_data
  dates <- as.Date(weekly$week_start)
  obs <- weekly$cases

  draws <- rstan::extract(fit, pars = c("logit_q", "alpha_R", "phi_obs", "expected_reported_cases", "C_pred", "S_prop", "immune_prop", "R0_t", "R_eff_t"), permuted = TRUE)
  n_draws <- nrow(draws$expected_reported_cases)

  early_idx <- which(dates >= early_bounds[1] & dates <= early_bounds[2])
  late_idx <- which(dates >= late_bounds[1] & dates <= late_bounds[2])

  loglik_per_draw <- vapply(seq_len(n_draws), function(i)
    sum(dnbinom(obs, mu = draws$expected_reported_cases[i, ], size = draws$phi_obs[i], log = TRUE)), numeric(1))
  early_loglik_per_draw <- vapply(seq_len(n_draws), function(i)
    sum(dnbinom(obs[early_idx], mu = draws$expected_reported_cases[i, early_idx], size = draws$phi_obs[i], log = TRUE)), numeric(1))
  late_loglik_per_draw <- vapply(seq_len(n_draws), function(i)
    sum(dnbinom(obs[late_idx], mu = draws$expected_reported_cases[i, late_idx], size = draws$phi_obs[i], log = TRUE)), numeric(1))

  ppc_lo <- apply(draws$C_pred, 2, quantile, .025); ppc_hi <- apply(draws$C_pred, 2, quantile, .975)
  coverage <- mean(obs >= ppc_lo & obs <= ppc_hi)
  rmse <- sqrt(mean((obs - apply(draws$C_pred, 2, median))^2))

  years <- as.integer(format(dates, "%Y"))
  annual_pred <- t(apply(draws$C_pred, 1, function(row) tapply(row, years, sum)))
  annual_obs <- tapply(obs, years, sum)
  annual_lo <- apply(annual_pred, 2, quantile, .025); annual_hi <- apply(annual_pred, 2, quantile, .975)
  annual_coverage <- mean(as.numeric(annual_obs) >= annual_lo & as.numeric(annual_obs) <= annual_hi)

  peak_early_idx <- early_idx[which.max(obs[early_idx])]
  peak_late_idx <- late_idx[which.max(obs[late_idx])]
  peak_early_covered <- obs[peak_early_idx] >= quantile(draws$C_pred[, peak_early_idx], .025) && obs[peak_early_idx] <= quantile(draws$C_pred[, peak_early_idx], .975)
  peak_late_covered <- obs[peak_late_idx] >= quantile(draws$C_pred[, peak_late_idx], .025) && obs[peak_late_idx] <= quantile(draws$C_pred[, peak_late_idx], .975)
  pred_peak_week_early <- dates[early_idx][round(median(apply(draws$C_pred[, early_idx, drop = FALSE], 1, which.max)))]
  pred_peak_week_late <- dates[late_idx][round(median(apply(draws$C_pred[, late_idx, drop = FALSE], 1, which.max)))]
  peak_timing_err_early_weeks <- as.numeric(pred_peak_week_early - dates[peak_early_idx]) / 7
  peak_timing_err_late_weeks <- as.numeric(pred_peak_week_late - dates[peak_late_idx]) / 7

  checkpoint_idx <- sapply(checkpoint_dates, function(d) which.min(abs(dates - d)))
  S_at <- setNames(apply(draws$S_prop[, checkpoint_idx, drop = FALSE], 2, median), paste0("S_", format(checkpoint_dates, "%Y")))
  U_at <- setNames(apply(draws$immune_prop[, checkpoint_idx, drop = FALSE], 2, median), paste0("immune_", format(checkpoint_dates, "%Y")))
  min_S <- apply(draws$S_prop, 1, min)
  max_R0_draw <- apply(draws$R0_t, 1, max)
  max_Reff_draw <- apply(draws$R_eff_t, 1, max)

  list(
    q_target = b$q_target, hmc = b$hmc,
    loglik_per_draw = loglik_per_draw, early_loglik_per_draw = early_loglik_per_draw, late_loglik_per_draw = late_loglik_per_draw,
    coverage = coverage, rmse = rmse, annual_coverage = annual_coverage,
    peak_early_covered = peak_early_covered, peak_late_covered = peak_late_covered,
    peak_timing_err_early_weeks = peak_timing_err_early_weeks, peak_timing_err_late_weeks = peak_timing_err_late_weeks,
    alpha_R_draws = draws$alpha_R, S_at = S_at, U_at = U_at, min_S = min_S, max_R0 = max_R0_draw, max_Reff = max_Reff_draw,
    n_weeks = length(obs)
  )
})
names(per_q) <- names(bundles)

# ================================================================
# Summary table
# ================================================================
total_loglik_mean <- sapply(per_q, function(x) mean(x$loglik_per_draw))
best_q_name <- names(which.max(total_loglik_mean))
message(sprintf("\nBest total log-lik at q=%.3f (delta=0 reference)", per_q[[best_q_name]]$q_target))

summary_tbl <- bind_rows(lapply(names(per_q), function(nm) {
  x <- per_q[[nm]]
  tibble(
    q = x$q_target,
    HMC_status = if (x$hmc$hmc_pass) "PASS" else "FAIL",
    divergences = x$hmc$divergences, max_treedepth_hits = x$hmc$max_treedepth_hits, max_rhat = x$hmc$maximum_rhat, min_bulk_ess = x$hmc$minimum_bulk_ess,
    total_loglik = mean(x$loglik_per_draw),
    delta_loglik = mean(x$loglik_per_draw) - total_loglik_mean[[best_q_name]],
    early_loglik = mean(x$early_loglik_per_draw), late_loglik = mean(x$late_loglik_per_draw),
    weekly_PPC_coverage = x$coverage, RMSE = x$rmse, annual_PPC_coverage = x$annual_coverage,
    peak_early_covered = x$peak_early_covered, peak_late_covered = x$peak_late_covered,
    peak_timing_err_early_weeks = x$peak_timing_err_early_weeks, peak_timing_err_late_weeks = x$peak_timing_err_late_weeks,
    alpha_R_median = median(x$alpha_R_draws), alpha_R_lo95 = quantile(x$alpha_R_draws, .025), alpha_R_hi95 = quantile(x$alpha_R_draws, .975),
    max_R0_median = median(x$max_R0), max_Reff_median = median(x$max_Reff),
    immune_2018 = x$U_at[["immune_2018"]], immune_2020 = x$U_at[["immune_2020"]], immune_2022 = x$U_at[["immune_2022"]], immune_2025 = x$U_at[["immune_2025"]],
    S_2018 = x$S_at[["S_2018"]], S_2020 = x$S_at[["S_2020"]], S_2022 = x$S_at[["S_2022"]], S_2025 = x$S_at[["S_2025"]],
    min_S_prop = median(x$min_S)
  )
}))
print(as.data.frame(summary_tbl), digits = 4)
write_csv(summary_tbl, file.path(table_dir, "MT_FIXED_Q_PROFILE_SUMMARY.csv"))
message("\n[saved] ", file.path(table_dir, "MT_FIXED_Q_PROFILE_SUMMARY.csv"))

# ================================================================
# Figure 1: case log-likelihood profile
# ================================================================
p_loglik <- ggplot(summary_tbl, aes(q, delta_loglik)) +
  geom_line(colour = "#315A7D") + geom_point(size = 2.5, colour = "#315A7D") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_x_continuous(breaks = Q_GRID) +
  labs(title = "MT fixed-q profile: relative case log-likelihood", subtitle = sprintf("Best q = %.3f set to Delta logLik = 0", per_q[[best_q_name]]$q_target),
       x = "Fixed q", y = "Delta total case log-likelihood") + theme_v4
ggsave(file.path(figure_dir, "MT_fixed_q_likelihood_profile.png"), p_loglik, width = 160, height = 110, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "MT_fixed_q_likelihood_profile.png"))

# ================================================================
# Figure 2: episode-specific log-likelihood
# ================================================================
episode_df <- bind_rows(
  summary_tbl |> transmute(q, era = "Early (2016-2020)", loglik = early_loglik - max(early_loglik)),
  summary_tbl |> transmute(q, era = "Late (2024-2025)", loglik = late_loglik - max(late_loglik))
)
p_episode <- ggplot(episode_df, aes(q, loglik, colour = era)) +
  geom_line() + geom_point(size = 2.5) + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("Early (2016-2020)" = "#D55E00", "Late (2024-2025)" = "#009E73")) +
  scale_x_continuous(breaks = Q_GRID) +
  labs(title = "MT fixed-q profile: episode-specific relative log-likelihood", subtitle = "Each era's best q set to Delta logLik = 0 (separately)",
       x = "Fixed q", y = "Delta episode case log-likelihood", colour = NULL) + theme_v4
ggsave(file.path(figure_dir, "MT_fixed_q_episode_loglik.png"), p_episode, width = 170, height = 110, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "MT_fixed_q_episode_loglik.png"))

# ================================================================
# Figure 3: case PPC across q (weekly, faceted)
# ================================================================
ppc_all <- bind_rows(lapply(names(bundles), function(nm) {
  b <- bundles[[nm]]; fit <- b$fit; weekly <- b$weekly_data; dates <- as.Date(weekly$week_start)
  C_pred <- rstan::extract(fit, pars = "C_pred")$C_pred
  tibble(q = b$q_target, week_start = dates, observed = weekly$cases,
         median = apply(C_pred, 2, median), lo95 = apply(C_pred, 2, quantile, .025), hi95 = apply(C_pred, 2, quantile, .975))
}))
p_ppc <- ggplot(ppc_all) +
  geom_ribbon(aes(week_start, ymin = lo95, ymax = hi95), fill = "#56B4E9", alpha = 0.25) +
  geom_line(aes(week_start, median), colour = "#315A7D", linewidth = 0.4) +
  geom_point(aes(week_start, observed), colour = "black", size = 0.25, alpha = 0.5) +
  facet_wrap(~ sprintf("q = %.3f", q), ncol = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "MT fixed-q profile: weekly case PPC across the q grid", x = NULL, y = "Weekly reported cases") + theme_v4
ggsave(file.path(figure_dir, "MT_fixed_q_case_PPC.png"), p_ppc, width = 200, height = 260, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "MT_fixed_q_case_PPC.png"))

# ================================================================
# Figure 4: susceptibility profile across q
# ================================================================
susceptibility_df <- bind_rows(lapply(names(bundles), function(nm) {
  b <- bundles[[nm]]; fit <- b$fit; weekly <- b$weekly_data; dates <- as.Date(weekly$week_start)
  S_draws <- rstan::extract(fit, pars = "S_prop")$S_prop
  tibble(q = b$q_target, week_start = dates, median = apply(S_draws, 2, median),
         lo95 = apply(S_draws, 2, quantile, .025), hi95 = apply(S_draws, 2, quantile, .975))
}))
p_susceptibility <- ggplot(susceptibility_df, aes(week_start, median, colour = factor(q), fill = factor(q))) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.08, colour = NA) +
  geom_line(linewidth = 0.6) +
  scale_colour_viridis_d(name = "Fixed q") + scale_fill_viridis_d(name = "Fixed q") +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "MT fixed-q profile: S/N susceptibility consequence of each fixed q", x = NULL, y = "S / N") + theme_v4
ggsave(file.path(figure_dir, "MT_fixed_q_susceptibility_profile.png"), p_susceptibility, width = 220, height = 120, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "MT_fixed_q_susceptibility_profile.png"))

# ================================================================
# Figure 5: q vs alpha_R / max R0 / immune_2025 trade-off
# ================================================================
p_q_alpha <- ggplot(summary_tbl, aes(q, alpha_R_median)) +
  geom_ribbon(aes(ymin = alpha_R_lo95, ymax = alpha_R_hi95), fill = "#D55E00", alpha = 0.2) +
  geom_line(colour = "#D55E00") + geom_point(size = 2.5, colour = "#D55E00") +
  scale_x_continuous(breaks = Q_GRID) +
  labs(title = "q vs alpha_R (transmission-level intercept)", x = "Fixed q", y = "alpha_R (median, 95% CrI)") + theme_v4

p_q_R0 <- ggplot(summary_tbl, aes(q, max_R0_median)) +
  geom_line(colour = "#315A7D") + geom_point(size = 2.5, colour = "#315A7D") +
  scale_x_continuous(breaks = Q_GRID) +
  labs(title = "q vs max R0(t)", x = "Fixed q", y = "max R0 (median)") + theme_v4

p_q_immune <- ggplot(summary_tbl, aes(q, immune_2025)) +
  geom_line(colour = "#76558F") + geom_point(size = 2.5, colour = "#76558F") +
  scale_x_continuous(breaks = Q_GRID) + scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(title = "q vs cumulative infection (immune_2025)", x = "Fixed q", y = "Immune fraction, end-2025") + theme_v4

figure_tradeoff <- p_q_alpha / p_q_R0 / p_q_immune +
  patchwork::plot_annotation(title = "Mato Grosso fixed-q profile: q-transmission trade-off",
                              subtitle = "Lower q -> more latent infections -> higher required transmission/depletion; higher q -> the reverse")
ggsave(file.path(figure_dir, "MT_fixed_q_q_alphaR_tradeoff.png"), figure_tradeoff, width = 160, height = 220, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "MT_fixed_q_q_alphaR_tradeoff.png"))

# Monotonicity check
message("\n=== Trade-off monotonicity check ===")
message(sprintf("q vs alpha_R: %s (Spearman rho=%.3f)",
                 if (all(diff(summary_tbl$alpha_R_median) > 0) || all(diff(summary_tbl$alpha_R_median) < 0)) "monotonic" else "NOT strictly monotonic",
                 cor(summary_tbl$q, summary_tbl$alpha_R_median, method = "spearman")))
message(sprintf("q vs immune_2025: %s (Spearman rho=%.3f)",
                 if (all(diff(summary_tbl$immune_2025) > 0) || all(diff(summary_tbl$immune_2025) < 0)) "monotonic" else "NOT strictly monotonic",
                 cor(summary_tbl$q, summary_tbl$immune_2025, method = "spearman")))

message("\n=== MT fixed-q profile analysis complete ===")
