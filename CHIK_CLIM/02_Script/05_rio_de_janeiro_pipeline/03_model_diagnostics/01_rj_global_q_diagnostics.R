# Rio de Janeiro Phase 1 -- Sections 10-13: HMC gate, q/S identification
# diagnostics, case PPC, latent trajectory figure.

required_packages <- c("rstan", "dplyr", "readr", "tibble", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr); library(tibble); library(ggplot2); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "03_Output/model_fits/rio_de_janeiro/v4_9_replication")
table_dir <- file.path(root, "03_Output/tables/rio_de_janeiro_v4_9_replication")
figure_dir <- file.path(root, "03_Output/figures/rio_de_janeiro_v4_9_replication")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

b <- readRDS(file.path(base_dir, "outputs/caseonly/rj_global_q_case_only.rds"))
fit <- b$fit
weekly <- b$weekly_data
dates <- as.Date(weekly$week_start)

message("=== Section 10: standard HMC gate ===")
sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
hmc_by_chain <- bind_rows(lapply(seq_along(sp), function(i) {
  tibble(chain = i, divergent = sum(sp[[i]][, "divergent__"]),
         max_treedepth_hits = sum(sp[[i]][, "treedepth__"] >= b$config$max_treedepth),
         mean_stepsize = mean(sp[[i]][, "stepsize__"]), mean_leapfrog = mean(sp[[i]][, "n_leapfrog__"]),
         mean_accept = mean(sp[[i]][, "accept_stat__"]))
}))
hmc_by_chain$bfmi <- rstan::get_bfmi(fit)
print(as.data.frame(hmc_by_chain), digits = 4)
mon <- rstan::monitor(fit, print = FALSE)
min_tail_ess <- min(mon[, "Tail_ESS"], na.rm = TRUE)
message(sprintf("\nOverall: divergences=%d, max_treedepth_hits=%d, max Rhat=%.4f, min bulk ESS=%.1f, min tail ESS=%.1f, hmc_pass=%s",
                 b$hmc$divergences, b$hmc$max_treedepth_hits, b$hmc$maximum_rhat, b$hmc$minimum_bulk_ess, min_tail_ess, b$hmc$hmc_pass))

hmc_classification <- if (b$hmc$hmc_pass) "HMC PASS" else if (b$hmc$divergences <= 5 && b$hmc$max_treedepth_hits == 0 && b$hmc$maximum_rhat <= 1.05) "HMC NEAR-PASS" else "HMC FAIL"
message("Classification: ", hmc_classification)
write_csv(hmc_by_chain, file.path(table_dir, "RJ_global_q_HMC.csv"))
message("[saved] ", file.path(table_dir, "RJ_global_q_HMC.csv"))

message("\n=== Section 11: q / S identification diagnostics ===")
draws <- rstan::extract(fit, pars = c("logit_q", "alpha_R", "S_prop", "immune_prop", "R0_t", "expected_reported_cases", "C_pred"), permuted = TRUE)
q_draws <- plogis(draws$logit_q)
q_prior_draws <- plogis(as.vector(rstan::extract(fit, "logit_q_prior_draw")$logit_q_prior_draw))

cor_q_alpha <- cor(draws$logit_q, draws$alpha_R)
immune_2025 <- draws$immune_prop[, ncol(draws$immune_prop)]
min_S <- apply(draws$S_prop, 1, min)
max_R0 <- apply(draws$R0_t, 1, max)
cor_q_immune <- cor(q_draws, immune_2025)
cor_q_minS <- cor(q_draws, min_S)
cor_q_maxR0 <- cor(q_draws, max_R0)

q_summary <- tibble(
  q_prior_median = median(q_prior_draws), q_prior_lo95 = quantile(q_prior_draws, .025), q_prior_hi95 = quantile(q_prior_draws, .975),
  q_posterior_median = median(q_draws), q_posterior_lo95 = quantile(q_draws, .025), q_posterior_hi95 = quantile(q_draws, .975),
  cor_logitq_alphaR = cor_q_alpha, cor_q_immune2025 = cor_q_immune, cor_q_minS = cor_q_minS, cor_q_maxR0 = cor_q_maxR0
)
print(as.data.frame(q_summary), digits = 4)
write_csv(q_summary, file.path(table_dir, "RJ_global_q_posterior_summary.csv"))
message("[saved] ", file.path(table_dir, "RJ_global_q_posterior_summary.csv"))

prior_width <- diff(quantile(q_prior_draws, c(.025, .975)))
posterior_width <- diff(quantile(q_draws, c(.025, .975)))
message(sprintf("\nPrior 95%% width=%.4f vs posterior 95%% width=%.4f (ratio=%.3f) -- %s",
                 prior_width, posterior_width, posterior_width / prior_width,
                 if (posterior_width / prior_width < 0.5) "posterior meaningfully narrower than prior" else "posterior NOT much narrower than prior -- weak identification signal"))

# ---- Checkpoints ----
checkpoint_dates <- as.Date(c("2016-12-31", "2018-12-31", "2020-12-31", "2022-12-31", "2025-12-21"))
checkpoint_idx <- sapply(checkpoint_dates, function(d) which.min(abs(dates - d)))
checkpoint_tbl <- bind_rows(lapply(seq_along(checkpoint_dates), function(i) {
  t_idx <- checkpoint_idx[i]
  tibble(checkpoint = checkpoint_dates[i],
         S_prop_median = median(draws$S_prop[, t_idx]), S_prop_lo95 = quantile(draws$S_prop[, t_idx], .025), S_prop_hi95 = quantile(draws$S_prop[, t_idx], .975),
         immune_prop_median = median(draws$immune_prop[, t_idx]), immune_prop_lo95 = quantile(draws$immune_prop[, t_idx], .025), immune_prop_hi95 = quantile(draws$immune_prop[, t_idx], .975),
         R0_median = median(draws$R0_t[, t_idx]))
}))
message("\n=== Checkpoints: S/N, cumulative immune fraction, R0 ===")
print(as.data.frame(checkpoint_tbl), digits = 3)
write_csv(checkpoint_tbl, file.path(table_dir, "RJ_global_q_checkpoints.csv"))
message("[saved] ", file.path(table_dir, "RJ_global_q_checkpoints.csv"))

# ================================================================
# Section 12: case PPC (weekly + annual)
# ================================================================
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
summarise_traj <- function(m) tibble(median = apply(m, 2, median), lo95 = apply(m, 2, quantile, .025), hi95 = apply(m, 2, quantile, .975))

pred_summary <- summarise_traj(draws$C_pred); pred_summary$week_start <- dates
exp_summary <- summarise_traj(draws$expected_reported_cases); exp_summary$week_start <- dates
obs_df <- tibble(week_start = dates, observed = weekly$cases)

p_weekly_ppc <- ggplot() +
  geom_ribbon(data = pred_summary, aes(week_start, ymin = lo95, ymax = hi95), fill = "#56B4E9", alpha = 0.2) +
  geom_line(data = exp_summary, aes(week_start, median), colour = "#315A7D", linewidth = 0.5) +
  geom_point(data = obs_df, aes(week_start, observed), colour = "black", size = 0.4, alpha = 0.6) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "RJ global-q case-only model: weekly case PPC",
       subtitle = "Black points = observed; blue ribbon = 95% posterior predictive; navy line = model expectation",
       x = NULL, y = "Weekly reported cases") + theme_v4
ggsave(file.path(figure_dir, "RJ_global_q_case_PPC.png"), p_weekly_ppc, width = 220, height = 110, units = "mm", dpi = 300, bg = "white")
message("\n[saved] ", file.path(figure_dir, "RJ_global_q_case_PPC.png"))

years <- as.integer(format(dates, "%Y"))
annual_pred <- t(apply(draws$C_pred, 1, function(row) tapply(row, years, sum)))
annual_pred_summary <- tibble(year = as.integer(colnames(annual_pred)), median = apply(annual_pred, 2, median),
                               lo95 = apply(annual_pred, 2, quantile, .025), hi95 = apply(annual_pred, 2, quantile, .975))
annual_obs <- tapply(weekly$cases, years, sum)
annual_pred_summary$observed <- as.numeric(annual_obs[as.character(annual_pred_summary$year)])

p_annual_ppc <- ggplot(annual_pred_summary, aes(year)) +
  geom_pointrange(aes(y = median, ymin = lo95, ymax = hi95), colour = "#315A7D") +
  geom_point(aes(y = observed), colour = "#D55E00", size = 2.5, shape = 18) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(title = "RJ global-q case-only model: annual case totals, observed vs predicted",
       subtitle = "Orange diamond = observed; blue point+range = posterior predictive median/95% CrI",
       x = NULL, y = "Annual total cases") + theme_v4
ggsave(file.path(figure_dir, "RJ_global_q_annual_PPC.png"), p_annual_ppc, width = 180, height = 110, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "RJ_global_q_annual_PPC.png"))
write_csv(annual_pred_summary, file.path(table_dir, "RJ_global_q_annual_PPC.csv"))

# Peak timing/magnitude check per major wave
wave_census <- read_csv(file.path(root, "03_Output/tables/national_wave_census_v2/brazil_chik_wave_census_v2.csv"), show_col_types = FALSE)
rj_waves <- wave_census |> dplyr::filter(state == "RJ", major_epidemic_primary == TRUE) |>
  transmute(wave_id, start_week = as.Date(onset_week, format = "%m/%d/%Y"), end_week = as.Date(end_week, format = "%m/%d/%Y"))
peak_check <- bind_rows(lapply(seq_len(nrow(rj_waves)), function(i) {
  w <- rj_waves[i, ]
  idx <- which(dates >= w$start_week & dates <= w$end_week)
  obs_peak_idx <- idx[which.max(weekly$cases[idx])]
  pred_peak_idx_per_draw <- idx[apply(draws$C_pred[, idx, drop = FALSE], 1, which.max)]
  tibble(wave_id = w$wave_id, observed_peak_week = dates[obs_peak_idx], observed_peak_cases = weekly$cases[obs_peak_idx],
         observed_total = sum(weekly$cases[idx]), predicted_total_median = median(rowSums(draws$C_pred[, idx, drop = FALSE])),
         predicted_peak_week_median = dates[round(median(pred_peak_idx_per_draw))])
}))
message("\n=== Per-wave peak timing/magnitude check ===")
print(as.data.frame(peak_check), digits = 3)
write_csv(peak_check, file.path(table_dir, "RJ_global_q_wave_peak_check.csv"))

# ================================================================
# Section 13: latent trajectory figure
# ================================================================
S_summary <- summarise_traj(draws$S_prop); S_summary$week_start <- dates
U_summary <- summarise_traj(draws$immune_prop); U_summary$week_start <- dates
R0_summary <- summarise_traj(draws$R0_t); R0_summary$week_start <- dates

p_cases <- ggplot() +
  geom_ribbon(data = exp_summary, aes(week_start, ymin = lo95, ymax = hi95), fill = "#315A7D", alpha = 0.2) +
  geom_line(data = exp_summary, aes(week_start, median), colour = "#315A7D", linewidth = 0.5) +
  geom_point(data = obs_df, aes(week_start, observed), colour = "black", size = 0.4, alpha = 0.6) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Observed + predicted cases", x = NULL, y = "Weekly cases") + theme_v4

p_S <- ggplot(S_summary, aes(week_start, median)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#009E73", alpha = 0.2) +
  geom_line(colour = "#009E73", linewidth = 0.6) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
  labs(title = "S/N (case-only model-implied)", x = NULL, y = "S / N") + theme_v4

p_U <- ggplot(U_summary, aes(week_start, median)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#76558F", alpha = 0.2) +
  geom_line(colour = "#76558F", linewidth = 0.6) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(title = "Cumulative immune fraction U/N (case-only model-implied)", x = NULL, y = "U / N") + theme_v4

p_R0 <- ggplot(R0_summary, aes(week_start, median)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#D55E00", alpha = 0.15) +
  geom_line(colour = "#D55E00", linewidth = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  labs(title = "R0(t)", x = NULL, y = "R0") + theme_v4

q_band <- tibble(week_start = dates, lo95 = quantile(q_draws, .025), hi95 = quantile(q_draws, .975), median = median(q_draws))
p_q <- ggplot(q_band, aes(week_start)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#E69F00", alpha = 0.25) +
  geom_line(aes(y = median), colour = "#E69F00", linewidth = 0.9) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA)) +
  labs(title = sprintf("q_RJ posterior: median=%.4f, 95%% CrI=[%.4f, %.4f]", median(q_draws), quantile(q_draws,.025), quantile(q_draws,.975)),
       x = NULL, y = "q") + theme_v4

figure <- (p_cases | p_S) / (p_U | p_R0) / p_q +
  patchwork::plot_annotation(title = "Rio de Janeiro: case-only model-implied latent trajectories (Phase 1, NOT validated absolute susceptibility)",
                              subtitle = sprintf("Global-q, frozen v4.9, %s", hmc_classification))
ggsave(file.path(figure_dir, "RJ_global_q_trajectories.png"), figure, width = 210, height = 240, units = "mm", dpi = 300, bg = "white")
message("\n[saved] ", file.path(figure_dir, "RJ_global_q_trajectories.png"))
