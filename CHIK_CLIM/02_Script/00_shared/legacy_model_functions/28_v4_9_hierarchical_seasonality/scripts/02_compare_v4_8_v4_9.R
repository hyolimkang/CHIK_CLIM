# Full v4.8 vs v4.9 comparison across the shared q grid: HMC, phi_obs,
# 2018 immune fraction, susceptibility before/after the 2022 seed, and
# 2017/2022 burden + peak diagnostics (50% and 95% intervals, NOT just
# 95% coverage). Produces the comparison CSV and the required figures.

required_packages <- c("rstan", "dplyr", "tibble", "readr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(readr); library(ggplot2); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
v4_8_dir <- file.path(root, "03_Output/model_fits/ceara/v4_8_seeded_recurrence")
v4_9_dir <- file.path(root, "03_Output/model_fits/ceara/v4_9_hierarchical_seasonality")
Q_GRID <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

load_bundle <- function(dir, version, q) {
  tag <- paste0("q", sprintf("%.2f", q))
  readRDS(file.path(dir, "outputs", tag, paste0("renewal_ceara_", version, "_fit_", tag, ".rds")))
}

# --- Per-draw peak magnitude/week within a calendar year (excludes seed weeks if any) ---
peak_stats <- function(C_pred, week_dates, year, exclude_idx = integer(0)) {
  idx <- which(format(week_dates, "%Y") == as.character(year))
  idx <- setdiff(idx, exclude_idx)
  sub <- C_pred[, idx, drop = FALSE]
  per_draw_peak <- apply(sub, 1, max)
  per_draw_peak_week_idx <- idx[apply(sub, 1, which.max)]
  list(
    peak_median = median(per_draw_peak), peak_q25 = quantile(per_draw_peak, .25), peak_q75 = quantile(per_draw_peak, .75),
    peak_lo = quantile(per_draw_peak, .025), peak_hi = quantile(per_draw_peak, .975),
    peak_week_median_date = week_dates[round(median(per_draw_peak_week_idx))],
    median_trajectory_peak = max(apply(sub, 2, median)),
    median_trajectory_peak_date = week_dates[idx[which.max(apply(sub, 2, median))]]
  )
}

burden_stats <- function(C_pred, week_dates, year, exclude_idx = integer(0)) {
  idx <- which(format(week_dates, "%Y") == as.character(year))
  idx <- setdiff(idx, exclude_idx)
  totals <- rowSums(C_pred[, idx, drop = FALSE])
  list(median = median(totals), q25 = quantile(totals, .25), q75 = quantile(totals, .75),
       lo = quantile(totals, .025), hi = quantile(totals, .975))
}

summarise_one <- function(bundle, version, q) {
  fit <- bundle$fit
  w <- bundle$weekly_data; dates <- as.Date(w$week_start)
  seed_idx <- if (!is.null(bundle$seed_idx)) bundle$seed_idx else integer(0)
  sero_idx <- bundle$sero_week$index

  core_pars <- if (version == "v4_8") c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs") else
    c("alpha_R", "z_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2", "sigma_season_year", "z_season_year", "phi_obs")
  bfmi <- rstan::get_bfmi(fit)

  C_pred <- rstan::extract(fit, "C_pred")$C_pred
  immune_prop <- rstan::extract(fit, "immune_prop")$immune_prop
  S_prop <- rstan::extract(fit, "S_prop")$S_prop
  phi_obs <- rstan::extract(fit, "phi_obs")$phi_obs

  pre_seed_idx <- if (length(seed_idx)) min(seed_idx) - 1L else which(format(dates, "%Y") == "2021" & format(dates, "%m") == "12")[4]
  post_2022_idx <- max(which(format(dates, "%Y") == "2022"))

  b2017 <- burden_stats(C_pred, dates, 2017)
  b2022 <- burden_stats(C_pred, dates, 2022, exclude_idx = seed_idx)
  p2017 <- peak_stats(C_pred, dates, 2017)
  p2022 <- peak_stats(C_pred, dates, 2022, exclude_idx = seed_idx)

  obs_2017 <- sum(w$cases[format(dates, "%Y") == "2017"])
  obs_2022_idx <- setdiff(which(format(dates, "%Y") == "2022"), seed_idx)
  obs_2022 <- sum(w$cases[obs_2022_idx])
  obs_peak_2017 <- max(w$cases[format(dates, "%Y") == "2017"])
  obs_peak_2017_date <- dates[format(dates, "%Y") == "2017"][which.max(w$cases[format(dates, "%Y") == "2017"])]
  obs_peak_2022 <- max(w$cases[obs_2022_idx])
  obs_peak_2022_date <- dates[obs_2022_idx][which.max(w$cases[obs_2022_idx])]

  tibble(
    model = version, q = q,
    divergences = bundle$hmc$divergences, max_treedepth_hits = bundle$hmc$max_treedepth_hits,
    max_rhat = bundle$hmc$maximum_rhat, min_bulk_ess = bundle$hmc$minimum_bulk_ess,
    min_bfmi = min(bfmi), hmc_pass = bundle$hmc$hmc_pass,
    phi_obs_median = median(phi_obs), phi_obs_lo = quantile(phi_obs, .025), phi_obs_hi = quantile(phi_obs, .975),
    immune_2018_median = median(immune_prop[, sero_idx]),
    S_prop_pre_2022_median = median(S_prop[, pre_seed_idx]),
    S_prop_post_2022_median = median(S_prop[, post_2022_idx]),
    obs_2017 = obs_2017, pred_2017_median = b2017$median, pred_2017_q25 = b2017$q25, pred_2017_q75 = b2017$q75,
    pred_2017_lo = b2017$lo, pred_2017_hi = b2017$hi,
    obs_2022 = obs_2022, pred_2022_median = b2022$median, pred_2022_q25 = b2022$q25, pred_2022_q75 = b2022$q75,
    pred_2022_lo = b2022$lo, pred_2022_hi = b2022$hi,
    obs_peak_2017 = obs_peak_2017, obs_peak_2017_date = as.character(obs_peak_2017_date),
    pred_peak_2017_median = p2017$peak_median, pred_peak_2017_lo = p2017$peak_lo, pred_peak_2017_hi = p2017$peak_hi,
    pred_peak_2017_date = as.character(p2017$median_trajectory_peak_date),
    peak_ratio_2017 = p2017$peak_median / obs_peak_2017,
    obs_peak_2022 = obs_peak_2022, obs_peak_2022_date = as.character(obs_peak_2022_date),
    pred_peak_2022_median = p2022$peak_median, pred_peak_2022_lo = p2022$peak_lo, pred_peak_2022_hi = p2022$peak_hi,
    pred_peak_2022_date = as.character(p2022$median_trajectory_peak_date),
    peak_ratio_2022 = p2022$peak_median / obs_peak_2022
  )
}

rows <- list()
for (q in Q_GRID) {
  rows[[length(rows) + 1]] <- summarise_one(load_bundle(v4_8_dir, "v4_8", q), "v4_8", q)
  rows[[length(rows) + 1]] <- summarise_one(load_bundle(v4_9_dir, "v4_9", q), "v4_9", q)
}
comparison <- bind_rows(rows)
write_csv(comparison, file.path(v4_9_dir, "v4_8_vs_v4_9_comparison.csv"))
message("[compare] saved: ", file.path(v4_9_dir, "v4_8_vs_v4_9_comparison.csv"))
print(as.data.frame(comparison))

theme_v4 <- theme_classic(base_size = 8.5) + theme(panel.grid.major.y = element_line(colour = "grey90"))
comparison$model_label <- ifelse(comparison$model == "v4_8", "v4.8 (single first harmonic)", "v4.9 (hierarchical seasonality)")

# --- Figure 2/3/4/7: burden + peak + phi_obs comparison dashboard ------------
p_burden_2017 <- ggplot(comparison, aes(q, pred_2017_median, colour = model_label)) +
  geom_ribbon(aes(ymin = pred_2017_lo, ymax = pred_2017_hi, fill = model_label), alpha = .15, colour = NA) +
  geom_ribbon(aes(ymin = pred_2017_q25, ymax = pred_2017_q75, fill = model_label), alpha = .3, colour = NA) +
  geom_line() + geom_point() +
  geom_hline(aes(yintercept = obs_2017), linetype = 2, colour = "black") +
  labs(title = "2017 total-case PPC: v4.8 vs v4.9", x = "q", y = "2017 total cases") + theme_v4 + theme(legend.position = "bottom")

p_burden_2022 <- ggplot(comparison, aes(q, pred_2022_median, colour = model_label)) +
  geom_ribbon(aes(ymin = pred_2022_lo, ymax = pred_2022_hi, fill = model_label), alpha = .15, colour = NA) +
  geom_ribbon(aes(ymin = pred_2022_q25, ymax = pred_2022_q75, fill = model_label), alpha = .3, colour = NA) +
  geom_line() + geom_point() +
  geom_hline(aes(yintercept = obs_2022), linetype = 2, colour = "black") +
  labs(title = "2022 total-case PPC: v4.8 vs v4.9 (post-seed)", x = "q", y = "2022 total cases") + theme_v4 + theme(legend.position = "bottom")

p_peak_2017 <- ggplot(comparison, aes(q, pred_peak_2017_median, colour = model_label)) +
  geom_ribbon(aes(ymin = pred_peak_2017_lo, ymax = pred_peak_2017_hi, fill = model_label), alpha = .15, colour = NA) +
  geom_line() + geom_point() +
  geom_hline(aes(yintercept = obs_peak_2017), linetype = 2, colour = "black") +
  labs(title = "2017 weekly peak PPC: v4.8 vs v4.9", x = "q", y = "peak weekly cases") + theme_v4 + theme(legend.position = "bottom")

p_peak_2022 <- ggplot(comparison, aes(q, pred_peak_2022_median, colour = model_label)) +
  geom_ribbon(aes(ymin = pred_peak_2022_lo, ymax = pred_peak_2022_hi, fill = model_label), alpha = .15, colour = NA) +
  geom_line() + geom_point() +
  geom_hline(aes(yintercept = obs_peak_2022), linetype = 2, colour = "black") +
  labs(title = "2022 weekly peak PPC: v4.8 vs v4.9", x = "q", y = "peak weekly cases") + theme_v4 + theme(legend.position = "bottom")

p_phi <- ggplot(comparison, aes(q, phi_obs_median, colour = model_label)) +
  geom_pointrange(aes(ymin = phi_obs_lo, ymax = phi_obs_hi), position = position_dodge(width = 0.01)) +
  labs(title = "phi_obs posterior by q (diagnostic only; prior unchanged)", x = "q", y = "phi_obs") +
  theme_v4 + theme(legend.position = "bottom")

# --- Susceptibility + R0/Reff trajectories at representative q=0.10 ----------
b9_ref <- load_bundle(v4_9_dir, "v4_9", 0.10)
w_ref <- b9_ref$weekly_data; dates_ref <- as.Date(w_ref$week_start)
S_prop_ref <- rstan::extract(b9_ref$fit, "S_prop")$S_prop
R0_ref <- rstan::extract(b9_ref$fit, "R0_t")$R0_t
Reff_ref <- rstan::extract(b9_ref$fit, "R_eff_t")$R_eff_t
traj_df <- tibble(
  week_start = dates_ref,
  S_median = apply(S_prop_ref, 2, median), S_lo = apply(S_prop_ref, 2, quantile, .025), S_hi = apply(S_prop_ref, 2, quantile, .975),
  R0_median = apply(R0_ref, 2, median), Reff_median = apply(Reff_ref, 2, median)
)
p_S <- ggplot(traj_df, aes(week_start)) +
  geom_ribbon(aes(ymin = S_lo, ymax = S_hi), fill = "seagreen", alpha = .25) +
  geom_line(aes(y = S_median), colour = "seagreen") +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
  labs(title = "v4.9 susceptible proportion (q=0.10, representative)", x = NULL, y = "S proportion") + theme_v4

p_R <- ggplot(traj_df, aes(week_start)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_line(aes(y = R0_median, colour = "R0(t)")) +
  geom_line(aes(y = Reff_median, colour = "Reff(t)")) +
  scale_colour_manual(values = c("R0(t)" = "#315A7D", "Reff(t)" = "#D55E00")) +
  labs(title = "v4.9 R0(t)/Reff(t) posterior median (q=0.10, representative)", x = NULL, y = "reproduction number", colour = NULL) +
  theme_v4 + theme(legend.position = "bottom")

figure2 <- (p_burden_2017 | p_burden_2022) / (p_peak_2017 | p_peak_2022) / (p_phi | p_S | p_R)
figure_dir <- file.path(root, "03_Output/figures/renewal_v4_9_hierarchical_seasonality")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(figure_dir, "v4_8_vs_v4_9_diagnostics.png"), figure2, width = 260, height = 300, units = "mm", dpi = 300)
message("[compare] figure saved: ", file.path(figure_dir, "v4_8_vs_v4_9_diagnostics.png"))

# --- Figure 1: full 2015-2025 weekly PPC by q, median clearly visible -------
weekly_ppc <- bind_rows(lapply(Q_GRID, function(q) {
  b <- load_bundle(v4_9_dir, "v4_9", q)
  C_pred <- rstan::extract(b$fit, "C_pred")$C_pred
  tibble(q = q, week_start = as.Date(b$weekly_data$week_start), observed = b$weekly_data$cases,
         pred_median = apply(C_pred, 2, median), pred_q25 = apply(C_pred, 2, quantile, .25), pred_q75 = apply(C_pred, 2, quantile, .75),
         pred_lo = apply(C_pred, 2, quantile, .025), pred_hi = apply(C_pred, 2, quantile, .975))
}))
weekly_ppc$q_label <- factor(sprintf("q = %.2f", weekly_ppc$q), levels = sprintf("q = %.2f", Q_GRID))
seed_dates <- b9_ref$seed_dates

p_full_ppc <- ggplot(weekly_ppc, aes(week_start)) +
  annotate("rect", xmin = min(seed_dates), xmax = max(seed_dates) + 7, ymin = -Inf, ymax = Inf, fill = "grey40", alpha = .18) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), fill = "steelblue", alpha = .2) +
  geom_ribbon(aes(ymin = pred_q25, ymax = pred_q75), fill = "steelblue", alpha = .35) +
  geom_line(aes(y = pred_median), colour = "#08306b", linewidth = .45) +
  geom_line(aes(y = observed), colour = "black", linewidth = .3) +
  facet_wrap(~q_label, ncol = 2) +
  labs(title = "v4.9 full 2015-2025 weekly case PPC by fixed q",
       subtitle = "Black = observed; dark blue = posterior median; bands = 50%/95% CrI; grey = 2022 seed (not scored)",
       x = NULL, y = "weekly cases") +
  theme_v4 + theme(strip.background = element_blank())
ggsave(file.path(figure_dir, "v4_9_full_period_ppc_by_q.png"), p_full_ppc, width = 220, height = 260, units = "mm", dpi = 300)
message("[compare] figure saved: ", file.path(figure_dir, "v4_9_full_period_ppc_by_q.png"))
