# v4.8 seeded-recurrence sweep, FULL 2015-2025 period: summary table +
# figures across the same q grid as v4.7, reproducing the v4.7 summary
# figure format (panels A-E) and adding the required v4.7-vs-v4.8
# comparison (F) and pre-2022 susceptible fraction (G). Seed-window weeks
# are excluded from all "predictive" PPC summaries and clearly marked in
# the weekly trajectory panel, per the v4.8 design spec.

required_packages <- c("rstan", "dplyr", "tibble", "readr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(readr); library(ggplot2); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
v4_7_dir <- file.path(root, "02_Script/40_renewal_model/26_v4_7_fixed_q_full_period")
v4_8_dir <- file.path(root, "02_Script/40_renewal_model/27_v4_8_seeded_recurrence")
Q_GRID <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

load_bundle <- function(dir, version, q) {
  tag <- paste0("q", sprintf("%.2f", q))
  readRDS(file.path(dir, "outputs", tag, paste0("renewal_ceara_", version, "_fit_", tag, ".rds")))
}
bundles_v8 <- lapply(Q_GRID, load_bundle, dir = v4_8_dir, version = "v4_8")
names(bundles_v8) <- sprintf("%.2f", Q_GRID)

# --- Per-q summary table (v4.8) ---------------------------------------------
summary_rows <- lapply(seq_along(Q_GRID), function(i) {
  q <- Q_GRID[i]; b <- bundles_v8[[i]]; fit <- b$fit
  sero_idx <- b$sero_week$index
  seed_idx <- b$seed_idx
  post_seed_2022_idx <- which(format(as.Date(b$weekly_data$week_start), "%Y") == "2022") |> setdiff(seed_idx)
  all_2022_idx <- which(format(as.Date(b$weekly_data$week_start), "%Y") == "2022")
  idx_2017 <- which(format(as.Date(b$weekly_data$week_start), "%Y") == "2017")
  pre_seed_idx <- min(seed_idx) - 1L

  immune_prop <- rstan::extract(fit, "immune_prop")$immune_prop
  S_prop <- rstan::extract(fit, "S_prop")$S_prop
  sero_pred <- rstan::extract(fit, "sero_pred")$sero_pred
  C_pred <- rstan::extract(fit, "C_pred")$C_pred
  observed <- b$weekly_data$cases

  pred_2017 <- rowSums(C_pred[, idx_2017])
  pred_2022_all <- rowSums(C_pred[, all_2022_idx])
  pred_2022_postseed <- rowSums(C_pred[, post_seed_2022_idx])
  obs_2022_postseed <- sum(observed[post_seed_2022_idx])

  # 2022 weekly peak timing/magnitude (observed vs. posterior median), restricted to non-seed weeks
  peak_obs_idx <- all_2022_idx[which.max(observed[all_2022_idx])]
  pred_2022_weekly_median <- apply(C_pred[, all_2022_idx], 2, median)
  peak_pred_idx <- all_2022_idx[which.max(pred_2022_weekly_median)]

  tibble(
    q = q, hmc_pass = b$hmc$hmc_pass, divergences = b$hmc$divergences,
    max_rhat = b$hmc$maximum_rhat, min_bulk_ess = b$hmc$minimum_bulk_ess,
    immune_2018_median = median(immune_prop[, sero_idx]),
    immune_2018_lo = quantile(immune_prop[, sero_idx], .025), immune_2018_hi = quantile(immune_prop[, sero_idx], .975),
    sero_pred_median = median(sero_pred), sero_pred_lo = quantile(sero_pred, .025), sero_pred_hi = quantile(sero_pred, .975),
    S_prop_pre_seed_median = median(S_prop[, pre_seed_idx]), S_prop_pre_seed_lo = quantile(S_prop[, pre_seed_idx], .025), S_prop_pre_seed_hi = quantile(S_prop[, pre_seed_idx], .975),
    obs_2017 = sum(observed[idx_2017]), pred_2017_median = median(pred_2017), pred_2017_lo = quantile(pred_2017, .025), pred_2017_hi = quantile(pred_2017, .975),
    obs_2022_all = sum(observed[all_2022_idx]), pred_2022_all_median = median(pred_2022_all), pred_2022_all_lo = quantile(pred_2022_all, .025), pred_2022_all_hi = quantile(pred_2022_all, .975),
    obs_2022_postseed = obs_2022_postseed, pred_2022_postseed_median = median(pred_2022_postseed),
    pred_2022_postseed_lo = quantile(pred_2022_postseed, .025), pred_2022_postseed_hi = quantile(pred_2022_postseed, .975),
    peak_week_observed = as.character(b$weekly_data$week_start[peak_obs_idx]), peak_magnitude_observed = observed[peak_obs_idx],
    peak_week_predicted_median = as.character(b$weekly_data$week_start[peak_pred_idx]), peak_magnitude_predicted_median = pred_2022_weekly_median[which.max(pred_2022_weekly_median)]
  )
})
summary_df <- bind_rows(summary_rows)
write_csv(summary_df, file.path(v4_8_dir, "v4_8_seeded_recurrence_comparison.csv"))
message("[v4.8] summary table saved."); print(as.data.frame(summary_df))

# --- Ascertainment-multiplier diagnostic (descriptive only, v4.7) -----------
v4_7_summary <- read_csv(file.path(v4_7_dir, "fixed_q_sweep_full_period_comparison.csv"), show_col_types = FALSE)
ascertainment_diagnostic <- v4_7_summary |>
  transmute(q, obs_2022 = obs_2022, pred_2022_v4_7_median = pred_2022_median,
            required_ascertainment_multiplier = obs_2022 / pred_2022_median)
write_csv(ascertainment_diagnostic, file.path(v4_8_dir, "v4_7_ascertainment_multiplier_diagnostic.csv"))
message("[v4.8] descriptive-only ascertainment multiplier (v4.7, NOT used in any likelihood):")
print(as.data.frame(ascertainment_diagnostic))

theme_v4 <- theme_classic(base_size = 8.5) + theme(panel.grid.major.y = element_line(colour = "grey90"))

# --- Panel A: weekly case PPC, seed window shaded ----------------------------
weekly_ppc <- bind_rows(lapply(seq_along(Q_GRID), function(i) {
  q <- Q_GRID[i]; b <- bundles_v8[[i]]; fit <- b$fit
  C_pred <- rstan::extract(fit, "C_pred")$C_pred
  tibble(q = q, week_start = as.Date(b$weekly_data$week_start), observed = b$weekly_data$cases,
         pred_median = apply(C_pred, 2, median), pred_lo = apply(C_pred, 2, quantile, .025), pred_hi = apply(C_pred, 2, quantile, .975))
}))
weekly_ppc$q_label <- factor(sprintf("q = %.2f", weekly_ppc$q), levels = sprintf("q = %.2f", Q_GRID))
seed_rect <- tibble(xmin = min(bundles_v8[[1]]$seed_dates), xmax = max(bundles_v8[[1]]$seed_dates) + 7)

pA <- ggplot(weekly_ppc, aes(week_start)) +
  annotate("rect", xmin = seed_rect$xmin, xmax = seed_rect$xmax, ymin = -Inf, ymax = Inf, fill = "grey40", alpha = .18) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), fill = "steelblue", alpha = .25) +
  geom_line(aes(y = pred_median), colour = "steelblue", linewidth = .3) +
  geom_line(aes(y = observed), colour = "black", linewidth = .3) +
  facet_wrap(~q_label, ncol = 2) +
  labs(title = "A. v4.8 weekly case PPC by fixed q (grey band = 2022 seed window, not scored)",
       x = NULL, y = "weekly cases") +
  theme_v4 + theme(strip.background = element_blank())

pB <- ggplot(summary_df, aes(q)) +
  geom_ribbon(aes(ymin = immune_2018_lo, ymax = immune_2018_hi), fill = "darkorange", alpha = .25) +
  geom_line(aes(y = immune_2018_median), colour = "darkorange") + geom_point(aes(y = immune_2018_median), colour = "darkorange") +
  geom_hline(yintercept = 103 / 404, linetype = 2, colour = "black") +
  labs(title = "B. Model-implied 2018 immune fraction vs q", x = "q", y = "immune fraction") + theme_v4

pC <- ggplot(summary_df, aes(q)) +
  geom_ribbon(aes(ymin = pred_2017_lo, ymax = pred_2017_hi), fill = "seagreen", alpha = .25) +
  geom_line(aes(y = pred_2017_median), colour = "seagreen") + geom_point(aes(y = pred_2017_median), colour = "seagreen") +
  geom_hline(aes(yintercept = obs_2017), linetype = 2, colour = "black") +
  labs(title = "C. 2017 case-count PPC vs q", x = "q", y = "2017 total cases") + theme_v4

pD <- ggplot(summary_df, aes(q)) +
  geom_ribbon(aes(ymin = pred_2022_postseed_lo, ymax = pred_2022_postseed_hi), fill = "firebrick", alpha = .25) +
  geom_line(aes(y = pred_2022_postseed_median), colour = "firebrick") + geom_point(aes(y = pred_2022_postseed_median), colour = "firebrick") +
  geom_hline(aes(yintercept = obs_2022_postseed), linetype = 2, colour = "black") +
  labs(title = "D. v4.8 2022 case-count PPC vs q (post-seed weeks only)", x = "q", y = "2022 total cases (excl. seed)") + theme_v4

pE <- ggplot(summary_df, aes(q, max_rhat, colour = hmc_pass)) +
  geom_hline(yintercept = 1.01, linetype = 2, colour = "grey50") +
  geom_point(size = 2) + geom_line(aes(group = 1), colour = "grey60") +
  labs(title = "E. HMC gate (max Rhat) across the q grid", x = "q", y = "max Rhat", colour = "HMC pass") +
  theme_v4 + theme(legend.position = "bottom")

# --- Panel F: v4.7 vs v4.8 predicted 2022 total cases by q -------------------
compare_2022 <- bind_rows(
  v4_7_summary |> transmute(q, model = "v4.7 (continuous)", pred_median = pred_2022_median, pred_lo = pred_2022_lo, pred_hi = pred_2022_hi, obs = obs_2022),
  summary_df |> transmute(q, model = "v4.8 (seeded, post-seed only)", pred_median = pred_2022_postseed_median, pred_lo = pred_2022_postseed_lo, pred_hi = pred_2022_postseed_hi, obs = obs_2022_postseed)
)
pF <- ggplot(compare_2022, aes(q, pred_median, colour = model)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi, fill = model), alpha = .18, colour = NA) +
  geom_line() + geom_point() +
  geom_hline(data = compare_2022 |> filter(model == "v4.8 (seeded, post-seed only)"), aes(yintercept = obs), linetype = 2, colour = "black") +
  labs(title = "F. v4.7 vs v4.8 predicted 2022 total cases by q (dashed = observed, post-seed)", x = "q", y = "2022 total cases") +
  theme_v4 + theme(legend.position = "bottom")

pG <- ggplot(summary_df, aes(q)) +
  geom_ribbon(aes(ymin = S_prop_pre_seed_lo, ymax = S_prop_pre_seed_hi), fill = "seagreen", alpha = .25) +
  geom_line(aes(y = S_prop_pre_seed_median), colour = "seagreen") + geom_point(aes(y = S_prop_pre_seed_median), colour = "seagreen") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
  labs(title = "G. Susceptible fraction immediately before the 2022 seed, by q", x = "q", y = "S proportion") + theme_v4

figure <- (pA) / (pB | pC | pD) / (pE | pF | pG) + plot_layout(heights = c(2.4, 1, 1))
figure_dir <- file.path(root, "03_Output/figures/renewal_v4_8_seeded_recurrence")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(figure_dir, "v4_8_seeded_recurrence_summary.png")
ggsave(out_path, figure, width = 240, height = 320, units = "mm", dpi = 300)
message("[v4.8] figure saved: ", out_path)
