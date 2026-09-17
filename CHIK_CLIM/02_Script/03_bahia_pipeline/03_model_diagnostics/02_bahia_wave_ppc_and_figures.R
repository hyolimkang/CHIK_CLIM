# Bahia wave-level PPC, susceptibility trajectories, HMC diagnostics table,
# and required figures 1-6 (BAHIA_V4_9_FROZEN_MODEL_ASSESSMENT.md Section
# 8/9/13). Uses the 9 objectively pre-defined Bahia major-epidemic waves
# from 03_Output/07_national_pipeline/tables/national_wave_analysis/brazil_chik_wave_analysis_master.csv
# (NOT re-derived here -- the existing wave-census pipeline's own onset/end
# dates are used as-is, per instruction not to redefine wave boundaries).
# All figures saved INSIDE the Bahia subfolder (not 03_Output/figures/).

required_packages <- c("rstan", "dplyr", "tibble", "readr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(readr); library(ggplot2); library(patchwork) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
base_dir <- file.path(root, "03_Output/03_bahia_pipeline/model_fits/v4_9_replication")
Q_GRID <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

waves <- read_csv(file.path(root, "03_Output/07_national_pipeline/tables/national_wave_analysis/brazil_chik_wave_analysis_master.csv"), show_col_types = FALSE) |>
  dplyr::filter(state == "BA") |>
  transmute(wave_id, wave_order, onset_week = as.Date(onset_week), end_week = as.Date(end_week), total_cases_census = total_cases)
message("[bahia] ", nrow(waves), " major-epidemic waves loaded from the existing wave-census pipeline (not re-derived).")
print(as.data.frame(waves))

load_bahia_bundle <- function(q) {
  tag <- paste0("q", sprintf("%.2f", q))
  readRDS(file.path(base_dir, "outputs", tag, paste0("renewal_bahia_v4_9_fit_", tag, ".rds")))
}

epidemic_width_weeks <- function(cases_vec) {
  cum <- cumsum(cases_vec) / sum(cases_vec)
  lo <- which(cum >= 0.10)[1]; hi <- which(cum >= 0.90)[1]
  hi - lo + 1
}

analyse_one_q <- function(q) {
  b <- load_bahia_bundle(q)
  fit <- b$fit
  w <- b$weekly_data; dates <- as.Date(w$week_start)
  C_pred <- rstan::extract(fit, "C_pred")$C_pred
  S_prop <- rstan::extract(fit, "S_prop")$S_prop

  wave_rows <- lapply(seq_len(nrow(waves)), function(i) {
    onset <- waves$onset_week[i]; end <- pmin(waves$end_week[i], max(dates)) # clip censored trailing wave to data end
    idx <- which(dates >= onset & dates <= end)
    if (length(idx) < 2) return(NULL)
    obs <- w$cases[idx]
    pred_totals <- rowSums(C_pred[, idx, drop = FALSE])
    pred_median_traj <- apply(C_pred[, idx, drop = FALSE], 2, median)
    per_draw_peak <- apply(C_pred[, idx, drop = FALSE], 1, max)

    pre_idx <- max(1, min(idx) - 1)
    post_idx <- min(length(dates), max(idx) + 1)

    tibble(
      state = "BA", q = q, wave_id = waves$wave_id[i], wave_onset = as.character(onset),
      observed_total = sum(obs), predicted_total_median = median(pred_totals),
      predicted_total_lo = quantile(pred_totals, .025), predicted_total_hi = quantile(pred_totals, .975),
      total_ratio = median(pred_totals) / sum(obs),
      observed_peak = max(obs), predicted_peak_median = median(per_draw_peak),
      peak_ratio = median(per_draw_peak) / max(obs),
      observed_peak_week = as.character(dates[idx][which.max(obs)]),
      predicted_peak_week = as.character(dates[idx][which.max(pred_median_traj)]),
      epidemic_width_weeks_observed = epidemic_width_weeks(obs),
      epidemic_width_weeks_predicted_median = epidemic_width_weeks(pmax(pred_median_traj, 1e-9)),
      pre_wave_susceptible = median(S_prop[, pre_idx]),
      post_wave_susceptible = median(S_prop[, post_idx]),
      max_rhat = b$hmc$maximum_rhat, divergences = b$hmc$divergences,
      min_bulk_ess = b$hmc$minimum_bulk_ess, min_tail_ess = NA_real_, bfmi_min = min(b$hmc$bfmi)
    )
  })
  bind_rows(wave_rows)
}

# min_tail_ess needs a direct posterior:: call since compute_hmc_gate() only stores bulk ESS
add_tail_ess <- function(df, q) {
  b <- load_bahia_bundle(q)
  core_pars <- c("alpha_R", "z_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2", "sigma_season_year", "z_season_year", "phi_obs")
  draws <- posterior::as_draws_array(rstan::extract(b$fit, pars = core_pars, permuted = FALSE))
  summ <- posterior::summarise_draws(draws, "ess_tail")
  df$min_tail_ess <- min(summ$ess_tail, na.rm = TRUE)
  df
}

message("[bahia] computing wave-level PPC across the q grid...")
all_rows <- bind_rows(lapply(Q_GRID, function(q) add_tail_ess(analyse_one_q(q), q)))
table_dir <- file.path(root, "03_Output/03_bahia_pipeline/tables/renewal_bahia_v4_9_replication")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(all_rows, file.path(table_dir, "bahia_wave_level_ppc_summary.csv"))
message("[bahia] saved: ", file.path(table_dir, "bahia_wave_level_ppc_summary.csv"))
print(as.data.frame(all_rows))

# --- HMC diagnostics table --------------------------------------------------
hmc_table <- bind_rows(lapply(Q_GRID, function(q) {
  b <- load_bahia_bundle(q)
  tibble(q = q, divergences = b$hmc$divergences, max_treedepth_hits = b$hmc$max_treedepth_hits,
         max_rhat = b$hmc$maximum_rhat, min_bulk_ess = b$hmc$minimum_bulk_ess,
         bfmi_min = min(b$hmc$bfmi), hmc_pass = b$hmc$hmc_pass, elapsed_seconds = b$config$elapsed_seconds)
}))
write_csv(hmc_table, file.path(table_dir, "bahia_hmc_diagnostics_by_q.csv"))
message("[bahia] HMC diagnostics:"); print(as.data.frame(hmc_table))

theme_v4 <- theme_classic(base_size = 8.5) + theme(panel.grid.major.y = element_line(colour = "grey90"))

# --- FIGURE 1: full weekly observed vs posterior median + 50/95% PPC, all q -
weekly_ppc <- bind_rows(lapply(Q_GRID, function(q) {
  b <- load_bahia_bundle(q)
  C_pred <- rstan::extract(b$fit, "C_pred")$C_pred
  tibble(q = q, week_start = as.Date(b$weekly_data$week_start), observed = b$weekly_data$cases,
         pred_median = apply(C_pred, 2, median), pred_q25 = apply(C_pred, 2, quantile, .25), pred_q75 = apply(C_pred, 2, quantile, .75),
         pred_lo = apply(C_pred, 2, quantile, .025), pred_hi = apply(C_pred, 2, quantile, .975))
}))
weekly_ppc$q_label <- factor(sprintf("q = %.2f", weekly_ppc$q), levels = sprintf("q = %.2f", Q_GRID))
fig1 <- ggplot(weekly_ppc, aes(week_start)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), fill = "steelblue", alpha = .2) +
  geom_ribbon(aes(ymin = pred_q25, ymax = pred_q75), fill = "steelblue", alpha = .35) +
  geom_line(aes(y = pred_median), colour = "#08306b", linewidth = .4) +
  geom_line(aes(y = observed), colour = "black", linewidth = .3) +
  facet_wrap(~q_label, ncol = 2) +
  labs(title = "FIGURE 1. Bahia weekly case PPC by fixed q", subtitle = "Black = observed; dark blue = posterior median; bands = 50%/95% CrI. No seed window (none used).",
       x = NULL, y = "weekly cases") + theme_v4 + theme(strip.background = element_blank())

# --- FIGURE 2: susceptibility/immunity trajectories by q --------------------
susc_df <- bind_rows(lapply(Q_GRID, function(q) {
  b <- load_bahia_bundle(q)
  S_prop <- rstan::extract(b$fit, "S_prop")$S_prop
  tibble(q = q, week_start = as.Date(b$weekly_data$week_start), S_median = apply(S_prop, 2, median),
         S_lo = apply(S_prop, 2, quantile, .025), S_hi = apply(S_prop, 2, quantile, .975))
}))
fig2 <- ggplot(susc_df, aes(week_start, S_median, colour = factor(q), fill = factor(q))) +
  geom_ribbon(aes(ymin = S_lo, ymax = S_hi), alpha = .12, colour = NA) + geom_line() +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
  labs(title = "FIGURE 2. Bahia susceptible proportion by q", x = NULL, y = "S proportion", colour = "q", fill = "q") +
  theme_v4 + theme(legend.position = "bottom")

# --- FIGURE 3: R0(t)/Reff(t) for representative q scenarios -----------------
repr_q <- c(0.05, 0.15, 0.30)
r_df <- bind_rows(lapply(repr_q, function(q) {
  b <- load_bahia_bundle(q)
  R0 <- rstan::extract(b$fit, "R0_t")$R0_t; Reff <- rstan::extract(b$fit, "R_eff_t")$R_eff_t
  tibble(q = q, week_start = as.Date(b$weekly_data$week_start), R0_median = apply(R0, 2, median), Reff_median = apply(Reff, 2, median))
}))
fig3 <- ggplot(r_df, aes(week_start)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_line(aes(y = R0_median, colour = "R0(t)")) + geom_line(aes(y = Reff_median, colour = "Reff(t)")) +
  facet_wrap(~sprintf("q = %.2f", q), ncol = 1) +
  scale_colour_manual(values = c("R0(t)" = "#315A7D", "Reff(t)" = "#D55E00")) +
  labs(title = "FIGURE 3. Bahia R0(t)/Reff(t), representative q", x = NULL, y = "reproduction number", colour = NULL) +
  theme_v4 + theme(legend.position = "bottom", strip.background = element_blank())

# --- FIGURE 4/5: wave-level total-case and peak PPC by q --------------------
fig4 <- ggplot(all_rows, aes(q, predicted_total_median, colour = wave_id)) +
  geom_ribbon(aes(ymin = predicted_total_lo, ymax = predicted_total_hi, fill = wave_id), alpha = .08, colour = NA) +
  geom_line() + geom_point(size = .8) +
  geom_hline(aes(yintercept = observed_total, colour = wave_id), linetype = 3, linewidth = .3) +
  facet_wrap(~wave_id, scales = "free_y") +
  labs(title = "FIGURE 4. Bahia wave-level total-case PPC by q", subtitle = "Dotted = observed total", x = "q", y = "wave total cases") +
  theme_v4 + theme(legend.position = "none", strip.text = element_text(size = 6.5))

fig5 <- ggplot(all_rows, aes(q, predicted_peak_median, colour = wave_id)) +
  geom_line() + geom_point(size = .8) +
  geom_hline(aes(yintercept = observed_peak, colour = wave_id), linetype = 3, linewidth = .3) +
  facet_wrap(~wave_id, scales = "free_y") +
  labs(title = "FIGURE 5. Bahia wave-level peak PPC by q", subtitle = "Dotted = observed peak", x = "q", y = "peak weekly cases") +
  theme_v4 + theme(legend.position = "none", strip.text = element_text(size = 6.5))

# --- FIGURE 6: HMC diagnostics across q -------------------------------------
fig6 <- ggplot(hmc_table, aes(q, max_rhat, colour = hmc_pass)) +
  geom_hline(yintercept = 1.01, linetype = 2, colour = "grey50") +
  geom_point(size = 2) + geom_line(aes(group = 1), colour = "grey60") +
  labs(title = "FIGURE 6. Bahia HMC gate (max Rhat) across q", x = "q", y = "max Rhat", colour = "HMC pass") +
  theme_v4 + theme(legend.position = "bottom")

figure_dir <- file.path(root, "03_Output/03_bahia_pipeline/figures/renewal_bahia_v4_9_replication")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(figure_dir, "bahia_figure1_weekly_ppc_by_q.png"), fig1, width = 220, height = 260, units = "mm", dpi = 300)
ggsave(file.path(figure_dir, "bahia_figure2_susceptibility_by_q.png"), fig2, width = 200, height = 140, units = "mm", dpi = 300)
ggsave(file.path(figure_dir, "bahia_figure3_R0_Reff_representative_q.png"), fig3, width = 200, height = 220, units = "mm", dpi = 300)
ggsave(file.path(figure_dir, "bahia_figure4_wave_total_ppc_by_q.png"), fig4, width = 240, height = 220, units = "mm", dpi = 300)
ggsave(file.path(figure_dir, "bahia_figure5_wave_peak_ppc_by_q.png"), fig5, width = 240, height = 220, units = "mm", dpi = 300)
ggsave(file.path(figure_dir, "bahia_figure6_hmc_diagnostics.png"), fig6, width = 160, height = 120, units = "mm", dpi = 300)
message("[bahia] figures 1-6 saved under: ", figure_dir)
