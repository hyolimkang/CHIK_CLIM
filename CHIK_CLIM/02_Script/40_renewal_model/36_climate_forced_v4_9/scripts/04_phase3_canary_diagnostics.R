# Climate-forced v4.9 canary -- Phase 3 diagnostics: HMC (classic + rank-
# normalised Rhat), reconstruction stability (baseline v4.9 vs climate
# model), climate identifiability/confounding, PPC comparison, and the
# climate effect itself.

required_packages <- c("here", "rstan", "posterior", "dplyr", "readr", "tibble", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(posterior); library(dplyr); library(readr); library(tibble); library(ggplot2); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "02_Script/40_renewal_model/36_climate_forced_v4_9")
table_dir <- file.path(root, "03_Output/tables/climate_forced_v4_9")
figure_dir <- file.path(root, "03_Output/figures/climate_forced_v4_9")
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

baseline <- readRDS(file.path(root, "02_Script/40_renewal_model/28_v4_9_hierarchical_seasonality/outputs/q0.05/renewal_ceara_v4_9_fit_q0.05.rds"))
climate <- readRDS(file.path(base_dir, "outputs/ce_canary/ce_climate_forced_canary_q0.05.rds"))
dates <- as.Date(climate$weekly_data$week_start)
stopifnot(identical(dates, as.Date(baseline$weekly_data$week_start)))

gate_pars_base <- c("alpha_R", "z_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2", "sigma_season_year", "z_season_year", "phi_obs")
gate_pars_clim <- c(gate_pars_base, "beta_T", "beta_P")

# ================================================================
# A. HMC -- classic (rstan) AND rank-normalised (posterior) Rhat
# ================================================================
rank_normalised_rhat <- function(fit, pars) {
  draws <- rstan::extract(fit, pars = pars, permuted = FALSE, inc_warmup = FALSE)
  sd <- posterior::summarise_draws(posterior::as_draws_array(draws), rhat = posterior::rhat, ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail)
  list(max_rhat = max(sd$rhat, na.rm = TRUE), min_ess_bulk = min(sd$ess_bulk, na.rm = TRUE), min_ess_tail = min(sd$ess_tail, na.rm = TRUE))
}
rn_base <- rank_normalised_rhat(baseline$fit, gate_pars_base)
rn_clim <- rank_normalised_rhat(climate$fit, gate_pars_clim)

hmc_compare <- tibble(
  model = c("baseline v4.9 (q=0.05)", "climate-forced canary (q=0.05)"),
  divergences = c(baseline$hmc$divergences, climate$hmc$divergences),
  max_treedepth_hits = c(baseline$hmc$max_treedepth_hits, climate$hmc$max_treedepth_hits),
  max_rhat_classic = c(baseline$hmc$maximum_rhat, climate$hmc$maximum_rhat),
  max_rhat_rank_normalised = c(rn_base$max_rhat, rn_clim$max_rhat),
  min_ess_bulk_rank_normalised = c(rn_base$min_ess_bulk, rn_clim$min_ess_bulk),
  min_ess_tail_rank_normalised = c(rn_base$min_ess_tail, rn_clim$min_ess_tail),
  hmc_pass_classic_gate = c(baseline$hmc$hmc_pass, climate$hmc$hmc_pass)
)
message("=== A. HMC comparison (classic vs rank-normalised Rhat) ===")
print(as.data.frame(hmc_compare), digits = 4)
write_csv(hmc_compare, file.path(table_dir, "CE_canary_HMC_comparison.csv"))

# ================================================================
# B. Reconstruction stability: S_prop, cumulative infection, alpha_R, etc.
# ================================================================
draws_base <- rstan::extract(baseline$fit, pars = c("S_prop", "immune_prop", "R0_t", "alpha_R", "beta_sin1", "beta_cos1", "sigma_season_year", "phi_obs"))
draws_clim <- rstan::extract(climate$fit, pars = c("S_prop", "immune_prop", "R0_t", "alpha_R", "beta_sin1", "beta_cos1", "sigma_season_year", "phi_obs", "beta_T", "beta_P"))

checkpoint_dates <- as.Date(c("2016-12-31", "2018-12-31", "2020-12-31", "2022-12-31", "2025-12-21"))
checkpoint_idx <- sapply(checkpoint_dates, function(d) which.min(abs(dates - d)))

checkpoint_compare <- bind_rows(lapply(seq_along(checkpoint_dates), function(i) {
  t_idx <- checkpoint_idx[i]
  tibble(checkpoint = checkpoint_dates[i],
         S_base = median(draws_base$S_prop[, t_idx]), S_climate = median(draws_clim$S_prop[, t_idx]),
         immune_base = median(draws_base$immune_prop[, t_idx]), immune_climate = median(draws_clim$immune_prop[, t_idx]),
         R0_base = median(draws_base$R0_t[, t_idx]), R0_climate = median(draws_clim$R0_t[, t_idx])) |>
    mutate(delta_S = S_climate - S_base, delta_immune = immune_climate - immune_base, delta_R0 = R0_climate - R0_base)
}))
message("\n=== B. Checkpoint comparison: baseline vs climate-forced (S_prop, immune_prop, R0) ===")
print(as.data.frame(checkpoint_compare), digits = 4)
write_csv(checkpoint_compare, file.path(table_dir, "CE_canary_checkpoint_comparison.csv"))

param_compare <- tibble(
  parameter = c("alpha_R", "beta_sin1", "beta_cos1", "sigma_season_year", "phi_obs", "min_S_prop", "cumulative_infection_2025"),
  baseline_median = c(median(draws_base$alpha_R), median(draws_base$beta_sin1), median(draws_base$beta_cos1),
                       median(draws_base$sigma_season_year), median(draws_base$phi_obs),
                       median(apply(draws_base$S_prop, 1, min)), median(draws_base$immune_prop[, ncol(draws_base$immune_prop)])),
  climate_median = c(median(draws_clim$alpha_R), median(draws_clim$beta_sin1), median(draws_clim$beta_cos1),
                      median(draws_clim$sigma_season_year), median(draws_clim$phi_obs),
                      median(apply(draws_clim$S_prop, 1, min)), median(draws_clim$immune_prop[, ncol(draws_clim$immune_prop)]))
) |> mutate(pct_change = 100 * (climate_median - baseline_median) / baseline_median)
message("\n=== B. Key parameter comparison ===")
print(as.data.frame(param_compare), digits = 4)
write_csv(param_compare, file.path(table_dir, "CE_canary_parameter_comparison.csv"))

max_abs_pct_change_S <- max(abs(checkpoint_compare$delta_S / checkpoint_compare$S_base)) * 100
message(sprintf("\nMax relative change in S_prop at any checkpoint: %.2f%% -- %s",
                 max_abs_pct_change_S, if (max_abs_pct_change_S < 5) "SMALL: climate does not radically change inferred susceptibility" else "LARGE: potential identifiability problem, flag before proceeding"))

# ================================================================
# C. Climate identifiability / confounding: posterior correlations
# ================================================================
cor_tbl <- tibble(
  with = c("alpha_R", "beta_sin1", "beta_cos1", "sigma_season_year", "phi_obs", "min_S_prop_draw", "final_immune_draw"),
  cor_beta_T = c(cor(draws_clim$beta_T, draws_clim$alpha_R), cor(draws_clim$beta_T, draws_clim$beta_sin1), cor(draws_clim$beta_T, draws_clim$beta_cos1),
                 cor(draws_clim$beta_T, draws_clim$sigma_season_year), cor(draws_clim$beta_T, draws_clim$phi_obs),
                 cor(draws_clim$beta_T, apply(draws_clim$S_prop, 1, min)), cor(draws_clim$beta_T, draws_clim$immune_prop[, ncol(draws_clim$immune_prop)])),
  cor_beta_P = c(cor(draws_clim$beta_P, draws_clim$alpha_R), cor(draws_clim$beta_P, draws_clim$beta_sin1), cor(draws_clim$beta_P, draws_clim$beta_cos1),
                 cor(draws_clim$beta_P, draws_clim$sigma_season_year), cor(draws_clim$beta_P, draws_clim$phi_obs),
                 cor(draws_clim$beta_P, apply(draws_clim$S_prop, 1, min)), cor(draws_clim$beta_P, draws_clim$immune_prop[, ncol(draws_clim$immune_prop)]))
)
cor_beta_T_beta_P <- cor(draws_clim$beta_T, draws_clim$beta_P)
message("\n=== C. Posterior correlations: beta_T / beta_P with other parameters ===")
print(as.data.frame(cor_tbl), digits = 3)
message(sprintf("\ncor(beta_T, beta_P) = %.4f (prior/data anomaly correlation was -0.389)", cor_beta_T_beta_P))
write_csv(cor_tbl, file.path(table_dir, "CE_canary_confounding_correlations.csv"))

# ================================================================
# D. PPC comparison
# ================================================================
pred_base <- rstan::extract(baseline$fit, pars = "C_pred")$C_pred
pred_clim <- rstan::extract(climate$fit, pars = "C_pred")$C_pred
obs <- climate$weekly_data$cases

coverage_base <- mean(obs >= apply(pred_base, 2, quantile, .025) & obs <= apply(pred_base, 2, quantile, .975))
coverage_clim <- mean(obs >= apply(pred_clim, 2, quantile, .025) & obs <= apply(pred_clim, 2, quantile, .975))
rmse_base <- sqrt(mean((obs - apply(pred_base, 2, median))^2))
rmse_clim <- sqrt(mean((obs - apply(pred_clim, 2, median))^2))

years <- as.integer(format(dates, "%Y"))
annual_obs <- tapply(obs, years, sum)
annual_pred_base <- apply(t(apply(pred_base, 1, function(r) tapply(r, years, sum))), 2, median)
annual_pred_clim <- apply(t(apply(pred_clim, 1, function(r) tapply(r, years, sum))), 2, median)
annual_compare <- tibble(year = as.integer(names(annual_obs)), observed = as.numeric(annual_obs),
                          pred_base = annual_pred_base, pred_climate = annual_pred_clim,
                          ratio_base = annual_pred_base / as.numeric(annual_obs), ratio_climate = annual_pred_clim / as.numeric(annual_obs))
message("\n=== D. PPC comparison ===")
message(sprintf("Weekly 95%% coverage: baseline=%.1f%%, climate=%.1f%%", 100*coverage_base, 100*coverage_clim))
message(sprintf("Weekly RMSE: baseline=%.1f, climate=%.1f", rmse_base, rmse_clim))
print(as.data.frame(annual_compare), digits = 3)
write_csv(annual_compare, file.path(table_dir, "CE_canary_annual_PPC_comparison.csv"))

ppc_summary <- tibble(model = c("baseline", "climate"), weekly_coverage_pct = c(100*coverage_base, 100*coverage_clim), weekly_RMSE = c(rmse_base, rmse_clim))
write_csv(ppc_summary, file.path(table_dir, "CE_canary_PPC_summary.csv"))

# ================================================================
# E. Climate effect itself
# ================================================================
beta_T_summary <- tibble(parameter = "beta_T", median = median(draws_clim$beta_T), lo95 = quantile(draws_clim$beta_T, .025), hi95 = quantile(draws_clim$beta_T, .975),
                          prob_positive = mean(draws_clim$beta_T > 0))
beta_P_summary <- tibble(parameter = "beta_P", median = median(draws_clim$beta_P), lo95 = quantile(draws_clim$beta_P, .025), hi95 = quantile(draws_clim$beta_P, .975),
                          prob_positive = mean(draws_clim$beta_P > 0))
beta_summary <- bind_rows(beta_T_summary, beta_P_summary)
message("\n=== E. Climate effect (beta_T, beta_P posteriors) ===")
print(as.data.frame(beta_summary), digits = 4)
write_csv(beta_summary, file.path(table_dir, "CE_canary_climate_effect_summary.csv"))

mult_draws <- rstan::extract(climate$fit, pars = "climate_multiplier")$climate_multiplier
mult_summary <- tibble(week_start = dates, median = apply(mult_draws, 2, median), lo95 = apply(mult_draws, 2, quantile, .025), hi95 = apply(mult_draws, 2, quantile, .975))
message(sprintf("climate_multiplier(t) range (posterior median): [%.3f, %.3f]", min(mult_summary$median), max(mult_summary$median)))
write_csv(mult_summary, file.path(table_dir, "CE_canary_climate_multiplier_posterior.csv"))

p_beta <- ggplot(bind_rows(tibble(param = "beta_T", value = draws_clim$beta_T), tibble(param = "beta_P", value = draws_clim$beta_P)), aes(value, fill = param)) +
  geom_density(alpha = 0.5) + geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "CE canary: beta_T, beta_P posteriors (prior: N(0, 0.15))", x = "Coefficient value", y = "Density") + theme_v4
ggsave(file.path(figure_dir, "CE_canary_beta_posteriors.png"), p_beta, width = 160, height = 100, units = "mm", dpi = 300, bg = "white")

p_mult <- ggplot(mult_summary, aes(week_start)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#E69F00", alpha = 0.25) + geom_line(aes(y = median), colour = "#E69F00") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  labs(title = "CE canary: posterior climate_multiplier(t)", x = NULL, y = "climate_multiplier") + theme_v4
ggsave(file.path(figure_dir, "CE_canary_climate_multiplier_posterior.png"), p_mult, width = 200, height = 100, units = "mm", dpi = 300, bg = "white")

p_S_compare <- ggplot() +
  geom_line(data = tibble(week_start = dates, S = apply(draws_base$S_prop, 2, median)), aes(week_start, S, colour = "baseline v4.9"), linewidth = 0.5) +
  geom_line(data = tibble(week_start = dates, S = apply(draws_clim$S_prop, 2, median)), aes(week_start, S, colour = "climate-forced canary"), linewidth = 0.5) +
  scale_colour_manual(name = NULL, values = c("baseline v4.9" = "#0072B2", "climate-forced canary" = "#D55E00")) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(title = "CE: S/N reconstruction, baseline vs. climate-forced", x = NULL, y = "S / N") + theme_v4
ggsave(file.path(figure_dir, "CE_canary_S_prop_comparison.png"), p_S_compare, width = 200, height = 100, units = "mm", dpi = 300, bg = "white")

message("\n[saved] All Phase 3 tables and figures to ", table_dir, " / ", figure_dir)
