# Diagnostic comparison: v4.9 (harmonic-only) vs v5.0a (harmonic + AR(1)
# dynamic-R residual), q=0.05. 2017/2022 weekly zoom, phi_obs, A_year vs
# sigma_R_dynamic/rho_R, susceptible trajectory -- same layout as
# 28_v4_9.../04_fit_v4_9_serology_ablation.R's diagnostic companion plot.

required_packages <- c("rstan", "dplyr", "ggplot2", "patchwork", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(ggplot2); library(patchwork); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
b9 <- readRDS(file.path(root, "03_Output/model_fits/ceara/v4_9_hierarchical_seasonality/outputs/q0.05/renewal_ceara_v4_9_fit_q0.05.rds"))
b50 <- readRDS(file.path(root, "02_Script/90_development_archive/02_historical_renewal_versions/29_v5_0_dynamic_R/outputs/v5_0a_q0.05/renewal_ceara_v5_0a_fit_q0.05.rds"))
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

zoom_year <- function(bundle, label, year, exclude_idx = integer(0)) {
  w <- bundle$weekly_data; w$week_start <- as.Date(w$week_start)
  idx <- which(format(w$week_start, "%Y") == as.character(year))
  idx <- setdiff(idx, exclude_idx)
  C_pred <- rstan::extract(bundle$fit, "C_pred")$C_pred[, idx]
  tibble(model = label, week_start = w$week_start[idx], observed = w$cases[idx],
         pred_median = apply(C_pred, 2, median), pred_lo = apply(C_pred, 2, quantile, .025), pred_hi = apply(C_pred, 2, quantile, .975))
}
z2017 <- bind_rows(zoom_year(b9, "v4.9", 2017), zoom_year(b50, "v5.0a", 2017))
p_2017 <- ggplot(z2017, aes(week_start)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi, fill = model), alpha = .2) +
  geom_line(aes(y = pred_median, colour = model), linewidth = .8) +
  geom_point(aes(y = observed), colour = "black", size = 1) +
  labs(title = "2017 weekly zoom (q=0.05): observed (black) vs posterior median", x = NULL, y = "weekly cases") +
  theme_v4 + theme(legend.position = "bottom")

z2022 <- bind_rows(zoom_year(b9, "v4.9", 2022, b9$seed_idx), zoom_year(b50, "v5.0a", 2022, b50$seed_idx))
p_2022 <- ggplot(z2022, aes(week_start)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi, fill = model), alpha = .2) +
  geom_line(aes(y = pred_median, colour = model), linewidth = .8) +
  geom_point(aes(y = observed), colour = "black", size = 1) +
  labs(title = "2022 weekly zoom (q=0.05): observed (black) vs posterior median", x = NULL, y = "weekly cases") +
  theme_v4 + theme(legend.position = "bottom")

phi_df <- bind_rows(
  tibble(model = "v4.9", phi_obs = rstan::extract(b9$fit, "phi_obs")$phi_obs),
  tibble(model = "v5.0a", phi_obs = rstan::extract(b50$fit, "phi_obs")$phi_obs)
)
p_phi <- ggplot(phi_df, aes(phi_obs, fill = model)) + geom_histogram(position = "identity", alpha = .5, bins = 40) +
  labs(title = "phi_obs posterior (q=0.05): v4.9 vs v5.0a", x = "phi_obs", y = "count") + theme_v4 + theme(legend.position = "bottom")

rho_R <- rstan::extract(b50$fit, "rho_R")$rho_R
sigma_R_dynamic <- rstan::extract(b50$fit, "sigma_R_dynamic")$sigma_R_dynamic
dyn_df <- tibble(rho_R = rho_R, sigma_R_dynamic = sigma_R_dynamic)
p_dyn <- ggplot(dyn_df, aes(rho_R, sigma_R_dynamic)) + geom_point(alpha = .15, size = .6) +
  labs(title = sprintf("v5.0a dynamic-R posterior: rho_R=%.2f [%.2f,%.2f], sigma_R_dynamic=%.3f [%.3f,%.3f]",
                        median(rho_R), quantile(rho_R,.025), quantile(rho_R,.975),
                        median(sigma_R_dynamic), quantile(sigma_R_dynamic,.025), quantile(sigma_R_dynamic,.975)),
       x = "rho_R", y = "sigma_R_dynamic") + theme_v4

A_year9 <- rstan::extract(b9$fit, "A_year")$A_year; A_year50 <- rstan::extract(b50$fit, "A_year")$A_year
years <- b9$config$years
A_df <- bind_rows(
  tibble(model = "v4.9", year = years, median = apply(A_year9, 2, median), lo = apply(A_year9, 2, quantile, .025), hi = apply(A_year9, 2, quantile, .975)),
  tibble(model = "v5.0a", year = years, median = apply(A_year50, 2, median), lo = apply(A_year50, 2, quantile, .025), hi = apply(A_year50, 2, quantile, .975))
)
p_A <- ggplot(A_df, aes(factor(year), median, colour = model)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = .4)) +
  labs(title = "A_year: v4.9 vs v5.0a (does delta_R absorb A_year's role?)", x = "year", y = "A_year") +
  theme_v4 + theme(legend.position = "bottom")

S9 <- rstan::extract(b9$fit, "S_prop")$S_prop; S50 <- rstan::extract(b50$fit, "S_prop")$S_prop
w9 <- as.Date(b9$weekly_data$week_start)
S_df <- bind_rows(
  tibble(model = "v4.9", week_start = w9, median = apply(S9, 2, median), lo = apply(S9, 2, quantile, .025), hi = apply(S9, 2, quantile, .975)),
  tibble(model = "v5.0a", week_start = w9, median = apply(S50, 2, median), lo = apply(S50, 2, quantile, .025), hi = apply(S50, 2, quantile, .975))
)
p_S <- ggplot(S_df, aes(week_start, median, colour = model, fill = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, colour = NA) + geom_line() +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
  labs(title = "Susceptible proportion: v4.9 vs v5.0a (q=0.05)", x = NULL, y = "S proportion") + theme_v4 + theme(legend.position = "bottom")

figure <- (p_2017 | p_2022) / (p_phi | p_dyn) / (p_A | p_S)
figure_dir <- file.path(root, "03_Output/figures/renewal_v5_0_dynamic_R")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(figure_dir, "v4_9_vs_v5_0a_diagnostic.png"), figure, width = 300, height = 300, units = "mm", dpi = 300)
message("saved: ", file.path(figure_dir, "v4_9_vs_v5_0a_diagnostic.png"))
