# Descriptive comparison of INDEPENDENT Ceara and Bahia posteriors
# (BAHIA_V4_9_FROZEN_MODEL_ASSESSMENT.md Section 11). NO pooling, NO
# hierarchical model -- just side-by-side posterior summaries at matching q,
# to see which parameters look similar/heterogeneous/weakly identified
# across states.

required_packages <- c("rstan", "dplyr", "tibble", "readr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(readr); library(ggplot2) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "02_Script/40_renewal_model/30_bahia_v4_9_replication")
Q_GRID <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
Q_GRID_CEARA_AVAILABLE <- c(0.05, 0.10, 0.15, 0.20) # Ceara v4.9 was only run at q=0.05 (privileged as reference) plus a partial grid up to 0.20; 0.25/0.30 were never completed for Ceara per user instruction

load_ceara <- function(q) readRDS(file.path(root, "02_Script/40_renewal_model/28_v4_9_hierarchical_seasonality/outputs",
                                              paste0("q", sprintf("%.2f", q)), paste0("renewal_ceara_v4_9_fit_q", sprintf("%.2f", q), ".rds")))
load_bahia <- function(q) readRDS(file.path(base_dir, "outputs", paste0("q", sprintf("%.2f", q)), paste0("renewal_bahia_v4_9_fit_q", sprintf("%.2f", q), ".rds")))

summarise_state <- function(bundle, state, q) {
  fit <- bundle$fit
  A_year <- rstan::extract(fit, "A_year")$A_year
  tibble(
    state = state, q = q,
    alpha_R_median = median(rstan::extract(fit, "alpha_R")$alpha_R),
    alpha_R_lo = quantile(rstan::extract(fit, "alpha_R")$alpha_R, .025),
    alpha_R_hi = quantile(rstan::extract(fit, "alpha_R")$alpha_R, .975),
    beta_sin1_median = median(rstan::extract(fit, "beta_sin1")$beta_sin1),
    beta_cos1_median = median(rstan::extract(fit, "beta_cos1")$beta_cos1),
    beta_sin2_median = median(rstan::extract(fit, "beta_sin2")$beta_sin2),
    beta_cos2_median = median(rstan::extract(fit, "beta_cos2")$beta_cos2),
    sigma_season_year_median = median(rstan::extract(fit, "sigma_season_year")$sigma_season_year),
    A_year_median_across_years = median(apply(A_year, 2, median)),
    A_year_sd_across_years = sd(apply(A_year, 2, median)),
    phi_obs_median = median(rstan::extract(fit, "phi_obs")$phi_obs),
    phi_obs_lo = quantile(rstan::extract(fit, "phi_obs")$phi_obs, .025),
    phi_obs_hi = quantile(rstan::extract(fit, "phi_obs")$phi_obs, .975)
  )
}

comparison <- bind_rows(lapply(Q_GRID, function(q) {
  rows <- list(summarise_state(load_bahia(q), "Bahia (frozen v4.9 structure, case-only, no seed)", q))
  if (q %in% Q_GRID_CEARA_AVAILABLE) {
    rows <- c(list(summarise_state(load_ceara(q), "Ceara (v4.9, w/ serology, w/ 2022 seed)", q)), rows)
  }
  bind_rows(rows)
}))
table_dir <- file.path(root, "03_Output/tables/renewal_bahia_v4_9_replication")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(comparison, file.path(table_dir, "ceara_vs_bahia_posterior_comparison.csv"))
message("[compare] saved: ", file.path(table_dir, "ceara_vs_bahia_posterior_comparison.csv"))
print(as.data.frame(comparison))

theme_v4 <- theme_classic(base_size = 8.5) + theme(panel.grid.major.y = element_line(colour = "grey90"))
p_alpha <- ggplot(comparison, aes(q, alpha_R_median, colour = state)) +
  geom_ribbon(aes(ymin = alpha_R_lo, ymax = alpha_R_hi, fill = state), alpha = .15, colour = NA) +
  geom_line() + geom_point() + labs(title = "alpha_R: Ceara vs Bahia (independent fits, no pooling)", x = "q", y = "alpha_R") +
  theme_v4 + theme(legend.position = "bottom")
p_phi <- ggplot(comparison, aes(q, phi_obs_median, colour = state)) +
  geom_ribbon(aes(ymin = phi_obs_lo, ymax = phi_obs_hi, fill = state), alpha = .15, colour = NA) +
  geom_line() + geom_point() + labs(title = "phi_obs: Ceara vs Bahia", x = "q", y = "phi_obs") +
  theme_v4 + theme(legend.position = "bottom")
p_Ayear <- ggplot(comparison, aes(q, A_year_sd_across_years, colour = state)) +
  geom_line() + geom_point() + labs(title = "SD of A_year across years: Ceara vs Bahia", x = "q", y = "SD(A_year)") +
  theme_v4 + theme(legend.position = "bottom")

figure_dir <- file.path(root, "03_Output/figures/renewal_bahia_v4_9_replication")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
library(patchwork)
ggsave(file.path(figure_dir, "ceara_vs_bahia_parameter_comparison.png"), (p_alpha | p_phi) / p_Ayear, width = 220, height = 200, units = "mm", dpi = 300)
message("[compare] figure saved: ", file.path(figure_dir, "ceara_vs_bahia_parameter_comparison.png"))
