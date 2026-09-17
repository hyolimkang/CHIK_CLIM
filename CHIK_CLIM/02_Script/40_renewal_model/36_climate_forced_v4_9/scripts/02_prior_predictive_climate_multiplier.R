# Climate-forced v4.9 extension -- Phase 2 preparatory step: prior-predictive
# check of the proposed climate_multiplier[t] = exp(beta_T*z_T_anom[t] +
# beta_P*z_P_anom[t]) BEFORE any Stan file is compiled or fit.
#
# This is a pure R simulation against the REAL observed z_T_anom/z_P_anom
# for the Phase-2 canary state (Ceara, primary lag spec) produced by
# 01_build_climate_anomaly_covariates.R. It does not touch Stan, does not
# compile anything, and does not fit anything -- it only asks: "given the
# proposed regularising priors beta_T, beta_P ~ normal(0, 0.15), how extreme
# can the resulting multiplicative shift on R0(t) plausibly get, over the
# range of anomalies actually observed?"
#
# Does NOT modify any existing file.

required_packages <- c("dplyr", "readr", "ggplot2", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(tidyr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir  <- file.path(root, "03_Output/tables/climate_forced_v4_9")
figure_dir <- file.path(root, "03_Output/figures/climate_forced_v4_9")

set.seed(20260916)
N_PRIOR_DRAWS <- 4000
BETA_SD <- 0.15  # proposed prior: beta_T, beta_P ~ normal(0, BETA_SD)

anomalies <- read_csv(file.path(table_dir, "climate_anomaly_covariates.csv"), show_col_types = FALSE)

ce_primary <- anomalies |>
  filter(state == "CE", lag_spec == "primary", !is.na(z_T_anom), !is.na(z_P_anom))

message(sprintf("CE primary spec: %d complete weeks (of %d total) used for the prior-predictive check",
                 nrow(ce_primary), sum(anomalies$state == "CE" & anomalies$lag_spec == "primary")))

beta_T_draws <- rnorm(N_PRIOR_DRAWS, 0, BETA_SD)
beta_P_draws <- rnorm(N_PRIOR_DRAWS, 0, BETA_SD)

# ---- 1. Multiplier at every OBSERVED (z_T_anom, z_P_anom) week, crossed with every prior draw ----
# This is the honest prior-predictive distribution: it uses the actual joint
# (zT, zP) pairs that occurred historically in Ceara (they are correlated --
# e.g. hot+dry weeks co-occur), not independently-extreme corners.
grid <- expand.grid(draw = seq_len(N_PRIOR_DRAWS), week_row = seq_len(nrow(ce_primary)))
grid$beta_T <- beta_T_draws[grid$draw]
grid$beta_P <- beta_P_draws[grid$draw]
grid$z_T <- ce_primary$z_T_anom[grid$week_row]
grid$z_P <- ce_primary$z_P_anom[grid$week_row]
grid$climate_multiplier <- exp(grid$beta_T * grid$z_T + grid$beta_P * grid$z_P)

summ_over_weeks <- quantile(grid$climate_multiplier, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))
message("\n=== Prior-predictive climate_multiplier, CE primary spec, over ALL observed (zT,zP) weeks x prior draws ===")
print(round(summ_over_weeks, 3))

# ---- 2. Multiplier at illustrative fixed (zT, zP) corners (percentiles of the OBSERVED CE range) ----
z_percentiles <- c(0.01, 0.25, 0.5, 0.75, 0.99)
zT_grid <- quantile(ce_primary$z_T_anom, z_percentiles)
zP_grid <- quantile(ce_primary$z_P_anom, z_percentiles)

corner_summary <- expand.grid(zT_pctl = names(zT_grid), zP_pctl = names(zP_grid)) |>
  mutate(z_T = zT_grid[zT_pctl], z_P = zP_grid[zP_pctl])

corner_summary <- bind_rows(lapply(seq_len(nrow(corner_summary)), function(i) {
  m <- exp(beta_T_draws * corner_summary$z_T[i] + beta_P_draws * corner_summary$z_P[i])
  tibble(zT_pctl = corner_summary$zT_pctl[i], zP_pctl = corner_summary$zP_pctl[i],
         z_T = corner_summary$z_T[i], z_P = corner_summary$z_P[i],
         mult_median = median(m), mult_p05 = quantile(m, 0.05), mult_p95 = quantile(m, 0.95))
}))

message("\n=== Prior-predictive climate_multiplier at fixed (zT, zP) percentile corners (CE observed range) ===")
print(as.data.frame(corner_summary), digits = 3)

# ---- 3. Worst-case realistic corner: 1st/99th percentile of BOTH zT and zP simultaneously ----
worst_case <- tibble(
  scenario = c("hottest+driest (p99 zT, p01 zP)", "coolest+wettest (p01 zT, p99 zP)"),
  z_T = c(zT_grid["99%"], zT_grid["1%"]),
  z_P = c(zP_grid["1%"], zP_grid["99%"])
) |>
  rowwise() |>
  mutate(mult_median = median(exp(beta_T_draws * z_T + beta_P_draws * z_P)),
         mult_p05 = quantile(exp(beta_T_draws * z_T + beta_P_draws * z_P), 0.05),
         mult_p95 = quantile(exp(beta_T_draws * z_T + beta_P_draws * z_P), 0.95)) |>
  ungroup()

message("\n=== Worst-case realistic corners (CE observed extremes, both anomalies extreme simultaneously) ===")
print(as.data.frame(worst_case), digits = 3)

write_csv(bind_rows(
  corner_summary |> mutate(scenario = paste0("zT=", zT_pctl, ", zP=", zP_pctl)) |>
    select(scenario, z_T, z_P, mult_median, mult_p05, mult_p95),
  worst_case |> select(scenario, z_T, z_P, mult_median, mult_p05, mult_p95)
), file.path(table_dir, "prior_predictive_climate_multiplier_corners.csv"))

# ---- QC figure: prior-predictive climate_multiplier distribution ----
p <- ggplot(grid, aes(climate_multiplier)) +
  geom_histogram(bins = 100, fill = "#0072B2", alpha = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30") +
  scale_x_continuous(limits = c(0, 3)) +
  labs(title = "Prior-predictive climate_multiplier = exp(beta_T*z_T_anom + beta_P*z_P_anom)",
       subtitle = sprintf("beta_T, beta_P ~ normal(0, %.2f); evaluated over all observed CE (primary spec) weekly (zT,zP) pairs x %d prior draws",
                           BETA_SD, N_PRIOR_DRAWS),
       x = "climate_multiplier (1 = no climate effect)", y = "count") +
  theme_classic(base_size = 10)
ggsave(file.path(figure_dir, "QC_E_prior_predictive_climate_multiplier.png"), p, width = 200, height = 130, units = "mm", dpi = 300, bg = "white")

message("\n[saved] prior_predictive_climate_multiplier_corners.csv, QC_E_prior_predictive_climate_multiplier.png")
message("\n=== Prior-predictive check complete. No Stan file compiled or fit. ===")
