# Climate-forced v4.9 extension -- Phase 1: state-week climate anomaly
# covariates for the three PRE-SPECIFIED lag definitions already used in
# the first-epidemic climate GAM work (no new lag search). Computed for
# all 27 UFs from the existing, unmodified national climate table; the
# state-pooled standardisation is fit on the 4 states with an accepted
# long-term v4.9 reconstruction (CE, BA, RJ, MT) -- the pooled set this
# analysis is actually for -- and documented explicitly as a choice.
#
# Does NOT modify any existing file. Read-only against
# brazil_chik_uf_weekly_climate.csv.

required_packages <- c("here", "dplyr", "readr", "tibble", "lubridate", "tidyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(lubridate); library(tidyr); library(ggplot2) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
table_dir <- file.path(root, "03_Output/tables/climate_forced_v4_9")
figure_dir <- file.path(root, "03_Output/figures/climate_forced_v4_9")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

PRIMARY_STATES <- c("CE", "BA", "RJ", "MT") # the pooled standardisation set -- the 4 accepted reconstructions
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

LAG_SPECS <- list(
  primary       = list(temp_offsets = -2:0, precip_offsets = -6:-2),
  sensitivity_A = list(temp_offsets = -1:0, precip_offsets = -4:-1),
  sensitivity_B = list(temp_offsets = -3:-1, precip_offsets = -8:-3)
)

climate <- read_csv(file.path(root, "03_Output/tables/national_wave_analysis/brazil_chik_uf_weekly_climate.csv"), show_col_types = FALSE) |>
  mutate(week_start = as.Date(week_start), week_of_year = lubridate::isoweek(week_start)) |>
  arrange(state, week_start)

if (anyNA(climate$Tmean) || anyNA(climate$PRCP)) stop("STOP: NA in national climate file -- resolve before proceeding.")

# ---- Rolling-window statistic: stat(x[t+offsets]) per state, NA if any offset is out of range ----
lagged_window_stat <- function(x, offsets, stat = c("mean", "sum")) {
  stat <- match.arg(stat)
  n <- length(x)
  out <- rep(NA_real_, n)
  for (t in seq_len(n)) {
    idx <- t + offsets
    if (all(idx >= 1L & idx <= n)) {
      out[t] <- if (stat == "mean") mean(x[idx]) else sum(x[idx])
    }
  }
  out
}

# ================================================================
# Build lagged climate + anomalies, per state per lag spec
# ================================================================
build_one_spec <- function(spec_name, spec) {
  bind_rows(lapply(split(climate, climate$state), function(d) {
    d <- d |> arrange(week_start)
    lagged_temp <- lagged_window_stat(d$Tmean, spec$temp_offsets, "mean")
    lagged_precip <- lagged_window_stat(d$PRCP, spec$precip_offsets, "sum")
    tibble(state = d$state[1], week_start = d$week_start, week_of_year = d$week_of_year,
           lag_spec = spec_name, lagged_temp = lagged_temp, lagged_precip = lagged_precip,
           log_precip = log1p(lagged_precip))
  }))
}

lagged_all <- bind_rows(lapply(names(LAG_SPECS), function(nm) build_one_spec(nm, LAG_SPECS[[nm]])))
message(sprintf("Built lagged climate for %d states x %d lag specs x %d weeks = %d rows (%d with complete lag window)",
                 n_distinct(lagged_all$state), length(LAG_SPECS), n_distinct(lagged_all$week_start), nrow(lagged_all),
                 sum(!is.na(lagged_all$lagged_temp) & !is.na(lagged_all$log_precip))))

# ---- Climatology: mean lagged_temp / log_precip by (state, week_of_year, lag_spec), across all available years ----
climatology <- lagged_all |> dplyr::filter(!is.na(lagged_temp), !is.na(log_precip)) |>
  group_by(state, week_of_year, lag_spec) |>
  summarise(T_clim = mean(lagged_temp), logP_clim = mean(log_precip), n_years = n(), .groups = "drop")

anomalies <- lagged_all |>
  left_join(climatology, by = c("state", "week_of_year", "lag_spec")) |>
  mutate(T_anom = lagged_temp - T_clim, P_anom = log_precip - logP_clim)

# ---- Standardisation: ONE common mean/SD per lag spec, pooled across the 4 primary states ----
scaling <- anomalies |> dplyr::filter(state %in% PRIMARY_STATES, !is.na(T_anom), !is.na(P_anom)) |>
  group_by(lag_spec) |>
  summarise(states_pooled = paste(PRIMARY_STATES, collapse = ";"),
            T_anom_mean = mean(T_anom), T_anom_sd = sd(T_anom),
            P_anom_mean = mean(P_anom), P_anom_sd = sd(P_anom), n_obs_pooled = n(), .groups = "drop")

anomalies <- anomalies |> left_join(scaling |> select(lag_spec, T_anom_mean, T_anom_sd, P_anom_mean, P_anom_sd), by = "lag_spec") |>
  mutate(z_T_anom = (T_anom - T_anom_mean) / T_anom_sd, z_P_anom = (P_anom - P_anom_mean) / P_anom_sd)

message("\n=== Standardisation (pooled across CE/BA/RJ/MT, per lag spec) ===")
print(as.data.frame(scaling), digits = 4)

# ---- Save everything ----
write_csv(lagged_all, file.path(table_dir, "climate_lagged_raw.csv"))
write_csv(climatology, file.path(table_dir, "climate_anomaly_climatology.csv"))
write_csv(anomalies, file.path(table_dir, "climate_anomaly_covariates.csv"))
write_csv(scaling, file.path(table_dir, "climate_anomaly_scaling_metadata.csv"))
message("\n[saved] climate_lagged_raw.csv, climate_anomaly_climatology.csv, climate_anomaly_covariates.csv, climate_anomaly_scaling_metadata.csv")

# ================================================================
# QC checks
# ================================================================
message("\n=== QC: mean/SD of z-scored anomalies (pooled states, should be ~0 / ~1 by construction) ===")
qc_pooled <- anomalies |> dplyr::filter(state %in% PRIMARY_STATES) |> group_by(lag_spec) |>
  summarise(mean_zT = mean(z_T_anom, na.rm = TRUE), sd_zT = sd(z_T_anom, na.rm = TRUE),
            mean_zP = mean(z_P_anom, na.rm = TRUE), sd_zP = sd(z_P_anom, na.rm = TRUE), .groups = "drop")
print(as.data.frame(qc_pooled), digits = 3)

message("\n=== QC: residual seasonal pattern in anomalies (correlation of anomaly with week_of_year, primary spec, primary states) ===")
qc_seasonal <- anomalies |> dplyr::filter(lag_spec == "primary", state %in% PRIMARY_STATES, !is.na(T_anom), !is.na(P_anom)) |>
  group_by(state) |>
  summarise(cor_Tanom_woy = cor(T_anom, week_of_year, use = "complete.obs"),
            cor_Panom_woy = cor(P_anom, week_of_year, use = "complete.obs"), .groups = "drop")
print(as.data.frame(qc_seasonal), digits = 3)
write_csv(qc_seasonal, file.path(table_dir, "climate_anomaly_seasonal_residual_check.csv"))

# ================================================================
# QC figures (primary lag spec, the 4 primary states)
# ================================================================
qc_data <- anomalies |> dplyr::filter(lag_spec == "primary", state %in% PRIMARY_STATES)

p_raw_vs_anom_T <- ggplot(qc_data, aes(week_start)) +
  geom_line(aes(y = lagged_temp, colour = "Raw (lagged mean T)"), linewidth = 0.3) +
  geom_line(aes(y = T_anom + mean(qc_data$lagged_temp, na.rm = TRUE), colour = "Anomaly (recentred for display)"), linewidth = 0.3) +
  facet_wrap(~ state, ncol = 1, scales = "free_y") +
  scale_colour_manual(name = NULL, values = c("Raw (lagged mean T)" = "#D55E00", "Anomaly (recentred for display)" = "#0072B2")) +
  labs(title = "QC A: raw lagged temperature vs. anomaly (primary lag spec)", x = NULL, y = "deg C") + theme_v4 + theme(legend.position = "top")
ggsave(file.path(figure_dir, "QC_A_raw_vs_anomaly_temperature.png"), p_raw_vs_anom_T, width = 200, height = 200, units = "mm", dpi = 300, bg = "white")

p_raw_vs_anom_P <- ggplot(qc_data, aes(week_start)) +
  geom_line(aes(y = lagged_precip, colour = "Raw (lagged sum PRCP)"), linewidth = 0.3) +
  geom_line(aes(y = exp(P_anom + mean(qc_data$logP_clim, na.rm = TRUE)) - 1, colour = "Anomaly (back-transformed for display)"), linewidth = 0.3) +
  facet_wrap(~ state, ncol = 1, scales = "free_y") +
  scale_colour_manual(name = NULL, values = c("Raw (lagged sum PRCP)" = "#D55E00", "Anomaly (back-transformed for display)" = "#0072B2")) +
  labs(title = "QC B: raw lagged precipitation vs. anomaly (primary lag spec)", x = NULL, y = "mm") + theme_v4 + theme(legend.position = "top")
ggsave(file.path(figure_dir, "QC_B_raw_vs_anomaly_precipitation.png"), p_raw_vs_anom_P, width = 200, height = 200, units = "mm", dpi = 300, bg = "white")

p_dist <- ggplot(qc_data |> pivot_longer(c(z_T_anom, z_P_anom), names_to = "variable", values_to = "value"), aes(value, fill = state)) +
  geom_density(alpha = 0.4) + facet_wrap(~ variable, scales = "free") +
  labs(title = "QC C: standardised anomaly distributions by state (primary lag spec)", x = NULL, y = "Density") + theme_v4
ggsave(file.path(figure_dir, "QC_C_anomaly_distributions_by_state.png"), p_dist, width = 200, height = 110, units = "mm", dpi = 300, bg = "white")

p_seasonal_check <- ggplot(qc_data, aes(week_of_year, T_anom, colour = state)) +
  geom_point(alpha = 0.3, size = 0.6) + geom_smooth(se = FALSE, linewidth = 0.6, method = "loess", span = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(title = "QC D: temperature anomaly vs. epidemiological week (should show no residual seasonal pattern)", x = "ISO week of year", y = "T_anom (deg C)") + theme_v4
ggsave(file.path(figure_dir, "QC_D_anomaly_vs_epiweek.png"), p_seasonal_check, width = 180, height = 120, units = "mm", dpi = 300, bg = "white")

message("\n[saved] 4 QC figures to ", figure_dir)
message("\n=== Phase 1 complete. STOP for review before Phase 2. ===")
