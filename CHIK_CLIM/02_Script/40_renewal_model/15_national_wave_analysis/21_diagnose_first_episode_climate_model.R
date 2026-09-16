# Brazil first-epidemic climate-transmission pilot -- Section 9: model
# diagnostics -- eligible-UF/observation counts, lag-specification
# sensitivity (primary vs. sensitivity A/B, fixed grids per instruction),
# and leave-one-UF-out sensitivity (does a small number of UFs dominate
# the climate-response curves?).

required_packages <- c("here", "dplyr", "readr", "mgcv", "tibble", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(mgcv); library(tibble); library(ggplot2) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
source(file.path(root, "02_Script/40_renewal_model/15_national_wave_analysis/00_national_wave_analysis_helpers.R"))
paths <- ensure_national_output_dirs()
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
K_SMOOTH <- 4L

all_specs <- read_csv(file.path(paths$table, "FIRST_EPISODE_CLIMATE_ALL_SPECS.csv"), show_col_types = FALSE) |> mutate(UF = factor(UF))
training <- all_specs |> filter(lag_spec == "primary")

fit_gam <- function(d) mgcv::gam(log(Re_median) ~ s(temperature, k = K_SMOOTH) + s(precipitation, k = K_SMOOTH) + s(UF, bs = "re"), data = d, method = "REML")

# ================================================================
# Lag-specification sensitivity
# ================================================================
message("=== Lag-specification sensitivity (primary vs. sensitivity A/B) ===")
lag_results <- bind_rows(lapply(unique(all_specs$lag_spec), function(spec_name) {
  d <- all_specs |> filter(lag_spec == spec_name)
  fit <- fit_gam(d)
  sm <- summary(fit)
  tibble(lag_spec = spec_name, n = nrow(d),
         edf_temperature = sm$s.table["s(temperature)", "edf"], p_temperature = sm$s.table["s(temperature)", "p-value"],
         edf_precipitation = sm$s.table["s(precipitation)", "edf"], p_precipitation = sm$s.table["s(precipitation)", "p-value"],
         dev_explained = sm$dev.expl, adj_r2 = sm$r.sq)
}))
print(as.data.frame(lag_results), digits = 3)
write_csv(lag_results, file.path(paths$table, "first_epidemic_climate_lag_sensitivity.csv"))
message("[saved] ", file.path(paths$table, "first_epidemic_climate_lag_sensitivity.csv"))

# ================================================================
# Leave-one-UF-out sensitivity
# ================================================================
message("\n=== Leave-one-UF-out sensitivity ===")
ufs <- levels(training$UF)
temp_grid_1d <- seq(min(training$temperature), max(training$temperature), length.out = 50)
precip_med <- median(training$precipitation)

loo_curves <- bind_rows(lapply(ufs, function(excluded_uf) {
  d <- training |> filter(UF != excluded_uf) |> mutate(UF = droplevels(UF))
  fit <- tryCatch(fit_gam(d), error = function(e) NULL)
  if (is.null(fit)) return(tibble())
  pred <- predict(fit, newdata = data.frame(temperature = temp_grid_1d, precipitation = precip_med), exclude = "s(UF)", newdata.guaranteed = TRUE)
  tibble(excluded_uf = excluded_uf, temperature = temp_grid_1d, R0_climate_hat = exp(as.numeric(pred)))
}))
full_fit <- fit_gam(training)
full_curve <- tibble(excluded_uf = "(none - full model)",
                      temperature = temp_grid_1d,
                      R0_climate_hat = exp(as.numeric(predict(full_fit, newdata = data.frame(temperature = temp_grid_1d, precipitation = precip_med), exclude = "s(UF)", newdata.guaranteed = TRUE))))

p_loo <- ggplot(loo_curves, aes(temperature, R0_climate_hat, group = excluded_uf)) +
  geom_line(colour = "grey70", alpha = 0.7, linewidth = 0.4) +
  geom_line(data = full_curve, aes(temperature, R0_climate_hat), colour = "#D55E00", linewidth = 1.1, inherit.aes = FALSE) +
  labs(title = "Leave-one-UF-out sensitivity: temperature response curve", subtitle = "Grey = one UF excluded (18 refits); orange = full model. Precipitation held at observed median.",
       x = "Temperature (deg C)", y = expression(R[0]^{climate})) + theme_v4
ggsave(file.path(paths$figure, "FIRST_EPISODE_LOO_UF_SENSITIVITY.png"), p_loo, width = 180, height = 120, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(paths$figure, "FIRST_EPISODE_LOO_UF_SENSITIVITY.png"))

# Range of LOO curves relative to full-model curve, at each temperature -- how much does excluding one UF move the prediction?
loo_range <- loo_curves |> group_by(temperature) |> summarise(loo_min = min(R0_climate_hat), loo_max = max(R0_climate_hat), .groups = "drop") |>
  left_join(full_curve |> select(temperature, R0_climate_hat_full = R0_climate_hat), by = "temperature") |>
  mutate(pct_range = 100 * (loo_max - loo_min) / R0_climate_hat_full)
message(sprintf("LOO curve spread: median %.1f%% of full-model prediction, max %.1f%%", median(loo_range$pct_range), max(loo_range$pct_range)))

# Which single UF's exclusion moves the temperature-effect p-value/edf the most?
loo_gam_stats <- bind_rows(lapply(ufs, function(excluded_uf) {
  d <- training |> filter(UF != excluded_uf) |> mutate(UF = droplevels(UF))
  fit <- tryCatch(fit_gam(d), error = function(e) NULL)
  if (is.null(fit)) return(tibble())
  sm <- summary(fit)
  tibble(excluded_uf = excluded_uf, p_temperature = sm$s.table["s(temperature)", "p-value"], p_precipitation = sm$s.table["s(precipitation)", "p-value"], dev_explained = sm$dev.expl)
}))
message("\n=== Leave-one-UF-out: GAM significance stability ===")
print(as.data.frame(loo_gam_stats), digits = 3)
write_csv(loo_gam_stats, file.path(paths$table, "first_epidemic_climate_loo_uf_sensitivity.csv"))
message("[saved] ", file.path(paths$table, "first_epidemic_climate_loo_uf_sensitivity.csv"))

n_flip_precip <- sum((loo_gam_stats$p_precipitation < 0.05) != (summary(full_fit)$s.table["s(precipitation)", "p-value"] < 0.05))
message(sprintf("\nPrecipitation significance (p<0.05) flips sign of conclusion in %d / %d leave-one-out refits", n_flip_precip, nrow(loo_gam_stats)))

# ================================================================
# Combined observation/eligibility summary (Section 9 header stats)
# ================================================================
overview <- tibble(
  n_eligible_ufs = n_distinct(training$UF), n_weekly_re_observations = nrow(training),
  obs_per_uf = nrow(training) / n_distinct(training$UF),
  temperature_min = min(training$temperature), temperature_max = max(training$temperature),
  precipitation_min = min(training$precipitation), precipitation_max = max(training$precipitation),
  loo_median_pct_spread = median(loo_range$pct_range), loo_max_pct_spread = max(loo_range$pct_range),
  n_loo_precip_significance_flips = n_flip_precip
)
message("\n=== Section 9 overview ===")
print(as.data.frame(overview), digits = 4)
write_csv(overview, file.path(paths$table, "first_epidemic_climate_diagnostics_overview.csv"))
message("[saved] ", file.path(paths$table, "first_epidemic_climate_diagnostics_overview.csv"))
