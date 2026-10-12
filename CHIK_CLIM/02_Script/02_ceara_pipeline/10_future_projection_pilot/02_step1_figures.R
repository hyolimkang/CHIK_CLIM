# STEP 1 FIGURES -- FIG01-04.
#
# Reads the tables/results written by 01_step1_climate_posterior_reanalysis.R
# (run that script first) plus the two accepted fits directly for PPC.
# Colour convention: reuses this repo's existing Okabe-Ito-derived
# colourblind-safe qualitative palette (e.g. #0072B2/#D55E00/#E69F00), the
# same one used throughout every other figure in 02_ceara_pipeline, for
# visual consistency and because it is a well-documented CVD-safe set.
# (Node.js, needed for this session's automated palette validator, is not
# installed on this machine -- Okabe-Ito's CVD safety is long-established
# in the literature and is this repo's existing convention, used here
# unchanged rather than substituted.)

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "ggplot2", "patchwork", "scales", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(ggplot2); library(patchwork) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/10_future_projection_pilot/00_config_future_pilot.R"))

OKABE <- c(blue = "#0072B2", orange = "#D55E00", yellow = "#E69F00", green = "#009E73", purple = "#CC79A7", grey = "#56B4E9")
theme_pub <- theme_classic(base_size = 9) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey30"), legend.position = "top")

save_fig <- function(plot, name, width_mm = 180, height_mm = 120) {
  ggsave(file.path(DIR_FIGURES, paste0(name, ".pdf")), plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(file.path(DIR_FIGURES, paste0(name, ".png")), plot, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
  ggsave(file.path(DIR_FIGURES, paste0(name, ".svg")), plot, width = width_mm, height = height_mm, units = "mm")
  message("[saved] ", name, " (.pdf/.png@600dpi/.svg)")
}

NOTE_CAPTION <- NULL # Step 1 figures describe the already-fitted historical model; the Step-2 scenario note is not applicable here.

# ================================================================
# FIG01 -- Climate coefficient inference: posterior vs prior, beta_T & beta_P
# ================================================================
beta_draws <- readRDS(file.path(DIR_RESULTS, "step1_beta_and_multiplier_draws.rds"))
beta_summary <- read_csv(file.path(DIR_TABLES, "STEP1_beta_posterior_vs_prior_summary.csv"), show_col_types = FALSE) |>
  rename(param = parameter)

prior_x <- seq(-0.6, 0.6, length.out = 400)
prior_df <- tibble(param = rep(c("beta_T", "beta_P"), each = length(prior_x)),
                    x = rep(prior_x, 2), density = dnorm(rep(prior_x, 2), 0, beta_draws$prior_sd))
post_df <- bind_rows(tibble(param = "beta_T", value = beta_draws$beta_T), tibble(param = "beta_P", value = beta_draws$beta_P))
param_labels <- c(beta_T = "beta[T]~(temperature)", beta_P = "beta[P]~(precipitation)")

fig01 <- ggplot() +
  geom_area(data = prior_df, aes(x, density), fill = "grey70", alpha = 0.35) +
  geom_line(data = prior_df, aes(x, density), colour = "grey45", linetype = "dashed", linewidth = 0.4) +
  geom_density(data = post_df, aes(value, fill = param, colour = param), alpha = 0.45, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "black", linewidth = 0.4) +
  geom_segment(data = beta_summary, aes(x = lo95, xend = hi95, y = -0.4, yend = -0.4), colour = "black", linewidth = 0.8) +
  geom_point(data = beta_summary, aes(x = posterior_median, y = -0.4), size = 1.6) +
  facet_wrap(~ param, scales = "free", labeller = as_labeller(param_labels, label_parsed)) +
  scale_fill_manual(values = c(beta_T = OKABE[["orange"]], beta_P = OKABE[["blue"]]), guide = "none") +
  scale_colour_manual(values = c(beta_T = OKABE[["orange"]], beta_P = OKABE[["blue"]]), guide = "none") +
  geom_text(data = beta_summary, aes(x = posterior_median, y = Inf,
            label = sprintf("median=%.3f\n95%% CrI [%.3f, %.3f]\nP(>0)=%.0f%%", posterior_median, lo95, hi95, 100 * prob_positive)),
            vjust = 1.1, size = 2.6, colour = "grey20") +
  labs(title = "FIG01: Climate coefficient inference (Ceara v4.9 climate-forced canary, q=0.05)",
       subtitle = "Grey: prior N(0, 0.15). Coloured: posterior density. Bar: posterior median and 95% credible interval. Dotted line: zero effect.",
       x = "Coefficient value (log R0 scale)", y = "Density") +
  theme_pub
save_fig(fig01, "FIG01_climate_coefficient_inference", width_mm = 180, height_mm = 95)

# ================================================================
# FIG02 -- Historical climate contribution: anomalies + posterior multiplier
# ================================================================
climate <- readRDS(ACCEPTED_CLIMATE_FIT)
dates <- as.Date(climate$weekly_data$week_start)
anom_panel <- tibble(week_start = dates, z_T_anom = climate$stan_data$z_T_anom, z_P_anom = climate$stan_data$z_P_anom) |>
  tidyr::pivot_longer(c(z_T_anom, z_P_anom), names_to = "variable", values_to = "value") |>
  mutate(variable = recode(variable, z_T_anom = "Temperature anomaly (z-score)", z_P_anom = "Precipitation anomaly (z-score)"))
mult_summary <- read_csv(file.path(DIR_TABLES, "STEP1_climate_multiplier_weekly.csv"), show_col_types = FALSE)

p_anom <- ggplot(anom_panel, aes(week_start, value, colour = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.3) +
  scale_colour_manual(name = NULL, values = c("Temperature anomaly (z-score)" = OKABE[["orange"]], "Precipitation anomaly (z-score)" = OKABE[["blue"]])) +
  labs(subtitle = "A. Standardised climate anomalies (primary lag spec) used by the fit", x = NULL, y = "Standardised anomaly (z)") +
  theme_pub
p_mult <- ggplot(mult_summary, aes(week_start)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = OKABE[["yellow"]], alpha = 0.3) +
  geom_line(aes(y = median), colour = OKABE[["yellow"]], linewidth = 0.6) +
  labs(subtitle = "B. Posterior climate multiplier exp(beta_T*z_T + beta_P*z_P), with 95% credible interval", x = NULL, y = "Multiplier on R0(t)") +
  theme_pub
fig02 <- p_anom / p_mult + plot_annotation(
  title = "FIG02: Historical climate contribution, Ceara 2015-2025",
  subtitle = "Temporal alignment: both panels share the same weekly x-axis.")
save_fig(fig02, "FIG02_historical_climate_contribution", width_mm = 180, height_mm = 130)

# ================================================================
# FIG03 -- Baseline vs climate-forced model: PPC, R0(t), S_prop(t)
# ================================================================
baseline <- readRDS(ACCEPTED_BASELINE_FIT)
pred_base <- rstan::extract(baseline$fit, pars = "C_pred")$C_pred
pred_clim <- rstan::extract(climate$fit, pars = "C_pred")$C_pred
obs <- climate$weekly_data$cases
years_vec <- as.integer(format(dates, "%Y"))
annual_obs <- tapply(obs, years_vec, sum)
annual_pred_base <- apply(t(apply(pred_base, 1, function(r) tapply(r, years_vec, sum))), 2, median)
annual_pred_clim <- apply(t(apply(pred_clim, 1, function(r) tapply(r, years_vec, sum))), 2, median)
ppc_df <- tibble(year = as.integer(names(annual_obs)), observed = as.numeric(annual_obs),
                  `baseline v4.9` = annual_pred_base, `climate-forced` = annual_pred_clim) |>
  tidyr::pivot_longer(c(`baseline v4.9`, `climate-forced`), names_to = "model", values_to = "predicted_median")

p_ppc <- ggplot(ppc_df, aes(year)) +
  geom_col(aes(y = observed), fill = "grey80", width = 0.6) +
  geom_point(aes(y = predicted_median, colour = model), size = 1.8) +
  geom_line(aes(y = predicted_median, colour = model, group = model), linewidth = 0.4) +
  scale_colour_manual(name = NULL, values = c("baseline v4.9" = OKABE[["blue"]], "climate-forced" = OKABE[["orange"]])) +
  scale_x_continuous(breaks = unique(ppc_df$year)) +
  labs(subtitle = "A. Annual posterior-predictive median vs. observed reported cases (bars)", x = NULL, y = "Reported cases / year") +
  theme_pub

draws_base_r0 <- rstan::extract(baseline$fit, pars = "R0_t")$R0_t
draws_clim_r0 <- rstan::extract(climate$fit, pars = "R0_t")$R0_t
r0_df <- tibble(week_start = dates, `baseline v4.9` = apply(draws_base_r0, 2, median), `climate-forced` = apply(draws_clim_r0, 2, median)) |>
  tidyr::pivot_longer(-week_start, names_to = "model", values_to = "R0")
p_r0 <- ggplot(r0_df, aes(week_start, R0, colour = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.3, alpha = 0.85) +
  scale_colour_manual(name = NULL, values = c("baseline v4.9" = OKABE[["blue"]], "climate-forced" = OKABE[["orange"]])) +
  labs(subtitle = "B. Posterior median R0(t)", x = NULL, y = "R0(t)") + theme_pub

s_compare <- read_csv(file.path(DIR_TABLES, "STEP1_S_prop_full_series_comparison.csv"), show_col_types = FALSE)
p_s <- ggplot(s_compare, aes(week_start)) +
  geom_line(aes(y = S_prop_baseline_median, colour = "baseline v4.9"), linewidth = 0.4) +
  geom_line(aes(y = S_prop_climate_median, colour = "climate-forced"), linewidth = 0.4) +
  scale_colour_manual(name = NULL, values = c("baseline v4.9" = OKABE[["blue"]], "climate-forced" = OKABE[["orange"]])) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(subtitle = sprintf("C. Susceptible proportion S(t)/N(t) (max relative difference across all %d weeks: %.2f%%)", nrow(s_compare), max(s_compare$rel_diff_pct)),
       x = NULL, y = "S / N") + theme_pub

fig03 <- p_ppc / p_r0 / p_s + plot_annotation(title = "FIG03: Baseline (non-climate) vs. climate-forced v4.9 -- reconstruction comparison")
save_fig(fig03, "FIG03_baseline_vs_climate_model", width_mm = 180, height_mm = 190)

# ================================================================
# FIG04 -- Transmission decomposition of log_R0(t)
# ================================================================
components <- read_csv(file.path(DIR_TABLES, "STEP1_log_R0_weekly_component_decomposition.csv"), show_col_types = FALSE)
decomp_long <- components |>
  select(week_start, Seasonal = seasonal, `Year effect` = year_effect_component, Climate = climate_component) |>
  tidyr::pivot_longer(-week_start, names_to = "component", values_to = "value") |>
  mutate(component = factor(component, levels = c("Seasonal", "Year effect", "Climate")))

decomp_summary <- read_csv(file.path(DIR_TABLES, "STEP1_variance_decomposition_summary.csv"), show_col_types = FALSE)
share_lab <- sprintf("Posterior median variance share: seasonal=%.0f%%, year=%.0f%%, climate=%.0f%%",
                      100 * decomp_summary$share_seasonal_median, 100 * decomp_summary$share_year_median, 100 * decomp_summary$share_climate_median)

p_decomp <- ggplot(decomp_long, aes(week_start, value, colour = component)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.3) +
  facet_wrap(~ component, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(Seasonal = OKABE[["blue"]], `Year effect` = OKABE[["green"]], Climate = OKABE[["orange"]]), guide = "none") +
  labs(title = "FIG04: Descriptive decomposition of log R0(t) - alpha_R (posterior-median-parameter draw)",
       subtitle = paste0(share_lab, "\nDescriptive mathematical decomposition only -- NOT a causal attribution of transmission drivers."),
       x = NULL, y = "Contribution to log R0(t)") +
  theme_pub
save_fig(p_decomp, "FIG04_transmission_decomposition", width_mm = 180, height_mm = 150)

message("\n[Step1 figures] FIG01-04 complete -> ", DIR_FIGURES)
