# STEP 2 FIGURES -- FIG05-09.
#
# Reads the reference-scenario trajectories/summaries and the sensitivity
# grid produced by 05_run_simulation_ensemble.R (run that first). Also
# computes a lightweight R0(t)-only sweep directly (cheap: no renewal
# simulation needed) for FIG05.
#
# MANDATORY per task instruction: every figure below carries the caption
# "Illustrative scenario simulation; no SSP forcing or vaccination applied."

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr", "ggplot2", "patchwork", "scales", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr); library(ggplot2); library(patchwork) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/10_future_projection_pilot/00_config_future_pilot.R"))
source(file.path(ROOT, "02_Script/02_ceara_pipeline/10_future_projection_pilot/03_future_renewal_functions.R"))

OKABE <- c(blue = "#0072B2", orange = "#D55E00", yellow = "#E69F00", green = "#009E73", purple = "#CC79A7", grey = "#56B4E9")
theme_pub <- theme_classic(base_size = 9) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey30"), legend.position = "top",
        plot.caption = element_text(size = 7, colour = "grey40", face = "italic"))
SCENARIO_NOTE <- "Illustrative scenario simulation; no SSP forcing or vaccination applied."

save_fig <- function(plot, name, width_mm = 180, height_mm = 120) {
  ggsave(file.path(DIR_FIGURES, paste0(name, ".pdf")), plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(file.path(DIR_FIGURES, paste0(name, ".png")), plot, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
  ggsave(file.path(DIR_FIGURES, paste0(name, ".svg")), plot, width = width_mm, height = height_mm, units = "mm")
  message("[saved] ", name, " (.pdf/.png@600dpi/.svg)")
}

message("[Step2 figures] Loading accepted fit + reference ensemble...")
climate <- readRDS(ACCEPTED_CLIMATE_FIT)
sd_clim <- climate$stan_data
dates_hist <- as.Date(climate$weekly_data$week_start)
N_hist <- length(dates_hist)
draws <- rstan::extract(climate$fit, pars = c("alpha_R", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2", "sigma_season_year", "beta_T", "beta_P"), permuted = TRUE)

future_dates <- build_future_dates()
n_future_weeks <- length(future_dates)
t_index_future <- (N_hist + 1):(N_hist + n_future_weeks)
years_future <- as.integer(format(future_dates, "%Y"))
future_year_labels <- sort(unique(years_future))
n_future_years <- length(future_year_labels)
year_id_future <- match(years_future, future_year_labels)
weeks_per_future_year <- as.integer(table(factor(years_future, levels = future_year_labels)))
hist_years <- as.integer(format(dates_hist, "%Y"))
hist_climate_by_year <- split(tibble(z_T_anom = sd_clim$z_T_anom, z_P_anom = sd_clim$z_P_anom), hist_years)

# ================================================================
# FIG05 -- Future transmission potential: R0(t) ensembles by year-effect x climate arm
# (R0(t) depends only on the posterior draw + arm construction, not on the
# full renewal simulation, so this is computed directly and cheaply here)
# ================================================================
set.seed(BASE_SEED + 100L)
n_draws_fig05 <- 30L
draw_ids_fig05 <- sample.int(length(draws$alpha_R), n_draws_fig05)
fig05_rows <- vector("list", 0)
for (year_arm in YEAR_EFFECT_ARMS) for (climate_arm in c("C0_stationary_neutral", "C1_resampled_weather_variability")) {
  for (d in draw_ids_fig05) {
    set.seed(BASE_SEED + d); year_struct <- sample_future_year_structure(year_arm, n_future_years, 0.40, draws$sigma_season_year[d], rng = function(n) rnorm(n))
    set.seed(BASE_SEED + d + 500L); climate_block <- sample_future_climate_block(climate_arm, future_year_labels, weeks_per_future_year, hist_climate_by_year, rng_discrete = function(n, ch) sample(ch, n, replace = TRUE))
    r0 <- compute_R0_future(t_index_future, year_id_future, draws$alpha_R[d], year_struct$year_effect, year_struct$A_year,
                             draws$beta_sin1[d], draws$beta_cos1[d], draws$beta_sin2[d], draws$beta_cos2[d], draws$beta_T[d], draws$beta_P[d], climate_block$z_T, climate_block$z_P)
    fig05_rows[[length(fig05_rows) + 1]] <- tibble(year_arm = year_arm, climate_arm = climate_arm, draw = d, week_start = future_dates, R0 = r0$R0)
  }
}
fig05_data <- bind_rows(fig05_rows) |>
  mutate(year_arm_label = recode(year_arm, Y0_reference = "Y0: reference (neutral) year effects", Y1_resampled_distributional = "Y1: resampled year effects"),
         climate_arm_label = recode(climate_arm, C0_stationary_neutral = "C0: climate-neutral", C1_resampled_weather_variability = "C1: resampled weather variability"))

fig05_summary <- fig05_data |> group_by(year_arm_label, climate_arm_label, week_start) |>
  summarise(median = median(R0), lo = quantile(R0, .1), hi = quantile(R0, .9), .groups = "drop")

fig05 <- ggplot() +
  geom_line(data = fig05_data, aes(week_start, R0, group = interaction(draw), colour = climate_arm_label), alpha = 0.08, linewidth = 0.2) +
  geom_line(data = fig05_summary, aes(week_start, median, colour = climate_arm_label), linewidth = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  facet_grid(year_arm_label ~ .) +
  scale_colour_manual(name = NULL, values = c("C0: climate-neutral" = OKABE[["blue"]], "C1: resampled weather variability" = OKABE[["orange"]])) +
  labs(title = "FIG05: Future transmission potential R0(t), 2026-2035", subtitle = "Thin lines: individual posterior draws (n=30). Thick line: ensemble median.",
       x = NULL, y = "R0(t)", caption = SCENARIO_NOTE) + theme_pub
save_fig(fig05, "FIG05_future_transmission_potential", width_mm = 180, height_mm = 140)

# ================================================================
# FIG06 -- Simulated epidemic trajectories (reference scenario, pilot ensemble)
# ================================================================
ref_traj <- readRDS(file.path(DIR_RESULTS, "STEP2_reference_scenario_trajectories.rds"))
weekly_pi <- ref_traj |> group_by(week_start) |>
  summarise(median = median(X_latent), lo25 = quantile(X_latent, .25), hi75 = quantile(X_latent, .75),
            lo025 = quantile(X_latent, .025), hi975 = quantile(X_latent, .975), .groups = "drop")
set.seed(BASE_SEED + 200L)
rep_keys <- ref_traj |> distinct(posterior_draw_id, stochastic_replicate_id) |> slice_sample(n = 18L)
rep_traj <- ref_traj |> inner_join(rep_keys, by = c("posterior_draw_id", "stochastic_replicate_id"))

fig06 <- ggplot() +
  geom_ribbon(data = weekly_pi, aes(week_start, ymin = lo025, ymax = hi975), fill = OKABE[["blue"]], alpha = 0.12) +
  geom_ribbon(data = weekly_pi, aes(week_start, ymin = lo25, ymax = hi75), fill = OKABE[["blue"]], alpha = 0.28) +
  geom_line(data = rep_traj, aes(week_start, X_latent, group = interaction(posterior_draw_id, stochastic_replicate_id)), colour = OKABE[["orange"]], alpha = 0.35, linewidth = 0.25) +
  geom_line(data = weekly_pi, aes(week_start, median), colour = OKABE[["blue"]], linewidth = 0.7) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1000), labels = scales::label_number(big.mark = ","), breaks = c(0, 1e3, 1e4, 1e5, 1e6)) +
  labs(title = "FIG06: Simulated weekly latent infections, 2026-2035",
       subtitle = sprintf("Reference scenario: Y1 x C1 x D1, lambda_import=%s. Orange: %d of %d simulated trajectories.\nBlue: ensemble median (line), 50%%/95%% prediction interval (bands). Pseudo-log y-axis.",
                           IMPORTATION_REFERENCE, nrow(rep_keys), n_distinct(ref_traj$posterior_draw_id) * n_distinct(ref_traj$stochastic_replicate_id)),
       x = NULL, y = "Weekly latent infections", caption = SCENARIO_NOTE) + theme_pub
save_fig(fig06, "FIG06_simulated_epidemic_trajectories", width_mm = 180, height_mm = 115)

# ================================================================
# FIG07 -- Epidemic mechanism: aligned R0/R_eff/S-fraction/infections for representative simulations
# ================================================================
cum_by_traj <- ref_traj |> group_by(posterior_draw_id, stochastic_replicate_id) |> summarise(cum10 = sum(X_latent), .groups = "drop") |> arrange(cum10)
pick <- cum_by_traj[round(c(0.1, 0.5, 0.9) * nrow(cum_by_traj)), ] |> mutate(label = c("Low-incidence realisation", "Median realisation", "High-incidence realisation"))
fig07_data <- ref_traj |> inner_join(pick, by = c("posterior_draw_id", "stochastic_replicate_id")) |>
  mutate(S_over_N = S / N, label = factor(label, levels = pick$label))

p_r0eff <- ggplot(fig07_data, aes(week_start)) +
  geom_line(aes(y = R0, colour = "R0(t)"), linewidth = 0.3) + geom_line(aes(y = R_eff, colour = "R_eff(t)"), linewidth = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  facet_wrap(~ label, ncol = 1) + scale_colour_manual(name = NULL, values = c("R0(t)" = OKABE[["grey"]], "R_eff(t)" = OKABE[["orange"]])) +
  labs(subtitle = "A. R0(t) vs. R_eff(t) = R0(t) x S(t)/N(t)", x = NULL, y = "Reproduction number") + theme_pub
p_s <- ggplot(fig07_data, aes(week_start, S_over_N)) + geom_line(colour = OKABE[["blue"]], linewidth = 0.4) +
  facet_wrap(~ label, ncol = 1) + scale_y_continuous(labels = scales::label_percent()) +
  labs(subtitle = "B. Susceptible fraction S(t)/N(t)", x = NULL, y = "S / N") + theme_pub
p_x <- ggplot(fig07_data, aes(week_start, X_latent)) + geom_line(colour = OKABE[["orange"]], linewidth = 0.3) +
  facet_wrap(~ label, ncol = 1) + scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(subtitle = "C. Weekly latent infections", x = NULL, y = "Infections / week") + theme_pub

fig07 <- p_r0eff | p_s | p_x
fig07 <- fig07 + plot_annotation(title = "FIG07: Epidemic mechanism in three representative future realisations",
                                  subtitle = "Reference scenario. Rows: low/median/high 10-year cumulative-incidence\nrealisation -- illustrates susceptible depletion and recurrence.",
                                  caption = SCENARIO_NOTE,
                                  theme = theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9, colour = "grey30"),
                                                plot.caption = element_text(size = 7, colour = "grey40", face = "italic")))
save_fig(fig07, "FIG07_epidemic_mechanism", width_mm = 190, height_mm = 150)

# ================================================================
# FIG08 -- Scenario sensitivity: cumulative infections by climate arm and importation pressure
# ================================================================
sens <- read_csv(file.path(DIR_RESULTS, "STEP2_sensitivity_grid_cumulative_summary.csv"), show_col_types = FALSE) |>
  mutate(climate_label = recode(climate_arm, C0_stationary_neutral = "C0: neutral", C1_resampled_weather_variability = "C1: resampled", C2_resampled_weather_zero_beta = "C2: resampled, beta=0"),
         year_label = recode(year_arm, Y0_reference = "Y0: reference", Y1_resampled_distributional = "Y1: resampled"))

p_climate <- sens |> dplyr::filter(lambda_import == IMPORTATION_REFERENCE) |>
  ggplot(aes(climate_label, cum_10y / 1e6, fill = year_label)) + geom_boxplot(outlier.size = 0.5, linewidth = 0.3) +
  scale_fill_manual(name = NULL, values = c("Y0: reference" = OKABE[["blue"]], "Y1: resampled" = OKABE[["orange"]])) +
  labs(subtitle = sprintf("A. By climate arm (importation fixed at reference lambda=%s)", IMPORTATION_REFERENCE), x = NULL, y = "10-year cumulative infections (millions)") + theme_pub
p_import <- sens |> dplyr::filter(climate_arm == "C1_resampled_weather_variability") |>
  ggplot(aes(factor(lambda_import), cum_10y / 1e6, fill = year_label)) + geom_boxplot(outlier.size = 0.5, linewidth = 0.3) +
  scale_fill_manual(name = NULL, values = c("Y0: reference" = OKABE[["blue"]], "Y1: resampled" = OKABE[["orange"]])) +
  labs(subtitle = "B. By weekly introduction pressure lambda_import (climate arm fixed at C1: resampled)", x = "lambda_import (weekly Poisson rate)", y = "10-year cumulative infections (millions)") + theme_pub

fig08 <- p_climate / p_import + plot_annotation(title = "FIG08: Scenario sensitivity -- 10-year cumulative infection distributions",
                                                 subtitle = "Each box: 10 posterior draws x 5 stochastic replicates (smoke-scale sensitivity grid).", caption = SCENARIO_NOTE)
save_fig(fig08, "FIG08_scenario_sensitivity", width_mm = 180, height_mm = 170)

# ================================================================
# FIG09 -- Five-year vs. ten-year cumulative incidence (reference scenario)
# ================================================================
ref_cum <- read_csv(file.path(DIR_RESULTS, "STEP2_reference_scenario_cumulative_summary.csv"), show_col_types = FALSE) |>
  select(posterior_draw_id, stochastic_replicate_id, `5-year` = incidence_per_100k_5y, `10-year` = incidence_per_100k_10y) |>
  pivot_longer(c(`5-year`, `10-year`), names_to = "horizon", values_to = "incidence_per_100k") |>
  mutate(horizon = factor(horizon, levels = c("5-year", "10-year")))

fig09 <- ggplot(ref_cum, aes(incidence_per_100k, fill = horizon)) +
  geom_density(alpha = 0.55, colour = NA) +
  scale_x_continuous(labels = scales::label_number(big.mark = ",")) +
  scale_fill_manual(name = "Horizon", values = c("5-year" = OKABE[["blue"]], "10-year" = OKABE[["orange"]])) +
  labs(title = "FIG09: 5-year vs. 10-year cumulative infection incidence (reference scenario, no vaccination)",
       subtitle = sprintf("Pilot ensemble: %d posterior draws x %d stochastic replicates. Structural scenario, not a forecast.",
                           n_distinct(ref_traj$posterior_draw_id), n_distinct(ref_traj$stochastic_replicate_id)),
       x = "Cumulative infection incidence (per 100,000 population)", y = "Density", caption = SCENARIO_NOTE) + theme_pub
save_fig(fig09, "FIG09_five_vs_ten_year_outcomes", width_mm = 170, height_mm = 100)

message("\n[Step2 figures] FIG05-09 complete -> ", DIR_FIGURES)
