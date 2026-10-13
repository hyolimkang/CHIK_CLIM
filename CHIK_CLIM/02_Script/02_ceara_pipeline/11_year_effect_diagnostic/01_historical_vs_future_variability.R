# Year-effect diagnostic -- Part 2: historical posterior vs. future Y1
# sampling distribution for year_effect and A_year.
#
# Re-analyses the ALREADY-FITTED accepted climate-forced v4.9 canary
# posterior. No Stan refit. Distinguishes POSTERIOR UNCERTAINTY (how
# uncertain are we about one specific historical year's value, across
# draws) from GENUINE INTERANNUAL VARIABILITY (how much do years differ
# from each other WITHIN a single posterior draw) -- these are different
# quantities and are reported separately throughout.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/11_year_effect_diagnostic/00_config_year_effect_diagnostic.R"))

message("[Diag 1] Loading accepted climate-forced fit (read-only; no refit)...")
climate <- readRDS(ACCEPTED_CLIMATE_FIT)
sd_clim <- climate$stan_data
dates_hist <- as.Date(climate$weekly_data$week_start)
N_hist <- length(dates_hist)
Y_hist <- sd_clim$Y
hist_years <- sort(unique(as.integer(format(dates_hist, "%Y"))))
year_id_hist <- sd_clim$year_id
stopifnot(length(hist_years) == Y_hist)

draws <- rstan::extract(climate$fit, pars = c(
  "alpha_R", "z_year", "year_effect", "A_year", "sigma_season_year", "z_season_year",
  "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2", "beta_T", "beta_P", "R0_t"
), permuted = TRUE)
n_draws_hist <- length(draws$alpha_R)
message(sprintf("[Diag 1] %d posterior draws extracted.", n_draws_hist))

YEAR_EFFECT_PRIOR_SD <- 0.40 # DATA value in the Stan model -- NOT estimated (see Section 3 below)

# ================================================================
# 2.A Historical POSTERIOR UNCERTAINTY per year (across draws, for each
# specific calendar year) -- "how uncertain are we about 2017's value?"
# ================================================================
message("[Diag 1] Per-year posterior summaries (year_effect, A_year, z_year)...")
per_year_summary <- tibble(
  year = hist_years,
  z_year_median = apply(draws$z_year, 2, median), z_year_lo95 = apply(draws$z_year, 2, quantile, .025), z_year_hi95 = apply(draws$z_year, 2, quantile, .975),
  z_year_posterior_sd = apply(draws$z_year, 2, sd), # posterior SD of z_year[y] across draws -- compare to prior SD=1
  year_effect_median = apply(draws$year_effect, 2, median), year_effect_lo95 = apply(draws$year_effect, 2, quantile, .025), year_effect_hi95 = apply(draws$year_effect, 2, quantile, .975),
  A_year_median = apply(draws$A_year, 2, median), A_year_lo95 = apply(draws$A_year, 2, quantile, .025), A_year_hi95 = apply(draws$A_year, 2, quantile, .975)
)
print(as.data.frame(per_year_summary), digits = 3)
write_csv(per_year_summary, file.path(DIR_TABLES_DIAG, "DIAG_historical_per_year_posterior_summary.csv"))

sigma_season_year_summary <- tibble(parameter = "sigma_season_year", median = median(draws$sigma_season_year),
                                     lo95 = quantile(draws$sigma_season_year, .025), hi95 = quantile(draws$sigma_season_year, .975),
                                     prior = "Normal(0, 0.20), truncated at 0 (lower=0 constraint) -- a GENUINE posterior-estimated parameter")
write_csv(sigma_season_year_summary, file.path(DIR_TABLES_DIAG, "DIAG_sigma_season_year_posterior.csv"))
message(sprintf("sigma_season_year posterior: median=%.4f, 95%% CrI [%.4f, %.4f] (prior SD=0.20)",
                 sigma_season_year_summary$median, sigma_season_year_summary$lo95, sigma_season_year_summary$hi95))

# ================================================================
# 2.B GENUINE INTERANNUAL VARIABILITY implied by the fitted model --
# within EACH draw, how much do the 11 years differ from one another?
# (NOT a single representative draw -- the full distribution across all
# n_draws_hist draws is reported.)
# ================================================================
message("[Diag 1] Within-draw (across-year) spread -- genuine interannual variability...")
within_draw_sd_year_effect <- apply(draws$year_effect, 1, sd)
within_draw_range_year_effect <- apply(draws$year_effect, 1, function(x) diff(range(x)))
within_draw_sd_A_year <- apply(draws$A_year, 1, sd)
within_draw_range_A_year <- apply(draws$A_year, 1, function(x) diff(range(x)))

# Theoretical prior-implied value for comparison: SD of (z - mean(z)) for
# z ~ iid N(0,1), Y_hist draws -- a finite-sample correction factor
# sqrt((Y-1)/Y) below 1, scaled by year_effect_prior_sd.
prior_implied_sd_year_effect <- YEAR_EFFECT_PRIOR_SD * sqrt((Y_hist - 1) / Y_hist)

historical_spread_summary <- tibble(
  quantity = c("year_effect (within-draw SD across 11 years)", "year_effect (within-draw range across 11 years)",
               "A_year (within-draw SD across 11 years)", "A_year (within-draw range across 11 years)"),
  median = c(median(within_draw_sd_year_effect), median(within_draw_range_year_effect), median(within_draw_sd_A_year), median(within_draw_range_A_year)),
  lo95 = c(quantile(within_draw_sd_year_effect, .025), quantile(within_draw_range_year_effect, .025), quantile(within_draw_sd_A_year, .025), quantile(within_draw_range_A_year, .025)),
  hi95 = c(quantile(within_draw_sd_year_effect, .975), quantile(within_draw_range_year_effect, .975), quantile(within_draw_sd_A_year, .975), quantile(within_draw_range_A_year, .975))
)
print(as.data.frame(historical_spread_summary), digits = 4)
write_csv(historical_spread_summary, file.path(DIR_TABLES_DIAG, "DIAG_historical_within_draw_spread.csv"))
message(sprintf("Prior-implied year_effect within-draw SD (if z_year were pure iid N(0,1)): %.4f. Fitted historical median: %.4f (ratio = %.2f).",
                 prior_implied_sd_year_effect, median(within_draw_sd_year_effect), median(within_draw_sd_year_effect) / prior_implied_sd_year_effect))

# ================================================================
# 2.C Historical annual mean/peak R0, per draw x year
# ================================================================
message("[Diag 1] Historical annual mean/peak R0 per draw x year...")
annual_R0_hist <- bind_rows(lapply(seq_len(n_draws_hist), function(d) {
  tibble(draw = d, year = hist_years, year_id = seq_along(hist_years)) |>
    mutate(mean_R0 = vapply(year_id, function(y) mean(draws$R0_t[d, year_id_hist == y]), numeric(1)),
           peak_R0 = vapply(year_id, function(y) max(draws$R0_t[d, year_id_hist == y]), numeric(1)))
}))
write_csv(annual_R0_hist |> sample_frac(min(1, 20000 / nrow(annual_R0_hist))), file.path(DIR_RESULTS_DIAG, "DIAG_historical_annual_R0_sample.csv"))
annual_R0_hist_summary <- annual_R0_hist |> group_by(year) |>
  summarise(mean_R0_median = median(mean_R0), peak_R0_median = median(peak_R0), peak_R0_lo95 = quantile(peak_R0, .025), peak_R0_hi95 = quantile(peak_R0, .975), .groups = "drop")
print(as.data.frame(annual_R0_hist_summary), digits = 3)
write_csv(annual_R0_hist_summary, file.path(DIR_TABLES_DIAG, "DIAG_historical_annual_R0_summary.csv"))
historical_peak_R0_overall <- quantile(annual_R0_hist$peak_R0, c(.50, .90, .95, .99, 1.0))
message("Historical annual peak R0 overall (all draws x years): ", paste(sprintf("%s=%.2f", names(historical_peak_R0_overall), historical_peak_R0_overall), collapse = ", "))

# ================================================================
# 2.D Future Y1 sampling distribution -- SAME spread metrics, computed
# from a large number of independent future-year-structure realisations
# (not simulated epidemics yet -- purely the R0-construction mechanism).
# One realisation = one posterior draw's sigma_season_year PLUS 10 fresh
# future years' z_year/z_season (the model's OWN documented distributional
# assumption, exactly as implemented in sample_future_year_structure()).
# ================================================================
message("[Diag 1] Future Y1 sampling distribution (10,000 realisations of the future-year mechanism)...")
n_future_years <- 10L
N_REALISATIONS <- 10000L
set.seed(DIAG_BASE_SEED)
draw_for_realisation <- sample.int(n_draws_hist, N_REALISATIONS, replace = TRUE)
future_spread <- bind_rows(lapply(seq_len(N_REALISATIONS), function(i) {
  d <- draw_for_realisation[i]
  z_year_f <- rnorm(n_future_years); z_season_f <- rnorm(n_future_years)
  year_effect_f <- YEAR_EFFECT_PRIOR_SD * (z_year_f - mean(z_year_f))
  A_year_f <- exp(draws$sigma_season_year[d] * z_season_f)
  tibble(realisation = i, draw = d,
         sd_year_effect = sd(year_effect_f), range_year_effect = diff(range(year_effect_f)),
         sd_A_year = sd(A_year_f), range_A_year = diff(range(A_year_f)),
         peak_year_effect = max(year_effect_f), peak_A_year = max(A_year_f))
}))
write_csv(future_spread, file.path(DIR_RESULTS_DIAG, "DIAG_future_Y1_spread_realisations.csv"))

future_spread_summary <- tibble(
  quantity = c("year_effect (within-realisation SD, 10 future years)", "year_effect (within-realisation range)",
               "A_year (within-realisation SD)", "A_year (within-realisation range)"),
  median = c(median(future_spread$sd_year_effect), median(future_spread$range_year_effect), median(future_spread$sd_A_year), median(future_spread$range_A_year)),
  lo95 = c(quantile(future_spread$sd_year_effect, .025), quantile(future_spread$range_year_effect, .025), quantile(future_spread$sd_A_year, .025), quantile(future_spread$range_A_year, .025)),
  hi95 = c(quantile(future_spread$sd_year_effect, .975), quantile(future_spread$range_year_effect, .975), quantile(future_spread$sd_A_year, .975), quantile(future_spread$range_A_year, .975))
)
print(as.data.frame(future_spread_summary), digits = 4)
write_csv(future_spread_summary, file.path(DIR_TABLES_DIAG, "DIAG_future_Y1_spread_summary.csv"))

ratio_table <- tibble(
  quantity = c("year_effect SD (future / historical)", "year_effect range (future / historical)",
               "A_year SD (future / historical)", "A_year range (future / historical)"),
  ratio_of_medians = c(median(future_spread$sd_year_effect) / median(within_draw_sd_year_effect),
                        median(future_spread$range_year_effect) / median(within_draw_range_year_effect),
                        median(future_spread$sd_A_year) / median(within_draw_sd_A_year),
                        median(future_spread$range_A_year) / median(within_draw_range_A_year))
)
print(as.data.frame(ratio_table), digits = 3)
write_csv(ratio_table, file.path(DIR_TABLES_DIAG, "DIAG_future_vs_historical_spread_ratios.csv"))

# ================================================================
# 2.E Frequency of extreme future R0 -- using the REAL posterior draw's
# alpha_R/seasonal coefficients plus the resampled year_effect/A_year,
# under C0 (climate-neutral, matching the isolation design used
# throughout this diagnostic) so the comparison is climate-free on BOTH
# sides as far as possible (historical peak R0 above still includes the
# fitted climate term, which Step 1 showed contributes ~1% of variance --
# noted explicitly as a minor asymmetry, not corrected away).
# ================================================================
message("[Diag 1] Frequency of extreme future peak R0 under Y1 (C0 climate-neutral)...")
future_dates <- build_future_dates()
t_index_future <- (N_hist + 1):(N_hist + length(future_dates))
years_future <- as.integer(format(future_dates, "%Y"))
future_year_labels <- sort(unique(years_future))
year_id_future <- match(years_future, future_year_labels)

set.seed(DIAG_BASE_SEED + 1L)
N_R0_REAL <- 2000L
draw_for_r0 <- sample.int(n_draws_hist, N_R0_REAL, replace = TRUE)
future_peak_R0 <- vapply(seq_len(N_R0_REAL), function(i) {
  d <- draw_for_r0[i]
  z_year_f <- rnorm(n_future_years); z_season_f <- rnorm(n_future_years)
  year_effect_f <- YEAR_EFFECT_PRIOR_SD * (z_year_f - mean(z_year_f))
  A_year_f <- exp(draws$sigma_season_year[d] * z_season_f)
  r0 <- compute_R0_future(t_index_future, year_id_future, draws$alpha_R[d], year_effect_f, A_year_f,
                           draws$beta_sin1[d], draws$beta_cos1[d], draws$beta_sin2[d], draws$beta_cos2[d],
                           0, 0, rep(0, length(future_dates)), rep(0, length(future_dates)))
  max(r0$R0)
}, numeric(1))

extreme_freq <- tibble(
  threshold_label = c("historical median annual peak", "historical 90th pct annual peak", "historical 95th pct annual peak", "historical 99th pct annual peak", "historical maximum observed/fitted peak"),
  threshold_value = as.numeric(historical_peak_R0_overall[c("50%", "90%", "95%", "99%", "100%")]),
  fraction_future_years_exceeding = vapply(as.numeric(historical_peak_R0_overall[c("50%", "90%", "95%", "99%", "100%")]), function(th) mean(future_peak_R0 > th), numeric(1))
)
print(as.data.frame(extreme_freq), digits = 3)
write_csv(extreme_freq, file.path(DIR_TABLES_DIAG, "DIAG_extreme_future_R0_frequency.csv"))
write_csv(tibble(realisation = seq_len(N_R0_REAL), draw = draw_for_r0, peak_R0 = future_peak_R0), file.path(DIR_RESULTS_DIAG, "DIAG_future_peak_R0_realisations.csv"))

message(sprintf("\n[Diag 1] Future Y1 peak-R0 realisations exceeding the historical maximum fitted annual peak (%.2f): %.1f%%",
                 as.numeric(historical_peak_R0_overall["100%"]), 100 * mean(future_peak_R0 > as.numeric(historical_peak_R0_overall["100%"]))))
message("[Diag 1] Historical vs. future variability analysis complete. Tables -> ", DIR_TABLES_DIAG)
