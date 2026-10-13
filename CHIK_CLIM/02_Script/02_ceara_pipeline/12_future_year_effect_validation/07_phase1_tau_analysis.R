# PHASE 1 (analysis) -- tau_year posterior summaries across P1-P3.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/12_future_year_effect_validation/00_config_year_effect_validation.R"))

message("[Phase1 analysis] Loading P1-P3 tau_year fits (read-only)...")
bundles <- lapply(names(TAU_PRIORS), function(p) readRDS(file.path(DIR_FITS_YV, sprintf("ce_tau_year_%s.rds", p))))
names(bundles) <- names(TAU_PRIORS)
dates_hist <- as.Date(bundles$P2$weekly_data$week_start)
N_hist <- length(dates_hist)
hist_years <- sort(unique(as.integer(format(dates_hist, "%Y"))))

draws_all <- lapply(names(bundles), function(p) {
  rstan::extract(bundles[[p]]$fit, pars = c("alpha_R", "z_year", "year_effect", "tau_year", "sigma_season_year", "A_year",
                                             "beta_T", "beta_P", "S_prop", "immune_prop", "R0_t", "phi_obs"), permuted = TRUE)
})
names(draws_all) <- names(bundles)

# ================================================================
# tau_year posterior vs. prior
# ================================================================
tau_summary <- bind_rows(lapply(names(bundles), function(p) {
  d <- draws_all[[p]]$tau_year
  tibble(prior_id = p, prior_scale = TAU_PRIORS[[p]], posterior_median = median(d), posterior_mean = mean(d), posterior_sd = sd(d),
         lo95 = quantile(d, .025), hi95 = quantile(d, .975), posterior_sd_over_prior_scale = sd(d) / TAU_PRIORS[[p]],
         prior_posterior_median_ratio = median(d) / TAU_PRIORS[[p]])
}))
print(as.data.frame(tau_summary), digits = 4)
write_csv(tau_summary, file.path(DIR_TABLES_YV, "PHASE1_tau_year_posterior_summary.csv"))

# Raw draws for FIG17
saveRDS(lapply(draws_all, `[[`, "tau_year"), file.path(DIR_RESULTS_YV, "PHASE1_tau_year_draws.rds"))

# ================================================================
# Prior-dependence diagnostic: if the posterior median scales ~linearly
# with the prior scale, tau_year is prior-dominated; if it is roughly
# CONSTANT across P1/P2/P3, the likelihood is informative.
# ================================================================
message(sprintf("\n[Phase1 analysis] tau_year posterior medians: P1=%.4f, P2=%.4f, P3=%.4f (prior scales 0.20/0.40/0.80)",
                 tau_summary$posterior_median[tau_summary$prior_id == "P1"], tau_summary$posterior_median[tau_summary$prior_id == "P2"], tau_summary$posterior_median[tau_summary$prior_id == "P3"]))
message("If these three medians are similar to each other (NOT similar to 0.20/0.40/0.80 respectively), tau_year is likelihood-informed, not prior-dependent.")

# ================================================================
# Posterior correlations: tau_year with alpha_R, sigma_season_year, beta_T, beta_P, terminal S_prop
# ================================================================
cor_summary <- bind_rows(lapply(names(bundles), function(p) {
  d <- draws_all[[p]]
  tibble(prior_id = p, with_quantity = c("alpha_R", "sigma_season_year", "beta_T", "beta_P", "terminal S_prop (2025-12-21)"),
         correlation = c(cor(d$tau_year, d$alpha_R), cor(d$tau_year, d$sigma_season_year), cor(d$tau_year, d$beta_T), cor(d$tau_year, d$beta_P), cor(d$tau_year, d$S_prop[, N_hist])))
}))
print(as.data.frame(cor_summary), digits = 3)
write_csv(cor_summary, file.path(DIR_TABLES_YV, "PHASE1_tau_year_correlations.csv"))

# ================================================================
# Implied historical year_effect[y] distributions under each prior
# ================================================================
year_effect_by_prior <- bind_rows(lapply(names(bundles), function(p) {
  d <- draws_all[[p]]$year_effect
  tibble(prior_id = p, year = hist_years, median = apply(d, 2, median), lo95 = apply(d, 2, quantile, .025), hi95 = apply(d, 2, quantile, .975))
}))
write_csv(year_effect_by_prior, file.path(DIR_TABLES_YV, "PHASE1_implied_year_effect_by_prior.csv"))

# within-draw spread (genuine interannual variability implied by each fit) -- directly comparable to the EARLIER
# year_effect_diagnostic's fixed-SD=0.40 finding (historical within-draw SD was 0.138 there)
within_draw_spread <- bind_rows(lapply(names(bundles), function(p) {
  sd_vals <- apply(draws_all[[p]]$year_effect, 1, sd)
  tibble(prior_id = p, within_draw_sd_median = median(sd_vals), within_draw_sd_lo95 = quantile(sd_vals, .025), within_draw_sd_hi95 = quantile(sd_vals, .975))
}))
print(as.data.frame(within_draw_spread), digits = 4)
write_csv(within_draw_spread, file.path(DIR_TABLES_YV, "PHASE1_within_draw_year_effect_spread.csv"))

# ================================================================
# Annual mean/peak R0(t) per prior + HMC
# ================================================================
year_id_hist <- bundles$P2$stan_data$year_id
annual_R0 <- bind_rows(lapply(names(bundles), function(p) {
  d <- draws_all[[p]]
  bind_rows(lapply(seq_len(min(400L, nrow(d$R0_t))), function(i) {
    tibble(prior_id = p, draw = i, year = hist_years, year_id = seq_along(hist_years)) |>
      mutate(mean_R0 = vapply(year_id, function(y) mean(d$R0_t[i, year_id_hist == y]), numeric(1)),
             peak_R0 = vapply(year_id, function(y) max(d$R0_t[i, year_id_hist == y]), numeric(1)))
  }))
}))
write_csv(annual_R0 |> sample_frac(min(1, 20000 / nrow(annual_R0))), file.path(DIR_RESULTS_YV, "PHASE1_annual_R0_sample_by_prior.csv"))
annual_R0_summary <- annual_R0 |> group_by(prior_id, year) |> summarise(peak_R0_median = median(peak_R0), peak_R0_lo95 = quantile(peak_R0, .025), peak_R0_hi95 = quantile(peak_R0, .975), .groups = "drop")
write_csv(annual_R0_summary, file.path(DIR_TABLES_YV, "PHASE1_annual_R0_summary_by_prior.csv"))

hmc_table <- read_csv(file.path(DIR_TABLES_YV, "PHASE1_hmc_diagnostics_P1_P3.csv"), show_col_types = FALSE)
message("\n=== HMC (from fitting step) ===")
print(as.data.frame(hmc_table), digits = 4)

message("\n[Phase1 analysis] Complete. Tables -> ", DIR_TABLES_YV)
