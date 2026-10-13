# PHASE 0 addendum -- Phase 1's lambda=1 (I2) arm for H1/G2 showed
# prob_major_epidemic=19.9%, roughly 2x Step 6's reported 10.55% for the
# SAME model and lambda, a larger gap than Phase 0's own entry/S/Reff
# checks (<6%). Investigate explicitly: is this within plausible
# between-posterior-draw-subsample Monte Carlo noise (independent
# sample.int() of 100 of ~1600 available draws, different seed stream
# from Step 6), or a genuine implementation difference?

suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble) })
filter <- dplyr::filter
source(file.path(ROOT, "02_Script/02_ceara_pipeline/16_step7_importation_persistence/00_config_step7.R"))
source(file.path(DIR_STEP6, "00_config_step6.R"))
max_run <- function(x) { rr <- rle(replace(x, is.na(x), FALSE)); if (!any(rr$values)) 0L else max(rr$lengths[rr$values]) }

phase0 <- readRDS(file.path(DIR_RESULTS_STEP7, "PHASE0_baseline_lambda1_simulations.rds"))
traj <- phase0$H1$trajectories |> mutate(year = as.integer(format(week_start, "%Y")),
                                          import_hazard = M_import / N, local_hazard = pmax(lambda_t - M_import / N, 0),
                                          favourable = R_eff > 1, locally_dominant = X_latent > 0 & local_hazard > import_hazard)
annual <- traj |> filter(year == 2022) |> group_by(posterior_draw_id, stochastic_replicate_id) |>
  summarise(annual_reported_cases = sum(Y_reported), favourable_pressure = any(favourable & M_import > 0),
            longest_locally_dominant_run = max_run(locally_dominant), population_mean = mean(N), .groups = "drop") |>
  mutate(established = favourable_pressure & longest_locally_dominant_run >= ESTABLISHMENT_RUN_WEEKS,
         annual_reported_incidence_per_100k = annual_reported_cases / population_mean * 1e5,
         major_epidemic = annual_reported_cases >= MAJOR_ABS_CASE_MINIMUM & annual_reported_incidence_per_100k >= MAJOR_INCIDENCE_MINIMUM_PER_100K)

# Between-draw heterogeneity: P(major) computed WITHIN each of the 100 sampled posterior draws (20 replicates each).
per_draw <- annual |> group_by(posterior_draw_id) |> summarise(p_major_within_draw = mean(major_epidemic), n = n(), .groups = "drop")
write_csv(per_draw, file.path(DIR_TABLES_STEP7, "PHASE0_addendum_prob_major_per_draw.csv"))

set.seed(1L)
boot <- replicate(4000L, mean(sample(per_draw$p_major_within_draw, nrow(per_draw), replace = TRUE)))
boot_ci <- tibble(statistic = "prob_major_epidemic_H1_G2_lambda1", point_estimate_new_run = mean(annual$major_epidemic),
                   step6_reported_value = 0.1055,
                   cluster_bootstrap_lo95 = as.numeric(quantile(boot, .025)), cluster_bootstrap_hi95 = as.numeric(quantile(boot, .975)),
                   n_posterior_draws = nrow(per_draw), sd_p_major_across_draws = sd(per_draw$p_major_within_draw))
write_csv(boot_ci, file.path(DIR_TABLES_STEP7, "PHASE0_addendum_prob_major_discrepancy_check.csv"))
message("=== Phase0 addendum: is the H1/G2 lambda=1 P(major)=19.9% vs Step6's 10.55% within MC noise? ===")
print(as.data.frame(boot_ci), digits = 4)
message(sprintf("SD of within-draw P(major) across the %d sampled posterior draws: %.3f -- a wide between-draw spread means WHICH 100 of ~1600 draws get sampled materially changes the aggregate probability.",
                 nrow(per_draw), sd(per_draw$p_major_within_draw)))
if (boot_ci$step6_reported_value[1] >= boot_ci$cluster_bootstrap_lo95[1] && boot_ci$step6_reported_value[1] <= boot_ci$cluster_bootstrap_hi95[1]) {
  message("[Phase0 addendum] Step 6's reported value falls WITHIN this run's cluster-bootstrap 95% interval -- consistent with between-posterior-draw-subsample Monte Carlo noise, not an implementation difference.")
} else {
  message("[Phase0 addendum WARNING] Step 6's reported value falls OUTSIDE this run's cluster-bootstrap 95% interval -- investigate further before treating this as MC noise alone.")
}
