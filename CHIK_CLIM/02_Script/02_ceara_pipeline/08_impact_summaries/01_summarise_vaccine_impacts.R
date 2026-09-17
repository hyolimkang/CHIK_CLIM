# Primary outcomes (Section 6) and Tables 1-3 (Section 10).

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(tidyr) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
script_dir <- file.path(root, "02_Script/02_ceara_pipeline/07_counterfactuals")
source(file.path(script_dir, "01_config.R"))

res <- readRDS(file.path(DIR_RESULTS, "three_arm_results.rds"))
totals <- res$totals
vax <- res$vax_bookkeeping
dates <- res$dates
agg_Inew <- res$agg_Inew
agg_S <- res$agg_S
agg_N <- res$agg_N
Q_ASCERTAINMENT <- 0.05 # accepted CE q=0.05, fixed -- for reported-cases-averted only

q95 <- function(x) quantile(x, c(.025, .5, .975), names = FALSE)

# ================================================================
# TABLE 1: historical counterfactual summary, by scenario
# ================================================================
by_draw <- totals |> group_by(draw) |>
  summarise(
    cum_inf_A = sum(X_A), cum_inf_B = sum(X_B), cum_inf_C = sum(X_C),
    peak_A = max(X_A), peak_B = max(X_B), peak_C = max(X_C),
    final_S_prop_A = S_A[which.max(t)] / N_A[which.max(t)], final_S_prop_B = S_B[which.max(t)] / N_B[which.max(t)],
    .groups = "drop"
  )
doses_by_draw <- vax |> group_by(draw) |> summarise(doses_B = sum(doses_B), .groups = "drop")
by_draw <- by_draw |> left_join(doses_by_draw, by = "draw") |>
  mutate(infections_averted_B = cum_inf_A - cum_inf_B, pct_reduction_B = 100 * infections_averted_B / cum_inf_A,
         reported_cases_A = Q_ASCERTAINMENT * cum_inf_A, reported_cases_B = Q_ASCERTAINMENT * cum_inf_B)

summarise_col <- function(x) { q <- q95(x); tibble(median = q[2], lo95 = q[1], hi95 = q[3]) }
table1 <- bind_rows(
  bind_cols(scenario = "A: No vaccine", summarise_col(by_draw$cum_inf_A)) |> rename_with(~paste0("cumulative_infections_", .), -scenario),
  bind_cols(scenario = "B: Routine age-12 (stress test)", summarise_col(by_draw$cum_inf_B)) |> rename_with(~paste0("cumulative_infections_", .), -scenario)
)
table1_full <- tibble(
  scenario = c("A: No vaccine", "B: Routine age-12 (stress test)"),
  cumulative_infections_median = c(median(by_draw$cum_inf_A), median(by_draw$cum_inf_B)),
  cumulative_infections_lo95 = c(quantile(by_draw$cum_inf_A, .025), quantile(by_draw$cum_inf_B, .025)),
  cumulative_infections_hi95 = c(quantile(by_draw$cum_inf_A, .975), quantile(by_draw$cum_inf_B, .975)),
  cumulative_reported_cases_median = c(median(by_draw$reported_cases_A), median(by_draw$reported_cases_B)),
  peak_weekly_infections_median = c(median(by_draw$peak_A), median(by_draw$peak_B)),
  final_S_prop_median = c(median(by_draw$final_S_prop_A), median(by_draw$final_S_prop_B)),
  doses_median = c(0, median(by_draw$doses_B)),
  infections_averted_median = c(NA, median(by_draw$infections_averted_B)),
  infections_averted_lo95 = c(NA, quantile(by_draw$infections_averted_B, .025)),
  infections_averted_hi95 = c(NA, quantile(by_draw$infections_averted_B, .975)),
  percent_reduction_median = c(NA, median(by_draw$pct_reduction_B))
)
message("=== TABLE 1: historical counterfactual summary ===")
print(as.data.frame(table1_full), digits = 4)
write_csv(table1_full, file.path(DIR_TABLES, "TABLE1_historical_counterfactual_summary.csv"))

# ================================================================
# TABLE 2: age-specific impact
# ================================================================
age_group_cum <- bind_rows(lapply(REPORTING_AGE_LABELS, function(g) {
  cum_A_by_draw <- rowSums(agg_Inew[, , g, "A"])
  cum_B_by_draw <- rowSums(agg_Inew[, , g, "B"])
  tibble(age_group = g, infections_no_vaccine_median = median(cum_A_by_draw), infections_no_vaccine_lo95 = quantile(cum_A_by_draw, .025), infections_no_vaccine_hi95 = quantile(cum_A_by_draw, .975),
         infections_vaccine_median = median(cum_B_by_draw), infections_vaccine_lo95 = quantile(cum_B_by_draw, .025), infections_vaccine_hi95 = quantile(cum_B_by_draw, .975),
         infections_averted_median = median(cum_A_by_draw - cum_B_by_draw),
         percent_reduction_median = median(100 * (cum_A_by_draw - cum_B_by_draw) / cum_A_by_draw))
})) |> mutate(age_group = factor(age_group, levels = REPORTING_AGE_LABELS))
message("\n=== TABLE 2: age-specific impact ===")
print(as.data.frame(age_group_cum), digits = 4)
write_csv(age_group_cum, file.path(DIR_TABLES, "TABLE2_age_specific_impact.csv"))

# ================================================================
# TABLE 3: effect decomposition (total / direct-only / indirect)
# ================================================================
decomp_by_draw <- totals |> group_by(draw) |> summarise(cum_X_A = sum(X_A), cum_X_B = sum(X_B), cum_X_C = sum(X_C), .groups = "drop") |>
  mutate(total_effect = cum_X_A - cum_X_B, direct_only = cum_X_A - cum_X_C, indirect = cum_X_C - cum_X_B, indirect_fraction = indirect / total_effect)
table3 <- tibble(
  quantity = c("Total infections averted", "Direct-only component (Arm A - Arm C)", "Transmission-mediated indirect component (Arm C - Arm B)", "Indirect fraction of total effect"),
  median = c(median(decomp_by_draw$total_effect), median(decomp_by_draw$direct_only), median(decomp_by_draw$indirect), median(decomp_by_draw$indirect_fraction)),
  lo95 = c(quantile(decomp_by_draw$total_effect, .025), quantile(decomp_by_draw$direct_only, .025), quantile(decomp_by_draw$indirect, .025), quantile(decomp_by_draw$indirect_fraction, .025)),
  hi95 = c(quantile(decomp_by_draw$total_effect, .975), quantile(decomp_by_draw$direct_only, .975), quantile(decomp_by_draw$indirect, .975), quantile(decomp_by_draw$indirect_fraction, .975))
)
message("\n=== TABLE 3: effect decomposition ===")
print(as.data.frame(table3), digits = 4)
write_csv(table3, file.path(DIR_TABLES, "TABLE3_effect_decomposition.csv"))
write_csv(decomp_by_draw, file.path(DIR_TABLES, "TABLE3_decomposition_by_draw.csv"))
write_csv(by_draw, file.path(DIR_TABLES, "table1_by_draw.csv"))

message(sprintf("\nSUMMARY: median total infections averted = %.0f (%.1f%% reduction); indirect fraction = %.1f%%",
                 median(decomp_by_draw$total_effect), median(by_draw$pct_reduction_B), 100 * median(decomp_by_draw$indirect_fraction)))
message("[saved] TABLE1_historical_counterfactual_summary.csv, TABLE2_age_specific_impact.csv, TABLE3_effect_decomposition.csv")
