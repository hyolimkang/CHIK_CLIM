# Number-Needed-to-Vaccinate (NNV) module.
# Extends the existing historical age-12 vaccine counterfactual analysis.
# Reads the already-saved three_arm_results.rds -- does NOT re-run the
# transmission model or the age-cohort simulator, and does NOT modify any
# existing results/tables/figures.
#
# NUMERATOR CONVENTION (per instruction): cumulative DOSES ADMINISTERED,
# never "effectively immunised" / "successfully protected".
# All NNV ratios are computed DRAW-BY-DRAW first, then summarised --
# never median(doses)/median(effect).

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(tidyr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
script_dir <- file.path(root, "02_Script/40_renewal_model/37_ceara_historical_age12_vaccine_counterfactual/scripts")
source(file.path(script_dir, "00_config.R"))

DIR_RESULTS_NNV <- file.path(DIR_RESULTS, "nnv"); DIR_TABLES_NNV <- file.path(DIR_TABLES, "nnv"); DIR_DIAG_NNV <- file.path(DIR_DIAGNOSTICS, "nnv")
for (d in c(DIR_RESULTS_NNV, DIR_TABLES_NNV, DIR_DIAG_NNV)) dir.create(d, recursive = TRUE, showWarnings = FALSE) # figures/nnv comes from 00_config.R (DIR_FIG_NNV, under 03_Output/)

res <- readRDS(file.path(DIR_RESULTS, "three_arm_results.rds"))
totals <- res$totals; vax <- res$vax_bookkeeping; dates <- res$dates; agg_Inew <- res$agg_Inew

Q_ASCERTAINMENT <- 0.05 # fixed CE q=0.05 canary -- NOT a posterior; constant-q observation mapping (documented, see NNV_VALIDATION.md item 6)
NNV_TOL_INFECTIONS <- 1 # denominator smaller than this (in infection count) is treated as ~0, NNV undefined
NNV_TOL_REPORTED_CASES <- Q_ASCERTAINMENT * NNV_TOL_INFECTIONS

stresstest <- SCENARIOS$stresstest

# ================================================================
# SECTION 2: draw-by-draw NNV (primary + decomposition)
# ================================================================
doses_check <- vax |> group_by(draw) |> summarise(doses_B = sum(doses_B), doses_C = sum(doses_C), .groups = "drop")
stopifnot(all(abs(doses_check$doses_B - doses_check$doses_C) < 1e-6)) # same programme in Arms B and C, by design

by_draw <- totals |> group_by(draw) |>
  summarise(infections_no_vaccine = sum(X_A), infections_vaccine = sum(X_B), infections_direct_only_arm = sum(X_C), .groups = "drop")

vax_by_draw <- vax |> group_by(draw) |>
  summarise(cumulative_eligible = sum(entrants_total_B), cumulative_doses = sum(doses_B), cumulative_effectively_protected = sum(effective_B), .groups = "drop")

nnv_draws <- by_draw |> left_join(vax_by_draw, by = "draw") |>
  mutate(
    scenario = "stresstest_coverage1_VE1",
    coverage = stresstest$coverage, VE_infection = stresstest$ve_infection,
    vaccination_start = WINDOW_START, vaccination_end = WINDOW_END, target_age = TARGET_AGE,

    infections_averted = infections_no_vaccine - infections_vaccine,
    reported_cases_no_vaccine = Q_ASCERTAINMENT * infections_no_vaccine,
    reported_cases_vaccine = Q_ASCERTAINMENT * infections_vaccine,
    reported_cases_averted = Q_ASCERTAINMENT * infections_averted,

    NNV_infection = ifelse(infections_averted > 0, cumulative_doses / infections_averted, NA_real_),
    NNV_reported_case = ifelse(reported_cases_averted > 0, cumulative_doses / reported_cases_averted, NA_real_),

    # Section 7: total / direct-only / indirect decomposition -- diagnostic ratios only.
    direct_only_effect = infections_no_vaccine - infections_direct_only_arm,       # A - C
    indirect_effect = infections_direct_only_arm - infections_vaccine,            # C - B
    NNV_total_infection = NNV_infection,                                          # primary programme-level NNV (A - B)
    NNV_direct_only_infection = ifelse(direct_only_effect > 0, cumulative_doses / direct_only_effect, NA_real_),
    NNV_indirect_infection = ifelse(indirect_effect > 0, cumulative_doses / indirect_effect, NA_real_)
  ) |>
  select(draw, scenario, coverage, VE_infection, vaccination_start, vaccination_end, target_age,
         cumulative_eligible, cumulative_doses, cumulative_effectively_protected,
         infections_no_vaccine, infections_vaccine, infections_averted,
         reported_cases_no_vaccine, reported_cases_vaccine, reported_cases_averted,
         NNV_infection, NNV_reported_case,
         infections_direct_only_arm, direct_only_effect, indirect_effect,
         NNV_total_infection, NNV_direct_only_infection, NNV_indirect_infection)

saveRDS(nnv_draws, file.path(DIR_RESULTS_NNV, "nnv_posterior_draws.rds"))
write_csv(nnv_draws, file.path(DIR_RESULTS_NNV, "nnv_posterior_draws.csv"))
message("[saved] results/nnv/nnv_posterior_draws.{rds,csv} -- ", nrow(nnv_draws), " draws")

# ================================================================
# SECTION 4: primary NNV summary table
# ================================================================
qsum <- function(x, na.rm = TRUE) {
  q <- quantile(x, c(.025, .25, .5, .75, .975), na.rm = na.rm, names = FALSE)
  tibble(p2_5 = q[1], p25 = q[2], median = q[3], p75 = q[4], p97_5 = q[5])
}
summary_block <- function(x, label) bind_cols(quantity = label, qsum(x))

table_nnv_primary <- bind_rows(
  summary_block(nnv_draws$cumulative_doses, "Cumulative doses administered"),
  summary_block(nnv_draws$infections_averted, "Infections averted (total, A-B)"),
  summary_block(nnv_draws$reported_cases_averted, "Reported cases averted (q=0.05 fixed)"),
  summary_block(nnv_draws$NNV_infection, "NNV per infection averted"),
  summary_block(nnv_draws$NNV_reported_case, "NNV per reported case averted")
) |>
  mutate(pr_benefit_gt0 = c(NA, mean(nnv_draws$infections_averted > 0), mean(nnv_draws$reported_cases_averted > 0), mean(nnv_draws$infections_averted > 0), mean(nnv_draws$reported_cases_averted > 0)))

pr_benefit <- tibble(
  quantity = c("Pr(infections_averted > 0)", "Pr(infections_averted <= 0)", "Pr(reported_cases_averted > 0)", "Pr(reported_cases_averted <= 0)"),
  probability = c(mean(nnv_draws$infections_averted > 0), mean(nnv_draws$infections_averted <= 0),
                  mean(nnv_draws$reported_cases_averted > 0), mean(nnv_draws$reported_cases_averted <= 0))
)

saveRDS(table_nnv_primary, file.path(DIR_TABLES_NNV, "table_nnv_primary.rds"))
write_csv(table_nnv_primary, file.path(DIR_TABLES_NNV, "table_nnv_primary.csv"))
write_csv(pr_benefit, file.path(DIR_TABLES_NNV, "table_nnv_probability_benefit.csv"))
message("\n=== TABLE: primary NNV summary ===")
print(as.data.frame(table_nnv_primary), digits = 4)
print(as.data.frame(pr_benefit), digits = 4)

# ================================================================
# SECTION 5: NNV by evaluation horizon
# ================================================================
# Horizons determined from the actual programme start (WINDOW_START =
# 2015-01-04) and available follow-up (through WINDOW_END = 2021-12-26).
# All horizons lie within available data -- none extrapolated.
horizon_dates <- list(
  "1yr" = as.Date("2016-01-04"), "2yr" = as.Date("2017-01-04"), "3yr" = as.Date("2018-01-04"),
  "5yr" = as.Date("2020-01-04"), "end_of_analysis_2021" = WINDOW_END
)
horizon_idx <- sapply(horizon_dates, function(dt) max(which(dates <= dt)))

cum_by_draw_week <- totals |> arrange(draw, t) |> group_by(draw) |>
  mutate(cum_X_A = cumsum(X_A), cum_X_B = cumsum(X_B), cum_X_C = cumsum(X_C)) |> ungroup() |>
  left_join(vax |> arrange(draw, t) |> group_by(draw) |> mutate(cum_doses = cumsum(doses_B)) |> ungroup() |> select(draw, t, cum_doses), by = c("draw", "t"))

nnv_by_horizon <- bind_rows(lapply(names(horizon_idx), function(h) {
  idx <- horizon_idx[[h]]
  sub <- cum_by_draw_week |> filter(t == idx) |>
    mutate(horizon = h, horizon_date = dates[idx],
           infections_averted_h = cum_X_A - cum_X_B,
           reported_cases_averted_h = Q_ASCERTAINMENT * (cum_X_A - cum_X_B),
           NNV_infection_h = ifelse(infections_averted_h > 0, cum_doses / infections_averted_h, NA_real_),
           NNV_reported_case_h = ifelse(reported_cases_averted_h > 0, cum_doses / reported_cases_averted_h, NA_real_))
  sub |> select(draw, horizon, horizon_date, cumulative_doses_h = cum_doses, infections_averted_h, reported_cases_averted_h, NNV_infection_h, NNV_reported_case_h)
})) |> mutate(horizon = factor(horizon, levels = names(horizon_idx)))

saveRDS(nnv_by_horizon, file.path(DIR_RESULTS_NNV, "nnv_by_horizon.rds"))

table_nnv_by_horizon <- nnv_by_horizon |> group_by(horizon, horizon_date) |>
  summarise(
    doses_median = median(cumulative_doses_h),
    infections_averted_median = median(infections_averted_h), infections_averted_lo95 = quantile(infections_averted_h, .025), infections_averted_hi95 = quantile(infections_averted_h, .975),
    reported_cases_averted_median = median(reported_cases_averted_h),
    NNV_infection_median = median(NNV_infection_h, na.rm = TRUE), NNV_infection_lo95 = quantile(NNV_infection_h, .025, na.rm = TRUE), NNV_infection_hi95 = quantile(NNV_infection_h, .975, na.rm = TRUE),
    NNV_reported_case_median = median(NNV_reported_case_h, na.rm = TRUE),
    pr_benefit_gt0 = mean(infections_averted_h > 0),
    .groups = "drop"
  )
write_csv(table_nnv_by_horizon, file.path(DIR_TABLES_NNV, "table_nnv_by_horizon.csv"))
message("\n=== TABLE: NNV by evaluation horizon ===")
print(as.data.frame(table_nnv_by_horizon), digits = 4)

# ================================================================
# REVISION (fig_nnv2): cumulative NNV at CALENDAR YEAR-END (not
# programme-start-relative horizons -- those remain above, unchanged, as a
# secondary table). Each year's value is cumulative from programme start
# (2015-01-04) through 31 December of that calendar year -- i.e. "how many
# doses, cumulatively, to avert one outcome cumulatively, by this point" --
# NEVER an annual (year-only) doses/effect ratio.
# ================================================================
analysis_years <- sort(unique(as.integer(format(dates, "%Y"))))
year_end_idx <- sapply(analysis_years, function(yr) max(which(as.integer(format(dates, "%Y")) == yr)))
names(year_end_idx) <- analysis_years

nnv_year_end_by_draw <- bind_rows(lapply(names(year_end_idx), function(yr) {
  idx <- year_end_idx[[yr]]
  cum_by_draw_week |> filter(t == idx) |>
    transmute(draw, year = as.integer(yr), year_end_date = dates[idx],
              cumulative_doses_y = cum_doses,
              cumulative_infections_averted_y = cum_X_A - cum_X_B,
              cumulative_reported_cases_averted_y = Q_ASCERTAINMENT * (cum_X_A - cum_X_B),
              NNV_infection_y = ifelse(cumulative_infections_averted_y > NNV_TOL_INFECTIONS, cumulative_doses_y / cumulative_infections_averted_y, NA_real_),
              NNV_reported_case_y = ifelse(cumulative_reported_cases_averted_y > NNV_TOL_REPORTED_CASES, cumulative_doses_y / cumulative_reported_cases_averted_y, NA_real_))
}))
saveRDS(nnv_year_end_by_draw, file.path(DIR_RESULTS_NNV, "nnv_year_end_cumulative_by_draw.rds"))

table_nnv_year_end <- nnv_year_end_by_draw |> group_by(year, year_end_date) |>
  summarise(
    doses_median = median(cumulative_doses_y),
    infections_averted_median = median(cumulative_infections_averted_y), infections_averted_lo95 = quantile(cumulative_infections_averted_y, .025), infections_averted_hi95 = quantile(cumulative_infections_averted_y, .975),
    pr_infections_averted_gt0 = mean(cumulative_infections_averted_y > NNV_TOL_INFECTIONS),
    reported_cases_averted_median = median(cumulative_reported_cases_averted_y), reported_cases_averted_lo95 = quantile(cumulative_reported_cases_averted_y, .025), reported_cases_averted_hi95 = quantile(cumulative_reported_cases_averted_y, .975),
    pr_reported_cases_averted_gt0 = mean(cumulative_reported_cases_averted_y > NNV_TOL_REPORTED_CASES),
    NNV_infection_median = median(NNV_infection_y, na.rm = TRUE), NNV_infection_lo95 = quantile(NNV_infection_y, .025, na.rm = TRUE), NNV_infection_hi95 = quantile(NNV_infection_y, .975, na.rm = TRUE),
    NNV_reported_case_median = median(NNV_reported_case_y, na.rm = TRUE), NNV_reported_case_lo95 = quantile(NNV_reported_case_y, .025, na.rm = TRUE), NNV_reported_case_hi95 = quantile(NNV_reported_case_y, .975, na.rm = TRUE),
    frac_draws_NNV_infection_defined = mean(!is.na(NNV_infection_y)), frac_draws_NNV_reported_case_defined = mean(!is.na(NNV_reported_case_y)),
    .groups = "drop"
  ) |>
  mutate(
    # "Defined" (denominator > tiny tolerance) is not the same as "stable":
    # e.g. 2015 has infections_averted > 0 in every draw but only ~30
    # infections accrued, against ~135,000 doses -- technically finite but
    # not a meaningfully precise cumulative NNV. Flag as unstable whenever
    # the median cumulative infections averted is below a materially
    # meaningful floor (chosen well below the smallest "mature" year,
    # ~150,000 infections averted by end of 2016) OR when <50% of draws
    # even clear the strict near-zero tolerance.
    nnv_infection_stable = frac_draws_NNV_infection_defined >= 0.5 & infections_averted_median >= 1000,
    stability_note = ifelse(nnv_infection_stable, "", "Insufficient accrued benefit for stable NNV")
  )
write_csv(table_nnv_year_end, file.path(DIR_TABLES_NNV, "table_nnv_year_end_cumulative.csv"))
message("\n=== TABLE: cumulative NNV at calendar year-end ===")
print(as.data.frame(table_nnv_year_end), digits = 4)
if (any(!table_nnv_year_end$nnv_infection_stable)) message("[note] Unstable/undefined-NNV years (<50% of draws have a positive denominator): ", paste(table_nnv_year_end$year[!table_nnv_year_end$nnv_infection_stable], collapse = ", "))

# ================================================================
# SECTION 6: NNV through calendar time (cumulative, weekly)
# ================================================================
nnv_weekly <- cum_by_draw_week |>
  mutate(cum_infections_averted = cum_X_A - cum_X_B,
         cum_reported_cases_averted = Q_ASCERTAINMENT * cum_infections_averted,
         NNV_infection_cum = ifelse(cum_infections_averted > NNV_TOL_INFECTIONS, cum_doses / cum_infections_averted, NA_real_),
         NNV_reported_case_cum = ifelse(cum_reported_cases_averted > NNV_TOL_REPORTED_CASES, cum_doses / cum_reported_cases_averted, NA_real_)) |>
  select(draw, t, week_start, cum_doses, cum_infections_averted, cum_reported_cases_averted, NNV_infection_cum, NNV_reported_case_cum)

saveRDS(nnv_weekly, file.path(DIR_RESULTS_NNV, "nnv_cumulative_by_week.rds"))

nnv_weekly_summary <- nnv_weekly |> group_by(t, week_start) |>
  summarise(doses_median = median(cum_doses), infections_averted_median = median(cum_infections_averted),
            NNV_infection_median = median(NNV_infection_cum, na.rm = TRUE), NNV_infection_lo95 = quantile(NNV_infection_cum, .025, na.rm = TRUE), NNV_infection_hi95 = quantile(NNV_infection_cum, .975, na.rm = TRUE),
            NNV_reported_case_median = median(NNV_reported_case_cum, na.rm = TRUE),
            frac_draws_denominator_above_tol = mean(cum_infections_averted > NNV_TOL_INFECTIONS), .groups = "drop")
write_csv(nnv_weekly_summary, file.path(DIR_TABLES_NNV, "table_nnv_by_week_summary.csv"))
message("[saved] results/nnv/nnv_cumulative_by_week.rds, tables/nnv/table_nnv_by_week_summary.csv")

# ================================================================
# SECTION 8 (REVISED): age-specific VACCINE BENEFIT -- infections averted and
# percent reduction ONLY. No age-group is assigned its own NNV: the age-12
# routine programme's doses are the numerator for exactly ONE, programme-
# level NNV (Sections 2-6 above); age groups other than 12 can receive
# substantial benefit (transmission-mediated) but that benefit is never
# divided by a per-group dose count, because no doses were administered to
# those groups. See reports/NNV_ASSESSMENT.md Section "Direct/indirect
# interpretation" for the full explanation.
# ================================================================
age_effect_by_draw <- bind_rows(lapply(REPORTING_AGE_LABELS, function(g) {
  tibble(draw = res$draw_idx, age_group = g,
         infections_no_vaccine_g = rowSums(agg_Inew[, , g, "A"]), infections_vaccine_g = rowSums(agg_Inew[, , g, "B"]))
})) |> mutate(infections_averted_g = infections_no_vaccine_g - infections_vaccine_g,
              percent_reduction_g = 100 * infections_averted_g / infections_no_vaccine_g)

table_age_benefit <- age_effect_by_draw |> group_by(age_group) |>
  summarise(infections_averted_median = median(infections_averted_g), infections_averted_lo95 = quantile(infections_averted_g, .025), infections_averted_hi95 = quantile(infections_averted_g, .975),
            percent_reduction_median = median(percent_reduction_g), percent_reduction_lo95 = quantile(percent_reduction_g, .025), percent_reduction_hi95 = quantile(percent_reduction_g, .975),
            .groups = "drop") |>
  mutate(age_group = factor(age_group, levels = REPORTING_AGE_LABELS),
         is_vaccine_targeted = age_group == "12") |> arrange(age_group)
write_csv(table_age_benefit, file.path(DIR_TABLES_NNV, "table_nnv_age_specific_benefit.csv"))
saveRDS(age_effect_by_draw, file.path(DIR_RESULTS_NNV, "nnv_age_specific_by_draw.rds"))
message("\n=== TABLE: age-specific vaccine BENEFIT (infections averted, % reduction) -- NOT age-specific NNV ===")
print(as.data.frame(table_age_benefit), digits = 4)

# Diagnostic-only ratio (doses of the AGE-12 PROGRAMME / infections averted
# IN A GIVEN AGE GROUP) -- kept ONLY as a non-main diagnostic file, per
# instruction ("REMOVE this quantity from the MAIN figures and tables").
# Explicitly NOT called NNV anywhere in this block.
age_dose_ratio_diagnostic <- age_effect_by_draw |>
  left_join(nnv_draws |> select(draw, cumulative_doses), by = "draw") |>
  mutate(programme_doses_per_infection_averted_in_group = ifelse(infections_averted_g > 0, cumulative_doses / infections_averted_g, NA_real_)) |>
  group_by(age_group) |>
  summarise(programme_doses_per_infection_averted_in_group_median = median(programme_doses_per_infection_averted_in_group, na.rm = TRUE),
            programme_doses_per_infection_averted_in_group_lo95 = quantile(programme_doses_per_infection_averted_in_group, .025, na.rm = TRUE),
            programme_doses_per_infection_averted_in_group_hi95 = quantile(programme_doses_per_infection_averted_in_group, .975, na.rm = TRUE),
            .groups = "drop") |>
  mutate(age_group = factor(age_group, levels = REPORTING_AGE_LABELS), note = "DIAGNOSTIC ONLY -- programme doses are administered to age 12, not to this group; this is NOT an age-specific NNV") |> arrange(age_group)
write_csv(age_dose_ratio_diagnostic, file.path(DIR_DIAG_NNV, "DIAGNOSTIC_doses_per_infection_by_group.csv")) # short filename -- long path names hit Windows MAX_PATH under this OneDrive-synced repo root
message("[relocated] doses/infections-averted-by-group ratio moved out of tables/nnv/ to diagnostics/nnv/ (diagnostic only, not NNV, not a main output)")

# ================================================================
# SECTION 9: baseline susceptibility vs NNV (descriptive only)
# ================================================================
agg_S <- res$agg_S; agg_N <- res$agg_N
baseline_S_age <- tibble( # positional: dim-1 index i of agg_S/agg_N corresponds to res$draw_idx[i], by construction in 06_execute_three_arms_all_draws.R
  draw = res$draw_idx,
  S_age12_N_baseline = agg_S[, 1, "12", "A"] / agg_N[, 1, "12", "A"],
  S_0_11_N_baseline = agg_S[, 1, "0-11", "A"] / agg_N[, 1, "0-11", "A"],
  S_65plus_N_baseline = agg_S[, 1, "65+", "A"] / agg_N[, 1, "65+", "A"]
)
baseline_S_total <- totals |> filter(t == 1) |> transmute(draw, S_total_N_baseline = S_A / N_A)
baseline_S <- baseline_S_total |> inner_join(baseline_S_age, by = "draw") # explicit join by draw id -- never relies on row order
nnv_susceptibility <- nnv_draws |> select(draw, NNV_infection, infections_averted) |> left_join(baseline_S, by = "draw")
saveRDS(nnv_susceptibility, file.path(DIR_RESULTS_NNV, "nnv_susceptibility_link.rds"))
write_csv(nnv_susceptibility, file.path(DIR_TABLES_NNV, "table_nnv_susceptibility_link.csv"))

spearman_cors <- nnv_susceptibility |> filter(is.finite(NNV_infection)) |>
  summarise(cor_S_total = suppressWarnings(cor(S_total_N_baseline, NNV_infection, method = "spearman")),
            cor_S_age12 = suppressWarnings(cor(S_age12_N_baseline, NNV_infection, method = "spearman")),
            cor_S_0_11 = suppressWarnings(cor(S_0_11_N_baseline, NNV_infection, method = "spearman")),
            cor_S_65plus = suppressWarnings(cor(S_65plus_N_baseline, NNV_infection, method = "spearman")))
write_csv(spearman_cors, file.path(DIR_TABLES_NNV, "table_nnv_susceptibility_spearman.csv"))
message("\n=== Descriptive Spearman correlation: baseline S/N vs NNV_infection (NOT causal) ===")
print(as.data.frame(spearman_cors), digits = 3)

message("\n[11_calculate_nnv] Done. All outputs under results/nnv/, tables/nnv/.")
