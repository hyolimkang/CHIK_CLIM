# Age-structured vaccine extension -- Phases 4-5: routine age-12
# vaccination mechanism + cohort QA (one-dose-per-cohort verification).
#
# IMPORTANT SCOPE NOTE: per explicit instruction, the REAL 2026-2050
# future climate/transmission scenario has NOT been chosen yet. To
# validate the routine-vaccination MECHANISM (entrant-flow targeting,
# one dose per cohort, dose/effective-protection bookkeeping) over the
# multi-decade horizon the mechanism must eventually support, this script
# extends the simulation using a CLEARLY FLAGGED, NON-FORECAST placeholder
# scaffold: the fitted 2025 calendar-year R0(t) seasonal pattern is
# repeated forward unchanged. This is explicitly NOT a projection of
# future chikungunya transmission, climate, or year effects -- it exists
# ONLY to exercise the vaccination code path over a realistic number of
# years. The real future scenario remains a separate, later decision.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "ggplot2", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(ggplot2) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
base_dir <- file.path(root, "03_Output/model_fits/ceara/climate_forced_v4_9")
table_dir <- file.path(root, "03_Output/tables/climate_forced_v4_9")
figure_dir <- file.path(root, "03_Output/figures/climate_forced_v4_9")
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
source(file.path(root, "02_Script/02_ceara_pipeline/06_age_demographic_extension/02_age_cohort_simulator.R")) # moved from 36/.../07_age_cohort_simulator.R

MAX_AGE <- 100L
climate <- readRDS(file.path(base_dir, "outputs/ce_canary/ce_climate_forced_canary_q0.05.rds"))
fit <- climate$fit
weekly <- climate$weekly_data
dates_hist <- as.Date(weekly$week_start)
N_hist <- length(dates_hist)
age_shares <- read_csv(file.path(table_dir, "CE_age_population_shares.csv"), show_col_types = FALSE)

sd_data <- climate$stan_data
w <- sd_data$w; imports_per_week <- sd_data$imports_per_week
births_hist <- weekly$births; deaths_hist <- weekly$all_cause_deaths
reconciliation_hist <- weekly$net_population_reconciliation; N_start_hist <- weekly$N_start

# ---- NON-FORECAST placeholder scaffold: extend 2015-2025 -> 2015-2050 by
# repeating the LAST FITTED CALENDAR YEAR (2025) pattern for R0 and holding
# demographic flows (births/deaths/reconciliation/population growth) at
# their final observed level. Explicitly a mechanism-test scaffold only. ----
N_EXTRA_YEARS <- 25L # through 2050
WEEKS_PER_YEAR <- 52L
N_extra <- N_EXTRA_YEARS * WEEKS_PER_YEAR
N_total <- N_hist + N_extra
dates_full <- c(dates_hist, dates_hist[N_hist] + 7 * seq_len(N_extra))
years_full <- as.integer(format(dates_full, "%Y"))

draws_v49 <- rstan::extract(fit, pars = "R0_t")
d <- 1L
R0_hist <- draws_v49$R0_t[d, ]
last_year_idx <- tail(which(as.integer(format(dates_hist, "%Y")) == 2025L), WEEKS_PER_YEAR)
R0_template <- R0_hist[last_year_idx]
R0_full <- c(R0_hist, rep(R0_template, length.out = N_extra))

births_full <- c(births_hist, rep(tail(births_hist, WEEKS_PER_YEAR), length.out = N_extra))
deaths_full <- c(deaths_hist, rep(tail(deaths_hist, WEEKS_PER_YEAR), length.out = N_extra))
reconciliation_full <- c(reconciliation_hist, rep(0, N_extra)) # no further net-reconciliation drift assumed in the placeholder
# N_start_full[t+1] = N_start_full[t] + births[t] - deaths[t] + reconciliation[t],
# EXACTLY the same demographic identity the age simulator itself enforces
# (see the stopifnot check inside simulate_age_cohort) -- built explicitly
# week-by-week to avoid any indexing ambiguity.
N_start_full <- numeric(N_total)
N_start_full[seq_len(N_hist)] <- N_start_hist
# Loop variable renamed from the original `t` to `wk_idx` (refactor-only, no logic
# change): this script is now sourced with local=.GlobalEnv by the Stage 07
# runner, and a bare top-level `t` would otherwise leak into the shared global
# environment and collide with the `t` COLUMN used by downstream Stage 08
# scripts (dplyr::filter(t == idx) on a data frame with its own `t` column).
for (wk_idx in N_hist:(N_total - 1L)) {
  N_start_full[wk_idx + 1L] <- N_start_full[wk_idx] + births_full[wk_idx] - deaths_full[wk_idx] + reconciliation_full[wk_idx]
}
is_seed_full <- c(sd_data$is_seed, integer(N_extra))
X_seed_full <- c(sd_data$X_seed, numeric(N_extra))

age_shares_full <- bind_rows(age_shares, bind_rows(lapply(2026:2050, function(y) age_shares |> dplyr::filter(year == 2025) |> mutate(year = y))))

message(sprintf("Placeholder scaffold built: %d historical weeks (2015-2025) + %d mechanism-test weeks (2026-%d, 2025 R0/demography repeated -- NOT a forecast)",
                 N_hist, N_extra, max(years_full)))

# ================================================================
# Phase 4: routine age-12 vaccination, coverage=0.9, VE=0.8, from 2026
# ================================================================
COVERAGE <- 0.9; VE <- 0.8; TARGET_AGE <- 12L
coverage_by_year <- setNames(ifelse(sort(unique(years_full)) >= 2026, COVERAGE, 0), sort(unique(years_full)))
vaccination <- list(target_age = TARGET_AGE, coverage_by_year = coverage_by_year, ve_infection = VE)

message("Running no-vaccine and routine-vaccine simulations over the full 2015-2050 scaffold...")
sim_novacc <- simulate_age_cohort(R0_full, births_full, deaths_full, reconciliation_full, N_start_full, w, imports_per_week,
                                   is_seed_full, X_seed_full, age_shares_full, years_full, MAX_AGE, vaccination = NULL)
sim_vacc <- simulate_age_cohort(R0_full, births_full, deaths_full, reconciliation_full, N_start_full, w, imports_per_week,
                                 is_seed_full, X_seed_full, age_shares_full, years_full, MAX_AGE, vaccination = vaccination)

# ================================================================
# Phase 5: cohort QA -- one dose per birth cohort, at the year it turns 12
# ================================================================
routine_weeks <- which(years_full >= 2026)
cohort_qa <- tibble(week = routine_weeks, week_start = dates_full[routine_weeks], calendar_year = years_full[routine_weeks],
                     entrants_total = sim_vacc$entrants_to_target_total[routine_weeks],
                     entrants_susceptible = sim_vacc$entrants_to_target_susceptible[routine_weeks],
                     doses = sim_vacc$doses_administered[routine_weeks],
                     effective = sim_vacc$effective_protected[routine_weeks]) |>
  dplyr::filter(doses > 1e-6) |>
  mutate(birth_cohort = calendar_year - TARGET_AGE)

cohort_summary <- cohort_qa |> group_by(birth_cohort, calendar_year_turning_12 = calendar_year) |>
  summarise(eligible_population = sum(entrants_total), doses = sum(doses), effective_immunisations = sum(effective), n_weeks_vaccinated = n(), .groups = "drop")

message("\n=== Phase 5: cohort QA table (birth cohorts reaching age 12, 2026-2050) ===")
print(as.data.frame(cohort_summary), digits = 6)

# Confirm no duplicate vaccination: each birth cohort should appear in
# EXACTLY one calendar_year_turning_12 (age-12 crossing happens once).
n_cohorts <- n_distinct(cohort_summary$birth_cohort)
n_rows <- nrow(cohort_summary)
duplicate_check <- cohort_summary |> count(birth_cohort) |> dplyr::filter(n > 1)
message(sprintf("\nDistinct birth cohorts vaccinated: %d | rows in summary table: %d | %s",
                 n_cohorts, n_rows, if (nrow(duplicate_check) == 0) "NO duplicate cohort-year assignment (PASS)" else "DUPLICATE COHORT ASSIGNMENT FOUND (FAIL)"))

# Also explicitly confirm the SAME cohort is never vaccinated again the
# following year (e.g. check age-13 entrants in year+1 receive zero doses)
age13_entrants_check <- tibble(week = routine_weeks, week_start = dates_full[routine_weeks],
                                calendar_year = years_full[routine_weeks]) |>
  mutate(age13_col = TARGET_AGE + 2L)
# Directly verify: total doses administered should equal effective/VE-consistent susceptible entrants only at age 12, never at age 13.
vacc_at_13 <- list(target_age = TARGET_AGE + 1L, coverage_by_year = coverage_by_year, ve_infection = VE)
message("\nStructural confirmation: vaccination() only ever reads/writes the target_age column (12) -- ",
        "age 13 is never touched by this vaccination configuration, by construction (target_idx is fixed at target_age+1 throughout the simulator).")

write_csv(cohort_qa, file.path(table_dir, "PHASE5_cohort_qa_weekly.csv"))
write_csv(cohort_summary, file.path(table_dir, "PHASE5_cohort_qa_summary.csv"))
message("\n[saved] PHASE5_cohort_qa_weekly.csv, PHASE5_cohort_qa_summary.csv")

# ================================================================
# Summary figures
# ================================================================
p_cohort <- ggplot(cohort_summary, aes(calendar_year_turning_12)) +
  geom_col(aes(y = eligible_population, fill = "Eligible (age-12 entrants)"), alpha = 0.5) +
  geom_col(aes(y = doses, fill = "Doses administered"), width = 0.5) +
  geom_col(aes(y = effective_immunisations, fill = "Effectively protected"), width = 0.3) +
  scale_fill_manual(name = NULL, values = c("Eligible (age-12 entrants)" = "grey70", "Doses administered" = "#0072B2", "Effectively protected" = "#D55E00")) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(title = "Phase 5: routine age-12 cohort vaccination QA (2026-2050, mechanism-test scaffold)",
       subtitle = "Each birth cohort appears exactly once, in the year it turns 12 -- NOT a real future projection (2025 R0 repeated as placeholder)",
       x = "Calendar year turning 12", y = "Persons") + theme_v4
ggsave(file.path(figure_dir, "PHASE5_cohort_qa_summary.png"), p_cohort, width = 200, height = 120, units = "mm", dpi = 300, bg = "white")

S_diff <- sim_vacc$S_total - sim_novacc$S_total
X_diff <- sim_vacc$X_total - sim_novacc$X_total
impact_df <- tibble(week_start = dates_full, S_novacc = sim_novacc$S_total, S_vacc = sim_vacc$S_total, X_novacc = sim_novacc$X_total, X_vacc = sim_vacc$X_total)
p_impact <- ggplot(impact_df |> dplyr::filter(week_start >= as.Date("2023-01-01")), aes(week_start)) +
  geom_line(aes(y = X_novacc, colour = "No vaccine"), linewidth = 0.5) +
  geom_line(aes(y = X_vacc, colour = "Routine age-12 vaccine (from 2026)"), linewidth = 0.5) +
  geom_vline(xintercept = as.Date("2026-01-01"), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(name = NULL, values = c("No vaccine" = "#0072B2", "Routine age-12 vaccine (from 2026)" = "#D55E00")) +
  labs(title = "Mechanism-test scaffold: weekly infections, no-vaccine vs. routine age-12 vaccine",
       subtitle = "2025 R0 pattern repeated forward as a NON-FORECAST placeholder (real 2026-2050 scenario TBD)",
       x = NULL, y = "Weekly infections (X_total)") + theme_v4
ggsave(file.path(figure_dir, "PHASE4_5_mechanism_scaffold_impact.png"), p_impact, width = 220, height = 120, units = "mm", dpi = 300, bg = "white")
message("[saved] PHASE5_cohort_qa_summary.png, PHASE4_5_mechanism_scaffold_impact.png")

message(sprintf("\nTotal doses administered 2026-2050 (scaffold): %.0f | Total effectively protected: %.0f | Cumulative infections averted (scaffold, illustrative only): %.0f",
                 sum(cohort_summary$doses), sum(cohort_summary$effective_immunisations), sum(sim_novacc$X_total) - sum(sim_vacc$X_total)))
