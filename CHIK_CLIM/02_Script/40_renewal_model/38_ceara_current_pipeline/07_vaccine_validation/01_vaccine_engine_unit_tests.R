# Age-structured vaccine extension -- Phase 3: vaccine-engine unit tests.
# Uses ONE representative posterior draw (mechanism tests do not require
# posterior uncertainty) over the SAME historical 2015-2025 window, purely
# to exercise the vaccination code path -- NOT a claim about real
# historical vaccination.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(ggplot2); library(patchwork) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
base_dir <- file.path(root, "02_Script/40_renewal_model/36_climate_forced_v4_9")
table_dir <- file.path(root, "03_Output/tables/climate_forced_v4_9")
figure_dir <- file.path(root, "03_Output/figures/climate_forced_v4_9")
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
source(file.path(root, "02_Script/40_renewal_model/38_ceara_current_pipeline/04_age_reconstruction/02_age_cohort_simulator.R")) # moved from 36/.../07_age_cohort_simulator.R

MAX_AGE <- 100L
climate <- readRDS(file.path(base_dir, "outputs/ce_canary/ce_climate_forced_canary_q0.05.rds"))
fit <- climate$fit
weekly <- climate$weekly_data
dates <- as.Date(weekly$week_start)
years <- as.integer(format(dates, "%Y"))
N <- length(dates)
age_shares <- read_csv(file.path(table_dir, "CE_age_population_shares.csv"), show_col_types = FALSE)

sd_data <- climate$stan_data
w <- sd_data$w; imports_per_week <- sd_data$imports_per_week
is_seed <- sd_data$is_seed; X_seed <- sd_data$X_seed
births <- weekly$births; deaths <- weekly$all_cause_deaths
reconciliation <- weekly$net_population_reconciliation; N_start <- weekly$N_start

draws_v49 <- rstan::extract(fit, pars = "R0_t")
d <- 1L
R0_draw <- draws_v49$R0_t[d, ]

run <- function(vaccination) simulate_age_cohort(R0_draw, births, deaths, reconciliation, N_start, w, imports_per_week,
                                                  is_seed, X_seed, age_shares, years, MAX_AGE, vaccination = vaccination)

message("=== Baseline (no-vaccine, vaccination = NULL) ===")
sim_baseline <- run(NULL)

# ================================================================
# TEST A: coverage = 0 everywhere -- must be IDENTICAL to no-vaccine
# ================================================================
vacc_cov0 <- list(target_age = 12L, coverage_by_year = setNames(rep(0, length(unique(years))), unique(years)), ve_infection = 0.8)
sim_testA <- run(vacc_cov0)
diff_A_S <- max(abs(sim_testA$S_total - sim_baseline$S_total))
diff_A_X <- max(abs(sim_testA$X_total - sim_baseline$X_total))
message(sprintf("TEST A (coverage=0): max|diff S_total|=%.3e, max|diff X_total|=%.3e -- %s",
                 diff_A_S, diff_A_X, if (diff_A_S < 1e-6 && diff_A_X < 1e-6) "PASS (identical)" else "FAIL"))

# ================================================================
# TEST B: VE_infection = 0 -- epidemiologically IDENTICAL (doses counted, no protection)
# ================================================================
vacc_ve0 <- list(target_age = 12L, coverage_by_year = setNames(rep(0.9, length(unique(years))), unique(years)), ve_infection = 0)
sim_testB <- run(vacc_ve0)
diff_B_S <- max(abs(sim_testB$S_total - sim_baseline$S_total))
diff_B_X <- max(abs(sim_testB$X_total - sim_baseline$X_total))
doses_B <- sum(sim_testB$doses_administered)
effective_B <- sum(sim_testB$effective_protected)
message(sprintf("TEST B (VE=0, coverage=0.9): max|diff S_total|=%.3e, max|diff X_total|=%.3e, doses=%.0f, effective_protected=%.4f -- %s",
                 diff_B_S, diff_B_X, doses_B, effective_B, if (diff_B_S < 1e-6 && diff_B_X < 1e-6 && doses_B > 0 && effective_B < 1e-9) "PASS (epi-identical, doses counted)" else "FAIL"))

unit_test_summary <- tibble(
  test = c("A: coverage=0", "B: VE_infection=0"),
  max_abs_diff_S = c(diff_A_S, diff_B_S), max_abs_diff_X = c(diff_A_X, diff_B_X),
  doses_administered_total = c(sum(sim_testA$doses_administered), doses_B),
  effective_protected_total = c(sum(sim_testA$effective_protected), effective_B),
  result = c(if (diff_A_S < 1e-6 && diff_A_X < 1e-6) "PASS" else "FAIL",
             if (diff_B_S < 1e-6 && diff_B_X < 1e-6 && doses_B > 0 && effective_B < 1e-9) "PASS" else "FAIL")
)
write_csv(unit_test_summary, file.path(table_dir, "PHASE3_unit_tests_A_B_summary.csv"))

# ================================================================
# TEST C: one-time effective pulse at age 12, ONE week only -- direct + indirect effect
# ================================================================
PULSE_WEEK <- which(dates == as.Date("2018-06-01")) # a quiet inter-epidemic week, well before the 2022 seed
if (length(PULSE_WEEK) == 0) PULSE_WEEK <- 180L
coverage_by_week <- rep(0, N); coverage_by_week[PULSE_WEEK] <- 1.0 # 100% coverage, ONE week only, at age 12 entrants
vacc_pulse <- list(target_age = 12L, coverage_by_week = coverage_by_week, ve_infection = 0.8)
sim_pulse <- run(vacc_pulse)

message(sprintf("\n=== TEST C: one-time pulse at age 12, week %d (%s) ===", PULSE_WEEK, dates[PULSE_WEEK]))
# NOTE on indexing: coverage_by_week[PULSE_WEEK] triggers vaccination of the
# age11->12 flow COMPUTED while processing week (PULSE_WEEK-1) -> arrives
# and is recorded at index PULSE_WEEK (doses_administered[t+1] where
# t+1 = PULSE_WEEK), matching where the S_age12 drop is actually observed.
message(sprintf("Entrants to age 12 that week: total=%.2f, susceptible=%.2f, doses=%.2f, effective_protected=%.2f",
                 sim_pulse$entrants_to_target_total[PULSE_WEEK], sim_pulse$entrants_to_target_susceptible[PULSE_WEEK],
                 sim_pulse$doses_administered[PULSE_WEEK], sim_pulse$effective_protected[PULSE_WEEK]))

# same-week force of infection must be UNCHANGED by vaccination (built from PAST infections only)
lambda_diff_same_week <- sim_pulse$lambda[PULSE_WEEK] - sim_baseline$lambda[PULSE_WEEK]
message(sprintf("Same-week lambda(t) difference at pulse week (must be ~0, generated from PAST infections): %.3e", lambda_diff_same_week))

# age-12 (age index 13) susceptibles should drop starting the week AFTER the pulse (since vaccination acts on the age11->12 inflow landing at t+1)
age12_col <- 12L + 1L
S_age12_diff <- sim_pulse$S_age[, age12_col] - sim_baseline$S_age[, age12_col]
message(sprintf("Age-12 susceptibles: max drop vs baseline = %.2f, first week of drop = %s",
                 min(S_age12_diff), dates[which(abs(S_age12_diff) > 1e-6)[1]]))

# X_total (and hence infectiousness/lambda in FOLLOWING weeks) should decline after the pulse
X_diff <- sim_pulse$X_total - sim_baseline$X_total
first_X_decline <- which(X_diff < -1e-9 & seq_len(N) > PULSE_WEEK)[1]
message(sprintf("First week X_total declines vs baseline (indirect effect via reduced infectiousness): t=%d (%s), %d weeks after the pulse",
                 first_X_decline, dates[first_X_decline], first_X_decline - PULSE_WEEK))

age0_11_cols <- 1:12; age65plus_cols <- 66:101
Inew_0_11_diff <- rowSums(sim_pulse$Inew_age[, age0_11_cols, drop = FALSE]) - rowSums(sim_baseline$Inew_age[, age0_11_cols, drop = FALSE])
Inew_65plus_diff <- rowSums(sim_pulse$Inew_age[, age65plus_cols, drop = FALSE]) - rowSums(sim_baseline$Inew_age[, age65plus_cols, drop = FALSE])
Inew_age12_diff <- sim_pulse$Inew_age[, age12_col] - sim_baseline$Inew_age[, age12_col]

test_C_summary <- tibble(week_start = dates, t = seq_len(N),
                          Inew_age12_diff = Inew_age12_diff, Inew_0_11_diff = Inew_0_11_diff, Inew_65plus_diff = Inew_65plus_diff, X_total_diff = X_diff)
write_csv(test_C_summary, file.path(table_dir, "PHASE3_testC_pulse_diagnostic.csv"))

message(sprintf("Cumulative infections averted, age 12 (direct): %.2f", -sum(Inew_age12_diff)))
message(sprintf("Cumulative infections averted, ages 0-11 (indirect, NEVER vaccinated): %.2f", -sum(Inew_0_11_diff)))
message(sprintf("Cumulative infections averted, ages 65+ (indirect, NEVER vaccinated): %.2f", -sum(Inew_65plus_diff)))

# ---- Figure: direct effect (age 12) + indirect effects (0-11, 65+) ----
plot_window <- which(dates >= PULSE_WEEK - 10 & seq_len(N) >= max(1, PULSE_WEEK - 10) & seq_len(N) <= min(N, PULSE_WEEK + 150))
plot_df <- test_C_summary |> filter(t >= max(1, PULSE_WEEK - 10), t <= min(N, PULSE_WEEK + 150))

p_direct <- ggplot(plot_df, aes(week_start, -Inew_age12_diff)) +
  geom_col(fill = "#D55E00", width = 5) + geom_vline(xintercept = dates[PULSE_WEEK], linetype = "dashed") +
  labs(title = "Direct effect: infections averted, age 12 (vaccinated cohort)", x = NULL, y = "Infections averted/week") + theme_v4

p_indirect_young <- ggplot(plot_df, aes(week_start, -Inew_0_11_diff)) +
  geom_col(fill = "#0072B2", width = 5) + geom_vline(xintercept = dates[PULSE_WEEK], linetype = "dashed") +
  labs(title = "Indirect effect: infections averted, ages 0-11 (never vaccinated)", x = NULL, y = "Infections averted/week") + theme_v4

p_indirect_old <- ggplot(plot_df, aes(week_start, -Inew_65plus_diff)) +
  geom_col(fill = "#009E73", width = 5) + geom_vline(xintercept = dates[PULSE_WEEK], linetype = "dashed") +
  labs(title = "Indirect effect: infections averted, ages 65+ (never vaccinated)", x = NULL, y = "Infections averted/week") + theme_v4

p_X_total <- ggplot(plot_df, aes(week_start, -X_total_diff)) +
  geom_col(fill = "grey40", width = 5) + geom_vline(xintercept = dates[PULSE_WEEK], linetype = "dashed") +
  labs(title = "Total infections averted, all ages", subtitle = "Dashed line: pulse week", x = NULL, y = "Infections averted/week") + theme_v4

figure_C <- (p_direct | p_indirect_young) / (p_indirect_old | p_X_total) +
  patchwork::plot_annotation(title = "Phase 3 TEST C: one-time age-12 vaccination pulse -- direct + indirect (herd) effects",
                              subtitle = sprintf("100%% coverage, VE=0.8, ONE week only (%s). Ages 0-11 and 65+ are NEVER vaccinated -- any reduction there is indirect.", dates[PULSE_WEEK]))
ggsave(file.path(figure_dir, "PHASE3_testC_pulse_direct_indirect.png"), figure_C, width = 220, height = 170, units = "mm", dpi = 300, bg = "white")
message("\n[saved] PHASE3_testC_pulse_direct_indirect.png")
message("[saved] PHASE3_unit_tests_A_B_summary.csv, PHASE3_testC_pulse_diagnostic.csv")
