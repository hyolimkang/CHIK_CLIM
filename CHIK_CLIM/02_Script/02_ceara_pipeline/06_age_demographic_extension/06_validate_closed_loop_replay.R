# Age-structured vaccine extension -- Phases 1-2: CLOSED-LOOP no-vaccine
# forward replay validation.
#
# IMPORTANT CLARIFICATION: the age-cohort simulator (07_age_cohort_simulator.R)
# was ALREADY closed-loop when the earlier Phase-4 identity test was run --
# its `infectiousness` term is built from the simulator's OWN recursively
# generated X_total history (see the `for (g in ...) infectiousness <-
# infectiousness + w[g] * X_total[t-g]` loop), never from the original
# fit's X(t). The function signature never even accepts the fitted X(t)
# as an input. This script re-confirms that closed-loop property with the
# STRICTER validation battery now requested: infectiousness[t],
# force_of_infection[t] (lambda[t]) and R_eff[t] comparisons, era-specific
# error breakdowns, and an explicit first-divergence-week search.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(ggplot2); library(patchwork) })

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
dates <- as.Date(weekly$week_start)
years <- as.integer(format(dates, "%Y"))
N <- length(dates)
age_shares <- read_csv(file.path(table_dir, "CE_age_population_shares.csv"), show_col_types = FALSE)

sd_data <- climate$stan_data
w <- sd_data$w; imports_per_week <- sd_data$imports_per_week
is_seed <- sd_data$is_seed; X_seed <- sd_data$X_seed
births <- weekly$births; deaths <- weekly$all_cause_deaths
reconciliation <- weekly$net_population_reconciliation; N_start <- weekly$N_start

CE_2017_ONSET <- as.Date("2017-02-12"); CE_2017_END <- as.Date("2018-12-16")
CE_2022_SEED_START <- min(climate$seed_dates); CE_2022_SEED_END <- max(climate$seed_dates)
CE_2022_WAVE_END <- as.Date("2023-12-17")

era <- case_when(
  dates < CE_2017_ONSET ~ "pre-2017",
  dates >= CE_2017_ONSET & dates <= CE_2017_END ~ "2017 epidemic window",
  dates > CE_2017_END & dates < CE_2022_SEED_START ~ "inter-epidemic (2019-2021)",
  dates >= CE_2022_SEED_START & dates <= CE_2022_SEED_END ~ "2022 seed window",
  dates > CE_2022_SEED_END & dates <= CE_2022_WAVE_END ~ "2022 epidemic (post-seed)",
  dates > CE_2022_WAVE_END ~ "post-2022"
)

draws_v49 <- rstan::extract(fit, pars = c("R0_t", "S", "U", "X", "force_of_infection", "R_eff_t"), permuted = TRUE)
n_draws_total <- nrow(draws_v49$R0_t)
N_DRAWS <- 60L
set.seed(20260919L)
test_draw_idx <- sample.int(n_draws_total, N_DRAWS)
message(sprintf("Running closed-loop replay validation for %d / %d posterior draws...", N_DRAWS, n_draws_total))

results <- vector("list", N_DRAWS)
for (i in seq_along(test_draw_idx)) {
  d <- test_draw_idx[i]
  sim <- simulate_age_cohort(draws_v49$R0_t[d, ], births, deaths, reconciliation, N_start, w, imports_per_week,
                              is_seed, X_seed, age_shares, years, MAX_AGE, vaccination = NULL)

  # Reconstruct the ORIGINAL fit's infectiousness (not stored directly) from
  # its saved force_of_infection: infectiousness = (foi - imports/N) * N / R0.
  N_t_v49 <- N_start
  infectiousness_v49 <- (draws_v49$force_of_infection[d, ] - imports_per_week / N_t_v49) * N_t_v49 / draws_v49$R0_t[d, ]

  results[[i]] <- tibble(
    draw = d, t = seq_len(N), week_start = dates, era = era,
    S_v49 = draws_v49$S[d, ], U_v49 = draws_v49$U[d, ], X_v49 = draws_v49$X[d, ],
    infectiousness_v49 = infectiousness_v49, foi_v49 = draws_v49$force_of_infection[d, ], Reff_v49 = draws_v49$R_eff_t[d, ],
    S_age = sim$S_total, U_age = sim$Uinf_total, X_age = sim$X_total,
    infectiousness_age = sim$infectiousness, foi_age = sim$lambda, Reff_age = sim$R_eff
  )
  if (i %% 15 == 0) message(sprintf("  ...%d / %d draws done", i, N_DRAWS))
}
comparison <- bind_rows(results) |>
  mutate(
    abs_diff_S = abs(S_age - S_v49), rel_diff_S = abs_diff_S / pmax(S_v49, 1),
    abs_diff_U = abs(U_age - U_v49), rel_diff_U = abs_diff_U / pmax(U_v49, 1),
    abs_diff_X = abs(X_age - X_v49), rel_diff_X = abs_diff_X / pmax(X_v49, 1),
    abs_diff_infectiousness = abs(infectiousness_age - infectiousness_v49), rel_diff_infectiousness = abs_diff_infectiousness / pmax(abs(infectiousness_v49), 1e-9),
    abs_diff_foi = abs(foi_age - foi_v49), rel_diff_foi = abs_diff_foi / pmax(abs(foi_v49), 1e-12),
    abs_diff_Reff = abs(Reff_age - Reff_v49), rel_diff_Reff = abs_diff_Reff / pmax(abs(Reff_v49), 1e-9)
  )

# ================================================================
# Overall summary
# ================================================================
overall <- tibble(
  quantity = c("S", "U", "X", "infectiousness", "force_of_infection (lambda)", "R_eff"),
  max_abs_diff = c(max(comparison$abs_diff_S), max(comparison$abs_diff_U), max(comparison$abs_diff_X),
                    max(comparison$abs_diff_infectiousness), max(comparison$abs_diff_foi), max(comparison$abs_diff_Reff)),
  max_rel_diff = c(max(comparison$rel_diff_S), max(comparison$rel_diff_U), max(comparison$rel_diff_X),
                    max(comparison$rel_diff_infectiousness), max(comparison$rel_diff_foi), max(comparison$rel_diff_Reff))
)
message("\n=== Overall closed-loop replay error summary (across all draws x weeks) ===")
print(as.data.frame(overall), digits = 6)
write_csv(overall, file.path(table_dir, "PHASE1_2_closed_loop_overall_summary.csv"))

# ================================================================
# Era-specific breakdown
# ================================================================
era_summary <- comparison |> group_by(era) |>
  summarise(max_rel_diff_S = max(rel_diff_S), max_rel_diff_U = max(rel_diff_U), max_rel_diff_X = max(rel_diff_X),
            max_rel_diff_foi = max(rel_diff_foi), max_rel_diff_Reff = max(rel_diff_Reff), n_weeks = n_distinct(t), .groups = "drop")
message("\n=== Era-specific error breakdown (max relative difference within era) ===")
print(as.data.frame(era_summary), digits = 4)
write_csv(era_summary, file.path(table_dir, "PHASE1_2_closed_loop_era_summary.csv"))

# ================================================================
# Cumulative infections + final susceptible fraction
# ================================================================
cum_final <- comparison |> group_by(draw) |>
  summarise(cum_X_v49 = sum(X_v49), cum_X_age = sum(X_age),
            final_S_v49 = S_v49[which.max(t)], final_S_age = S_age[which.max(t)],
            final_N = N_start[length(N_start)], .groups = "drop") |>
  mutate(cum_rel_diff = (cum_X_age - cum_X_v49) / cum_X_v49,
         final_S_frac_v49 = final_S_v49 / final_N, final_S_frac_age = final_S_age / final_N,
         final_S_frac_diff = final_S_frac_age - final_S_frac_v49)
message(sprintf("\nCumulative infections: max |relative diff| across draws = %.2e", max(abs(cum_final$cum_rel_diff))))
message(sprintf("Final susceptible fraction: max |abs diff| across draws = %.2e", max(abs(cum_final$final_S_frac_diff))))
write_csv(cum_final, file.path(table_dir, "PHASE1_2_closed_loop_cumulative_final_summary.csv"))

# ================================================================
# First-divergence-week search (threshold: relative error > 1e-6)
# ================================================================
DIVERGENCE_THRESHOLD <- 1e-6
first_divergence <- comparison |> dplyr::filter(rel_diff_X > DIVERGENCE_THRESHOLD | rel_diff_S > DIVERGENCE_THRESHOLD | rel_diff_U > DIVERGENCE_THRESHOLD) |>
  arrange(t) |> slice(1)
if (nrow(first_divergence) == 0) {
  message(sprintf("\nNo week across any of the %d test draws exceeds the %.0e relative-error threshold for S/U/X.", N_DRAWS, DIVERGENCE_THRESHOLD))
  message("All observed error is floating-point accumulation noise -- NOT systematic drift.")
} else {
  message(sprintf("\nFirst week exceeding %.0e relative error: t=%d (%s), draw=%d", DIVERGENCE_THRESHOLD, first_divergence$t, first_divergence$week_start, first_divergence$draw))
  print(as.data.frame(first_divergence))
}
write_csv(comparison |> dplyr::filter(rel_diff_X > DIVERGENCE_THRESHOLD | rel_diff_S > DIVERGENCE_THRESHOLD | rel_diff_U > DIVERGENCE_THRESHOLD),
          file.path(table_dir, "PHASE1_2_divergence_weeks_above_threshold.csv"))

identity_pass <- all(overall$max_rel_diff < 1e-6)
message(sprintf("\n=== PHASE 1-2 CLOSED-LOOP REPLAY: %s (all quantities max relative diff < 1e-6) ===", if (identity_pass) "PASS" else "FAIL"))

# ================================================================
# Figures
# ================================================================
rep_draw <- test_draw_idx[1]
rep_data <- comparison |> dplyr::filter(draw == rep_draw)

p_X <- ggplot(rep_data, aes(week_start)) +
  geom_line(aes(y = X_v49, colour = "closed-loop age-model"), linewidth = 0.4, linetype = "dashed") +
  geom_line(aes(y = X_age, colour = "original v4.9 (fitted)"), linewidth = 0.8, alpha = 0.6) +
  scale_colour_manual(name = NULL, values = c("original v4.9 (fitted)" = "#0072B2", "closed-loop age-model" = "#D55E00")) +
  labs(title = "X(t): closed-loop replay vs. original fitted v4.9", x = NULL, y = "Weekly infections") + theme_v4

p_foi <- ggplot(rep_data, aes(week_start)) +
  geom_line(aes(y = foi_v49, colour = "closed-loop age-model"), linewidth = 0.4, linetype = "dashed") +
  geom_line(aes(y = foi_age, colour = "original v4.9 (fitted)"), linewidth = 0.8, alpha = 0.6) +
  scale_colour_manual(name = NULL, values = c("original v4.9 (fitted)" = "#0072B2", "closed-loop age-model" = "#D55E00")) +
  scale_y_log10() +
  labs(title = "Force of infection lambda(t) (log scale)", x = NULL, y = "lambda(t)") + theme_v4

p_Reff <- ggplot(rep_data, aes(week_start)) +
  geom_line(aes(y = Reff_v49, colour = "closed-loop age-model"), linewidth = 0.4, linetype = "dashed") +
  geom_line(aes(y = Reff_age, colour = "original v4.9 (fitted)"), linewidth = 0.8, alpha = 0.6) +
  geom_hline(yintercept = 1, linetype = "22", colour = "grey50") +
  scale_colour_manual(name = NULL, values = c("original v4.9 (fitted)" = "#0072B2", "closed-loop age-model" = "#D55E00")) +
  labs(title = "R_eff(t) = R0(t) x S(t)/N(t)", x = NULL, y = "R_eff") + theme_v4

p_rel_all <- ggplot(comparison, aes(week_start, rel_diff_X, group = draw)) +
  geom_line(alpha = 0.12, colour = "#D55E00") +
  scale_y_log10(labels = scales::label_scientific()) +
  labs(title = "Relative difference in X(t), all test draws (log scale)", subtitle = sprintf("Threshold line at %.0e", DIVERGENCE_THRESHOLD),
       x = NULL, y = "|X_age - X_v49| / X_v49") +
  geom_hline(yintercept = DIVERGENCE_THRESHOLD, linetype = "dashed", colour = "black") + theme_v4

figure <- (p_X | p_foi) / (p_Reff | p_rel_all) +
  patchwork::plot_annotation(title = "Phase 1-2: closed-loop no-vaccine forward replay validation",
                              subtitle = sprintf("%s across %d draws", if (identity_pass) "PASS" else "FAIL", N_DRAWS))
ggsave(file.path(figure_dir, "PHASE1_2_closed_loop_replay_validation.png"), figure, width = 240, height = 180, units = "mm", dpi = 300, bg = "white")
message("\n[saved] PHASE1_2_closed_loop_replay_validation.png")
message("[saved] PHASE1_2_closed_loop_overall_summary.csv, PHASE1_2_closed_loop_era_summary.csv, PHASE1_2_closed_loop_cumulative_final_summary.csv")
