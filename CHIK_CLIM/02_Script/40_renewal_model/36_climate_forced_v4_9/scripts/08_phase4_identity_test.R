# Age-structured vaccine extension -- Phase 4: THE critical identity test.
# Run the age-cohort simulator (vaccination = NULL, i.e. OFF) for many
# posterior draws of the accepted climate-forced v4.9 CE fit, and compare
# the aggregated S_total(t)/Uinf_total(t)/X_total(t) to the SAME draws'
# original S[t]/U[t]/X[t]. No refitting; this is pure bookkeeping
# verification.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(ggplot2); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "02_Script/40_renewal_model/36_climate_forced_v4_9")
table_dir <- file.path(root, "03_Output/tables/climate_forced_v4_9")
figure_dir <- file.path(root, "03_Output/figures/climate_forced_v4_9")
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
source(file.path(base_dir, "scripts/07_age_cohort_simulator.R"))

climate <- readRDS(file.path(base_dir, "outputs/ce_canary/ce_climate_forced_canary_q0.05.rds"))
fit <- climate$fit
weekly <- climate$weekly_data
dates <- as.Date(weekly$week_start)
years <- as.integer(format(dates, "%Y"))
N <- length(dates)

age_shares <- read_csv(file.path(table_dir, "CE_age_population_shares.csv"), show_col_types = FALSE)
MAX_AGE <- 100L

sd_data <- climate$stan_data
w <- sd_data$w
imports_per_week <- sd_data$imports_per_week
is_seed <- sd_data$is_seed
X_seed <- sd_data$X_seed
births <- weekly$births
deaths <- weekly$all_cause_deaths
reconciliation <- weekly$net_population_reconciliation
N_start <- weekly$N_start

draws <- rstan::extract(fit, pars = c("R0_t", "S", "U", "X"), permuted = TRUE)
n_draws_total <- nrow(draws$R0_t)
N_TEST_DRAWS <- 60L
set.seed(20260917L)
test_draw_idx <- sample.int(n_draws_total, N_TEST_DRAWS)
message(sprintf("Running age-cohort identity test for %d / %d posterior draws...", N_TEST_DRAWS, n_draws_total))

results <- vector("list", N_TEST_DRAWS)
for (i in seq_along(test_draw_idx)) {
  d <- test_draw_idx[i]
  sim <- simulate_age_cohort(
    R0_t_draw = draws$R0_t[d, ], births = births, deaths = deaths, reconciliation = reconciliation,
    N_start = N_start, w = w, imports_per_week = imports_per_week, is_seed = is_seed, X_seed = X_seed,
    age_shares = age_shares, years = years, max_age = MAX_AGE, vaccination = NULL
  )
  results[[i]] <- tibble(
    draw = d, t = seq_len(N), week_start = dates,
    S_v49 = draws$S[d, ], U_v49 = draws$U[d, ], X_v49 = draws$X[d, ],
    S_age = sim$S_total, U_age = sim$Uinf_total, X_age = sim$X_total
  )
  if (i %% 10 == 0) message(sprintf("  ...%d / %d draws done", i, N_TEST_DRAWS))
}
comparison <- bind_rows(results) |>
  mutate(abs_diff_S = abs(S_age - S_v49), abs_diff_U = abs(U_age - U_v49), abs_diff_X = abs(X_age - X_v49),
         rel_diff_S = abs_diff_S / pmax(S_v49, 1), rel_diff_U = abs_diff_U / pmax(U_v49, 1), rel_diff_X = abs_diff_X / pmax(X_v49, 1))

message("\n=== Phase 4 identity test: overall difference summary (across all draws x weeks) ===")
overall <- tibble(
  quantity = c("S", "U", "X"),
  max_abs_diff = c(max(comparison$abs_diff_S), max(comparison$abs_diff_U), max(comparison$abs_diff_X)),
  max_rel_diff = c(max(comparison$rel_diff_S), max(comparison$rel_diff_U), max(comparison$rel_diff_X)),
  mean_abs_diff = c(mean(comparison$abs_diff_S), mean(comparison$abs_diff_U), mean(comparison$abs_diff_X))
)
print(as.data.frame(overall), digits = 6)
write_csv(overall, file.path(table_dir, "PHASE4_identity_test_overall_summary.csv"))

# ---- Yearly checkpoints (posterior median across test draws, at end of each year) ----
checkpoint_dates <- as.Date(paste0(2015:2025, "-12-31"))
checkpoint_idx <- sapply(checkpoint_dates, function(cd) which.min(abs(dates - cd)))
checkpoint_tbl <- bind_rows(lapply(seq_along(checkpoint_dates), function(i) {
  t_idx <- checkpoint_idx[i]
  sub <- comparison |> filter(t == t_idx)
  tibble(year = (2015:2025)[i], week_start = dates[t_idx],
         median_abs_diff_S = median(sub$abs_diff_S), median_rel_diff_S = median(sub$rel_diff_S),
         median_abs_diff_U = median(sub$abs_diff_U), median_rel_diff_U = median(sub$rel_diff_U))
}))
message("\n=== Yearly checkpoint differences (median across test draws) ===")
print(as.data.frame(checkpoint_tbl), digits = 4)
write_csv(checkpoint_tbl, file.path(table_dir, "PHASE4_identity_test_yearly_checkpoints.csv"))

# ---- Cumulative infection difference (final week, X summed over time) ----
cum_diff <- comparison |> group_by(draw) |>
  summarise(cum_X_v49 = sum(X_v49), cum_X_age = sum(X_age), .groups = "drop") |>
  mutate(cum_diff = cum_X_age - cum_X_v49, cum_rel_diff = cum_diff / cum_X_v49)
message(sprintf("\nCumulative infection (sum of X over 2015-2025): max |relative difference| across draws = %.2e", max(abs(cum_diff$cum_rel_diff))))
write_csv(cum_diff, file.path(table_dir, "PHASE4_identity_test_cumulative_infection_diff.csv"))

pass_threshold <- 1e-6
identity_pass <- overall$max_rel_diff[overall$quantity == "S"] < pass_threshold &&
  overall$max_rel_diff[overall$quantity == "U"] < pass_threshold &&
  overall$max_rel_diff[overall$quantity == "X"] < pass_threshold
message(sprintf("\n=== PHASE 4 IDENTITY TEST: %s (threshold: max relative difference < %.0e) ===",
                 if (identity_pass) "PASS" else "FAIL", pass_threshold))

# ---- Comparison plot (one representative draw + envelope across all test draws) ----
rep_draw <- test_draw_idx[1]
rep_data <- comparison |> filter(draw == rep_draw)
p_S <- ggplot(rep_data, aes(week_start)) +
  geom_line(aes(y = S_v49, colour = "v4.9 (original)"), linewidth = 0.8) +
  geom_line(aes(y = S_age, colour = "age-expanded (vax OFF)"), linewidth = 0.4, linetype = "dashed") +
  scale_colour_manual(name = NULL, values = c("v4.9 (original)" = "#0072B2", "age-expanded (vax OFF)" = "#D55E00")) +
  labs(title = "Phase 4: S(t) identity check (1 representative draw)", x = NULL, y = "S (absolute count)") + theme_v4
p_U <- ggplot(rep_data, aes(week_start)) +
  geom_line(aes(y = U_v49, colour = "v4.9 (original)"), linewidth = 0.8) +
  geom_line(aes(y = U_age, colour = "age-expanded (vax OFF)"), linewidth = 0.4, linetype = "dashed") +
  scale_colour_manual(name = NULL, values = c("v4.9 (original)" = "#0072B2", "age-expanded (vax OFF)" = "#D55E00")) +
  labs(title = "U(t) identity check", x = NULL, y = "U (absolute count)") + theme_v4
p_X <- ggplot(rep_data, aes(week_start)) +
  geom_line(aes(y = X_v49, colour = "v4.9 (original)"), linewidth = 0.8) +
  geom_line(aes(y = X_age, colour = "age-expanded (vax OFF)"), linewidth = 0.4, linetype = "dashed") +
  scale_colour_manual(name = NULL, values = c("v4.9 (original)" = "#0072B2", "age-expanded (vax OFF)" = "#D55E00")) +
  labs(title = "X(t) (latent infections) identity check", x = NULL, y = "X (weekly infections)") + theme_v4
p_reldiff <- ggplot(comparison, aes(week_start, rel_diff_X)) +
  geom_line(aes(group = draw), alpha = 0.15, colour = "#D55E00") +
  scale_y_log10(labels = scales::label_number()) +
  labs(title = "Relative difference in X(t) across all test draws (log scale)", x = NULL, y = "|X_age - X_v49| / X_v49") + theme_v4

figure <- (p_S | p_U) / (p_X | p_reldiff) +
  patchwork::plot_annotation(title = "Phase 4 identity test: age-expanded (vaccination OFF) vs. accepted climate-v4.9",
                              subtitle = sprintf("%s across %d test draws (threshold 1e-6)", if (identity_pass) "PASS" else "FAIL", N_TEST_DRAWS))
ggsave(file.path(figure_dir, "PHASE4_identity_test_comparison.png"), figure, width = 220, height = 180, units = "mm", dpi = 300, bg = "white")
message("\n[saved] ", file.path(figure_dir, "PHASE4_identity_test_comparison.png"))
message("[saved] PHASE4_identity_test_overall_summary.csv, PHASE4_identity_test_yearly_checkpoints.csv, PHASE4_identity_test_cumulative_infection_diff.csv")
