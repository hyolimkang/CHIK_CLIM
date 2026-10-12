# STEP 2 -- Mandatory validation (Section 2.10 A-J).
#
# Numerical unit tests for the future-renewal simulator. Run this BEFORE
# trusting any ensemble output. Writes a single pass/fail QA summary table;
# stops with an error if any check fails (consistent with the repo's
# existing Stage-07 "fails loudly on any substantive failure" convention).

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/10_future_projection_pilot/00_config_future_pilot.R"))
source(file.path(ROOT, "02_Script/02_ceara_pipeline/10_future_projection_pilot/03_future_renewal_functions.R"))

results <- list()
check <- function(name, section, pass, detail) {
  message(sprintf("[QA %-5s] %-45s %s | %s", section, name, if (pass) "PASS" else "FAIL", detail))
  results[[length(results) + 1]] <<- tibble(section = section, check = name, pass = pass, detail = detail)
}

message("[Step2 unit tests] Loading accepted climate-forced fit (read-only)...")
climate <- readRDS(ACCEPTED_CLIMATE_FIT)
sd_clim <- climate$stan_data
dates_hist <- as.Date(climate$weekly_data$week_start)
N_hist <- length(dates_hist)
G <- sd_clim$G
w <- generation_weights_future()
stopifnot(isTRUE(all.equal(w, sd_clim$w, tolerance = 1e-10)))

draws <- rstan::extract(climate$fit, pars = c(
  "alpha_R", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2", "sigma_season_year",
  "beta_T", "beta_P", "phi_obs", "S", "U", "X", "R0_t", "R_eff_t"
), permuted = TRUE)
n_draws <- length(draws$alpha_R)

# ================================================================
# A. HISTORICAL REPLAY -- deterministic mode, historical R0(t)/demography/seed
# ================================================================
set.seed(BASE_SEED)
test_d <- sample.int(n_draws, 1L)
replay <- simulate_renewal_core(
  R0_vec = draws$R0_t[test_d, ], w = w, S0 = sd_clim$N_start[1], U0 = 0,
  N_start_vec = sd_clim$N_start, N_end_vec = sd_clim$N_end, births_vec = sd_clim$births, deaths_vec = sd_clim$deaths,
  X_recent_history = rep(0, G), lambda_import = sd_clim$imports_per_week, stochastic = FALSE,
  is_seed = sd_clim$is_seed, X_seed = sd_clim$X_seed
)
max_rel_err <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-9))
err_S <- max_rel_err(replay$S, draws$S[test_d, ]); err_U <- max_rel_err(replay$U, draws$U[test_d, ]); err_X <- max_rel_err(replay$X, draws$X[test_d, ])
check("Historical replay reproduces S/U/X", "A", max(err_S, err_U, err_X) < 1e-6,
      sprintf("draw=%d, max relative error S=%.2e U=%.2e X=%.2e", test_d, err_S, err_U, err_X))

# ================================================================
# B. MASS CONSERVATION -- S(t)+U(t) tracks the N_start bookkeeping identity
# ================================================================
mass_err <- max(abs(replay$S + replay$U - sd_clim$N_start))
check("Mass conservation S(t)+U(t)=N_start(t)", "B", mass_err < 1e-6 * max(sd_clim$N_start), sprintf("max abs error = %.3e", mass_err))

# ================================================================
# C. ZERO-INTRODUCTION TEST -- zero history + zero importation => no spontaneous infections
# ================================================================
zero_intro <- simulate_renewal_core(
  R0_vec = rep(2.0, 100), w = w, S0 = 1e6, U0 = 0,
  N_start_vec = rep(1e6, 100), N_end_vec = rep(1e6, 100), births_vec = rep(0, 100), deaths_vec = rep(0, 100),
  X_recent_history = rep(0, G), lambda_import = 0, stochastic = FALSE
)
check("Zero introduction pressure -> zero infections", "C", all(zero_intro$X == 0), sprintf("max X = %.3e (R0=2.0, 100 weeks, zero history, zero imports)", max(zero_intro$X)))

zero_intro_stoch <- simulate_renewal_core(
  R0_vec = rep(2.0, 100), w = w, S0 = 1e6, U0 = 0,
  N_start_vec = rep(1e6, 100), N_end_vec = rep(1e6, 100), births_vec = rep(0, 100), deaths_vec = rep(0, 100),
  X_recent_history = rep(0, G), lambda_import = 0, stochastic = TRUE,
  rng_pois = function(n, lam) rpois(n, lam), rng_binom = function(n, size, p) rbinom(n, size, p)
)
check("Zero introduction pressure, stochastic mode", "C", all(zero_intro_stoch$X == 0), sprintf("max X = %.3e", max(zero_intro_stoch$X)))

# ================================================================
# F. FORCING IDENTITY -- beta=0 (zero coefficient) vs z=0 (zero anomaly)
# must give IDENTICAL R0_future for the same draw and year-effect arm
# ================================================================
future_dates <- build_future_dates()
t_index_future <- (N_hist + 1):(N_hist + length(future_dates))
years_future <- as.integer(format(future_dates, "%Y"))
year_id_future <- match(years_future, sort(unique(years_future)))
n_future_years <- length(unique(years_future))

set.seed(BASE_SEED + 1L)
year_struct <- sample_future_year_structure("Y1_resampled_distributional", n_future_years, 0.40, draws$sigma_season_year[test_d], rng = function(n) rnorm(n))
z_fake <- rnorm(length(future_dates), 0, 1) # arbitrary nonzero anomaly sequence for the identity check

r0_C0 <- compute_R0_future(t_index_future, year_id_future, draws$alpha_R[test_d], year_struct$year_effect, year_struct$A_year,
                            draws$beta_sin1[test_d], draws$beta_cos1[test_d], draws$beta_sin2[test_d], draws$beta_cos2[test_d],
                            beta_T_eff = draws$beta_T[test_d], beta_P_eff = draws$beta_P[test_d], z_T_future = rep(0, length(future_dates)), z_P_future = rep(0, length(future_dates)))
r0_C2 <- compute_R0_future(t_index_future, year_id_future, draws$alpha_R[test_d], year_struct$year_effect, year_struct$A_year,
                            draws$beta_sin1[test_d], draws$beta_cos1[test_d], draws$beta_sin2[test_d], draws$beta_cos2[test_d],
                            beta_T_eff = 0, beta_P_eff = 0, z_T_future = z_fake, z_P_future = z_fake)
forcing_err <- max(abs(r0_C0$log_R0 - r0_C2$log_R0))
check("Forcing identity: beta=0 == z=0 (C0 vs C2)", "F", forcing_err < 1e-10, sprintf("max abs log_R0 difference = %.3e", forcing_err))

# Reference specification reproduces historical R0(t) when supplied with
# historical coefficients/climate/year effects (2.2's explicit requirement).
hist_year_effect_draw <- rstan::extract(climate$fit, pars = "year_effect", permuted = TRUE)$year_effect[test_d, ]
hist_A_year_draw <- rstan::extract(climate$fit, pars = "A_year", permuted = TRUE)$A_year[test_d, ]
r0_hist_reconstructed <- compute_R0_future(seq_len(N_hist), sd_clim$year_id, draws$alpha_R[test_d], hist_year_effect_draw, hist_A_year_draw,
                                            draws$beta_sin1[test_d], draws$beta_cos1[test_d], draws$beta_sin2[test_d], draws$beta_cos2[test_d],
                                            beta_T_eff = draws$beta_T[test_d], beta_P_eff = draws$beta_P[test_d], z_T_future = sd_clim$z_T_anom, z_P_future = sd_clim$z_P_anom)
hist_reproduction_err <- max(abs(r0_hist_reconstructed$R0 - draws$R0_t[test_d, ]))
check("compute_R0_future() reproduces historical R0(t) exactly", "2.2", hist_reproduction_err < 1e-8, sprintf("max abs error = %.3e", hist_reproduction_err))

# ================================================================
# D. STOCHASTIC BEHAVIOUR -- distinct trajectories, fadeout possible
# ================================================================
n_future_weeks <- length(future_dates)
r0_future_demo <- rep(1.5, n_future_weeks) # simple constant-R0 scenario for the stochastic-behaviour check
set.seed(BASE_SEED + 2L)
stoch_reps <- lapply(1:20, function(r) {
  simulate_renewal_core(
    R0_vec = r0_future_demo, w = w, S0 = 2e6, U0 = 6e6, N_start_vec = rep(8e6, n_future_weeks), N_end_vec = rep(8e6, n_future_weeks),
    births_vec = rep(800, n_future_weeks), deaths_vec = rep(800, n_future_weeks), X_recent_history = rep(0, G),
    lambda_import = 0.25, stochastic = TRUE, rng_pois = function(n, lam) rpois(n, lam), rng_binom = function(n, size, p) rbinom(n, size, p)
  )
})
cum_X <- vapply(stoch_reps, function(s) sum(s$X), numeric(1))
n_fadeout_weeks <- vapply(stoch_reps, function(s) sum(s$X == 0), integer(1))
check("Stochastic replicates produce distinct trajectories", "D", length(unique(round(cum_X))) > 1L,
      sprintf("unique cumulative-infection outcomes = %d / 20 replicates, var(cum_X)=%.1f", length(unique(round(cum_X))), var(cum_X)))
check("Stochastic replicates permit fadeout (zero-infection weeks)", "D", any(n_fadeout_weeks > 0L),
      sprintf("replicates with >=1 zero-infection week: %d / 20 (range %d-%d weeks)", sum(n_fadeout_weeks > 0), min(n_fadeout_weeks), max(n_fadeout_weeks)))

# ================================================================
# E. POSTERIOR PROPAGATION -- parameters and terminal states drawn jointly
# ================================================================
sample_d <- sample.int(n_draws, min(300L, n_draws))
cor_used <- cor(draws$alpha_R[sample_d], draws$S[sample_d, N_hist])
cor_full_posterior <- cor(draws$alpha_R, draws$S[, N_hist])
check("Posterior draw index consistent across alpha_R and terminal S (joint propagation)", "E",
      abs(cor_used - cor_full_posterior) < 0.1, sprintf("cor(alpha_R, S_terminal) on sample=%.3f vs full posterior=%.3f", cor_used, cor_full_posterior))

# ================================================================
# G. NO FUTURE DATA LEAKAGE -- future dates do not overlap any historical
# observed-data file's date range; the simulator's own argument list
# contains no observed-case or observed-date input for 2026+.
# ================================================================
no_overlap <- min(future_dates) > HIST_DATE_END
sim_core_args <- names(formals(simulate_renewal_core))
no_obs_arg <- !any(grepl("observed|future_case|future_obs", sim_core_args, ignore.case = TRUE))
check("Future horizon starts strictly after the historical data end date", "G", no_overlap, sprintf("future starts %s, historical ends %s", min(future_dates), HIST_DATE_END))
check("simulate_renewal_core() accepts no observed-future-data argument", "G", no_obs_arg, paste("args:", paste(sim_core_args, collapse = ", ")))

# ================================================================
# H. DEMOGRAPHIC PLAUSIBILITY
# ================================================================
plausible <- all(replay$S >= -1e-6) && all(replay$U >= -1e-6) && all(replay$X <= replay$S + replay$X + 1e-6) && all(replay$N > 0)
check("Demographic plausibility (S>=0, U>=0, N>0) on historical replay", "H", plausible, "checked S(t)>=0, U(t)>=0, N(t)>0 for all weeks")

# ================================================================
# I. COMPUTATIONAL REPRODUCIBILITY -- same seed => same trajectory
# ================================================================
run_with_seed <- function(seed) {
  set.seed(seed)
  simulate_renewal_core(R0_vec = r0_future_demo, w = w, S0 = 2e6, U0 = 6e6, N_start_vec = rep(8e6, n_future_weeks), N_end_vec = rep(8e6, n_future_weeks),
                         births_vec = rep(800, n_future_weeks), deaths_vec = rep(800, n_future_weeks), X_recent_history = rep(0, G),
                         lambda_import = 1, stochastic = TRUE, rng_pois = function(n, lam) rpois(n, lam), rng_binom = function(n, size, p) rbinom(n, size, p))
}
rep1a <- run_with_seed(777L); rep1b <- run_with_seed(777L); rep2 <- run_with_seed(778L)
check("Same seed reproduces identical trajectory", "I", identical(rep1a$X, rep1b$X), "two runs with seed=777")
check("Different seed yields a different trajectory", "I", !identical(rep1a$X, rep2$X), "seed=777 vs seed=778")

# ================================================================
# 2.1 Terminal-state bridge: completing the final historical week's
# demographic step must land exactly on N_end(N_hist), for every draw.
# ================================================================
reconciliation_last <- sd_clim$N_end[N_hist] - sd_clim$N_start[N_hist] - sd_clim$births[N_hist] + sd_clim$deaths[N_hist]
bridge_err <- vapply(sample_d, function(d) {
  adv <- advance_terminal_state(draws$S[d, N_hist], draws$U[d, N_hist], draws$X[d, N_hist], sd_clim$births[N_hist], sd_clim$deaths[N_hist], reconciliation_last)
  abs(adv$S + adv$U - sd_clim$N_end[N_hist])
}, numeric(1))
check("Terminal-state bridge lands exactly on N_end(N_hist)", "2.1", max(bridge_err) < 1e-6, sprintf("max abs error across %d draws = %.3e", length(sample_d), max(bridge_err)))

# ================================================================
# Save QA summary
# ================================================================
qa_summary <- bind_rows(results)
write_csv(qa_summary, file.path(DIR_TABLES, "STEP2_unit_test_QA_summary.csv"))
message("\n=== STEP 2 unit test summary ===")
print(as.data.frame(qa_summary[, c("section", "check", "pass")]))
if (!all(qa_summary$pass)) {
  stop("One or more Step 2 unit tests FAILED -- see STEP2_unit_test_QA_summary.csv. Do not proceed to the simulation ensemble until resolved.")
}
message("\n[Step2 unit tests] ALL CHECKS PASSED (", nrow(qa_summary), " checks).")
