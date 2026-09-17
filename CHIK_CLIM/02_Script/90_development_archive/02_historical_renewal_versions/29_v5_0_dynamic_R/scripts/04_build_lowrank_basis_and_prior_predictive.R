# v5.0b low-rank cyclic dynamic-R basis: construction + prior-predictive
# check, BOTH done before any model fitting per the pre-registration
# requirement ("Choose the prior BEFORE examining the new model fit").
#
# K=8 equally-spaced knots over a cyclic day-of-year coordinate (0-52,
# smoothly wrapping) -> periodic piecewise-linear ("hat") basis, each row
# summing to 1. The basis is then PRE-CENTRED WITHIN EACH CALENDAR YEAR in
# R (each column's within-year mean subtracted), so that for ANY choice of
# year-specific coefficients theta[y,], the resulting weekly deviation has
# exactly zero mean within that year by construction -- this is what
# prevents the dynamic component from duplicating year_effect's role.

build_dynamic_basis <- function(dates, year_id, K = 8L) {
  n <- length(dates)
  day_of_year <- as.numeric(format(dates, "%j")) - 1
  days_in_year <- ifelse(as.integer(format(dates, "%Y")) %% 4 == 0 &
                            (as.integer(format(dates, "%Y")) %% 100 != 0 | as.integer(format(dates, "%Y")) %% 400 == 0),
                          366, 365)
  cyclic_pos <- day_of_year / days_in_year * 52 # continuous position on a 0-52 cyclic scale

  knot_pos <- (seq_len(K) - 1) / K * 52 # K equally spaced knots, e.g. 0, 6.5, 13, ..., 45.5
  knot_spacing <- 52 / K

  B <- matrix(0, n, K)
  for (i in seq_len(n)) {
    p <- cyclic_pos[i]
    # index of the knot at or before p (cyclic)
    k0 <- floor(p / knot_spacing) %% K + 1L
    k1 <- (k0 %% K) + 1L
    frac <- (p - knot_pos[k0]) / knot_spacing
    frac <- (frac %% 1 + 1) %% 1 # guard against tiny negative fractions from floating point
    B[i, k0] <- 1 - frac
    B[i, k1] <- B[i, k1] + frac
  }
  stopifnot(all(abs(rowSums(B) - 1) < 1e-8))

  # Pre-centre within each calendar year: subtract each column's within-year
  # mean from that year's rows, so any theta[y,] combination is automatically
  # zero-mean within year y.
  B_centred <- B
  for (y in sort(unique(year_id))) {
    idx <- which(year_id == y)
    col_means <- colMeans(B[idx, , drop = FALSE])
    B_centred[idx, ] <- sweep(B[idx, , drop = FALSE], 2, col_means, "-")
  }
  list(B_raw = B, B_centred = B_centred, knot_pos = knot_pos, cyclic_pos = cyclic_pos)
}

run_prior_predictive_check <- function(K = 8L, Y = 11L, tau_prior_sd = 0.15, n_draws = 20000L, seed = 20260912L) {
  set.seed(seed)
  tau_dynamic <- abs(rnorm(n_draws, 0, tau_prior_sd)) # half-normal
  # Typical basis entry at a knot after within-year centring is ~ (1 - mean_weight);
  # with ~52 weeks/year and K=8 knots, mean_weight per column ~ 52/8/52 = 1/8,
  # so peak centred-basis value is ~ 1 - 1/8 = 0.875 (matches build_dynamic_basis()).
  peak_basis_value <- 0.875
  z <- rnorm(n_draws)
  delta_R_at_knot <- tau_dynamic * z * peak_basis_value
  multiplicative_deviation <- exp(delta_R_at_knot)
  list(
    tau_dynamic_summary = quantile(tau_dynamic, c(.5, .8, .95, .99)),
    delta_R_at_knot_summary = quantile(delta_R_at_knot, c(.025, .25, .5, .75, .975)),
    multiplicative_deviation_summary = quantile(multiplicative_deviation, c(.025, .25, .5, .75, .975)),
    prob_deviation_exceeds_50pct = mean(abs(delta_R_at_knot) > log(1.5)),
    prob_deviation_exceeds_2x = mean(multiplicative_deviation > 2 | multiplicative_deviation < 0.5)
  )
}

if (sys.nframe() == 0L) {
  root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
  b9 <- readRDS(file.path(root, "03_Output/model_fits/ceara/v4_9_hierarchical_seasonality/outputs/q0.05/renewal_ceara_v4_9_fit_q0.05.rds"))
  dates <- as.Date(b9$weekly_data$week_start)
  year_id <- b9$stan_data$year_id
  basis <- build_dynamic_basis(dates, year_id, K = 8L)
  cat("Basis dims:", dim(basis$B_centred), "\n")
  cat("Row sums of RAW basis (should all be 1):", range(rowSums(basis$B_raw)), "\n")
  cat("Column means within year 1 (should be ~0 after centring):\n")
  print(round(colMeans(basis$B_centred[year_id == 1, ]), 6))
  cat("Range of centred basis entries:", range(basis$B_centred), "\n\n")

  cat("=== PRE-REGISTERED PRIOR-PREDICTIVE CHECK (tau_dynamic ~ half-normal(0, 0.15)) ===\n")
  pp <- run_prior_predictive_check(tau_prior_sd = 0.15)
  cat("tau_dynamic quantiles (50/80/95/99%):\n"); print(pp$tau_dynamic_summary)
  cat("\ndelta_R at a knot, quantiles:\n"); print(pp$delta_R_at_knot_summary)
  cat("\nimplied multiplicative R0 deviation at a knot, quantiles:\n"); print(pp$multiplicative_deviation_summary)
  cat(sprintf("\nP(|deviation| > 50%%) = %.3f\n", pp$prob_deviation_exceeds_50pct))
  cat(sprintf("P(deviation doubles or halves R0) = %.4f\n", pp$prob_deviation_exceeds_2x))
}
