# Mode audit Part 3: weekly posterior-predictive residual diagnostics for
# both modes, using the already-fit clean mode-targeted fits
# (modeB_lowclean, modeB_high). Uses randomized (probability-integral-
# transform) residuals for the negative-binomial observation model, which
# are approximately standard normal under correct specification and are the
# appropriate residual definition for discrete count models (Dunn & Smyth
# 1996) -- ordinary Pearson/standardized residuals are not well-behaved for
# NB counts with small means.

required_packages <- c("here", "rstan", "dplyr", "tibble", "readr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(readr); library(ggplot2); library(patchwork) })

acf_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

# Randomized quantile residuals for NB2(mu, phi): for each posterior draw,
# residual_t = qnorm(u_t) where u_t ~ Uniform(F(y_t - 1), F(y_t)) (F = NB
# CDF), then AVERAGE the resulting residual time series across draws (a
# standard posterior-predictive-residual summary).
randomized_residuals <- function(observed, mu_draws, phi_draws, seed = 20260912L) {
  set.seed(seed)
  n_draws <- length(phi_draws)
  n_t <- length(observed)
  resid <- matrix(NA_real_, n_draws, n_t)
  for (d in seq_len(n_draws)) {
    mu <- mu_draws[d, ]; phi <- phi_draws[d]
    lower <- ifelse(observed == 0, 0, pnbinom(observed - 1, size = phi, mu = mu))
    upper <- pnbinom(observed, size = phi, mu = mu)
    u <- runif(n_t, lower, upper)
    u <- pmin(pmax(u, 1e-9), 1 - 1e-9)
    resid[d, ] <- qnorm(u)
  }
  colMeans(resid)
}

run_weekly_residual_acf <- function() {
  root <- acf_root()
  base_dir <- file.path(root, "03_Output", "model_fits", "ceara", "v4_3_short_2015_2019_q_calibration")

  analyse_mode <- function(tag, label) {
    bundle <- readRDS(file.path(base_dir, "outputs", tag, paste0("renewal_ceara_v4_3_fit_", tag, ".rds")))
    fit <- bundle$fit
    mu <- rstan::extract(fit, pars = "expected_reported_cases")$expected_reported_cases
    phi <- rstan::extract(fit, pars = "phi_obs")$phi_obs
    observed <- bundle$weekly_data$cases
    dates <- as.Date(bundle$weekly_data$week_start)
    resid <- randomized_residuals(observed, mu, phi)
    acf_obj <- acf(resid, lag.max = 12, plot = FALSE)
    list(
      label = label, dates = dates, observed = observed, expected_median = apply(mu, 2, median), residual = resid,
      acf_table = tibble(mode = label, lag = as.integer(acf_obj$lag), acf = as.numeric(acf_obj$acf))
    )
  }

  low <- analyse_mode("modeB_lowclean", "low_q")
  high <- analyse_mode("modeB_high", "high_q")

  acf_table <- bind_rows(low$acf_table, high$acf_table) |> filter(lag >= 1, lag <= 12)
  write_csv(acf_table, file.path(base_dir, "weekly_residual_acf.csv"))
  message("[residual ACF] lag 1/2/4/8 by mode:")
  print(as.data.frame(acf_table |> filter(lag %in% c(1, 2, 4, 8))))

  theme_v4 <- theme_classic(base_size = 8.5) + theme(panel.grid.major.y = element_line(colour = "grey90"))
  p_acf <- ggplot(acf_table, aes(lag, acf, fill = mode)) +
    geom_col(position = "dodge") + geom_hline(yintercept = c(-1, 1) * 1.96 / sqrt(length(low$residual)), linetype = 2, colour = "grey40") +
    labs(title = "A. Weekly residual ACF (randomized quantile residuals)", x = "lag (weeks)", y = "ACF") + theme_v4

  resid_df <- bind_rows(
    tibble(mode = "low_q", week_start = low$dates, residual = low$residual, observed = low$observed, expected = low$expected_median),
    tibble(mode = "high_q", week_start = high$dates, residual = high$residual, observed = high$observed, expected = high$expected_median)
  )
  write_csv(resid_df, file.path(base_dir, "weekly_residual_timeseries.csv"))
  p_resid <- ggplot(resid_df, aes(week_start, residual, colour = mode)) + geom_line(alpha = .7) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = "B. Weekly residual time series", x = NULL, y = "randomized quantile residual") + theme_v4
  p_obs <- ggplot(resid_df, aes(week_start)) +
    geom_point(aes(y = observed), size = .3, colour = "grey30") +
    geom_line(aes(y = expected, colour = mode)) +
    labs(title = "C. Observed vs posterior-expected weekly cases", x = NULL, y = "cases") + theme_v4

  figure <- p_acf / p_resid / p_obs
  figure_dir <- file.path(root, "03_Output", "figures", "renewal_v4_3_short_2015_2019_q_calibration")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(figure_dir, "weekly_residual_acf.png"), figure, width = 200, height = 220, units = "mm", dpi = 300)
  message("[residual ACF] figure saved: ", file.path(figure_dir, "weekly_residual_acf.png"))

  # Flag epidemic-period residual patterns explicitly (Part 3 requirement).
  years <- as.integer(format(resid_df$week_start, "%Y"))
  period_summary <- resid_df |> mutate(year = as.integer(format(week_start, "%Y")),
    period = case_when(year == 2016 ~ "2016 epidemic", year == 2017 ~ "2017 epidemic",
                        year %in% c(2018, 2019) ~ "2018-2019", TRUE ~ "other/trough")) |>
    group_by(mode, period) |> summarise(mean_residual = mean(residual), sd_residual = sd(residual), .groups = "drop")
  write_csv(period_summary, file.path(base_dir, "weekly_residual_period_summary.csv"))
  message("[residual ACF] residual patterns by epidemic period:")
  print(as.data.frame(period_summary))

  invisible(list(acf_table = acf_table, resid_df = resid_df, period_summary = period_summary))
}

if (sys.nframe() == 0L) run_weekly_residual_acf()
