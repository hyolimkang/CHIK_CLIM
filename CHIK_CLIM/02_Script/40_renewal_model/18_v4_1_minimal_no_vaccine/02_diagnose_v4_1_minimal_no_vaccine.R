# Read-only diagnostics and validation for Ceará renewal v4.1.
#
# v4.1 changes exactly two things relative to the HMC-passing v4.0: (A) the
# annual-effect parameterisation (explicit orthonormal sum-to-zero basis,
# learned sigma_year) and (B) the reporting fraction (estimated constant q
# with an informative Beta prior, replacing v4.0's q_fixed = 0.10). This
# script therefore adds, beyond the v4.0 six-panel diagnostics: prior-vs-
# posterior overlays for q and sigma_year; the six prespecified posterior
# correlations (q vs alpha_R, q vs sigma_year, q vs S_2025, q vs R0_2017,
# q vs R0_2022); an explicit annual PPC table for 2016/2017/2021/2022/2023
# with expected/observed ratios; and an eight-panel figure.

required_packages <- c("here", "rstan", "posterior", "dplyr", "tibble", "ggplot2", "patchwork", "scales", "readr", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(rstan); library(posterior); library(dplyr); library(tibble)
  library(ggplot2); library(patchwork); library(scales); library(readr); library(tidyr)
})

v4_1_diagnostic_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

v4_1_diagnostic_paths <- function(root = v4_1_diagnostic_root()) {
  run_tag <- Sys.getenv("RENEWAL_V4_1_RUN_TAG", unset = "base")
  if (!grepl("^[A-Za-z0-9_-]+$", run_tag)) stop("RENEWAL_V4_1_RUN_TAG may use only letters, numbers, _ and -.")
  suffix <- if (identical(run_tag, "base")) "" else paste0("_", run_tag)
  list(
    fit = file.path(root, "02_Script", "stan", paste0("renewal_ceara_v4_1_minimal_no_vaccine_fit", suffix, ".rds")),
    table = file.path(root, "03_Output", "tables", paste0("renewal_v4_1_minimal_no_vaccine", suffix)),
    figure = file.path(root, "03_Output", "figures", paste0("renewal_v4_1_minimal_no_vaccine", suffix)),
    v4_0_base_table = file.path(root, "03_Output", "tables", "renewal_v4_0_minimal_no_vaccine_td14"),
    run_tag = run_tag
  )
}

summarise_matrix <- function(x, dates, variable) {
  tibble(
    week_start = dates, variable = variable,
    q025 = apply(x, 2, quantile, probs = .025), q25 = apply(x, 2, quantile, probs = .25),
    median = apply(x, 2, median), q75 = apply(x, 2, quantile, probs = .75), q975 = apply(x, 2, quantile, probs = .975)
  )
}

add_band <- function(plot, summary, colour) {
  plot +
    geom_ribbon(data = summary, aes(ymin = q025, ymax = q975), fill = colour, alpha = .13) +
    geom_ribbon(data = summary, aes(ymin = q25, ymax = q75), fill = colour, alpha = .24) +
    geom_line(data = summary, aes(y = median), colour = colour, linewidth = .65)
}

run_v4_1_diagnostics <- function() {
  paths <- v4_1_diagnostic_paths()
  if (!file.exists(paths$fit)) stop("Fit not found: ", paths$fit)
  dir.create(paths$table, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$figure, recursive = TRUE, showWarnings = FALSE)
  bundle <- readRDS(paths$fit)
  if (!inherits(bundle$fit, "stanfit") || !identical(bundle$config$model_version, "v4_1_minimal_no_vaccine")) {
    stop("Expected a v4.1 minimal no-vaccination fit bundle.")
  }

  fit <- bundle$fit
  dates <- as.Date(bundle$weekly_data$week_start)
  years <- as.integer(format(dates, "%Y"))
  core_pars <- c("alpha_R", "z_year_free", "sigma_year", "beta_sin", "beta_cos", "phi_obs", "q")
  raw <- rstan::extract(fit, pars = core_pars, permuted = FALSE, inc_warmup = FALSE)
  sampler_draws <- posterior::as_draws_array(raw)
  parameter_summary <- as.data.frame(posterior::summarise_draws(sampler_draws, posterior::rhat, posterior::ess_bulk, posterior::ess_tail))
  names(parameter_summary) <- sub("^posterior::", "", names(parameter_summary))
  names(parameter_summary)[names(parameter_summary) == ".variable"] <- "parameter"
  bfmi <- rstan::get_bfmi(fit)
  hmc <- tibble(
    metric = c("divergences", "max_treedepth_hits", "maximum_Rhat", "n_Rhat_gt_1.01", "minimum_bulk_ESS", "minimum_tail_ESS", "elapsed_seconds", paste0("BFMI_chain_", seq_along(bfmi))),
    value = c(
      rstan::get_num_divergent(fit), rstan::get_num_max_treedepth(fit),
      max(parameter_summary$rhat, na.rm = TRUE), sum(parameter_summary$rhat > 1.01, na.rm = TRUE),
      min(parameter_summary$ess_bulk, na.rm = TRUE), min(parameter_summary$ess_tail, na.rm = TRUE),
      bundle$config$elapsed_seconds, bfmi
    )
  )
  hmc_pass <- hmc$value[hmc$metric == "divergences"] == 0 &&
    hmc$value[hmc$metric == "max_treedepth_hits"] == 0 &&
    hmc$value[hmc$metric == "maximum_Rhat"] <= 1.01 &&
    hmc$value[hmc$metric == "minimum_bulk_ESS"] >= 100 && all(bfmi >= .3)
  if (!identical(hmc_pass, bundle$hmc$hmc_pass)) warning("Stored and recomputed HMC gate disagree.")
  write_csv(hmc, file.path(paths$table, "hmc_diagnostics.csv"))
  write_csv(arrange(parameter_summary, desc(rhat)), file.path(paths$table, "parameter_diagnostics.csv"))

  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  treedepth_all <- unlist(lapply(sampler_params, function(x) x[, "treedepth__"]))
  chain_diag <- tibble(
    chain = seq_along(sampler_params),
    treedepth_hits = vapply(sampler_params, function(x) sum(x[, "treedepth__"] >= bundle$config$max_treedepth), integer(1)),
    mean_treedepth = vapply(sampler_params, function(x) mean(x[, "treedepth__"]), numeric(1)),
    mean_accept_stat = vapply(sampler_params, function(x) mean(x[, "accept_stat__"]), numeric(1)),
    n_divergent = vapply(sampler_params, function(x) sum(x[, "divergent__"]), numeric(1)),
    bfmi = bfmi
  )
  write_csv(chain_diag, file.path(paths$table, "chain_specific_diagnostics.csv"))
  write_csv(
    tibble(quantile = c("50%", "90%", "95%", "99%", "max"), treedepth = c(quantile(treedepth_all, c(.5, .9, .95, .99)), max(treedepth_all))),
    file.path(paths$table, "treedepth_quantiles.csv")
  )

  # --- q and sigma_year: point summaries, prior-vs-posterior --------------
  q_draws <- rstan::extract(fit, pars = "q")$q
  sigma_year_draws <- rstan::extract(fit, pars = "sigma_year")$sigma_year
  q_prior_a <- bundle$stan_data$q_prior_a; q_prior_b <- bundle$stan_data$q_prior_b
  q_prior_mean <- q_prior_a / (q_prior_a + q_prior_b)
  q_summary <- tibble(
    parameter = "q", prior_mean = q_prior_mean, prior_a = q_prior_a, prior_b = q_prior_b,
    posterior_mean = mean(q_draws), posterior_median = median(q_draws), posterior_sd = sd(q_draws),
    posterior_q025 = quantile(q_draws, .025), posterior_q975 = quantile(q_draws, .975),
    posterior_minus_prior_mean_in_prior_sd = (mean(q_draws) - q_prior_mean) / sqrt(q_prior_a * q_prior_b / ((q_prior_a + q_prior_b)^2 * (q_prior_a + q_prior_b + 1)))
  )
  sigma_year_summary <- tibble(
    parameter = "sigma_year", prior_mean = 0.35 * sqrt(2 / pi), prior_sd_param = 0.35,
    posterior_mean = mean(sigma_year_draws), posterior_median = median(sigma_year_draws), posterior_sd = sd(sigma_year_draws),
    posterior_q025 = quantile(sigma_year_draws, .025), posterior_q975 = quantile(sigma_year_draws, .975)
  )
  write_csv(bind_rows(q_summary, sigma_year_summary), file.path(paths$table, "q_and_sigma_year_prior_vs_posterior.csv"))
  message(sprintf("[v4.1 diagnostics] q: prior mean=%.4f -> posterior mean=%.4f [%.4f, %.4f] (%.2f prior SD shift)",
                   q_prior_mean, mean(q_draws), quantile(q_draws, .025), quantile(q_draws, .975), q_summary$posterior_minus_prior_mean_in_prior_sd))

  # --- Derived quantities needed for correlations and PPC -----------------
  alpha_R_draws <- rstan::extract(fit, pars = "alpha_R")$alpha_R
  draws <- rstan::extract(fit, pars = c("C_pred", "X", "S_prop", "immune_prop", "R0_t", "R_eff_t", "expected_reported_cases", "year_effect", "S"))
  sorted_years <- sort(unique(years))
  year_first_index <- vapply(sorted_years, function(y) which(years == y)[1], integer(1))
  R0_by_year_draws <- sapply(seq_along(sorted_years), function(i) exp(alpha_R_draws + draws$year_effect[, i]))
  colnames(R0_by_year_draws) <- sorted_years
  S_2025_draws <- draws$S[, which(years == 2025L)[1]]
  N_start_2025 <- bundle$stan_data$N_start[which(years == 2025L)[1]]

  correlations <- tibble(
    pair = c("q_vs_alpha_R", "q_vs_sigma_year", "q_vs_S_2025", "q_vs_R0_2017", "q_vs_R0_2022"),
    pearson_r = c(
      cor(q_draws, alpha_R_draws), cor(q_draws, sigma_year_draws), cor(q_draws, S_2025_draws / N_start_2025),
      cor(q_draws, R0_by_year_draws[, "2017"]), cor(q_draws, R0_by_year_draws[, "2022"])
    )
  ) |> mutate(pathological_confounding_flag = abs(pearson_r) > 0.9)
  write_csv(correlations, file.path(paths$table, "q_posterior_correlations.csv"))
  message("[v4.1 diagnostics] q posterior correlations:")
  for (i in seq_len(nrow(correlations))) message(sprintf("  %s: r=%.3f%s", correlations$pair[i], correlations$pearson_r[i], if (correlations$pathological_confounding_flag[i]) " *** |r|>0.9 ***" else ""))

  # --- Trajectory summary, PPC, annual decomposition -----------------------
  summary <- bind_rows(
    summarise_matrix(draws$C_pred, dates, "posterior_predictive_cases"), summarise_matrix(draws$X, dates, "latent_infections"),
    summarise_matrix(draws$S_prop, dates, "susceptible_proportion"), summarise_matrix(draws$immune_prop, dates, "infection_derived_immune_proportion"),
    summarise_matrix(draws$R0_t, dates, "R0"), summarise_matrix(draws$R_eff_t, dates, "Reff"), summarise_matrix(draws$expected_reported_cases, dates, "expected_reported_cases")
  )
  write_csv(summary, file.path(paths$table, "trajectory_summary.csv"))
  get_summary <- function(variable_name) dplyr::filter(summary, .data$variable == variable_name)
  case_summary <- get_summary("posterior_predictive_cases")
  predictive_coverage <- mean(bundle$weekly_data$cases >= case_summary$q025 & bundle$weekly_data$cases <= case_summary$q975)
  write_csv(
    tibble(metric = c("weekly_case_95pct_predictive_coverage", "weekly_case_median_RMSE", "weekly_case_median_MAE"),
           value = c(predictive_coverage, sqrt(mean((bundle$weekly_data$cases - case_summary$median)^2)), mean(abs(bundle$weekly_data$cases - case_summary$median)))),
    file.path(paths$table, "posterior_predictive_checks.csv")
  )

  annual <- bind_rows(lapply(seq_along(sorted_years), function(i) {
    year <- sorted_years[i]; index <- which(years == year)
    infections <- rowSums(draws$X[, index, drop = FALSE]); expected_cases <- rowSums(draws$expected_reported_cases[, index, drop = FALSE])
    observed <- sum(bundle$weekly_data$cases[index])
    tibble(
      year = year, observed_cases = observed,
      expected_cases_q025 = quantile(expected_cases, .025), expected_cases_median = median(expected_cases), expected_cases_q975 = quantile(expected_cases, .975),
      infections_q025 = quantile(infections, .025), infections_median = median(infections), infections_q975 = quantile(infections, .975),
      S_start_median = median(draws$S_prop[, index[1]]), S_end_median = median(draws$S_prop[, tail(index, 1)]),
      annual_attack_median = median(infections / bundle$stan_data$N_start[index[1]]),
      relative_residual = (observed - median(expected_cases)) / pmax(median(expected_cases), 1),
      annual_R0_median = median(R0_by_year_draws[, i]), annual_R0_q025 = quantile(R0_by_year_draws[, i], .025), annual_R0_q975 = quantile(R0_by_year_draws[, i], .975),
      ratio_expected_over_observed = median(expected_cases) / observed
    )
  }))
  write_csv(annual, file.path(paths$table, "annual_attack_decomposition.csv"))

  highlight_years <- c(2016L, 2017L, 2021L, 2022L, 2023L)
  annual_ppc_highlight <- filter(annual, year %in% highlight_years) |>
    select(year, observed_cases, expected_cases_median, expected_cases_q025, expected_cases_q975, ratio_expected_over_observed)
  write_csv(annual_ppc_highlight, file.path(paths$table, "annual_ppc_highlight_years.csv"))
  message("[v4.1 diagnostics] annual PPC (2016/2017/2021/2022/2023):")
  for (i in seq_len(nrow(annual_ppc_highlight))) {
    r <- annual_ppc_highlight[i, ]
    message(sprintf("  %d: observed=%s expected_median=%.0f [%.0f, %.0f] ratio_exp/obs=%.3f",
                     r$year, format(r$observed_cases, big.mark = ","), r$expected_cases_median, r$expected_cases_q025, r$expected_cases_q975, r$ratio_expected_over_observed))
  }

  v4_0_annual_path <- file.path(paths$v4_0_base_table, "annual_attack_decomposition.csv")
  if (file.exists(v4_0_annual_path)) {
    v4_0_annual <- read_csv(v4_0_annual_path, show_col_types = FALSE) |>
      transmute(year, v4_0_expected_median = expected_cases_median, v4_0_relative_residual = (observed_cases - expected_cases_median) / pmax(expected_cases_median, 1))
    mismatch_comparison <- annual |>
      transmute(year, observed_cases, v4_1_expected_median = expected_cases_median, v4_1_relative_residual = relative_residual) |>
      left_join(v4_0_annual, by = "year") |>
      mutate(abs_relative_residual_improved = abs(v4_1_relative_residual) < abs(v4_0_relative_residual))
    write_csv(mismatch_comparison, file.path(paths$table, "annual_mismatch_comparison_v4_0_vs_v4_1.csv"))
    for (yr in c(2017L, 2022L)) {
      row <- filter(mismatch_comparison, year == yr)
      if (nrow(row)) message(sprintf("[v4.1 diagnostics] %d relative residual: v4.0(td14)=%.3f -> v4.1=%.3f (%s)", yr, row$v4_0_relative_residual, row$v4_1_relative_residual, if (row$abs_relative_residual_improved) "IMPROVED" else "NOT improved"))
    }
  } else {
    message("[v4.1 diagnostics] v4.0 td14 annual_attack_decomposition.csv not found -- skipping direct mismatch comparison.")
  }

  serology_windows <- tibble(
    site = c("Juazeiro do Norte", "Quixada"), observed_positive = c(103L, 289L), observed_n = c(404L, 409L),
    start = as.Date(c("2018-06-03", "2018-06-03")), end = as.Date(c("2018-12-30", "2019-12-29"))
  ) |> rowwise() |> mutate(
    state_immune_q025 = quantile(rowMeans(draws$immune_prop[, dates >= start & dates <= end, drop = FALSE]), .025),
    state_immune_median = median(rowMeans(draws$immune_prop[, dates >= start & dates <= end, drop = FALSE])),
    state_immune_q975 = quantile(rowMeans(draws$immune_prop[, dates >= start & dates <= end, drop = FALSE]), .975),
    observed_proportion = observed_positive / observed_n
  ) |> ungroup()
  write_csv(serology_windows, file.path(paths$table, "external_serology_consistency_check.csv"))

  # --- Figures --------------------------------------------------------------
  theme_v4 <- theme_classic(base_size = 8.5) + theme(legend.position = "top", legend.title = element_blank(), panel.grid.major.y = element_line(colour = "grey90"))
  xscale <- scale_x_date(breaks = seq(as.Date("2015-01-01"), as.Date("2026-01-01"), by = "2 years"), date_labels = "%Y")
  p1 <- add_band(ggplot(case_summary, aes(week_start)), case_summary, "#56B4E9") +
    geom_point(data = tibble(week_start = dates, observed_cases = bundle$weekly_data$cases), aes(week_start, observed_cases), size = .3) +
    xscale + labs(title = "A. Reported cases and posterior prediction", x = NULL, y = "Weekly cases") + theme_v4
  p2_data <- get_summary("latent_infections")
  p2 <- add_band(ggplot(p2_data, aes(week_start)), p2_data, "#D55E00") + xscale + labs(title = "B. Latent infection incidence", x = NULL, y = "Infections") + theme_v4
  p3_data <- bind_rows(mutate(get_summary("susceptible_proportion"), stock = "Susceptible"), mutate(get_summary("infection_derived_immune_proportion"), stock = "Infection-derived immune"))
  p3 <- ggplot(p3_data, aes(week_start, median, colour = stock, fill = stock)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .12, colour = NA) + geom_line() +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) + xscale + labs(title = "C. Susceptible and immune stocks", x = NULL, y = "Population proportion") + theme_v4
  p4_data <- bind_rows(mutate(get_summary("R0"), reproduction_number = "R0"), mutate(get_summary("Reff"), reproduction_number = "Reff"))
  p4 <- ggplot(p4_data, aes(week_start, median, colour = reproduction_number, fill = reproduction_number)) +
    geom_hline(yintercept = 1, linetype = 2) + geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .12, colour = NA) + geom_line() +
    xscale + labs(title = "D. Reproduction numbers", x = NULL, y = "R") + theme_v4
  annual_long <- bind_rows(
    transmute(annual, year, type = "Observed reported cases", estimate = observed_cases, lower = observed_cases, upper = observed_cases),
    transmute(annual, year, type = "Posterior expected reported cases", estimate = expected_cases_median, lower = expected_cases_q025, upper = expected_cases_q975)
  )
  p5 <- ggplot(annual_long, aes(year, estimate, colour = type)) + geom_line() + geom_point() +
    geom_errorbar(data = filter(annual_long, type != "Observed reported cases"), aes(ymin = lower, ymax = upper), width = .15) +
    labs(title = "E. Annual reported-case validation", x = NULL, y = "Cases") + theme_v4
  p6 <- ggplot(serology_windows, aes(site, state_immune_median)) +
    geom_errorbar(aes(ymin = state_immune_q025, ymax = state_immune_q975), width = .15, colour = "#0072B2") +
    geom_point(colour = "#0072B2") + geom_point(aes(y = observed_proportion), colour = "#D55E00", shape = 18, size = 2.5) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    labs(title = "F. External serology consistency", subtitle = "Blue: model; orange: survey; not fitted", x = NULL, y = "Immune proportion") + theme_v4
  annual_plot_data <- mutate(annual, highlight = year %in% c(2017, 2022), label = as.character(year))
  p7 <- ggplot(annual_plot_data, aes(year, relative_residual, fill = highlight)) +
    geom_col() + geom_hline(yintercept = 0, linetype = 2) +
    geom_text(data = filter(annual_plot_data, highlight), aes(label = label), vjust = -.4, size = 2.6) +
    scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey60"), guide = "none") +
    scale_y_continuous(labels = scales::label_percent()) +
    labs(title = "G. Annual mismatch by calendar year", subtitle = "(observed - expected median) / expected median", x = NULL, y = "Relative residual") + theme_v4

  q_prior_curve <- tibble(x = seq(0, 0.4, length.out = 400)) |> mutate(density = dbeta(x, q_prior_a, q_prior_b))
  p8a <- ggplot() +
    geom_histogram(data = tibble(q = q_draws), aes(q, after_stat(density)), bins = 60, fill = "#56B4E9", alpha = .6) +
    geom_line(data = q_prior_curve, aes(x, density), colour = "black", linewidth = .6, linetype = 2) +
    labs(title = "H1. q: prior (dashed) vs posterior", x = "q", y = "Density") + theme_v4
  sigma_prior_curve <- tibble(x = seq(0, 1.2, length.out = 400)) |> mutate(density = 2 * dnorm(x, 0, 0.35))
  p8b <- ggplot() +
    geom_histogram(data = tibble(sigma_year = sigma_year_draws), aes(sigma_year, after_stat(density)), bins = 60, fill = "#009E73", alpha = .6) +
    geom_line(data = sigma_prior_curve, aes(x, density), colour = "black", linewidth = .6, linetype = 2) +
    labs(title = "H2. sigma_year: prior (dashed) vs posterior", x = "sigma_year", y = "Density") + theme_v4

  subtitle <- if (hmc_pass) {
    sprintf("HMC gate passed; q posterior mean=%.3f [%.3f,%.3f] (prior mean %.3f); sigma_year median=%.3f",
            mean(q_draws), quantile(q_draws, .025), quantile(q_draws, .975), q_prior_mean, median(sigma_year_draws))
  } else "UNCONVERGED - do not interpret posterior trajectories"
  figure <- (p1 | p2) / (p3 | p4) / (p5 | p6) / (p7 | plot_spacer()) / (p8a | p8b) + plot_annotation(
    title = "Ceara v4.1 minimal all-age no-vaccination renewal model (estimated q + learned sigma_year)",
    subtitle = subtitle, tag_levels = "A",
    theme = theme(plot.subtitle = element_text(face = "bold", colour = if (hmc_pass) "grey30" else "#D55E00"))
  )
  figure_path <- file.path(paths$figure, paste0("renewal_v4_1_minimal_no_vaccine_panel", if (hmc_pass) "" else "_UNCONVERGED", ".pdf"))
  ggsave(figure_path, figure, width = 183, height = 340, units = "mm", device = cairo_pdf, limitsize = FALSE)
  message("[v4.1 diagnostics] HMC gate: ", if (hmc_pass) "PASS" else "FAIL", "; figure: ", figure_path)
  invisible(list(hmc = hmc, annual = annual, correlations = correlations, q_summary = q_summary, sigma_year_summary = sigma_year_summary))
}

if (sys.nframe() == 0L) run_v4_1_diagnostics()
