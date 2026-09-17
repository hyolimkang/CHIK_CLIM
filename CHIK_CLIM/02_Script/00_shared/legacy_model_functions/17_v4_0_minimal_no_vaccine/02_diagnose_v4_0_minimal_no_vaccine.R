# Read-only diagnostics and six-panel validation for Ceará renewal v4.0.

required_packages <- c("here", "rstan", "posterior", "dplyr", "tibble", "ggplot2", "patchwork", "scales", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(rstan)
  library(posterior)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(readr)
})

v4_diagnostic_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

v4_diagnostic_paths <- function(root = v4_diagnostic_root()) {
  run_tag <- Sys.getenv("RENEWAL_V4_RUN_TAG", unset = "base")
  if (!grepl("^[A-Za-z0-9_-]+$", run_tag)) stop("RENEWAL_V4_RUN_TAG may use only letters, numbers, _ and -.")
  suffix <- if (identical(run_tag, "base")) "" else paste0("_", run_tag)
  list(
    fit = file.path(root, "03_Output", "02_ceara_pipeline", "model_fits", "initialization", paste0("renewal_ceara_v4_0_minimal_no_vaccine_fit", suffix, ".rds")),
    table = file.path(root, "03_Output", "02_ceara_pipeline", "tables", paste0("renewal_v4_0_minimal_no_vaccine", suffix)),
    figure = file.path(root, "03_Output", "02_ceara_pipeline", "figures", "transmission", paste0("renewal_v4_0_minimal_no_vaccine", suffix)),
    run_tag = run_tag
  )
}

summarise_matrix <- function(x, dates, variable) {
  tibble(
    week_start = dates, variable = variable,
    q025 = apply(x, 2, quantile, probs = .025),
    q25 = apply(x, 2, quantile, probs = .25),
    median = apply(x, 2, median),
    q75 = apply(x, 2, quantile, probs = .75),
    q975 = apply(x, 2, quantile, probs = .975)
  )
}

add_band <- function(plot, summary, colour) {
  plot +
    geom_ribbon(data = summary, aes(ymin = q025, ymax = q975), fill = colour, alpha = .13) +
    geom_ribbon(data = summary, aes(ymin = q25, ymax = q75), fill = colour, alpha = .24) +
    geom_line(data = summary, aes(y = median), colour = colour, linewidth = .65)
}

make_run_record <- function(path, bundle, hmc, predictive_coverage) {
  lines <- c(
    "# Ceará v4.0 minimal no-vaccination renewal run",
    "",
    "## Model status",
    "",
    paste0("- HMC gate: **", if (isTRUE(bundle$hmc$hmc_pass)) "PASS" else "FAIL", "**"),
    paste0("- Divergences: ", hmc$value[hmc$metric == "divergences"]),
    paste0("- Maximum treedepth hits: ", hmc$value[hmc$metric == "max_treedepth_hits"]),
    paste0("- Maximum Rhat: ", signif(hmc$value[hmc$metric == "maximum_Rhat"], 5)),
    paste0("- Minimum bulk ESS: ", signif(hmc$value[hmc$metric == "minimum_bulk_ESS"], 5)),
    paste0("- Weekly observed-case 95% posterior-predictive coverage: ", scales::percent(predictive_coverage, accuracy = .1)),
    "",
    "## Fixed scenario inputs",
    "",
    paste0("- q_fixed: ", bundle$config$q_fixed, " (not estimated)"),
    paste0("- imports_per_week: ", bundle$config$imports_per_week, " (not estimated)"),
    paste0("- annual-effect prior SD: ", bundle$config$year_effect_prior_sd),
    "",
    "## Interpretation guardrail",
    "",
    "This benchmark is eligible for scientific interpretation only when the HMC gate passes. It is an all-age, no-vaccination model; it does not identify age-specific transmission, reporting, or local serology offsets."
  )
  writeLines(lines, path)
}

make_convergence_assessment <- function(path, fit, hmc) {
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  tree_depth <- unlist(lapply(sampler, function(x) x[, "treedepth__"]))
  chain_hits <- vapply(sampler, function(x) sum(x[, "treedepth__"] >= 12), numeric(1))
  lines <- c(
    "# Convergence assessment", "",
    "## What passed", "",
    paste0("- Divergences: ", hmc$value[hmc$metric == "divergences"], "."),
    paste0("- Maximum R-hat: ", signif(hmc$value[hmc$metric == "maximum_Rhat"], 5), "."),
    paste0("- Minimum bulk ESS: ", signif(hmc$value[hmc$metric == "minimum_bulk_ESS"], 6), "."),
    "- BFMI was above the prespecified 0.30 threshold in every chain.", "",
    "## Why this run is not accepted", "",
    paste0("- Maximum-treedepth hits: ", hmc$value[hmc$metric == "max_treedepth_hits"],
      " (", scales::percent(mean(tree_depth >= 12), accuracy = .01), " of retained transitions)."),
    paste0("- Hits by chain: ", paste(chain_hits, collapse = ", "), "."),
    paste0("- Tree-depth quantiles: median ", quantile(tree_depth, .5),
      ", 95th percentile ", quantile(tree_depth, .95),
      ", 99th percentile ", quantile(tree_depth, .99), "."),
    "", "## Interpretation", "",
    "This is not a catastrophic failure: there are no divergences, R-hat is below 1.01, and effective sample sizes are adequate. However, the sampler repeatedly reached its numerical trajectory cap, concentrated in one chain. Posterior intervals and the six-panel figure may be used for exploratory model development only; they are not accepted as the baseline for vaccine-impact inference.",
    "", "## Prespecified next steps", "",
    "1. Run the identical model with a larger maximum treedepth in a separate output folder. If all other diagnostics remain satisfactory and no trajectories hit the new cap, retain that numerical run.",
    "2. If cap hits remain, retain the epidemiological assumptions but reparameterize the annual effects as an explicit sum-to-zero basis, then reassess HMC. This removes a redundant annual-effect direction without adding biological complexity.",
    "3. If reparameterization still fails, reduce annual transmission flexibility before adding climate, reporting, serology, age, or vaccination components. Do not solve this by loosening convergence thresholds."
  )
  writeLines(lines, path)
}

run_v4_0_diagnostics <- function() {
  paths <- v4_diagnostic_paths()
  if (!file.exists(paths$fit)) stop("Fit not found: ", paths$fit)
  dir.create(paths$table, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$figure, recursive = TRUE, showWarnings = FALSE)
  bundle <- readRDS(paths$fit)
  if (!inherits(bundle$fit, "stanfit") || !identical(bundle$config$model_version, "v4_0_minimal_no_vaccine")) {
    stop("Expected a v4.0 minimal no-vaccination fit bundle.")
  }

  fit <- bundle$fit
  dates <- as.Date(bundle$weekly_data$week_start)
  years <- as.integer(format(dates, "%Y"))
  raw <- rstan::extract(fit, pars = c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs"), permuted = FALSE, inc_warmup = FALSE)
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

  draws <- rstan::extract(fit, pars = c("C_pred", "X", "S_prop", "immune_prop", "R0_t", "R_eff_t", "expected_reported_cases"))
  summary <- bind_rows(
    summarise_matrix(draws$C_pred, dates, "posterior_predictive_cases"),
    summarise_matrix(draws$X, dates, "latent_infections"),
    summarise_matrix(draws$S_prop, dates, "susceptible_proportion"),
    summarise_matrix(draws$immune_prop, dates, "infection_derived_immune_proportion"),
    summarise_matrix(draws$R0_t, dates, "R0"),
    summarise_matrix(draws$R_eff_t, dates, "Reff"),
    summarise_matrix(draws$expected_reported_cases, dates, "expected_reported_cases")
  )
  write_csv(summary, file.path(paths$table, "trajectory_summary.csv"))
  get_summary <- function(variable_name) {
    dplyr::filter(summary, .data$variable == variable_name)
  }
  case_summary <- get_summary("posterior_predictive_cases")
  predictive_coverage <- mean(bundle$weekly_data$cases >= case_summary$q025 & bundle$weekly_data$cases <= case_summary$q975)
  predictive_checks <- tibble(
    metric = c("weekly_case_95pct_predictive_coverage", "weekly_case_median_RMSE", "weekly_case_median_MAE"),
    value = c(
      predictive_coverage,
      sqrt(mean((bundle$weekly_data$cases - case_summary$median)^2)),
      mean(abs(bundle$weekly_data$cases - case_summary$median))
    )
  )
  write_csv(predictive_checks, file.path(paths$table, "posterior_predictive_checks.csv"))

  annual <- bind_rows(lapply(sort(unique(years)), function(year) {
    index <- which(years == year)
    infections <- rowSums(draws$X[, index, drop = FALSE])
    expected_cases <- rowSums(draws$expected_reported_cases[, index, drop = FALSE])
    tibble(
      year = year, observed_cases = sum(bundle$weekly_data$cases[index]),
      expected_cases_q025 = quantile(expected_cases, .025), expected_cases_median = median(expected_cases), expected_cases_q975 = quantile(expected_cases, .975),
      infections_q025 = quantile(infections, .025), infections_median = median(infections), infections_q975 = quantile(infections, .975),
      S_start_median = median(draws$S_prop[, index[1]]), S_end_median = median(draws$S_prop[, tail(index, 1)]),
      annual_attack_median = median(infections / bundle$stan_data$N_start[index[1]])
    )
  }))
  write_csv(annual, file.path(paths$table, "annual_attack_decomposition.csv"))

  serology_windows <- tibble(
    site = c("Juazeiro do Norte", "Quixada"), observed_positive = c(103L, 289L), observed_n = c(404L, 409L),
    start = as.Date(c("2018-06-03", "2018-06-03")), end = as.Date(c("2018-12-30", "2019-12-29"))
  ) |>
    rowwise() |>
    mutate(
      state_immune_q025 = quantile(rowMeans(draws$immune_prop[, dates >= start & dates <= end, drop = FALSE]), .025),
      state_immune_median = median(rowMeans(draws$immune_prop[, dates >= start & dates <= end, drop = FALSE])),
      state_immune_q975 = quantile(rowMeans(draws$immune_prop[, dates >= start & dates <= end, drop = FALSE]), .975),
      observed_proportion = observed_positive / observed_n
    ) |>
    ungroup()
  write_csv(serology_windows, file.path(paths$table, "external_serology_consistency_check.csv"))
  make_run_record(file.path(paths$table, "v4_0_run_validation.md"), bundle, hmc, predictive_coverage)
  make_convergence_assessment(file.path(paths$table, "convergence_assessment.md"), fit, hmc)

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
  p5 <- ggplot(annual_long, aes(year, estimate, colour = type)) +
    geom_line() + geom_point() + geom_errorbar(data = filter(annual_long, type != "Observed reported cases"), aes(ymin = lower, ymax = upper), width = .15) +
    labs(title = "E. Annual reported-case validation", x = NULL, y = "Cases") + theme_v4
  p6 <- ggplot(serology_windows, aes(site, state_immune_median)) +
    geom_errorbar(aes(ymin = state_immune_q025, ymax = state_immune_q975), width = .15, colour = "#0072B2") +
    geom_point(colour = "#0072B2") + geom_point(aes(y = observed_proportion), colour = "#D55E00", shape = 18, size = 2.5) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    labs(title = "F. External serology consistency", subtitle = "Blue: state model; orange: municipality survey; not fitted", x = NULL, y = "Immune proportion") + theme_v4
  subtitle <- if (hmc_pass) "HMC gate passed; fixed q/importation scenario" else "UNCONVERGED - do not interpret posterior trajectories"
  figure <- (p1 | p2) / (p3 | p4) / (p5 | p6) + plot_annotation(
    title = "Ceara v4.0 minimal all-age no-vaccination renewal model",
    subtitle = subtitle, tag_levels = "A",
    theme = theme(plot.subtitle = element_text(face = "bold", colour = if (hmc_pass) "grey30" else "#D55E00"))
  )
  figure_path <- file.path(paths$figure, paste0("renewal_v4_0_minimal_no_vaccine_six_panel", if (hmc_pass) "" else "_UNCONVERGED", ".pdf"))
  ggsave(figure_path, figure, width = 183, height = 225, units = "mm", device = cairo_pdf)
  message("[v4.0 diagnostics] HMC gate: ", if (hmc_pass) "PASS" else "FAIL", "; figure: ", figure_path)
  invisible(list(hmc = hmc, predictive_checks = predictive_checks, annual = annual, serology = serology_windows))
}

if (sys.nframe() == 0L) run_v4_0_diagnostics()
