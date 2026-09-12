# =============================================================================
# 05_diagnose_renewal_v2_1_full_period.R
#
# Diagnostics for the v2.1 2015-2025 stress-test fit. This script is read-only:
# it does not compile Stan, sample, or alter the model.
#
# Environment overrides:
#   RENEWAL_V2_1_FULL_DIAG_FIT, RENEWAL_V2_1_FULL_DIAG_TAG
# =============================================================================

required_packages <- c("here", "rstan", "ggplot2", "dplyr", "tidyr", "tibble", "patchwork", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
})

summarise_draw_matrix <- function(draw_matrix, dates, variable) {
  if (ncol(draw_matrix) != length(dates)) {
    stop("Draw matrix and dates have incompatible dimensions for ", variable)
  }
  quantiles <- t(vapply(seq_len(ncol(draw_matrix)), function(t) {
    stats::quantile(draw_matrix[, t], probs = c(0.025, 0.5, 0.975))
  }, numeric(3)))
  tibble::tibble(
    week_start = dates,
    variable = variable,
    q025 = quantiles[, 1],
    median = quantiles[, 2],
    q975 = quantiles[, 3]
  )
}

consecutive_runs <- function(condition, dates) {
  runs <- rle(condition)
  end_index <- cumsum(runs$lengths)
  start_index <- end_index - runs$lengths + 1L
  tibble::tibble(
    flagged = runs$values,
    weeks = runs$lengths,
    start = dates[start_index],
    end = dates[end_index]
  ) |>
    dplyr::filter(flagged) |>
    dplyr::select(-flagged)
}

summarise_sampler_diagnostics <- function(fit) {
  parameter_summary <- summary(fit)$summary
  bfmi <- rstan::get_bfmi(fit)
  tibble::tibble(
    max_rhat = max(parameter_summary[, "Rhat"], na.rm = TRUE),
    n_rhat_gt_1.01 = sum(parameter_summary[, "Rhat"] > 1.01, na.rm = TRUE),
    min_ess = min(parameter_summary[, "n_eff"], na.rm = TRUE),
    median_ess = stats::median(parameter_summary[, "n_eff"], na.rm = TRUE),
    divergences = rstan::get_num_divergent(fit),
    max_treedepth_hits = rstan::get_num_max_treedepth(fit),
    min_bfmi = min(bfmi),
    median_bfmi = stats::median(bfmi)
  )
}

if (sys.nframe() == 0) {
  default_fit <- here::here("02_Script/stan/renewal_ceara_v2_1_full_period_fit.rds")
  fit_path <- Sys.getenv("RENEWAL_V2_1_FULL_DIAG_FIT", unset = default_fit)
  if (!file.exists(fit_path)) stop("Stress-test fit bundle not found: ", fit_path)
  bundle <- readRDS(fit_path)
  if (!isTRUE(bundle$config$stress_test) || !inherits(bundle$fit, "stanfit")) {
    stop("Expected a v2.1 full-period stress-test fit bundle")
  }

  # Thresholds create screening flags, not scientific pass/fail criteria.
  LOW_OBSERVED_CASES <- 5L
  EXTREME_R0_MEDIAN <- 10
  WIDE_R0_95_WIDTH <- 10
  BRIDGING_INCIDENCE_PER_100K <- 1
  BRIDGING_RUN_WEEKS <- 8L
  SUSCEPTIBLE_REPLENISHMENT_PER_YEAR <- 0.02
  REPORTING_LOWER <- 0.01
  REPORTING_UPPER <- 0.95
  REPORTING_SPAN <- 0.50

  pars <- c(
    "C_pred", "X", "S", "immune_prop", "R0_t", "R_eff_t", "rho_sym_t",
    "overall_detection", "p_symp", "expected_reported_cases", "sigma_R", "phi_obs"
  )
  draws <- rstan::extract(bundle$fit, pars = pars)
  weekly <- bundle$weekly_data
  dates <- weekly$week_start
  N_pop_t <- bundle$stan_data$N_pop_t
  if (length(N_pop_t) != length(dates) || any(!is.finite(N_pop_t)) || any(N_pop_t <= 0)) {
    stop("Stress-test bundle has invalid weekly population data")
  }

  susceptible_prop <- sweep(draws$S, 2, N_pop_t, "/")
  summaries <- dplyr::bind_rows(
    summarise_draw_matrix(draws$C_pred, dates, "posterior_predictive_cases"),
    summarise_draw_matrix(draws$X, dates, "latent_infections"),
    summarise_draw_matrix(susceptible_prop, dates, "susceptible_proportion"),
    summarise_draw_matrix(draws$immune_prop, dates, "cumulative_attack_rate"),
    summarise_draw_matrix(draws$R0_t, dates, "R0"),
    summarise_draw_matrix(draws$R_eff_t, dates, "Reff"),
    summarise_draw_matrix(draws$rho_sym_t, dates, "symptomatic_reporting"),
    summarise_draw_matrix(draws$overall_detection, dates, "overall_detection"),
    summarise_draw_matrix(draws$expected_reported_cases, dates, "expected_reported_cases")
  )
  get_summary <- function(name) dplyr::filter(summaries, variable == name)

  cases <- get_summary("posterior_predictive_cases")
  p_cases <- ggplot(cases, aes(week_start, median)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#56B4E9", alpha = 0.25) +
    geom_line(colour = "#0072B2", linewidth = 0.6) +
    geom_point(data = weekly, aes(week_start, cases), inherit.aes = FALSE, size = 0.45) +
    scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    labs(title = "Observed versus posterior-predicted weekly cases", x = NULL, y = "Cases") +
    theme_classic(base_size = 10)

  p_latent <- ggplot(get_summary("latent_infections"), aes(week_start, median)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#D55E00", alpha = 0.20) +
    geom_line(colour = "#D55E00", linewidth = 0.6) +
    scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    labs(title = "Latent weekly infections", x = NULL, y = "Infections") +
    theme_classic(base_size = 10)

  state_plot_data <- dplyr::bind_rows(
    get_summary("susceptible_proportion") |> dplyr::mutate(state = "Susceptible"),
    get_summary("cumulative_attack_rate") |> dplyr::mutate(state = "Cumulative attack rate")
  )
  p_states <- ggplot(state_plot_data, aes(week_start, median, colour = state, fill = state)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.14, colour = NA) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = c("Susceptible" = "#009E73", "Cumulative attack rate" = "#CC79A7")) +
    scale_fill_manual(values = c("Susceptible" = "#009E73", "Cumulative attack rate" = "#CC79A7")) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
    labs(title = "Susceptibility and cumulative attack rate", x = NULL, y = "Proportion") +
    theme_classic(base_size = 10) +
    theme(legend.position = "top")

  reproduction <- dplyr::bind_rows(
    get_summary("R0") |> dplyr::mutate(number = "R0"),
    get_summary("Reff") |> dplyr::mutate(number = "Reff")
  ) |>
    dplyr::group_by(number) |>
    dplyr::mutate(t = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::filter(t > bundle$stan_data$seed_weeks)
  p_reproduction <- ggplot(reproduction, aes(week_start, median, colour = number, fill = number)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.14, colour = NA) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~number, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = c("R0" = "#0072B2", "Reff" = "#D55E00")) +
    scale_fill_manual(values = c("R0" = "#0072B2", "Reff" = "#D55E00")) +
    labs(title = "Basic and effective reproduction numbers", x = NULL, y = "Reproduction number") +
    theme_classic(base_size = 10) +
    theme(legend.position = "top")

  reporting <- dplyr::bind_rows(
    get_summary("symptomatic_reporting") |> dplyr::mutate(component = "Symptomatic reporting"),
    get_summary("overall_detection") |> dplyr::mutate(component = "Overall detection")
  )
  p_reporting <- ggplot(reporting, aes(week_start, median, colour = component, fill = component)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.14, colour = NA) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = c("Symptomatic reporting" = "#76558F", "Overall detection" = "#E69F00")) +
    scale_fill_manual(values = c("Symptomatic reporting" = "#76558F", "Overall detection" = "#E69F00")) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
    labs(title = "Long-term case ascertainment", x = NULL, y = "Probability") +
    theme_classic(base_size = 10) +
    theme(legend.position = "top")

  p_population <- ggplot(tibble::tibble(week_start = dates, population = N_pop_t), aes(week_start, population)) +
    geom_line(colour = "#4D4D4D", linewidth = 0.6) +
    scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    labs(title = "Observed state population supplied to v2.1", x = NULL, y = "Population") +
    theme_classic(base_size = 10)

  parameter_summary <- summary(bundle$fit)$summary
  parameter_diagnostics <- tibble::tibble(
    parameter = rownames(parameter_summary),
    ess = parameter_summary[, "n_eff"],
    rhat = parameter_summary[, "Rhat"]
  )
  sampler_diagnostics <- summarise_sampler_diagnostics(bundle$fit)
  diagnostic_plot_data <- sampler_diagnostics |>
    tidyr::pivot_longer(dplyr::everything(), names_to = "metric", values_to = "value")
  p_hmc <- ggplot(diagnostic_plot_data, aes(metric, value)) +
    geom_col(fill = "#4D4D4D") +
    geom_text(aes(label = format(round(value, 3), trim = TRUE)), vjust = -0.25, size = 2.8) +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    labs(title = "HMC diagnostics", subtitle = "Per-parameter Rhat and ESS are saved as a table", x = NULL, y = NULL) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  r0 <- get_summary("R0") |>
    dplyr::mutate(observed_cases = weekly$cases, r0_width = q975 - q025)
  r0_flags <- r0 |>
    dplyr::filter(
      observed_cases <= LOW_OBSERVED_CASES,
      median > EXTREME_R0_MEDIAN | r0_width > WIDE_R0_95_WIDTH
    ) |>
    dplyr::transmute(
      flag = "Extreme or weakly identified R0 in low-incidence week",
      week_start,
      value = median,
      threshold = ifelse(median > EXTREME_R0_MEDIAN, EXTREME_R0_MEDIAN, WIDE_R0_95_WIDTH),
      detail = paste0("observed cases=", observed_cases, "; 95% R0 width=", round(r0_width, 2))
    )

  latent <- get_summary("latent_infections") |>
    dplyr::mutate(
      observed_cases = weekly$cases,
      latent_per_100k = median / N_pop_t * 1e5,
      later_period = as.integer(format(week_start, "%Y")) >= 2020L
    )
  bridging_runs <- consecutive_runs(
    latent$later_period & latent$observed_cases <= LOW_OBSERVED_CASES &
      latent$latent_per_100k > BRIDGING_INCIDENCE_PER_100K,
    latent$week_start
  ) |>
    dplyr::filter(weeks >= BRIDGING_RUN_WEEKS)
  bridging_flags <- bridging_runs |>
    dplyr::transmute(
      flag = "Persistent latent incidence during low observed incidence",
      week_start = start,
      value = weeks,
      threshold = BRIDGING_RUN_WEEKS,
      detail = paste0("run through ", end, "; median latent incidence > ", BRIDGING_INCIDENCE_PER_100K, " per 100,000/week")
    )

  susceptible <- get_summary("susceptible_proportion") |>
    dplyr::mutate(change_over_52_weeks = median - dplyr::lag(median, 52L))
  replenishment_flags <- susceptible |>
    dplyr::filter(change_over_52_weeks > SUSCEPTIBLE_REPLENISHMENT_PER_YEAR) |>
    dplyr::transmute(
      flag = "Large one-year susceptible-proportion replenishment",
      week_start,
      value = change_over_52_weeks,
      threshold = SUSCEPTIBLE_REPLENISHMENT_PER_YEAR,
      detail = "Posterior median susceptible proportion increased over 52 weeks"
    )

  reporting_summary <- get_summary("symptomatic_reporting")
  reporting_flags <- reporting_summary |>
    dplyr::filter(median < REPORTING_LOWER | median > REPORTING_UPPER) |>
    dplyr::transmute(
      flag = "Extreme long-term symptomatic reporting probability",
      week_start,
      value = median,
      threshold = ifelse(median < REPORTING_LOWER, REPORTING_LOWER, REPORTING_UPPER),
      detail = "Posterior median outside pre-specified screening bounds"
    )
  reporting_span <- max(reporting_summary$median) - min(reporting_summary$median)
  if (reporting_span > REPORTING_SPAN) {
    reporting_flags <- dplyr::bind_rows(
      reporting_flags,
      tibble::tibble(
        flag = "Large long-term symptomatic-reporting trend",
        week_start = as.Date(NA),
        value = reporting_span,
        threshold = REPORTING_SPAN,
        detail = "Range of posterior median reporting probability across the stress-test period"
      )
    )
  }

  sampler_flags <- tibble::tibble(
    flag = c(
      "Sampler Rhat above 1.01",
      "Sampler divergences",
      "Sampler maximum-treedepth hits",
      "Sampler BFMI below 0.3"
    ),
    week_start = as.Date(NA),
    value = c(
      sampler_diagnostics$max_rhat,
      sampler_diagnostics$divergences,
      sampler_diagnostics$max_treedepth_hits,
      sampler_diagnostics$min_bfmi
    ),
    threshold = c(1.01, 0, 0, 0.3),
    detail = c(
      "Maximum parameter Rhat",
      "Total divergent transitions",
      "Total transitions at maximum treedepth",
      "Minimum chain E-BFMI"
    )
  ) |>
    dplyr::filter(
      (flag == "Sampler Rhat above 1.01" & value > threshold) |
        (flag %in% c("Sampler divergences", "Sampler maximum-treedepth hits") & value > threshold) |
        (flag == "Sampler BFMI below 0.3" & value < threshold)
    )

  baseline_path <- Sys.getenv(
    "RENEWAL_V2_1_BASELINE_DIAG_FIT",
    unset = here::here("02_Script/stan/renewal_ceara_v2_1_fit.rds")
  )
  sampler_worsening_flags <- tibble::tibble(
    flag = character(), week_start = as.Date(character()), value = numeric(),
    threshold = numeric(), detail = character()
  )
  make_sampler_worsening_flag <- function(flag, value, threshold, detail) {
    tibble::tibble(
      flag = flag,
      week_start = as.Date(NA),
      value = value,
      threshold = threshold,
      detail = detail
    )
  }
  baseline_sampler_diagnostics <- NULL
  if (file.exists(baseline_path)) {
    baseline_bundle <- readRDS(baseline_path)
    if (inherits(baseline_bundle$fit, "stanfit")) {
      baseline_sampler_diagnostics <- summarise_sampler_diagnostics(baseline_bundle$fit)
      sampler_worsening_flags <- dplyr::bind_rows(
        sampler_worsening_flags,
        if (sampler_diagnostics$max_rhat > baseline_sampler_diagnostics$max_rhat + 0.005) {
          make_sampler_worsening_flag("Worsened maximum Rhat versus v2.1 baseline", sampler_diagnostics$max_rhat, baseline_sampler_diagnostics$max_rhat, "Stress-test maximum Rhat exceeds baseline by more than 0.005")
        },
        if (sampler_diagnostics$min_ess < 0.5 * baseline_sampler_diagnostics$min_ess) {
          make_sampler_worsening_flag("Worsened minimum ESS versus v2.1 baseline", sampler_diagnostics$min_ess, 0.5 * baseline_sampler_diagnostics$min_ess, "Stress-test minimum ESS is less than half the baseline value")
        },
        if (sampler_diagnostics$divergences > baseline_sampler_diagnostics$divergences) {
          make_sampler_worsening_flag("More divergences than v2.1 baseline", sampler_diagnostics$divergences, baseline_sampler_diagnostics$divergences, "Stress-test divergent-transition count exceeds baseline")
        },
        if (sampler_diagnostics$max_treedepth_hits > baseline_sampler_diagnostics$max_treedepth_hits) {
          make_sampler_worsening_flag("More treedepth hits than v2.1 baseline", sampler_diagnostics$max_treedepth_hits, baseline_sampler_diagnostics$max_treedepth_hits, "Stress-test maximum-treedepth hit count exceeds baseline")
        },
        if (sampler_diagnostics$min_bfmi < baseline_sampler_diagnostics$min_bfmi - 0.05) {
          make_sampler_worsening_flag("Worsened minimum BFMI versus v2.1 baseline", sampler_diagnostics$min_bfmi, baseline_sampler_diagnostics$min_bfmi - 0.05, "Stress-test minimum E-BFMI is more than 0.05 below baseline")
        }
      )
    }
  }

  structural_flags <- dplyr::bind_rows(
    r0_flags, bridging_flags, replenishment_flags, reporting_flags, sampler_flags,
    sampler_worsening_flags
  )

  figure_dir <- here::here("03_Output/figures/renewal_v2_1/full_period_stress_test")
  table_dir <- here::here("03_Output/tables/renewal_v2_1/full_period_stress_test")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  tag <- Sys.getenv("RENEWAL_V2_1_FULL_DIAG_TAG", unset = "")
  suffix <- if (nzchar(tag)) paste0("_", tag) else ""
  output_base <- paste0("renewal_v2_1_full_period_stress_test", suffix)
  trajectory_figure <- (p_cases | p_latent) /
    (p_states | p_reproduction) /
    (p_reporting | p_population) +
    patchwork::plot_annotation(
      title = "Ceara renewal v2.1 full-period stress test, 2015-2025",
      subtitle = "Diagnostic extension only; not the final full-history model",
      tag_levels = "A"
    )
  trajectory_path <- file.path(
    figure_dir, paste0("renewal_v2_1_trajectories_ceara_full", suffix, ".png")
  )
  ggplot2::ggsave(
    trajectory_path, trajectory_figure,
    width = 15, height = 14,
    units = "in", dpi = 220, bg = "white"
  )
  message("[save] ", trajectory_path)
  save_plot <- function(name, plot, width, height) {
    path <- file.path(figure_dir, paste0(output_base, "_", name, ".png"))
    ggplot2::ggsave(path, plot, width = width, height = height, units = "in", dpi = 220, bg = "white")
    message("[save] ", path)
  }
  save_plot("cases", p_cases, 10, 4.5)
  save_plot("latent_infections", p_latent, 10, 4.5)
  save_plot("states", p_states, 10, 4.5)
  save_plot("reproduction", p_reproduction, 10, 7)
  save_plot("reporting", p_reporting, 10, 4.5)
  save_plot("population", p_population, 10, 4.5)
  save_plot("hmc", p_hmc, 8, 6)

  utils::write.csv(summaries, file.path(table_dir, paste0(output_base, "_trajectory_summary.csv")), row.names = FALSE)
  utils::write.csv(parameter_diagnostics, file.path(table_dir, paste0(output_base, "_parameter_diagnostics.csv")), row.names = FALSE)
  utils::write.csv(sampler_diagnostics, file.path(table_dir, paste0(output_base, "_sampler_diagnostics.csv")), row.names = FALSE)
  if (!is.null(baseline_sampler_diagnostics)) {
    utils::write.csv(baseline_sampler_diagnostics, file.path(table_dir, paste0(output_base, "_baseline_sampler_diagnostics.csv")), row.names = FALSE)
  }
  utils::write.csv(structural_flags, file.path(table_dir, paste0(output_base, "_structural_flags.csv")), row.names = FALSE)
  if (nrow(structural_flags)) {
    message("[flags] ", nrow(structural_flags), " stress-test warning(s); inspect the structural-flags table.")
  } else {
    message("[flags] No warnings crossed the pre-specified screening thresholds.")
  }
}
