# =============================================================================
# 03_plot_renewal_v2_1_fit.R
#
# Publication-quality trajectory figure for the Ceara renewal v2.1 model.
# Mirrors 03_plot_renewal_v2_fit.R's six-panel layout, adapted for v2.1's
# time-varying population (S is a raw count, not a fixed-N proportion),
# estimated p_symp, and two serosurvey sites each linked to a collection
# window rather than a single date.
# Set RENEWAL_V2_1_FULL_PERIOD=true to render the diagnostic 2015-2025 fit
# with the same layout and a `_full` output suffix.
# =============================================================================

required_packages <- c("here", "rstan", "ggplot2", "dplyr", "patchwork", "scales", "tibble")
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
  library(patchwork)
})

COL_INK <- "#202124"
COL_MUTED <- "#6B7280"
COL_GRID <- "#E5E7EB"
COL_NAVY <- "#315A7D"
COL_BLUE <- "#56B4E9"
COL_VERMILLION <- "#D55E00"
COL_ORANGE <- "#E69F00"
COL_GREEN <- "#009E73"
COL_PURPLE <- "#76558F"
COL_TEAL <- "#0072B2"

theme_publication <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_family = "Arial", base_size = base_size) +
    ggplot2::theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, color = COL_MUTED, margin = margin(b = 7)),
      axis.title = element_text(size = base_size, color = COL_INK),
      axis.text = element_text(size = base_size - 0.5, color = COL_MUTED),
      axis.line = element_line(color = COL_INK, linewidth = 0.35),
      axis.ticks = element_line(color = COL_INK, linewidth = 0.35),
      panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.width = grid::unit(13, "pt"),
      legend.key.height = grid::unit(7, "pt"),
      plot.margin = margin(7, 8, 7, 7)
    )
}

summarise_trajectory <- function(draw_matrix, prefix) {
  stopifnot(length(dim(draw_matrix)) == 2L)
  tibble::tibble(
    variable = prefix,
    t = seq_len(ncol(draw_matrix)),
    q025 = apply(draw_matrix, 2, stats::quantile, probs = 0.025),
    q25 = apply(draw_matrix, 2, stats::quantile, probs = 0.25),
    median = apply(draw_matrix, 2, stats::median),
    q75 = apply(draw_matrix, 2, stats::quantile, probs = 0.75),
    q975 = apply(draw_matrix, 2, stats::quantile, probs = 0.975)
  )
}

add_interval_layers <- function(plot, data, colour, fill = colour) {
  plot +
    geom_ribbon(
      data = data, aes(ymin = q025, ymax = q975),
      fill = fill, alpha = 0.13, colour = NA
    ) +
    geom_ribbon(
      data = data, aes(ymin = q25, ymax = q75),
      fill = fill, alpha = 0.24, colour = NA
    ) +
    geom_line(data = data, aes(y = median), colour = colour, linewidth = 0.65)
}

if (sys.nframe() == 0) {
  FULL_PERIOD <- tolower(Sys.getenv("RENEWAL_V2_1_FULL_PERIOD", unset = "false")) %in%
    c("1", "true", "yes")
  if (FULL_PERIOD) {
    FIT_TAG <- Sys.getenv("RENEWAL_V2_1_FULL_TAG", unset = "")
    fit_prefix <- "renewal_ceara_v2_1_full_period_fit"
  } else {
    FIT_TAG <- Sys.getenv("RENEWAL_V2_1_TAG", unset = "")
    fit_prefix <- "renewal_ceara_v2_1_fit"
  }
  suffix <- if (nzchar(FIT_TAG)) paste0("_", FIT_TAG) else ""
  fit_path <- here::here(
    "02_Script/stan", paste0(fit_prefix, suffix, ".rds")
  )
  if (!file.exists(fit_path)) stop("Fit bundle not found: ", fit_path)

  bundle <- readRDS(fit_path)
  fit <- bundle$fit
  weekly <- bundle$weekly_data
  dates <- weekly$week_start
  N_pop_t <- bundle$stan_data$N_pop_t

  pars <- c(
    "X", "S", "immune_prop", "R0_t", "R_eff_t",
    "rho_sym_t", "overall_detection", "expected_reported_cases", "C_pred",
    "p_symp", "state_attack_window", "p_site"
  )
  draws <- rstan::extract(fit, pars = pars, permuted = TRUE)
  S_prop <- sweep(draws$S, 2, N_pop_t, "/")

  summaries <- dplyr::bind_rows(
    summarise_trajectory(draws$X, "latent_infections"),
    summarise_trajectory(S_prop, "susceptible_prop"),
    summarise_trajectory(draws$immune_prop, "immune_prop"),
    summarise_trajectory(draws$R0_t, "R0_t"),
    summarise_trajectory(draws$R_eff_t, "R_eff_t"),
    summarise_trajectory(draws$rho_sym_t, "rho_sym_t"),
    summarise_trajectory(draws$overall_detection, "overall_detection"),
    summarise_trajectory(draws$expected_reported_cases, "expected_reported_cases"),
    summarise_trajectory(draws$C_pred, "posterior_predictive_cases")
  ) |>
    dplyr::mutate(week_start = dates[t])

  get_summary <- function(name) dplyr::filter(summaries, variable == name)
  date_breaks <- if (FULL_PERIOD) {
    seq(as.Date("2015-01-01"), as.Date("2026-01-01"), by = "2 years")
  } else {
    seq(as.Date("2015-01-01"), as.Date("2020-01-01"), by = "1 year")
  }
  x_scale <- scale_x_date(
    breaks = date_breaks,
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.015))
  )

  expected <- get_summary("expected_reported_cases")
  predictive <- get_summary("posterior_predictive_cases")
  cases_df <- tibble::tibble(week_start = dates, observed = weekly$cases)

  p_cases <- ggplot() +
    geom_ribbon(
      data = predictive,
      aes(week_start, ymin = q025, ymax = q975),
      fill = COL_BLUE, alpha = 0.15
    ) +
    geom_ribbon(
      data = expected,
      aes(week_start, ymin = q25, ymax = q75),
      fill = COL_NAVY, alpha = 0.24
    ) +
    geom_line(
      data = expected, aes(week_start, median, colour = "Model expectation"),
      linewidth = 0.65
    ) +
    geom_point(
      data = cases_df, aes(week_start, observed, colour = "Observed cases"),
      size = 0.55, alpha = 0.72
    ) +
    scale_colour_manual(values = c(
      "Observed cases" = COL_INK,
      "Model expectation" = COL_NAVY
    )) +
    scale_y_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = c(0, 0.06))) +
    x_scale +
    labs(
      title = "Reported cases and posterior prediction",
      subtitle = "Shading: 95% predictive interval; dark band: 50% CrI",
      x = NULL, y = "Weekly reported cases"
    ) +
    theme_publication()

  latent <- get_summary("latent_infections")
  p_infections <- add_interval_layers(
    ggplot(latent, aes(week_start)), latent, COL_VERMILLION
  ) +
    x_scale +
    scale_y_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = c(0, 0.06))) +
    labs(
      title = "Latent infection incidence",
      subtitle = "Median, 50% and 95% credible intervals",
      x = NULL, y = "True infections per week"
    ) +
    theme_publication()

  immunity <- dplyr::bind_rows(
    get_summary("susceptible_prop") |> mutate(state = "Susceptible"),
    get_summary("immune_prop") |> mutate(state = "Infection-derived immune")
  )
  p_immunity <- ggplot(immunity, aes(week_start, colour = state, fill = state)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.10, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.7) +
    scale_colour_manual(values = c(
      "Susceptible" = COL_GREEN,
      "Infection-derived immune" = COL_PURPLE
    )) +
    scale_fill_manual(values = c(
      "Susceptible" = COL_GREEN,
      "Infection-derived immune" = COL_PURPLE
    )) +
    scale_y_continuous(
      limits = c(0, 1), labels = scales::label_percent(accuracy = 1),
      breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.01))
    ) +
    x_scale +
    labs(
      title = "Population susceptibility and immunity",
      subtitle = "Weekly observed population; lifelong immunity",
      x = NULL, y = "Population proportion"
    ) +
    theme_publication()

  reproduction <- dplyr::bind_rows(
    get_summary("R0_t") |>
      filter(t > bundle$stan_data$seed_weeks) |>
      mutate(number = "R0(t)"),
    get_summary("R_eff_t") |>
      filter(t > bundle$stan_data$seed_weeks) |>
      mutate(number = "Reff(t)")
  )
  p_reproduction <- ggplot(reproduction, aes(week_start, colour = number, fill = number)) +
    geom_hline(yintercept = 1, colour = COL_MUTED, linewidth = 0.35, linetype = "22") +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.09, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.65) +
    scale_colour_manual(
      values = c("R0(t)" = COL_NAVY, "Reff(t)" = COL_VERMILLION),
      labels = c("R0(t)" = expression(R[0](t)), "Reff(t)" = expression(R[eff](t)))
    ) +
    scale_fill_manual(
      values = c("R0(t)" = COL_NAVY, "Reff(t)" = COL_VERMILLION),
      labels = c("R0(t)" = expression(R[0](t)), "Reff(t)" = expression(R[eff](t)))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.06))) +
    x_scale +
    labs(
      title = "Time-varying reproduction numbers",
      subtitle = expression(R[eff](t) == R[0](t) %*% S(t-1)/N(t-1)),
      x = NULL, y = "Reproduction number"
    ) +
    theme_publication()

  p_symp_summary <- stats::quantile(draws$p_symp, probs = c(0.025, 0.5, 0.975))
  reporting <- dplyr::bind_rows(
    get_summary("rho_sym_t") |> mutate(probability = "Symptomatic-case reporting"),
    get_summary("overall_detection") |> mutate(probability = "Infection-to-report")
  )
  p_reporting <- ggplot(reporting, aes(week_start, colour = probability, fill = probability)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.10, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.7) +
    scale_colour_manual(values = c(
      "Symptomatic-case reporting" = COL_PURPLE,
      "Infection-to-report" = COL_ORANGE
    )) +
    scale_fill_manual(values = c(
      "Symptomatic-case reporting" = COL_PURPLE,
      "Infection-to-report" = COL_ORANGE
    )) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA), expand = expansion(mult = c(0, 0.07))) +
    x_scale +
    labs(
      title = "Case ascertainment",
      subtitle = sprintf(
        "Estimated symptomatic proportion: %.2f [%.2f, %.2f]",
        p_symp_summary[2], p_symp_summary[1], p_symp_summary[3]
      ),
      x = NULL, y = "Probability"
    ) +
    theme_publication()

  attack <- get_summary("immune_prop")
  sero_sites <- bundle$config$serology_sites
  sero_window_start <- as.Date(bundle$config$serology_window_start)
  sero_window_end <- as.Date(bundle$config$serology_window_end)
  sero_mid <- sero_window_start + (sero_window_end - sero_window_start) / 2
  sero_positive <- bundle$config$serology_positive
  sero_n <- bundle$config$serology_n

  sero_observed_df <- dplyr::bind_rows(lapply(seq_along(sero_sites), function(j) {
    interval <- stats::binom.test(sero_positive[j], sero_n[j])$conf.int
    tibble::tibble(
      site = sero_sites[j], week_start = sero_mid[j],
      estimate = sero_positive[j] / sero_n[j], lower = interval[1], upper = interval[2]
    )
  }))
  site_summary <- dplyr::bind_rows(lapply(seq_along(sero_sites), function(j) {
    site_draws <- draws$p_site[, j]
    tibble::tibble(
      site = sero_sites[j], week_start = sero_mid[j] + 25,
      median = stats::median(site_draws),
      q025 = stats::quantile(site_draws, 0.025), q975 = stats::quantile(site_draws, 0.975)
    )
  }))
  window_rects <- tibble::tibble(
    site = sero_sites, xmin = sero_window_start, xmax = sero_window_end
  )

  p_attack <- ggplot(attack, aes(week_start)) +
    geom_rect(
      data = window_rects, inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = site),
      alpha = 0.07
    ) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = COL_GREEN, alpha = 0.14) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = COL_GREEN, alpha = 0.24) +
    geom_line(aes(y = median, colour = "Ceara attack rate"), linewidth = 0.7) +
    geom_errorbar(
      data = sero_observed_df,
      aes(week_start, ymin = lower, ymax = upper, colour = "Observed site seroprevalence"),
      width = 18, linewidth = 0.55
    ) +
    geom_point(
      data = sero_observed_df,
      aes(week_start, estimate, colour = "Observed site seroprevalence"),
      shape = 18, size = 2.3
    ) +
    geom_errorbar(
      data = site_summary,
      aes(week_start, ymin = q025, ymax = q975, colour = "Model-implied site attack rate"),
      width = 18, linewidth = 0.55
    ) +
    geom_point(
      data = site_summary,
      aes(week_start, median, colour = "Model-implied site attack rate"),
      shape = 16, size = 1.8
    ) +
    scale_colour_manual(values = c(
      "Ceara attack rate" = COL_GREEN,
      "Observed site seroprevalence" = COL_ORANGE,
      "Model-implied site attack rate" = COL_PURPLE
    ), labels = c(
      "Ceara attack rate" = "Ceara",
      "Observed site seroprevalence" = "Site observed",
      "Model-implied site attack rate" = "Site model"
    )) +
    scale_fill_manual(values = c("Juazeiro do Norte" = COL_ORANGE, "Quixada" = COL_TEAL), guide = "none") +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA), expand = expansion(mult = c(0, 0.07))) +
    x_scale +
    labs(
      title = "Cumulative infection and serological anchors",
      subtitle = "Shaded fields: Juazeiro do Norte (orange) and Quixadá (teal) survey windows",
      x = NULL, y = "Cumulative proportion infected"
    ) +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_publication(base_size = 8) +
    theme(legend.text = element_text(size = 7.2))

  # ---- Serology hierarchy: state pooling down to site-level validation ----
  # A three-stage slopegraph per site: the state-window mean (before any
  # geographic offset) -> the site-level posterior after the non-centered
  # geographic offset is applied (p_site) -> the raw observed seroprevalence.
  # This makes the partial-pooling structure legible in a way the time-series
  # attack-rate panel (p_attack) does not: how far the offset moves each site
  # away from the shared Ceara baseline, and whether that lands near the data.
  hierarchy_stages <- dplyr::bind_rows(lapply(seq_along(sero_sites), function(j) {
    dplyr::bind_rows(
      tibble::tibble(
        site = sero_sites[j], stage = "Ceara (state window mean)",
        median = stats::median(draws$state_attack_window[, j]),
        q025 = stats::quantile(draws$state_attack_window[, j], 0.025),
        q975 = stats::quantile(draws$state_attack_window[, j], 0.975)
      ),
      tibble::tibble(
        site = sero_sites[j], stage = "Site model (+ geographic offset)",
        median = stats::median(draws$p_site[, j]),
        q025 = stats::quantile(draws$p_site[, j], 0.025),
        q975 = stats::quantile(draws$p_site[, j], 0.975)
      )
    )
  }))
  hierarchy_observed <- tibble::tibble(
    site = sero_sites, stage = "Site observed",
    median = sero_positive / sero_n,
    q025 = vapply(seq_along(sero_sites), function(j) {
      stats::binom.test(sero_positive[j], sero_n[j])$conf.int[1]
    }, numeric(1)),
    q975 = vapply(seq_along(sero_sites), function(j) {
      stats::binom.test(sero_positive[j], sero_n[j])$conf.int[2]
    }, numeric(1))
  )
  hierarchy_data <- dplyr::bind_rows(hierarchy_stages, hierarchy_observed) |>
    dplyr::mutate(
      stage = factor(stage, levels = c(
        "Ceara (state window mean)", "Site model (+ geographic offset)", "Site observed"
      )),
      # The two sites' state-window means sit within a fraction of a percent
      # of each other, so a shared vjust makes their text labels collide;
      # placing one series above its point and the other below keeps every
      # stage legible, including where the site-level estimates diverge.
      label_vjust = dplyr::if_else(site == "Juazeiro do Norte", 1.9, -1.1)
    )

  p_hierarchy <- ggplot(hierarchy_data, aes(stage, median, colour = site, group = site)) +
    geom_line(position = position_dodge(width = 0.2), linewidth = 0.5, alpha = 0.6) +
    geom_pointrange(
      aes(ymin = q025, ymax = q975),
      position = position_dodge(width = 0.2), size = 0.45, linewidth = 0.6, fatten = 2.5
    ) +
    geom_text(
      aes(label = scales::label_percent(accuracy = 0.1)(median), vjust = label_vjust),
      position = position_dodge(width = 0.2), size = 2.9, show.legend = FALSE
    ) +
    scale_colour_manual(values = c("Juazeiro do Norte" = COL_ORANGE, "Quixada" = COL_TEAL)) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
    labs(
      title = "Serology hierarchy: state pooling to site-level validation",
      subtitle = "State window mean → geographic-offset-adjusted site estimate → observed seroprevalence",
      x = NULL, y = "Cumulative attack proportion"
    ) +
    theme_publication()

  figure <- (p_cases | p_infections) /
    (p_immunity | p_reproduction) /
    (p_reporting | p_attack) +
    patchwork::plot_annotation(
      title = "Chikungunya transmission and population susceptibility in Ceará, Brazil",
      subtitle = if (FULL_PERIOD) {
        "Bayesian renewal model v2.1 stress test, weekly observations from 2015 to 2025"
      } else {
        "Bayesian renewal model v2.1, weekly observations from 2015 to 2019"
      },
      tag_levels = "A",
      theme = theme(
        text = element_text(family = "Arial", colour = COL_INK),
        plot.title = element_text(size = 12.5, face = "bold", margin = margin(b = 3)),
        plot.subtitle = element_text(size = 9.5, colour = COL_MUTED, margin = margin(b = 8)),
        plot.tag = element_text(size = 11, face = "bold")
      )
    )

  figure_dir <- here::here("03_Output/figures")
  table_dir <- here::here("03_Output/tables")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  output_base <- paste0(
    if (FULL_PERIOD) "renewal_v2_1_trajectories_ceara_full" else "renewal_v2_1_trajectories_ceara",
    suffix
  )
  pdf_path <- file.path(figure_dir, paste0(output_base, ".pdf"))
  tiff_path <- file.path(figure_dir, paste0(output_base, ".tiff"))
  png_path <- file.path(figure_dir, paste0(output_base, ".png"))
  table_path <- file.path(
    table_dir,
    paste0(
      if (FULL_PERIOD) "renewal_v2_1_trajectory_summary_full" else "renewal_v2_1_trajectory_summary",
      suffix, ".csv"
    )
  )

  ggsave(pdf_path, figure, width = 183, height = 225, units = "mm", device = cairo_pdf)
  ggsave(tiff_path, figure, width = 183, height = 225, units = "mm", dpi = 600, compression = "lzw")
  ggsave(png_path, figure, width = 183, height = 225, units = "mm", dpi = 300, bg = "white")
  utils::write.csv(summaries, table_path, row.names = FALSE)

  message("[save] ", pdf_path)
  message("[save] ", tiff_path)
  message("[save] ", png_path)
  message("[save] ", table_path)

  hierarchy_base <- paste0(
    if (FULL_PERIOD) "renewal_v2_1_serology_hierarchy_full" else "renewal_v2_1_serology_hierarchy",
    suffix
  )
  hierarchy_png <- file.path(figure_dir, paste0(hierarchy_base, ".png"))
  hierarchy_pdf <- file.path(figure_dir, paste0(hierarchy_base, ".pdf"))
  hierarchy_table_path <- file.path(table_dir, paste0(hierarchy_base, ".csv"))

  ggsave(hierarchy_pdf, p_hierarchy, width = 140, height = 110, units = "mm", device = cairo_pdf)
  ggsave(hierarchy_png, p_hierarchy, width = 140, height = 110, units = "mm", dpi = 300, bg = "white")
  utils::write.csv(hierarchy_data, hierarchy_table_path, row.names = FALSE)

  message("[save] ", hierarchy_pdf)
  message("[save] ", hierarchy_png)
  message("[save] ", hierarchy_table_path)
}
