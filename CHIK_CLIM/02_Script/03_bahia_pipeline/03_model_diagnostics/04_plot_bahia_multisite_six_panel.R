# Six-panel trajectory figure for the Bahia multi-site geo-adjusted
# serology model (renewal_ceara_v4_9_bahia_multisite_serology.stan), same
# layout as 03_plot_bahia_six_panel.R. Panel F now shows the 6 primary
# survey sites: observed prevalence vs state-level model immunity
# (p_state_window, no geo offset) vs geographically-adjusted site
# prediction (p_site_window).
#
# Figures saved under 03_Output/figures/renewal_bahia_v4_9_multisite_serology/.
#
# Usage: QVAL=0.05 STAGE=pilot Rscript 08_plot_bahia_multisite_six_panel.R

required_packages <- c("here", "rstan", "ggplot2", "dplyr", "patchwork", "scales", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(ggplot2); library(dplyr); library(patchwork) })

COL_INK <- "#202124"; COL_MUTED <- "#6B7280"; COL_GRID <- "#E5E7EB"
COL_NAVY <- "#315A7D"; COL_BLUE <- "#56B4E9"; COL_VERMILLION <- "#D55E00"
COL_ORANGE <- "#E69F00"; COL_GREEN <- "#009E73"; COL_PURPLE <- "#76558F"

theme_publication <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, color = COL_MUTED, margin = margin(b = 7)),
      axis.title = element_text(size = base_size, color = COL_INK),
      axis.text = element_text(size = base_size - 0.5, color = COL_MUTED),
      axis.line = element_line(color = COL_INK, linewidth = 0.35),
      axis.ticks = element_line(color = COL_INK, linewidth = 0.35),
      panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "top", legend.justification = "left", legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.width = grid::unit(13, "pt"), legend.key.height = grid::unit(7, "pt"),
      plot.margin = margin(7, 8, 7, 7)
    )
}

summarise_trajectory <- function(draw_matrix, prefix) {
  stopifnot(length(dim(draw_matrix)) == 2L)
  tibble::tibble(
    variable = prefix, t = seq_len(ncol(draw_matrix)),
    q025 = apply(draw_matrix, 2, stats::quantile, probs = 0.025),
    q25 = apply(draw_matrix, 2, stats::quantile, probs = 0.25),
    median = apply(draw_matrix, 2, stats::median),
    q75 = apply(draw_matrix, 2, stats::quantile, probs = 0.75),
    q975 = apply(draw_matrix, 2, stats::quantile, probs = 0.975)
  )
}

add_interval_layers <- function(plot, data, colour, fill = colour) {
  plot +
    geom_ribbon(data = data, aes(ymin = q025, ymax = q975), fill = fill, alpha = 0.13, colour = NA) +
    geom_ribbon(data = data, aes(ymin = q25, ymax = q75), fill = fill, alpha = 0.24, colour = NA) +
    geom_line(data = data, aes(y = median), colour = colour, linewidth = 0.65)
}

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
base_dir <- file.path(root, "03_Output/model_fits/bahia/v4_9_replication")

plot_bahia_multisite_six_panel <- function(qval, stage = "pilot") {
  tag <- paste0("q", sprintf("%.2f", qval))
  bundle <- readRDS(file.path(base_dir, "outputs_multisite_serology", tag, paste0("renewal_bahia_multisite_fit_", tag, ".rds")))
  fit <- bundle$fit
  weekly <- bundle$weekly_data
  dates <- as.Date(weekly$week_start)

  pars <- c("X", "S_prop", "immune_prop", "R0_t", "R_eff_t", "expected_reported_cases", "C_pred")
  draws <- rstan::extract(fit, pars = pars, permuted = TRUE)

  summaries <- dplyr::bind_rows(
    summarise_trajectory(draws$X, "latent_infections"),
    summarise_trajectory(draws$S_prop, "susceptible_prop"),
    summarise_trajectory(draws$immune_prop, "immune_prop"),
    summarise_trajectory(draws$R0_t, "R0_t"),
    summarise_trajectory(draws$R_eff_t, "R_eff_t"),
    summarise_trajectory(draws$expected_reported_cases, "expected_reported_cases"),
    summarise_trajectory(draws$C_pred, "posterior_predictive_cases")
  ) |> dplyr::mutate(week_start = dates[t])

  get_summary <- function(name) dplyr::filter(summaries, variable == name)
  date_breaks <- seq(as.Date("2015-01-01"), as.Date("2026-01-01"), by = "1 year")
  x_scale <- scale_x_date(breaks = date_breaks, date_labels = "%Y", expand = expansion(mult = c(0.01, 0.015)))

  expected <- get_summary("expected_reported_cases")
  predictive <- get_summary("posterior_predictive_cases")
  cases_df <- tibble::tibble(week_start = dates, observed = weekly$cases)

  p_cases <- ggplot() +
    geom_ribbon(data = predictive, aes(week_start, ymin = q025, ymax = q975), fill = COL_BLUE, alpha = 0.15) +
    geom_ribbon(data = expected, aes(week_start, ymin = q25, ymax = q75), fill = COL_NAVY, alpha = 0.24) +
    geom_line(data = expected, aes(week_start, median, colour = "Model expectation"), linewidth = 0.65) +
    geom_point(data = cases_df, aes(week_start, observed, colour = "Observed cases"), size = 0.5, alpha = 0.65) +
    scale_colour_manual(values = c("Observed cases" = COL_INK, "Model expectation" = COL_NAVY)) +
    scale_y_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = c(0, 0.06))) +
    x_scale +
    labs(title = "Reported cases and posterior prediction", subtitle = "No seed mechanism used",
         x = NULL, y = "Weekly reported cases") +
    theme_publication()

  latent <- get_summary("latent_infections")
  p_infections <- add_interval_layers(ggplot(latent, aes(week_start)), latent, COL_VERMILLION) +
    x_scale + scale_y_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = c(0, 0.06))) +
    labs(title = "Latent infection incidence", subtitle = "Median, 50% and 95% credible intervals",
         x = NULL, y = "True infections per week") +
    theme_publication()

  immunity <- dplyr::bind_rows(
    get_summary("susceptible_prop") |> mutate(state = "Susceptible"),
    get_summary("immune_prop") |> mutate(state = "Infection-derived immune (p_state)")
  )
  p_immunity <- ggplot(immunity, aes(week_start, colour = state, fill = state)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.10, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.7) +
    scale_colour_manual(values = c("Susceptible" = COL_GREEN, "Infection-derived immune (p_state)" = COL_PURPLE)) +
    scale_fill_manual(values = c("Susceptible" = COL_GREEN, "Infection-derived immune (p_state)" = COL_PURPLE)) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.01))) +
    x_scale +
    labs(title = "Population susceptibility and immunity", subtitle = "S[t]/U[t] unchanged from frozen v4.9 (actual counts)",
         x = NULL, y = "Population proportion") +
    theme_publication()

  reproduction <- dplyr::bind_rows(
    get_summary("R0_t") |> mutate(number = "R0(t)"),
    get_summary("R_eff_t") |> mutate(number = "Reff(t)")
  )
  p_reproduction <- ggplot(reproduction, aes(week_start, colour = number, fill = number)) +
    geom_hline(yintercept = 1, colour = COL_MUTED, linewidth = 0.35, linetype = "22") +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.09, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.65) +
    scale_colour_manual(values = c("R0(t)" = COL_NAVY, "Reff(t)" = COL_VERMILLION),
                         labels = c("R0(t)" = expression(R[0](t)), "Reff(t)" = expression(R[eff](t)))) +
    scale_fill_manual(values = c("R0(t)" = COL_NAVY, "Reff(t)" = COL_VERMILLION),
                       labels = c("R0(t)" = expression(R[0](t)), "Reff(t)" = expression(R[eff](t)))) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.06))) +
    x_scale +
    labs(title = "Time-varying reproduction numbers", subtitle = "Frozen v4.9 structure (unchanged)",
         x = NULL, y = "Reproduction number") +
    theme_publication()

  ascertainment_df <- tibble::tibble(week_start = dates, q = qval)
  p_ascertainment <- ggplot(ascertainment_df, aes(week_start, q)) +
    geom_line(colour = COL_ORANGE, linewidth = 0.9) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, max(0.35, qval * 1.3)), expand = expansion(mult = c(0, 0.07))) +
    x_scale +
    labs(title = "Case ascertainment", subtitle = sprintf("q FIXED at %.2f (not estimated)", qval),
         x = NULL, y = "Ascertainment fraction (q)") +
    theme_publication()

  # Panel F: state-level cumulative infection (immune fraction) trajectory --
  # same standalone view as the case-only Bahia six-panel's panel F. The
  # detailed 6-survey serology comparison (state vs geo-adjusted vs
  # observed, across q) is produced separately and more fully by
  # 09_bahia_multisite_serology_diagnostics.R
  # (bahia_multisite_serology_comparison_by_q.png), rather than crammed
  # into this shared grid.
  attack <- get_summary("immune_prop")
  p_attack <- ggplot(attack, aes(week_start)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = COL_GREEN, alpha = 0.14) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = COL_GREEN, alpha = 0.24) +
    geom_line(aes(y = median), colour = COL_GREEN, linewidth = 0.7) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA), expand = expansion(mult = c(0, 0.07))) +
    x_scale +
    labs(title = "Cumulative infection (state-level immune fraction, p_state)",
         subtitle = "See bahia_multisite_serology_comparison_by_q.png for the 6-site geo-adjusted comparison",
         x = NULL, y = "Cumulative proportion infected") +
    theme_publication()

  # Panel G: 6-site serology comparison (observed vs state vs geo-adjusted site), this q only.
  p_state_window <- rstan::extract(fit, "p_state_window")$p_state_window
  p_site_window <- rstan::extract(fit, "p_site_window")$p_site_window
  audit <- bundle$sero_audit
  sero_df <- bind_rows(lapply(seq_len(nrow(audit)), function(j) {
    tibble::tibble(
      sero_id = audit$sero_id[j],
      quantity = c("Observed", "State-level (no geo offset)", "Site-adjusted (with geo offset)"),
      median = c(audit$observed_prevalence[j], median(p_state_window[, j]), median(p_site_window[, j])),
      lo = c(NA, quantile(p_state_window[, j], .025), quantile(p_site_window[, j], .025)),
      hi = c(NA, quantile(p_state_window[, j], .975), quantile(p_site_window[, j], .975))
    )
  }))
  p_sero <- ggplot(sero_df, aes(sero_id, median, colour = quantity)) +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = .5), fatten = 2) +
    scale_colour_manual(values = c("Observed" = COL_INK, "State-level (no geo offset)" = COL_GREEN, "Site-adjusted (with geo offset)" = COL_PURPLE)) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    labs(title = "6 primary Bahia serosurveys: observed vs model (this q)", x = "survey", y = "seroprevalence") +
    theme_publication(base_size = 8) + theme(legend.position = "bottom", legend.text = element_text(size = 6.5))

  figure <- (p_cases | p_infections) / (p_immunity | p_reproduction) / (p_ascertainment | p_attack) / p_sero +
    patchwork::plot_annotation(
      title = "Chikungunya transmission and multi-site serology in Bahia, Brazil (2015-2025)",
      subtitle = sprintf("Frozen v4.9 + geo-adjusted multi-site serology, q = %.2f (%s, %s)",
                          qval, toupper(stage), if (bundle$hmc$hmc_pass) "HMC PASS" else "HMC FAIL/marginal"),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(size = 12.5, face = "bold", margin = margin(b = 3)),
        plot.subtitle = element_text(size = 9.5, colour = COL_MUTED, margin = margin(b = 8)),
        plot.tag = element_text(size = 11, face = "bold")
      )
    )

  figure_dir <- file.path(root, "03_Output/figures/renewal_bahia_v4_9_multisite_serology")
  table_dir <- file.path(root, "03_Output/tables/renewal_bahia_v4_9_multisite_serology")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  output_base <- paste0("bahia_multisite_trajectories_", tag, "_", stage)
  png_path <- file.path(figure_dir, paste0(output_base, ".png"))
  pdf_path <- file.path(figure_dir, paste0(output_base, ".pdf"))
  table_path <- file.path(table_dir, paste0("bahia_multisite_trajectory_summary_", tag, "_", stage, ".csv"))

  ggsave(png_path, figure, width = 183, height = 290, units = "mm", dpi = 300, bg = "white")
  ggsave(pdf_path, figure, width = 183, height = 290, units = "mm", device = cairo_pdf)
  utils::write.csv(summaries, table_path, row.names = FALSE)
  message("[", tag, "] saved: ", png_path)
  invisible(figure)
}

if (sys.nframe() == 0L) {
  qval <- as.numeric(Sys.getenv("QVAL", "0.05"))
  stage <- Sys.getenv("STAGE", "pilot")
  plot_bahia_multisite_six_panel(qval, stage)
}
