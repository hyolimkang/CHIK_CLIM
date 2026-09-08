# =============================================================================
# 07_plot_renewal_v3_0_spatial_fit.R
#
# Publication-style figure for a v3.0 spatial renewal fit (debug pilot or a
# future all-Ceara fit), matching the v2.1/v2.2 six-panel layout as closely
# as v3.0's structure allows, plus a municipality-level susceptibility map
# with no v2.x equivalent, since v3.0's whole point is per-municipality S/U.
#
# Two structural differences from v2.1/v2.2 change what panels D and E show:
#   * v3.0 pools reporting into one scalar (p_symp * rho_sym_global), not a
#     time-varying trend, so panel E is a posterior density, not a trajectory.
#   * v3.0's R0(t) is a shared weekly random walk (mu_R) plus a *time-constant*
#     per-municipality offset, not a per-week municipality R0. Panel D shows
#     the shared mu_R trend and a population-weighted R_eff aggregate (built
#     from the per-municipality R_eff matrix already in the fit).
#
# Environment override (matches diagnose_renewal_v3_0_spatial.R):
#   RENEWAL_V3_0_FIT_PATH=path/to/fit_bundle.rds
# =============================================================================

required_packages <- c("here", "rstan", "ggplot2", "dplyr", "sf", "patchwork", "scales", "tibble")
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
  library(sf)
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
      plot.margin = margin(7, 8, 7, 7)
    )
}

posterior_interval <- function(draw_matrix) {
  stopifnot(length(dim(draw_matrix)) == 2L)
  tibble::tibble(
    t = seq_len(ncol(draw_matrix)),
    median = apply(draw_matrix, 2, stats::median),
    q025 = apply(draw_matrix, 2, stats::quantile, probs = 0.025),
    q975 = apply(draw_matrix, 2, stats::quantile, probs = 0.975)
  )
}

if (sys.nframe() == 0) {
  default_path <- here::here("02_Script/stan/renewal_ceara_v3_0_spatial_debug_pilot.rds")
  fit_path <- Sys.getenv("RENEWAL_V3_0_FIT_PATH", default_path)
  if (!file.exists(fit_path)) stop("v3.0 fit bundle not found: ", fit_path)

  bundle <- readRDS(fit_path)
  fit <- bundle$fit
  prepared <- bundle$prepared
  dates <- prepared$week_dates
  muni_ids <- prepared$municipality_ids
  N_start <- prepared$stan_data$N_start
  seed_weeks <- prepared$stan_data$seed_weeks
  is_debug_subset <- isTRUE(prepared$config$debug_subset)

  date_breaks <- seq(as.Date("2015-01-01"), as.Date("2020-01-01"), by = "1 year")
  x_scale <- scale_x_date(breaks = date_breaks, date_labels = "%Y", expand = expansion(mult = c(0.01, 0.015)))

  # ---- Panels A-C, F: Ceara-aggregate generated quantities (cheap) --------
  agg <- rstan::extract(
    fit,
    pars = c("C_ceara", "C_ceara_pred", "X_ceara_total", "S_ceara_prop", "immune_ceara_prop", "p_sero"),
    permuted = TRUE
  )
  observed_cases <- tibble::tibble(week_start = dates, observed = agg$C_ceara[1, ])

  cpred_summary <- posterior_interval(agg$C_ceara_pred) |> mutate(week_start = dates[t])
  p_cases <- ggplot() +
    geom_ribbon(data = cpred_summary, aes(week_start, ymin = q025, ymax = q975), fill = COL_BLUE, alpha = 0.2) +
    geom_line(data = cpred_summary, aes(week_start, median, colour = "Model expectation"), linewidth = 0.65) +
    geom_point(data = observed_cases, aes(week_start, observed, colour = "Observed cases"), size = 0.55, alpha = 0.72) +
    scale_colour_manual(values = c("Observed cases" = COL_INK, "Model expectation" = COL_NAVY)) +
    scale_y_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = c(0, 0.06))) +
    x_scale +
    labs(
      title = "Reported cases and posterior prediction",
      subtitle = if (is_debug_subset) sprintf("Debug pilot (%d/184 municipalities)", length(muni_ids)) else "All 184 municipalities",
      x = NULL, y = "Weekly reported cases"
    ) +
    theme_publication()

  latent_summary <- posterior_interval(agg$X_ceara_total) |> mutate(week_start = dates[t])
  p_infections <- ggplot(latent_summary, aes(week_start)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = COL_VERMILLION, alpha = 0.15) +
    geom_line(aes(y = median), colour = COL_VERMILLION, linewidth = 0.65) +
    scale_y_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = c(0, 0.06))) +
    x_scale +
    labs(
      title = "Latent infection incidence", subtitle = "Sum of municipality X[m,t]",
      x = NULL, y = "True infections per week"
    ) +
    theme_publication()

  immunity <- dplyr::bind_rows(
    posterior_interval(agg$S_ceara_prop) |> mutate(week_start = dates[t], state = "Susceptible"),
    posterior_interval(agg$immune_ceara_prop) |> mutate(week_start = dates[t], state = "Infection-derived immune")
  )
  p_immunity <- ggplot(immunity, aes(week_start, colour = state, fill = state)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.10, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.7) +
    scale_colour_manual(values = c("Susceptible" = COL_GREEN, "Infection-derived immune" = COL_PURPLE)) +
    scale_fill_manual(values = c("Susceptible" = COL_GREEN, "Infection-derived immune" = COL_PURPLE)) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.01))) +
    x_scale +
    labs(
      title = "Population susceptibility and immunity",
      subtitle = "Population-weighted S[m,t]/U[m,t]", x = NULL, y = "Population proportion"
    ) +
    theme_publication()

  # ---- Panel D: shared R0(t) random walk + population-weighted R_eff ------
  # R_eff[m,t] is an M x N matrix per draw -- kept as a separate, explicit
  # extraction since it is far larger than the Ceara-aggregate quantities.
  mu_r_draws <- rstan::extract(fit, pars = "mu_R", permuted = TRUE)$mu_R
  shared_r0 <- posterior_interval(exp(mu_r_draws)) |>
    mutate(week_start = dates[seed_weeks + t], number = "Shared R0(t)")

  r_eff_draws <- rstan::extract(fit, pars = "R_eff", permuted = TRUE)$R_eff  # [draws, M, N]
  weights <- N_start[, 1]  # municipality population share, fixed at series start
  r_eff_weighted <- apply(r_eff_draws, c(1, 3), function(x) stats::weighted.mean(x, w = weights))
  r_eff_summary <- posterior_interval(r_eff_weighted) |>
    mutate(week_start = dates[t], number = "Weighted Reff(t)") |>
    filter(t > seed_weeks)

  reproduction <- dplyr::bind_rows(shared_r0, r_eff_summary)
  p_reproduction <- ggplot(reproduction, aes(week_start, colour = number, fill = number)) +
    geom_hline(yintercept = 1, colour = COL_MUTED, linewidth = 0.35, linetype = "22") +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.12, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.65) +
    scale_colour_manual(values = c("Shared R0(t)" = COL_NAVY, "Weighted Reff(t)" = COL_VERMILLION)) +
    scale_fill_manual(values = c("Shared R0(t)" = COL_NAVY, "Weighted Reff(t)" = COL_VERMILLION)) +
    x_scale +
    labs(
      title = "Reproduction numbers",
      subtitle = "Municipality R0 offsets not shown",
      x = NULL, y = "Reproduction number"
    ) +
    theme_publication()

  # ---- Panel E: pooled reporting (no time trend in v3.0) -------------------
  reporting_draws <- rstan::extract(fit, pars = c("p_symp", "rho_sym_global", "overall_detection"), permuted = TRUE)
  reporting_labels <- c(p_symp = "p_symp", rho_sym_global = "rho_sym", overall_detection = "overall")
  reporting_df <- tibble::tibble(
    parameter = rep(unname(reporting_labels), each = length(reporting_draws$p_symp)),
    value = c(reporting_draws$p_symp, reporting_draws$rho_sym_global, reporting_draws$overall_detection)
  )
  p_reporting <- ggplot(reporting_df, aes(value, fill = parameter, colour = parameter)) +
    geom_density(alpha = 0.2, linewidth = 0.6) +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    scale_colour_manual(values = c(p_symp = COL_PURPLE, rho_sym = COL_ORANGE, overall = COL_INK)) +
    scale_fill_manual(values = c(p_symp = COL_PURPLE, rho_sym = COL_ORANGE, overall = COL_INK)) +
    labs(
      title = "Pooled case ascertainment",
      subtitle = "One scalar per parameter (no time trend)",
      x = "Probability", y = "Posterior density"
    ) +
    theme_publication()

  # ---- Panel F: serology anchors -------------------------------------------
  serology <- prepared$serology
  sero_observed <- dplyr::bind_rows(lapply(seq_len(nrow(serology)), function(j) {
    interval <- stats::binom.test(serology$positive[j], serology$n[j])$conf.int
    tibble::tibble(
      site = serology$site[j],
      week_start = serology$window_start[j] + (serology$window_end[j] - serology$window_start[j]) / 2,
      estimate = serology$positive[j] / serology$n[j], lower = interval[1], upper = interval[2]
    )
  }))
  sero_model <- dplyr::bind_rows(lapply(seq_len(nrow(serology)), function(j) {
    draws_j <- agg$p_sero[, j]
    tibble::tibble(
      site = serology$site[j],
      week_start = sero_observed$week_start[j] + 25,
      median = stats::median(draws_j), q025 = stats::quantile(draws_j, 0.025), q975 = stats::quantile(draws_j, 0.975)
    )
  }))
  immune_df <- posterior_interval(agg$immune_ceara_prop) |> mutate(week_start = dates[t])

  p_attack <- ggplot(immune_df, aes(week_start)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = COL_GREEN, alpha = 0.16) +
    geom_line(aes(y = median, colour = "Ceará (weighted)"), linewidth = 0.7) +
    geom_errorbar(data = sero_observed, aes(week_start, ymin = lower, ymax = upper, colour = "Site observed"), width = 18, linewidth = 0.55) +
    geom_point(data = sero_observed, aes(week_start, estimate, colour = "Site observed"), shape = 18, size = 2.3) +
    geom_errorbar(data = sero_model, aes(week_start, ymin = q025, ymax = q975, colour = "Site model (direct)"), width = 18, linewidth = 0.55) +
    geom_point(data = sero_model, aes(week_start, median, colour = "Site model (direct)"), shape = 16, size = 1.8) +
    scale_colour_manual(values = c(
      "Ceará (weighted)" = COL_GREEN, "Site observed" = COL_ORANGE, "Site model (direct)" = COL_PURPLE
    )) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA), expand = expansion(mult = c(0, 0.07))) +
    x_scale +
    labs(
      title = "Cumulative infection and serology",
      subtitle = "Direct municipality linkage (no geographic offset)",
      x = NULL, y = "Cumulative proportion infected"
    ) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_publication(base_size = 8) +
    theme(legend.text = element_text(size = 7.2))

  figure <- (p_cases | p_infections) /
    (p_immunity | p_reproduction) /
    (p_reporting | p_attack) +
    patchwork::plot_annotation(
      title = "Chikungunya transmission, municipality-level renewal model, Ceará, Brazil",
      subtitle = if (is_debug_subset) {
        sprintf("v3.0 spatial DEBUG PILOT (%d municipalities) -- computational smoke test only, not a scientific fit", length(muni_ids))
      } else {
        "Bayesian renewal model v3.0 (spatial), weekly observations from 2015 to 2019, all 184 Ceará municipalities"
      },
      tag_levels = "A",
      theme = theme(
        text = element_text(family = "Arial", colour = COL_INK),
        plot.title = element_text(size = 12.5, face = "bold", margin = margin(b = 3)),
        plot.subtitle = element_text(size = 9.5, colour = COL_MUTED, margin = margin(b = 8)),
        plot.tag = element_text(size = 11, face = "bold")
      )
    )

  # ---- New: municipality-level susceptibility map --------------------------
  # v3.0's whole contribution over v2.2 is per-municipality S/U; a state-
  # aggregate figure alone would hide it, so the final-week susceptible
  # proportion (and its posterior width) is mapped directly.
  s_prop_draws <- rstan::extract(fit, pars = "S_prop", permuted = TRUE)$S_prop  # [draws, M, N]
  final_week <- dim(s_prop_draws)[3]
  s_prop_final <- s_prop_draws[, , final_week]
  muni_summary <- tibble::tibble(
    muni6 = muni_ids,
    S_prop_median = apply(s_prop_final, 2, stats::median),
    S_prop_width_95 = apply(s_prop_final, 2, function(x) diff(stats::quantile(x, c(0.025, 0.975)))),
    total_cases = rowSums(prepared$stan_data$C)
  )

  polygons <- readRDS(here::here("01_Data/ibge_muni_polygons.rds"))
  map_data <- polygons |>
    dplyr::inner_join(muni_summary, by = "muni6")
  if (nrow(map_data) != nrow(muni_summary)) {
    warning(sprintf(
      "%d of %d fitted municipalities had no matching polygon and are omitted from the map",
      nrow(muni_summary) - nrow(map_data), nrow(muni_summary)
    ))
  }
  ceara_outline <- polygons |> dplyr::filter(substr(muni6, 1, 2) == "23")

  p_map_median <- ggplot() +
    geom_sf(data = ceara_outline, fill = "#F3F4F6", colour = "#D1D5DB", linewidth = 0.1) +
    geom_sf(data = map_data, aes(fill = S_prop_median), colour = "white", linewidth = 0.15) +
    scale_fill_viridis_c(option = "mako", direction = -1, labels = scales::label_percent(accuracy = 1), limits = c(0, 1), name = "Susceptible\n(median)") +
    labs(title = sprintf("Susceptible proportion, %s", format(dates[final_week]))) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = "right"
    )

  p_map_width <- ggplot() +
    geom_sf(data = ceara_outline, fill = "#F3F4F6", colour = "#D1D5DB", linewidth = 0.1) +
    geom_sf(data = map_data, aes(fill = S_prop_width_95), colour = "white", linewidth = 0.15) +
    scale_fill_viridis_c(option = "rocket", direction = -1, name = "95% CrI\nwidth") +
    labs(title = "Posterior uncertainty (95% CrI width)") +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = "right"
    )

  map_note <- if (is_debug_subset) {
    sprintf(
      "DEBUG PILOT: only %d of 184 municipalities were fitted; unfitted municipalities are shown as background outline only.",
      length(muni_ids)
    )
  } else {
    "All 184 Ceará municipalities fitted."
  }
  spatial_figure <- (p_map_median | p_map_width) +
    patchwork::plot_annotation(
      title = "Municipality-level susceptibility, v3.0 spatial renewal model",
      subtitle = map_note,
      theme = theme(
        plot.title = element_text(size = 12.5, face = "bold"),
        plot.subtitle = element_text(size = 9, colour = COL_MUTED)
      )
    )

  # ---- Save ------------------------------------------------------------------
  figure_dir <- here::here("03_Output/figures/renewal_v3_0")
  table_dir <- here::here("03_Output/tables/renewal_v3_0")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  suffix <- if (is_debug_subset) "_debug_pilot" else ""
  main_path <- file.path(figure_dir, paste0("renewal_v3_0_trajectories_ceara", suffix, ".png"))
  map_path <- file.path(figure_dir, paste0("renewal_v3_0_municipality_susceptibility", suffix, ".png"))
  table_path <- file.path(table_dir, paste0("renewal_v3_0_municipality_susceptibility_final_week", suffix, ".csv"))

  ggsave(main_path, figure, width = 183, height = 225, units = "mm", dpi = 300, bg = "white")
  ggsave(map_path, spatial_figure, width = 200, height = 120, units = "mm", dpi = 300, bg = "white")
  utils::write.csv(muni_summary, table_path, row.names = FALSE)

  message("[save] ", main_path)
  message("[save] ", map_path)
  message("[save] ", table_path)
}
