# Rio de Janeiro -- 6-panel trajectory figures for each COMPLETED targeted
# fixed-q fit, saved together under one fixedq/ figure subfolder. Reuses
# the plotting helpers from 07_plot_rj_six_panel.R (theme_publication(),
# summarise_trajectory(), add_interval_layers()) but writes to a shared
# subfolder with one filename per q, rather than one subfolder per tag.

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
source(file.path(root, "02_Script/05_rio_de_janeiro_pipeline/03_model_diagnostics/02_plot_rj_six_panel.R"))

base_dir <- file.path(root, "03_Output/model_fits/rio_de_janeiro/v4_9_replication")
fixedq_dir <- file.path(base_dir, "outputs/fixedq")
figure_dir <- file.path(root, "03_Output/figures/rio_de_janeiro_v4_9_replication/fixedq")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

plot_rj_fixedq_panel <- function(fit_path, q_target, show_u10 = FALSE) {
  bundle <- readRDS(fit_path)
  fit <- bundle$fit
  weekly <- bundle$weekly_data
  dates <- as.Date(weekly$week_start)

  pars <- c("X", "S_prop", "immune_prop", "R0_t", "R_eff_t", "expected_reported_cases", "C_pred")
  draws <- rstan::extract(fit, pars = pars, permuted = TRUE)
  q_draws <- as.vector(rstan::extract(fit, "q")$q)

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
    x_scale + labs(title = "Reported cases and posterior prediction", x = NULL, y = "Weekly reported cases") + theme_publication()

  latent <- get_summary("latent_infections")
  p_infections <- add_interval_layers(ggplot(latent, aes(week_start)), latent, COL_VERMILLION) +
    x_scale + scale_y_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = c(0, 0.06))) +
    labs(title = "Latent infection incidence", x = NULL, y = "True infections per week") + theme_publication()

  immunity <- dplyr::bind_rows(get_summary("susceptible_prop") |> mutate(state = "Susceptible"),
                                get_summary("immune_prop") |> mutate(state = "Infection-derived immune (p_state)"))
  p_immunity <- ggplot(immunity, aes(week_start, colour = state, fill = state)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.10, colour = NA) + geom_line(aes(y = median), linewidth = 0.7) +
    scale_colour_manual(values = c("Susceptible" = COL_GREEN, "Infection-derived immune (p_state)" = COL_PURPLE)) +
    scale_fill_manual(values = c("Susceptible" = COL_GREEN, "Infection-derived immune (p_state)" = COL_PURPLE)) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.01))) +
    x_scale + labs(title = "Population susceptibility and immunity", x = NULL, y = "Population proportion") + theme_publication()

  reproduction <- dplyr::bind_rows(get_summary("R0_t") |> mutate(number = "R0(t)"), get_summary("R_eff_t") |> mutate(number = "Reff(t)"))
  p_reproduction <- ggplot(reproduction, aes(week_start, colour = number, fill = number)) +
    geom_hline(yintercept = 1, colour = COL_MUTED, linewidth = 0.35, linetype = "22") +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.09, colour = NA) + geom_line(aes(y = median), linewidth = 0.65) +
    scale_colour_manual(values = c("R0(t)" = COL_NAVY, "Reff(t)" = COL_VERMILLION)) +
    scale_fill_manual(values = c("R0(t)" = COL_NAVY, "Reff(t)" = COL_VERMILLION)) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.06))) +
    x_scale + labs(title = "Time-varying reproduction numbers", x = NULL, y = "Reproduction number") + theme_publication()

  q_band <- tibble::tibble(week_start = dates, post_lo95 = quantile(q_draws, .025), post_hi95 = quantile(q_draws, .975),
                            post_lo50 = quantile(q_draws, .25), post_hi50 = quantile(q_draws, .75), post_median = median(q_draws))
  p_ascertainment <- ggplot(q_band, aes(week_start)) +
    geom_ribbon(aes(ymin = post_lo95, ymax = post_hi95), fill = COL_ORANGE, alpha = 0.18) +
    geom_ribbon(aes(ymin = post_lo50, ymax = post_hi50), fill = COL_ORANGE, alpha = 0.30) +
    geom_line(aes(y = post_median), colour = COL_ORANGE, linewidth = 0.9) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA), expand = expansion(mult = c(0, 0.07))) +
    x_scale + labs(title = "Case ascertainment (q FIXED, tight prior)", subtitle = sprintf("target q=%.4f (posterior median=%.4f)", q_target, median(q_draws)),
                    x = NULL, y = "Ascertainment fraction (q)") + theme_publication()

  attack <- get_summary("immune_prop")
  p_attack <- ggplot(attack, aes(week_start)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = COL_GREEN, alpha = 0.14) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = COL_GREEN, alpha = 0.24) + geom_line(aes(y = median), colour = COL_GREEN, linewidth = 0.7) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA), expand = expansion(mult = c(0, 0.07))) +
    x_scale + labs(title = "Cumulative infection (state-level immune fraction)", x = NULL, y = "Cumulative proportion infected") + theme_publication()

  if (show_u10) {
    # Observed U10 (Perisse 2020, Rio de Janeiro CITY) overlaid on the
    # STATEWIDE trajectory, same vertical-errorbar convention already
    # established for PE's U14 panel F (reported 95% CI, not raw n/N CI).
    # NOT a state=city equivalence claim -- see RJ_U10_CASEONLY_VALIDATION.md.
    u10_df <- tibble::tibble(
      week_start = as.Date("2018-07-01") + as.numeric(as.Date("2018-10-31") - as.Date("2018-07-01")) / 2,
      observed_prevalence = 0.180, observed_lo95 = 0.148, observed_hi95 = 0.212
    )
    p_attack <- p_attack +
      geom_errorbar(data = u10_df, aes(week_start, ymin = observed_lo95, ymax = observed_hi95),
                     width = 60, colour = COL_INK, linewidth = 0.4, inherit.aes = FALSE) +
      geom_point(data = u10_df, aes(week_start, observed_prevalence, shape = "Observed U10 (Rio de Janeiro CITY, weighted)"),
                 colour = COL_INK, size = 2.2, inherit.aes = FALSE) +
      scale_shape_manual(values = c("Observed U10 (Rio de Janeiro CITY, weighted)" = 18), name = NULL) +
      labs(subtitle = "Black diamond = observed U10 (Rio de Janeiro CITY, weighted 95% CI). NOT a state=city\nequivalence claim -- external plausibility check only (RJ_U10_CASEONLY_VALIDATION.md).") +
      theme(legend.position = "bottom")
  }

  hmc_txt <- if (bundle$hmc$hmc_pass) "HMC PASS" else sprintf("HMC FAIL (div=%d, td14=%d, maxRhat=%.3f)", bundle$hmc$divergences, bundle$hmc$max_treedepth_hits, bundle$hmc$maximum_rhat)
  figure <- (p_cases | p_infections) / (p_immunity | p_reproduction) / (p_ascertainment | p_attack) +
    patchwork::plot_annotation(title = sprintf("Rio de Janeiro: targeted fixed-q = %.4f", q_target),
                                subtitle = sprintf("Frozen v4.9, case-only, q pinned via tight prior (sd=0.02 on logit scale) -- %s", hmc_txt),
                                tag_levels = "A",
                                theme = theme(plot.title = element_text(size = 12.5, face = "bold", margin = margin(b = 3)),
                                              plot.subtitle = element_text(size = 9.5, colour = COL_MUTED, margin = margin(b = 8)),
                                              plot.tag = element_text(size = 11, face = "bold")))
  out_path <- file.path(figure_dir, sprintf("RJ_fixedq_q%.4f_panels.png", q_target))
  ggsave(out_path, figure, width = 183, height = 220, units = "mm", dpi = 300, bg = "white")
  message("[saved] ", out_path)
}

if (sys.nframe() == 0L) {
  rds_files <- list.files(fixedq_dir, pattern = "^rj_fixedq_q[0-9.]+\\.rds$", full.names = TRUE)
  message("Found ", length(rds_files), " completed fixed-q fit(s): ", paste(basename(rds_files), collapse = ", "))
  for (f in rds_files) {
    q_target <- as.numeric(sub("^rj_fixedq_q([0-9.]+)\\.rds$", "\\1", basename(f)))
    plot_rj_fixedq_panel(f, q_target)
  }
}
