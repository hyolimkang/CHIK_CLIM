# Diagnostics for the v3.2 coarse-regional pilot. This script never samples.
# Epidemiological plots are withheld unless the prespecified HMC gate passes.
required_packages <- c("here", "rstan", "ggplot2", "dplyr", "patchwork", "scales", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(here); library(rstan); library(ggplot2); library(dplyr); library(patchwork)
})

COL_INK <- "#202124"; COL_MUTED <- "#6B7280"; COL_GRID <- "#E5E7EB"
COL_NAVY <- "#315A7D"; COL_BLUE <- "#56B4E9"; COL_VERMILLION <- "#D55E00"
COL_ORANGE <- "#E69F00"; COL_GREEN <- "#009E73"; COL_PURPLE <- "#76558F"

theme_publication <- function(base_size = 8.5) {
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, colour = COL_MUTED, margin = margin(b = 7)),
      axis.title = element_text(size = base_size, colour = COL_INK),
      axis.text = element_text(size = base_size - .5, colour = COL_MUTED),
      axis.line = element_line(colour = COL_INK, linewidth = .35),
      axis.ticks = element_line(colour = COL_INK, linewidth = .35),
      panel.grid.major.y = element_line(colour = COL_GRID, linewidth = .3),
      panel.grid.minor = element_blank(), legend.position = "top", legend.justification = "left",
      legend.title = element_blank(), legend.text = element_text(size = base_size - .5),
      legend.key.width = grid::unit(13, "pt"), legend.key.height = grid::unit(7, "pt"),
      plot.margin = margin(7, 8, 7, 7))
}

summarise_trajectory <- function(draw_matrix, variable) {
  stopifnot(length(dim(draw_matrix)) == 2L)
  tibble::tibble(variable = variable, t = seq_len(ncol(draw_matrix)),
    q025 = apply(draw_matrix, 2, stats::quantile, probs = .025),
    q25 = apply(draw_matrix, 2, stats::quantile, probs = .25),
    median = apply(draw_matrix, 2, stats::median),
    q75 = apply(draw_matrix, 2, stats::quantile, probs = .75),
    q975 = apply(draw_matrix, 2, stats::quantile, probs = .975))
}
draws_to_summary <- function(a, name, dates) summarise_trajectory(a, name) |> mutate(week_start = dates[t])
interval_layers <- function(plot, data, colour, fill = colour) {
  plot + geom_ribbon(data = data, aes(ymin = q025, ymax = q975), fill = fill, alpha = .13, colour = NA) +
    geom_ribbon(data = data, aes(ymin = q25, ymax = q75), fill = fill, alpha = .24, colour = NA) +
    geom_line(data = data, aes(y = median), colour = colour, linewidth = .65)
}

run_v3_2_diagnostics <- function() {
  fit_path <- Sys.getenv("RENEWAL_V3_2_FIT_PATH", here::here("02_Script/stan/renewal_ceara_v3_2_coarse_regional_pilot.rds"))
  if (!file.exists(fit_path)) stop("No v3.2 pilot bundle; diagnostics do not sample: ", fit_path)
  bundle <- readRDS(fit_path)
  if (is.null(bundle$pilot_diagnostics) || is.null(bundle$prepared)) stop("Fit bundle lacks pilot diagnostics or prepared data")
  h <- bundle$pilot_diagnostics; prepared <- bundle$prepared; stan_data <- prepared$stan_data
  dates <- as.Date(prepared$week_dates)
  allow_unconverged <- tolower(Sys.getenv("RENEWAL_V3_2_ALLOW_UNCONVERGED_PLOTS", "false")) %in% c("1", "true", "yes")
  table_dir <- here::here("03_Output/tables/renewal_v3_2_coarse_regional")
  figure_dir <- here::here("03_Output/figures/renewal_v3_2_coarse_regional")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE); dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  gate <- isTRUE(h$divergences == 0) && isTRUE(h$treedepth_hits == 0) && isTRUE(h$max_rhat <= 1.01) &&
    isTRUE(h$min_bulk_ess >= 100) && isTRUE(h$min_tail_ess >= 100) && all(h$bfmi >= .3)
  gate_summary <- data.frame(metric = c("divergences", "treedepth_hits", "max_rhat", "rhat_gt_1.01", "min_bulk_ess", "min_tail_ess", "elapsed_seconds", "gate_pass"),
    value = c(h$divergences, h$treedepth_hits, h$max_rhat, h$rhat_gt_1_01, h$min_bulk_ess, h$min_tail_ess, bundle$elapsed_seconds, gate))
  utils::write.csv(gate_summary, file.path(table_dir, "v3_2_hmc_gate.csv"), row.names = FALSE)
  utils::write.csv(h$convergence, file.path(table_dir, "v3_2_convergence.csv"), row.names = FALSE)
  worst <- head(h$convergence[order(h$convergence$rhat, decreasing = TRUE, na.last = NA), ], 10)
  utils::write.csv(worst, file.path(table_dir, "v3_2_worst_mixing.csv"), row.names = FALSE)

  # Always save chain traces, even when the computational gate rejects the pilot.
  trace_pars <- c("alpha_R", "sigma_R", "sigma_region", "q", "phi_obs", "mu_seed", "sigma_seed")
  trace_draws <- rstan::extract(bundle$fit, pars = trace_pars, permuted = FALSE)
  trace_path <- file.path(figure_dir, "renewal_v3_2_coarse_regional_traceplots.pdf")
  grDevices::pdf(trace_path, width = 11, height = 8.5); par(mfrow = c(4, 2))
  for (p in dimnames(trace_draws)[[3]]) matplot(trace_draws[, , p], type = "l", lty = 1,
    col = seq_len(dim(trace_draws)[2]), xlab = "Post-warmup iteration", ylab = p, main = p)
  grDevices::dev.off()
  hmc_pdf_path <- file.path(figure_dir, "renewal_v3_2_coarse_regional_hmc_gate.pdf")
  grDevices::pdf(hmc_pdf_path, width = 7.2, height = 5.1)
  plot.new()
  text(.05, .92, "v3.2 coarse-regional pilot: HMC gate", adj = c(0, 1), cex = 1.45, font = 2)
  text(.05, .83, if (gate) "PASS: epidemiological figures may be inspected" else "FAIL: do not interpret posterior trajectories", adj = c(0, 1), cex = 1.1,
       col = if (gate) COL_GREEN else COL_VERMILLION, font = 2)
  gate_lines <- c(
    sprintf("Elapsed time: %.1f minutes", bundle$elapsed_seconds / 60),
    sprintf("Divergences: %d (threshold: 0)", h$divergences),
    sprintf("Maximum-treedepth hits: %d (threshold: 0)", h$treedepth_hits),
    sprintf("Maximum R-hat: %.3f (threshold: <= 1.01)", h$max_rhat),
    sprintf("Minimum bulk ESS: %.1f (threshold: >= 100)", h$min_bulk_ess),
    sprintf("Minimum tail ESS: %.1f (threshold: >= 100)", h$min_tail_ess),
    sprintf("BFMI by chain: %s (threshold: >= 0.30)", paste(sprintf("%.3f", h$bfmi), collapse = ", ")),
    "Interpretation: treedepth saturation and non-mixing invalidate this pilot."
  )
  text(.05, .69, paste(gate_lines, collapse = "\n\n"), adj = c(0, 1), cex = .95)
  grDevices::dev.off()
  if (!gate && !allow_unconverged) {
    message("[gate] FAILED; epidemiological trajectory figures intentionally not generated")
    return(invisible(list(gate = gate, gate_summary = gate_summary, trace_path = trace_path, hmc_pdf_path = hmc_pdf_path)))
  }
  if (!gate) message("[gate] FAILED; rendering explicitly labelled exploratory non-converged trajectories")

  draws <- rstan::extract(bundle$fit, pars = c("C_ceara", "C_ceara_pred", "X_ceara", "S_ceara_prop", "immune_ceara_prop", "mu_R", "R_eff_ceara", "q", "p_sero", "C_pred", "S_prop"), permuted = TRUE)
  ceara_summary <- bind_rows(
    draws_to_summary(draws$X_ceara, "latent_infections", dates),
    draws_to_summary(draws$S_ceara_prop, "susceptible_prop", dates),
    draws_to_summary(draws$immune_ceara_prop, "immune_prop", dates),
    draws_to_summary(exp(draws$mu_R), "shared_R0", dates),
    draws_to_summary(draws$R_eff_ceara, "Reff_ceara", dates))
  get_summary <- function(name) filter(ceara_summary, variable == name)
  x_scale <- scale_x_date(breaks = seq(as.Date("2015-01-01"), as.Date("2020-01-01"), by = "1 year"), date_labels = "%Y", expand = expansion(mult = c(.01, .015)))

  observed_cases <- tibble(week_start = dates, cases = draws$C_ceara[1, ])
  expected_cases <- draws_to_summary(apply(draws$C_pred, c(1, 3), sum), "expected_cases", dates)
  predictive_cases <- draws_to_summary(draws$C_ceara_pred, "predictive_cases", dates)
  p_cases <- ggplot() +
    geom_ribbon(data = predictive_cases, aes(week_start, ymin = q025, ymax = q975), fill = COL_BLUE, alpha = .15) +
    geom_ribbon(data = expected_cases, aes(week_start, ymin = q25, ymax = q75), fill = COL_NAVY, alpha = .24) +
    geom_line(data = expected_cases, aes(week_start, median, colour = "Model expectation"), linewidth = .65) +
    geom_point(data = observed_cases, aes(week_start, cases, colour = "Observed cases"), size = .55, alpha = .72) +
    scale_colour_manual(values = c("Observed cases" = COL_INK, "Model expectation" = COL_NAVY)) +
    scale_y_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = c(0, .06))) + x_scale +
    labs(title = "Reported cases and posterior prediction", subtitle = "Shading: 95% predictive interval; dark band: 50% CrI", x = NULL, y = "Weekly reported cases") + theme_publication()

  latent <- get_summary("latent_infections")
  p_latent <- interval_layers(ggplot(latent, aes(week_start)), latent, COL_VERMILLION) + x_scale +
    scale_y_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = c(0, .06))) +
    labs(title = "Latent infection incidence", subtitle = "Ceará total: median, 50% and 95% credible intervals", x = NULL, y = "True infections per week") + theme_publication()

  immunity <- bind_rows(get_summary("susceptible_prop") |> mutate(state = "Susceptible"), get_summary("immune_prop") |> mutate(state = "Infection-derived immune"))
  p_immunity <- ggplot(immunity, aes(week_start, colour = state, fill = state)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .10, colour = NA) + geom_line(aes(y = median), linewidth = .7) +
    scale_colour_manual(values = c("Susceptible" = COL_GREEN, "Infection-derived immune" = COL_PURPLE)) +
    scale_fill_manual(values = c("Susceptible" = COL_GREEN, "Infection-derived immune" = COL_PURPLE)) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1), breaks = seq(0, 1, .25), expand = expansion(mult = c(0, .01))) + x_scale +
    labs(title = "Population susceptibility and immunity", subtitle = "Population-weighted Ceará total", x = NULL, y = "Population proportion") + theme_publication()

  reproduction <- bind_rows(get_summary("shared_R0") |> filter(t > stan_data$seed_weeks) |> mutate(number = "Shared R0(t)"), get_summary("Reff_ceara") |> filter(t > stan_data$seed_weeks) |> mutate(number = "Ceará Reff(t)"))
  p_reproduction <- ggplot(reproduction, aes(week_start, colour = number, fill = number)) +
    geom_hline(yintercept = 1, colour = COL_MUTED, linewidth = .35, linetype = "22") + geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .09, colour = NA) + geom_line(aes(y = median), linewidth = .65) +
    scale_colour_manual(values = c("Shared R0(t)" = COL_NAVY, "Ceará Reff(t)" = COL_VERMILLION)) + scale_fill_manual(values = c("Shared R0(t)" = COL_NAVY, "Ceará Reff(t)" = COL_VERMILLION)) +
    scale_y_continuous(expand = expansion(mult = c(.03, .06))) + x_scale +
    labs(title = "Transmission potential", subtitle = "Shared R0 trajectory and population-weighted Ceará Reff", x = NULL, y = "Reproduction number") + theme_publication()

  q_ci <- stats::quantile(draws$q, c(.025, .5, .975))
  p_reporting <- ggplot(tibble(q = draws$q), aes(y = q)) + geom_boxplot(fill = COL_ORANGE, alpha = .65, width = .25, outlier.alpha = .12) +
    scale_y_continuous(labels = scales::label_percent(accuracy = .1), limits = c(0, NA)) +
    labs(title = "Pooled infection-to-report probability", subtitle = sprintf("q = %.3f [%.3f, %.3f]", q_ci[2], q_ci[1], q_ci[3]), x = NULL, y = "Probability") + theme_publication()

  serology <- prepared$serology
  sero_obs <- tibble(site = serology$site, observed = stan_data$sero_positive / stan_data$sero_n,
    lower = vapply(seq_len(stan_data$J), function(j) stats::binom.test(stan_data$sero_positive[j], stan_data$sero_n[j])$conf.int[1], numeric(1)),
    upper = vapply(seq_len(stan_data$J), function(j) stats::binom.test(stan_data$sero_positive[j], stan_data$sero_n[j])$conf.int[2], numeric(1)))
  sero_model <- bind_rows(lapply(seq_len(stan_data$J), function(j) tibble(site = serology$site[j], median = median(draws$p_sero[, j]), q025 = quantile(draws$p_sero[, j], .025), q975 = quantile(draws$p_sero[, j], .975))))
  sero_plot <- left_join(sero_obs, sero_model, by = "site")
  p_serology <- ggplot(sero_plot, aes(site, median)) +
    geom_errorbar(aes(ymin = q025, ymax = q975, colour = "Model site seroprevalence"), width = .15, linewidth = .6) + geom_point(aes(colour = "Model site seroprevalence"), size = 2) +
    geom_errorbar(aes(ymin = lower, ymax = upper, colour = "Observed seroprevalence"), width = .15, linewidth = .6, position = position_nudge(x = .12)) + geom_point(aes(y = observed, colour = "Observed seroprevalence"), size = 2, shape = 18, position = position_nudge(x = .12)) +
    scale_colour_manual(values = c("Model site seroprevalence" = COL_PURPLE, "Observed seroprevalence" = COL_ORANGE)) + scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
    labs(title = "Direct local serology validation", subtitle = "Model: regional survey-window mean; bars: 95% intervals", x = NULL, y = "Cumulative proportion infected") + theme_publication()

  figure_subtitle <- if (gate) {
    "v3.2 coarse-regional renewal pilot: statewide population-weighted trajectory"
  } else {
    "EXPLORATORY ONLY — UNCONVERGED CHAINS; DO NOT INTERPRET AS POSTERIOR INFERENCE"
  }
  ceara_page <- (p_cases | p_latent) / (p_immunity | p_reproduction) / (p_reporting | p_serology) +
    plot_annotation(title = "Chikungunya transmission and susceptibility in Ceará, Brazil", subtitle = figure_subtitle, tag_levels = "A", theme = theme(text = element_text(family = "Arial", colour = COL_INK), plot.title = element_text(size = 12.5, face = "bold"), plot.subtitle = element_text(size = 9.5, colour = if (gate) COL_MUTED else COL_VERMILLION, face = if (gate) "plain" else "bold"), plot.tag = element_text(size = 11, face = "bold")))

  # S(t)/N(t) for every fitted spatial unit, on a common 0--1 scale.
  regional_s <- bind_rows(lapply(seq_len(stan_data$K), function(k) {
    draws_to_summary(draws$S_prop[, k, ], "susceptible_prop", dates) |> mutate(unit = prepared$unit_labels[k])
  }))
  regional_page <- ggplot(regional_s, aes(week_start)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = COL_GREEN, alpha = .13, colour = NA) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = COL_GREEN, alpha = .24, colour = NA) + geom_line(aes(y = median), colour = COL_GREEN, linewidth = .6) +
    facet_wrap(~unit, ncol = 2) + scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1), breaks = seq(0, 1, .25), expand = expansion(mult = c(0, .01))) + x_scale +
    labs(title = "Sub-spatial susceptible trajectories", subtitle = "S(t)/N(t) at the start of each week; median, 50% and 95% credible intervals", x = NULL, y = "Susceptible proportion") +
    theme_publication() + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", colour = COL_INK), legend.position = "none")

  output_suffix <- if (gate) "" else "_UNCONVERGED_DO_NOT_INTERPRET"
  summary_path <- file.path(table_dir, paste0("renewal_v3_2_coarse_regional_trajectory_summary", output_suffix, ".csv"))
  utils::write.csv(bind_rows(ceara_summary, regional_s), summary_path, row.names = FALSE)
  pdf_path <- file.path(figure_dir, paste0("renewal_v3_2_coarse_regional_diagnostics", output_suffix, ".pdf"))
  grDevices::cairo_pdf(pdf_path, width = 7.2, height = 8.86, onefile = TRUE)
  print(ceara_page); print(regional_page); grDevices::dev.off()
  message("[save] ", pdf_path); message("[save] ", summary_path)
  invisible(list(gate = gate, gate_summary = gate_summary, trace_path = trace_path, pdf_path = pdf_path, summary_path = summary_path))
}
if (sys.nframe() == 0) run_v3_2_diagnostics()
