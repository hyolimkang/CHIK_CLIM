# Diagnostics and outputs for the Ceara episode-sequential susceptibility pilot.
#
# This script is deliberately read-only with respect to the Stan model and raw
# data. It labels the filtered reconstruction as primary and the full-data
# smoothed reconstruction as retrospective sensitivity only.

required_packages <- c("here", "rstan", "posterior", "dplyr", "tibble", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(rstan)
  library(posterior)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

project_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM R-project root.")
}

project_path <- function(...) file.path(project_root(), ...)

output_paths <- function() {
  list(
    table = project_path("03_Output", "tables", "episode_susceptibility", "ceara_pilot_v1"),
    figure = project_path("03_Output", "figures", "episode_susceptibility", "ceara_pilot_v1"),
    fit = project_path("02_Script", "stan", "ceara_episode_susceptibility_pilot_v1_fit.rds"),
    documentation = normalizePath(file.path(project_root(), "..", "docs"), mustWork = TRUE)
  )
}

publication_theme <- function() {
  theme_classic(base_size = 8.5) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90"),
      strip.background = element_rect(fill = "grey95", colour = NA)
    )
}

interval_summary <- function(x) {
  c(q025 = quantile(x, 0.025), q25 = quantile(x, 0.25), median = median(x), q75 = quantile(x, 0.75), q975 = quantile(x, 0.975))
}

summarise_hmc <- function(fit, label, elapsed_seconds) {
  sampler_draws <- posterior::as_draws_array(rstan::extract(
    fit,
    pars = c("log_lambda", "alpha_q", "beta_q", "log_phi"),
    permuted = FALSE,
    inc_warmup = FALSE
  ))
  parameter_summary <- as.data.frame(posterior::summarise_draws(
    sampler_draws,
    posterior::rhat,
    posterior::ess_bulk,
    posterior::ess_tail
  ))
  names(parameter_summary) <- sub("^posterior::", "", names(parameter_summary))
  bfmi <- rstan::get_bfmi(fit)
  divergences <- rstan::get_num_divergent(fit)
  max_treedepth_hits <- rstan::get_num_max_treedepth(fit)
  maximum_rhat <- max(parameter_summary$rhat, na.rm = TRUE)
  n_rhat_gt_1_01 <- sum(parameter_summary$rhat > 1.01, na.rm = TRUE)
  minimum_bulk_ess <- min(parameter_summary$ess_bulk, na.rm = TRUE)
  minimum_tail_ess <- min(parameter_summary$ess_tail, na.rm = TRUE)
  minimum_bfmi <- min(bfmi)
  tibble(
    fit = label,
    elapsed_seconds = elapsed_seconds,
    divergences = divergences,
    max_treedepth_hits = max_treedepth_hits,
    maximum_rhat = maximum_rhat,
    n_rhat_gt_1_01 = n_rhat_gt_1_01,
    minimum_bulk_ess = minimum_bulk_ess,
    minimum_tail_ess = minimum_tail_ess,
    minimum_bfmi = minimum_bfmi,
    bfmi_by_chain = paste(round(bfmi, 3), collapse = "; "),
    gate_pass =
      divergences == 0 &&
      max_treedepth_hits == 0 &&
      maximum_rhat <= 1.01 &&
      minimum_bulk_ess >= 100 &&
      minimum_tail_ess >= 100 &&
      minimum_bfmi >= 0.3
  )
}

extract_period_draws <- function(fit_result, reconstruction, target_episode = NULL) {
  variables <- c(
    "lambda", "attack_prob", "latent_infections", "S_start", "S_end",
    "q", "C_rep", "S_prop_start", "S_prop_end", "accounting_error"
  )
  draws <- rstan::extract(fit_result$fit, pars = variables)
  periods <- fit_result$periods
  n_draw <- nrow(draws$lambda)

  bind_rows(lapply(seq_len(nrow(periods)), function(period_index) {
    tibble(
      draw = seq_len(n_draw),
      state = periods$state[period_index],
      reconstruction = reconstruction,
      target_episode = if (is.null(target_episode)) NA_character_ else target_episode,
      period_id = periods$period_id[period_index],
      period_type = periods$period_type[period_index],
      episode_id = periods$episode_id[period_index],
      period_start = periods$period_start[period_index],
      period_end = periods$period_end[period_index],
      reported_cases = periods$reported_cases_period[period_index],
      N_start = periods$N_start[period_index],
      N_end = periods$N_end[period_index],
      lambda = draws$lambda[, period_index],
      attack_prob = draws$attack_prob[, period_index],
      latent_infections = draws$latent_infections[, period_index],
      S_start = draws$S_start[, period_index],
      S_end = draws$S_end[, period_index],
      S_prop_start = draws$S_prop_start[, period_index],
      S_prop_end = draws$S_prop_end[, period_index],
      q = draws$q[, period_index],
      C_rep = draws$C_rep[, period_index],
      accounting_error = draws$accounting_error[, period_index]
    )
  }))
}

summarise_period_draws <- function(period_draws) {
  numeric_columns <- c(
    "lambda", "attack_prob", "latent_infections", "S_start", "S_end",
    "S_prop_start", "S_prop_end", "q", "C_rep", "accounting_error"
  )
  period_draws |>
    group_by(
      state, reconstruction, period_id, period_type, episode_id, period_start,
      period_end, reported_cases, N_start, N_end
    ) |>
    summarise(
      across(
        all_of(numeric_columns),
        list(q025 = ~ quantile(.x, 0.025), q25 = ~ quantile(.x, 0.25), median = median, q75 = ~ quantile(.x, 0.75), q975 = ~ quantile(.x, 0.975), max_abs = ~ max(abs(.x))),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
}

episode_summary <- function(period_draws, reconstruction, episode_lookup, filtered = FALSE) {
  if (filtered) {
    # The terminal period in each cut fit ends one week before the target onset.
    onset_draws <- period_draws |>
      group_by(target_episode, draw) |>
      filter(period_end == max(period_end)) |>
      ungroup() |>
      transmute(
        episode_id = target_episode,
        draw,
        S_at_onset = S_end,
        S_prop_at_onset = S_prop_end
      )
  } else {
    onset_draws <- period_draws |>
      filter(period_type == "epidemic") |>
      transmute(
        episode_id,
        draw,
        S_at_onset = S_start,
        S_prop_at_onset = S_prop_start
      )
  }

  end_draws <- if (filtered) {
    # Filtered analyses deliberately do not condition on the target episode,
    # so an episode-end posterior is unavailable by design.
    tibble(episode_id = character(), draw = integer(), S_prop_end = numeric())
  } else {
    period_draws |>
      filter(period_type == "epidemic") |>
      transmute(episode_id, draw, S_prop_end)
  }

  onset_summary <- onset_draws |>
    group_by(episode_id) |>
    summarise(
      across(
        c(S_at_onset, S_prop_at_onset),
        list(q025 = ~ quantile(.x, 0.025), q25 = ~ quantile(.x, 0.25), median = median, q75 = ~ quantile(.x, 0.75), q975 = ~ quantile(.x, 0.975)),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
  if (nrow(end_draws) > 0L) {
    onset_summary <- onset_summary |>
      left_join(
        end_draws |>
          group_by(episode_id) |>
          summarise(
            across(S_prop_end, list(q025 = ~ quantile(.x, .025), q25 = ~ quantile(.x, .25), median = median, q75 = ~ quantile(.x, .75), q975 = ~ quantile(.x, .975)), .names = "{.col}_{.fn}"),
            .groups = "drop"
          ),
        by = "episode_id"
      )
  }
  episode_lookup |>
    left_join(onset_summary, by = "episode_id") |>
    mutate(reconstruction = reconstruction)
}

make_pilot_episode_table <- function(smoothed_draws, episode_lookup) {
  smoothed_draws |>
    filter(period_type == "epidemic") |>
    group_by(episode_id) |>
    summarise(
      across(
        c(latent_infections, S_prop_start, S_prop_end),
        list(q025 = ~ quantile(.x, .025), median = median, q975 = ~ quantile(.x, .975)),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    ) |>
    left_join(episode_lookup, by = "episode_id") |>
    transmute(
      episode_id,
      onset_date = episode_start,
      episode_end,
      reported_cases = reported_cases_episode_audit,
      latent_infections_q025,
      latent_infections_median,
      latent_infections_q975,
      S_prop_at_onset_q025 = S_prop_start_q025,
      S_prop_at_onset_median = S_prop_start_median,
      S_prop_at_onset_q975 = S_prop_start_q975,
      S_prop_end_q025,
      S_prop_end_median,
      S_prop_end_q975
    )
}

make_prior_draws <- function(periods, priors, n = 10000L) {
  set.seed(12012026L)
  bind_rows(lapply(seq_len(nrow(periods)), function(index) {
    lambda <- rlnorm(
      n,
      meanlog = if (periods$period_type[index] == "epidemic") priors$log_lambda_epi_prior_mean else priors$log_lambda_inter_prior_mean,
      sdlog = if (periods$period_type[index] == "epidemic") priors$log_lambda_epi_prior_sd else priors$log_lambda_inter_prior_sd
    )
    tibble(
      period_id = periods$period_id[index],
      period_type = periods$period_type[index],
      lambda = lambda,
      attack_prob = 1 - exp(-lambda),
      q = plogis(rnorm(n, priors$q_alpha_prior_mean + priors$q_beta_prior_mean * periods$year_centered[index], priors$q_alpha_prior_sd))
    )
  }))
}

write_pilot_documentation <- function(paths, bundle, hmc, smoothed_draws, filtered_episode, smoothed_episode) {
  if (!all(hmc$gate_pass)) return(FALSE)

  total_infections <- smoothed_draws |>
    group_by(draw) |>
    summarise(total = sum(latent_infections), .groups = "drop") |>
    pull(total)
  foi_path <- project_path("03_Output", "tables", "ce_foi_summary_brazil_ceara.csv")
  burden_path <- project_path("03_Output", "tables", "ce_burden_infection_summary.csv")
  foi_note <- if (file.exists(foi_path)) paste(readLines(foi_path, n = 2L), collapse = " / ") else "Not available."
  burden_note <- if (file.exists(burden_path)) paste(readLines(burden_path, n = 2L), collapse = " / ") else "Not available."
  v22_short_fit <- project_path("02_Script", "stan", "renewal_ceara_v2_2_fit.rds")

  document <- c(
    "# Ceara episode-sequential susceptibility pilot",
    "",
    "## Status",
    "",
    "The computational diagnostic gate passed for the smoothed fit and every filtered cut fit. These results are a period-level susceptibility reconstruction, not a weekly transmission model. No weekly R0, Reff, introductions, or climate covariates were estimated.",
    "",
    "## Data and episode definitions",
    "",
    paste0("- Surveillance definition: ", bundle$prepared$config$case_definition),
    "- Analysis period: 2015-01-04 through 2024-12-29.",
    "- Episodes: seven existing Ceara audit-defined episodes. The timeline is fully partitioned into 15 non-overlapping epidemic/inter-episode periods.",
    "- Population: period sums of the corrected weekly demographic series; births enter S and deaths/reconciliation are allocated proportionally after infections, matching v2.2.",
    "- Serology: the two available surveys are municipality-scale, not defensible as state-representative primary likelihood inputs. They are stored as metadata and excluded from this primary pilot.",
    "",
    "## Identification and priors",
    "",
    "Surveillance alone identifies the product of detection and infections much more strongly than either component alone. Therefore the detection process is a strongly regularized state-level intercept plus linear calendar-time trend in this one-state pilot. Results should be read with the saved prior-vs-posterior plots and not as independently serology-calibrated infection totals.",
    "",
    sprintf("The smoothed full-timeline posterior cumulative infections have median %.0f (95%% CrI %.0f to %.0f).", median(total_infections), quantile(total_infections, .025), quantile(total_infections, .975)),
    "",
    "## Filtered versus smoothed",
    "",
    "Filtered S_at_onset is the primary recurrence quantity: each episode uses a separate fit ending immediately before its onset. Smoothed S_at_onset conditions on the full history and is supplied only as a retrospective sensitivity. The filtered output intentionally has no posterior S at the end of the target episode, because doing so would require target-episode outcomes.",
    "",
    "## Comparison with existing work",
    "",
    "Existing FOI/burden files provide cross-sectional overall infection summaries rather than an annual/episode-resolved susceptibility trajectory; they are useful only for broad magnitude comparison.",
    paste0("- Existing FOI file preview: ", foi_note),
    paste0("- Existing burden file preview: ", burden_note),
    if (file.exists(v22_short_fit)) "- A v2.2 2015-2019 posterior fit was found and should be compared externally." else "- No valid saved v2.2 2015-2019 posterior fit was found. The available v2.2 2015-2025 stress fit failed its HMC gate, so it is not used for posterior trajectory comparison.",
    "",
    "Episode-level timing can differ from annual susceptibility summaries because this model updates S after every epidemic and inter-episode period rather than averaging the whole year.",
    "",
    "## Outputs",
    "",
    "- `filtered_episode_susceptibility.csv`: primary onset susceptibility summaries.",
    "- `smoothed_episode_susceptibility.csv`: full-history sensitivity summaries.",
    "- `smoothed_state_period_draws.rds` and `filtered_state_period_draws.rds`: every posterior draw by state-period.",
    "- Figures and diagnostic CSV files are in the matching `ceara_pilot_v1` output folders."
  )
  writeLines(document, file.path(paths$documentation, "episode_susceptibility_pilot_results.md"))
  TRUE
}

run_ceara_episode_susceptibility_diagnostics <- function() {
  paths <- output_paths()
  if (!file.exists(paths$fit)) stop("Fit bundle not found: ", paths$fit)
  dir.create(paths$table, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$figure, recursive = TRUE, showWarnings = FALSE)
  bundle <- readRDS(paths$fit)

  smoothed_draws <- extract_period_draws(bundle$smoothed, "smoothed")
  filtered_draws <- bind_rows(lapply(bundle$filtered, function(x) {
    extract_period_draws(x, "filtered", x$target_episode$episode_id)
  }))
  smoothed_summary <- summarise_period_draws(smoothed_draws)
  filtered_summary <- summarise_period_draws(filtered_draws)
  filtered_episode <- episode_summary(filtered_draws, "filtered", bundle$prepared$episodes, filtered = TRUE)
  smoothed_episode <- episode_summary(smoothed_draws, "smoothed", bundle$prepared$episodes)
  pilot_episode <- make_pilot_episode_table(smoothed_draws, bundle$prepared$episodes)

  hmc <- bind_rows(
    summarise_hmc(bundle$smoothed$fit, "smoothed", bundle$smoothed$elapsed_seconds),
    bind_rows(lapply(bundle$filtered, function(x) {
      summarise_hmc(x$fit, paste0("filtered_", x$target_episode$episode_id), x$elapsed_seconds)
    }))
  )
  hmc_gate_label <- if (all(hmc$gate_pass)) "HMC diagnostic gate passed" else "UNCONVERGED: do not interpret posterior reconstruction"

  saveRDS(smoothed_draws, file.path(paths$table, "smoothed_state_period_draws.rds"))
  saveRDS(filtered_draws, file.path(paths$table, "filtered_state_period_draws.rds"))
  write.csv(smoothed_summary, file.path(paths$table, "smoothed_period_susceptibility_summary.csv"), row.names = FALSE)
  write.csv(filtered_summary, file.path(paths$table, "filtered_period_susceptibility_summary.csv"), row.names = FALSE)
  write.csv(filtered_episode, file.path(paths$table, "filtered_episode_susceptibility.csv"), row.names = FALSE)
  write.csv(smoothed_episode, file.path(paths$table, "smoothed_episode_susceptibility.csv"), row.names = FALSE)
  write.csv(pilot_episode, file.path(paths$table, "ceara_episode_pilot_summary.csv"), row.names = FALSE)
  write.csv(hmc, file.path(paths$table, "hmc_diagnostics.csv"), row.names = FALSE)

  ppc <- smoothed_draws |>
    group_by(period_id, period_type, episode_id, period_start, period_end, reported_cases) |>
    summarise(
      C_rep_q025 = quantile(C_rep, .025), C_rep_median = median(C_rep), C_rep_q975 = quantile(C_rep, .975),
      .groups = "drop"
    )
  write.csv(ppc, file.path(paths$table, "smoothed_period_ppc.csv"), row.names = FALSE)

  timeline_s <- smoothed_summary |>
    transmute(date = period_start, S_prop_q025 = S_prop_start_q025, S_prop_median = S_prop_start_median, S_prop_q975 = S_prop_start_q975) |>
    bind_rows(smoothed_summary |>
      transmute(date = period_end + 7, S_prop_q025 = S_prop_end_q025, S_prop_median = S_prop_end_median, S_prop_q975 = S_prop_end_q975)) |>
    arrange(date)
  onset_points <- filtered_episode |>
    transmute(episode_id, date = episode_start, S_prop_q025 = S_prop_at_onset_q025, S_prop_median = S_prop_at_onset_median, S_prop_q975 = S_prop_at_onset_q975)

  case_plot <- ggplot(bundle$prepared$weekly_cases, aes(week_start, cases)) +
    geom_rect(
      data = bundle$prepared$episodes,
      aes(xmin = episode_start, xmax = episode_end + 7, ymin = 0, ymax = Inf),
      inherit.aes = FALSE, fill = "#E69F00", alpha = .18
    ) +
    geom_line(colour = "#D55E00", linewidth = .35) +
    scale_y_continuous(trans = scales::pseudo_log_trans(10), breaks = c(0, 1, 10, 100, 1000, 10000), labels = scales::label_number(big.mark = ",")) +
    labs(title = "Ceara reported chikungunya cases", subtitle = "Gold bands are audit-defined epidemic periods", x = NULL, y = "Weekly notified cases") +
    publication_theme()
  susceptibility_plot <- ggplot(timeline_s, aes(date, S_prop_median)) +
    geom_ribbon(aes(ymin = S_prop_q025, ymax = S_prop_q975), fill = "#56B4E9", alpha = .25) +
    geom_step(colour = "#0072B2", linewidth = .6) +
    geom_errorbar(data = onset_points, aes(x = date, ymin = S_prop_q025, ymax = S_prop_q975), inherit.aes = FALSE, width = 12, colour = "#009E73") +
    geom_point(data = onset_points, aes(x = date, y = S_prop_median), inherit.aes = FALSE, colour = "#009E73", size = 1.7) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
    labs(title = "Sequential susceptible fraction", subtitle = "Blue: smoothed full-history reconstruction; green: primary filtered onset estimates", x = NULL, y = "S / N") +
    publication_theme()
  ggsave(file.path(paths$figure, "ceara_episode_susceptibility_timeline.pdf"), case_plot / susceptibility_plot + plot_annotation(title = hmc_gate_label), width = 183, height = 180, units = "mm", device = cairo_pdf)

  ppc_plot <- ggplot(ppc, aes(reported_cases, C_rep_median, colour = period_type)) +
    geom_errorbar(aes(ymin = C_rep_q025, ymax = C_rep_q975), alpha = .65) +
    geom_point(size = 1.7) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
    scale_x_continuous(trans = scales::pseudo_log_trans(10), labels = scales::label_number(big.mark = ",")) +
    scale_y_continuous(trans = scales::pseudo_log_trans(10), labels = scales::label_number(big.mark = ",")) +
    labs(title = "Period-level posterior predictive check", subtitle = hmc_gate_label, x = "Observed notified cases", y = "Predicted notified cases") +
    publication_theme()
  ggsave(file.path(paths$figure, "ceara_episode_susceptibility_ppc.pdf"), ppc_plot, width = 183, height = 115, units = "mm", device = cairo_pdf)

  prior_draws <- make_prior_draws(bundle$prepared$periods, bundle$priors)
  posterior_parameter_draws <- smoothed_draws |>
    select(period_id, period_type, lambda, attack_prob, q) |>
    mutate(source = "Posterior")
  prior_parameter_draws <- prior_draws |>
    mutate(source = "Prior")
  prior_posterior <- bind_rows(prior_parameter_draws, posterior_parameter_draws)
  prior_plot <- ggplot(prior_posterior, aes(attack_prob, colour = source, fill = source)) +
    geom_density(alpha = .12, linewidth = .45) +
    facet_wrap(~ period_type, scales = "free_y") +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    labs(title = "Prior versus posterior attack probability", x = "Period attack probability", y = "Density") +
    publication_theme()
  detection_plot <- ggplot(prior_posterior, aes(q, colour = source, fill = source)) +
    geom_density(alpha = .12, linewidth = .45) +
    scale_x_continuous(labels = scales::label_percent(accuracy = .1)) +
    labs(title = "Prior versus posterior detection probability", x = "q", y = "Density") +
    publication_theme()
  ggsave(file.path(paths$figure, "ceara_episode_susceptibility_prior_posterior.pdf"), prior_plot / detection_plot, width = 183, height = 180, units = "mm", device = cairo_pdf)

  serology_plot <- ggplot(bundle$prepared$serology_metadata, aes(study_id, n_positive / n_tested)) +
    geom_point(size = 2, colour = "#CC79A7") +
    geom_errorbar(aes(ymin = pmax(0, n_positive / n_tested - 1.96 * sqrt((n_positive / n_tested) * (1 - n_positive / n_tested) / n_tested)), ymax = pmin(1, n_positive / n_tested + 1.96 * sqrt((n_positive / n_tested) * (1 - n_positive / n_tested) / n_tested))), width = .1, colour = "#CC79A7") +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
    labs(title = "Available serology: metadata only", subtitle = "Municipality-scale surveys are excluded from the primary Ceara state likelihood", x = NULL, y = "Observed seropositive proportion") +
    publication_theme()
  ggsave(file.path(paths$figure, "ceara_episode_susceptibility_serology_status.pdf"), serology_plot, width = 183, height = 95, units = "mm", device = cairo_pdf)

  documentation_written <- write_pilot_documentation(paths, bundle, hmc, smoothed_draws, filtered_episode, smoothed_episode)
  message("[gate] ", hmc_gate_label)
  message("[documentation] ", if (documentation_written) "written" else "not written because the HMC gate failed")
  invisible(list(hmc = hmc, documentation_written = documentation_written))
}


if (sys.nframe() == 0L) {
  run_ceara_episode_susceptibility_diagnostics()
}
