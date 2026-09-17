# v2.3_state: read-only diagnostics and six-panel trajectory figure.
#
# This file reads an existing fit only. It writes HMC diagnostics, annual
# attack-rate summaries, weak-identification flags, and the six-panel PDF.

## 1. Packages, paths, and reusable helpers ------------------------------------

required_packages <- c(
  "here", "rstan", "posterior", "ggplot2", "dplyr", "tibble",
  "patchwork", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(posterior)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(patchwork)
})

summarise_weekly_draws <- function(draws, dates, variable_name) {
  tibble(
    week_start = dates,
    variable = variable_name,
    q025 = apply(draws, 2, quantile, probs = 0.025),
    q25 = apply(draws, 2, quantile, probs = 0.25),
    median = apply(draws, 2, median),
    q75 = apply(draws, 2, quantile, probs = 0.75),
    q975 = apply(draws, 2, quantile, probs = 0.975)
  )
}

diagnostic_theme <- function() {
  theme_classic(base_size = 8.5) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90")
    )
}

add_interval_bands <- function(plot, summary_data, colour) {
  plot +
    geom_ribbon(
      data = summary_data, aes(ymin = q025, ymax = q975),
      fill = colour, alpha = 0.13
    ) +
    geom_ribbon(
      data = summary_data, aes(ymin = q25, ymax = q75),
      fill = colour, alpha = 0.24
    ) +
    geom_line(
      data = summary_data, aes(y = median), colour = colour, linewidth = 0.65
    )
}
## 2. Main diagnostic workflow --------------------------------------------------

run_v2_3_diagnostics <- function() {
  ## 2. Read fit and create output directories ---------------------------------
  fit_path <- Sys.getenv(
    "RENEWAL_V2_3_FIT",
    here::here("02_Script", "stan", "renewal_ceara_v2_3_state_fit.rds")
  )
  if (!file.exists(fit_path)) stop("Fit not found: ", fit_path)

  fit_bundle <- readRDS(fit_path)
  fit <- fit_bundle$fit
  dates <- as.Date(fit_bundle$weekly_data$week_start)
  year <- as.integer(format(dates, "%Y"))

  table_directory <- here::here("03_Output", "tables", "renewal_v2_3")
  figure_directory <- here::here("03_Output", "figures", "renewal_v2_3")
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
  ## 3. HMC diagnostics and computational gate ----------------------------------
  # Generated quantities are excluded: this gate is based on sampler parameters.
  sampler_parameters <- c(
    "alpha_R", "beta_sin", "beta_cos", "phi_R", "sigma_R", "z_R",
    "log_seed_hazard", "q", "sigma_geo", "z_site", "phi_obs"
  )
  sampler_draws <- posterior::as_draws_array(rstan::extract(
    fit,
    pars = sampler_parameters,
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
  if (".variable" %in% names(parameter_summary)) {
    names(parameter_summary)[names(parameter_summary) == ".variable"] <- "parameter"
  }

  bfmi <- rstan::get_bfmi(fit)
  hmc_summary <- tibble(
    metric = c(
      "divergences", "max_treedepth_hits", "maximum_Rhat",
      "n_Rhat_gt_1.01", "minimum_bulk_ESS", "minimum_tail_ESS",
      "elapsed_seconds", paste0("BFMI_chain_", seq_along(bfmi))
    ),
    value = c(
      rstan::get_num_divergent(fit),
      rstan::get_num_max_treedepth(fit),
      max(parameter_summary$rhat, na.rm = TRUE),
      sum(parameter_summary$rhat > 1.01, na.rm = TRUE),
      min(parameter_summary$ess_bulk, na.rm = TRUE),
      min(parameter_summary$ess_tail, na.rm = TRUE),
      fit_bundle$config$elapsed_seconds,
      bfmi
    )
  )
  passed_gate <-
    hmc_summary$value[hmc_summary$metric == "divergences"] == 0 &&
      hmc_summary$value[hmc_summary$metric == "max_treedepth_hits"] == 0 &&
      hmc_summary$value[hmc_summary$metric == "maximum_Rhat"] <= 1.01 &&
      hmc_summary$value[hmc_summary$metric == "minimum_bulk_ESS"] >= 100 &&
      hmc_summary$value[hmc_summary$metric == "minimum_tail_ESS"] >= 100 &&
      all(bfmi >= 0.3)

  write.csv(hmc_summary, file.path(table_directory, "hmc_diagnostics.csv"), row.names = FALSE)
  write.csv(
    head(arrange(parameter_summary, desc(rhat)), 30),
    file.path(table_directory, "worst_30_rhat.csv"),
    row.names = FALSE
  )
  write.csv(
    head(arrange(parameter_summary, ess_bulk), 30),
    file.path(table_directory, "worst_30_bulk_ess.csv"),
    row.names = FALSE
  )
  write.csv(
    head(arrange(parameter_summary, ess_tail), 30),
    file.path(table_directory, "worst_30_tail_ess.csv"),
    row.names = FALSE
  )
  ## 4. Posterior trajectories and annual attack-rate decomposition -------------
  draws <- rstan::extract(
    fit,
    pars = c(
      "C_pred", "X", "S_prop", "immune_prop", "R0_t", "R_eff_t",
      "expected_reported_cases", "q", "p_site", "state_attack_window"
    )
  )
  trajectory_summary <- bind_rows(
    summarise_weekly_draws(draws$C_pred, dates, "C_pred"),
    summarise_weekly_draws(draws$X, dates, "X"),
    summarise_weekly_draws(draws$S_prop, dates, "S_prop"),
    summarise_weekly_draws(draws$immune_prop, dates, "attack"),
    summarise_weekly_draws(draws$R0_t, dates, "R0"),
    summarise_weekly_draws(draws$R_eff_t, dates, "Reff"),
    summarise_weekly_draws(draws$expected_reported_cases, dates, "expected_cases")
  )
  get_trajectory <- function(variable_name) {
    filter(trajectory_summary, variable == variable_name)
  }

  annual_attack <- bind_rows(lapply(sort(unique(year)), function(current_year) {
    index <- which(year == current_year)
    annual_infections <- rowSums(draws$X[, index, drop = FALSE])
    attack_increment <- annual_infections / fit_bundle$stan_data$N_start[index[1]]

    data.frame(
      year = current_year,
      observed_cases = sum(fit_bundle$weekly_data$cases[index]),
      posterior_infections_q025 = quantile(annual_infections, 0.025),
      posterior_infections_median = median(annual_infections),
      posterior_infections_q975 = quantile(annual_infections, 0.975),
      overall_detection_q025 = quantile(draws$q, 0.025),
      overall_detection_median = median(draws$q),
      overall_detection_q975 = quantile(draws$q, 0.975),
      S_start_q025 = quantile(draws$S_prop[, index[1]], 0.025),
      S_start_median = median(draws$S_prop[, index[1]]),
      S_start_q975 = quantile(draws$S_prop[, index[1]], 0.975),
      S_end_q025 = quantile(draws$S_prop[, tail(index, 1)], 0.025),
      S_end_median = median(draws$S_prop[, tail(index, 1)]),
      S_end_q975 = quantile(draws$S_prop[, tail(index, 1)], 0.975),
      attack_increment_q025 = quantile(attack_increment, 0.025),
      attack_increment_median = median(attack_increment),
      attack_increment_q975 = quantile(attack_increment, 0.975)
    )
  }))
  write.csv(
    annual_attack,
    file.path(table_directory, "annual_attack_decomposition.csv"),
    row.names = FALSE
  )
  write.csv(
    trajectory_summary,
    file.path(table_directory, "trajectory_summary.csv"),
    row.names = FALSE
  )

  low_incidence_r0_flags <- get_trajectory("R0") |>
    mutate(
      observed_cases = fit_bundle$weekly_data$cases,
      width = q975 - q025
    ) |>
    filter(observed_cases <= 5 & (median > 10 | width > 10))
  write.csv(
    low_incidence_r0_flags,
    file.path(table_directory, "low_incidence_extreme_or_weak_R0_weeks.csv"),
    row.names = FALSE
  )

  # Backward-compatible aliases for the original six-panel plotting block.
  # The quantities are named explicitly above; these aliases prevent any
  # change in plotted values while that block remains behaviour-identical.
  b <- fit_bundle
  d <- dates
  z <- draws
  get <- get_trajectory
  gate <- passed_gate
  fd <- figure_directory
  theme2 <- diagnostic_theme
  band <- add_interval_bands

  ## 5. Six-panel figure --------------------------------------------------------
  xs <- scale_x_date(breaks = seq(as.Date("2015-01-01"), as.Date("2026-01-01"), by = "2 years"), date_labels = "%Y")
  p1 <- band(ggplot(get("C_pred"), aes(week_start)), get("C_pred"), "#56B4E9") + geom_point(data = tibble(week_start = d, cases = b$weekly_data$cases), aes(week_start, cases), size = .3) + xs + labs(title = "Reported cases and posterior prediction", x = NULL, y = "Cases") + theme2()
  p2 <- band(ggplot(get("X"), aes(week_start)), get("X"), "#D55E00") + xs + labs(title = "Latent infection incidence", x = NULL, y = "Infections") + theme2()
  st <- bind_rows(mutate(get("S_prop"), state = "Susceptible"), mutate(get("attack"), state = "Infection-derived immune"))
  p3 <- ggplot(st, aes(week_start, median, colour = state, fill = state)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .12, colour = NA) +
    geom_line() +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    xs +
    labs(title = "Population susceptibility and immunity", x = NULL, y = "Proportion") +
    theme2()
  rr <- bind_rows(mutate(get("R0"), n = "R0"), mutate(get("Reff"), n = "Reff"))
  p4 <- ggplot(rr, aes(week_start, median, colour = n, fill = n)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .12, colour = NA) +
    geom_line() +
    xs +
    labs(title = "Reproduction numbers", x = NULL, y = "R") +
    theme2()
  p5 <- ggplot(tibble(q = z$q), aes(y = q)) +
    geom_boxplot(fill = "#E69F00") +
    scale_y_continuous(labels = scales::label_percent(), limits = c(0, NA)) +
    labs(title = "Overall infection-to-report probability", x = NULL, y = "q") +
    theme2()
  obs <- tibble(site = c("Juazeiro do Norte", "Quixada"), v = c(103 / 404, 289 / 409))
  mod <- bind_rows(lapply(1:2, function(j) tibble(site = obs$site[j], q025 = quantile(z$p_site[, j], .025), median = median(z$p_site[, j]), q975 = quantile(z$p_site[, j], .975))))
  p6 <- ggplot(left_join(obs, mod, by = "site"), aes(site, median)) +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = .15) +
    geom_point() +
    geom_point(aes(y = v), shape = 18, colour = "#D55E00") +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    labs(title = "Serology posterior fit", x = NULL, y = "Proportion") +
    theme2()
  sub <- if (gate) "v2.3_state full-period fit" else "UNCONVERGED <U+2014> DO NOT INTERPRET POSTERIOR TRAJECTORIES"
  fig <- (p1 | p2) / (p3 | p4) / (p5 | p6) + plot_annotation(title = "Cear<U+00E1> renewal v2.3 state model", subtitle = sub, tag_levels = "A", theme = theme(plot.subtitle = element_text(colour = if (gate) "grey30" else "#D55E00", face = "bold")))
  ggsave(file.path(fd, paste0("renewal_v2_3_state_trajectories", if (gate) "" else "_UNCONVERGED", ".pdf")), fig, width = 183, height = 225, units = "mm", device = cairo_pdf)
  message("[gate] ", if (gate) "PASS" else "FAIL")
}
if (sys.nframe() == 0) run_v2_3_diagnostics()
