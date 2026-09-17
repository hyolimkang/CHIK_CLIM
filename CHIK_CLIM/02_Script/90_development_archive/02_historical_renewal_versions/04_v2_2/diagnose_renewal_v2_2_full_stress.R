# Read-only diagnostics for the unchanged v2.2 full-period stress test.
required_packages <- c("here", "rstan", "posterior", "ggplot2", "dplyr", "tidyr", "tibble", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(posterior)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
})
summ <- function(x, dates, variable) tibble(week_start = dates, variable = variable, q025 = apply(x, 2, quantile, .025), q25 = apply(x, 2, quantile, .25), median = apply(x, 2, median), q75 = apply(x, 2, quantile, .75), q975 = apply(x, 2, quantile, .975))
iv <- function(x) c(q025 = quantile(x, .025), median = median(x), q975 = quantile(x, .975))
theme_v2 <- function() theme_classic(base_size = 8.5) + theme(legend.position = "top", legend.title = element_blank(), panel.grid.major.y = element_line(colour = "grey90"), plot.title = element_text(face = "bold"))
add_band <- function(p, d, colour) p + geom_ribbon(data = d, aes(ymin = q025, ymax = q975), fill = colour, alpha = .13) + geom_ribbon(data = d, aes(ymin = q25, ymax = q75), fill = colour, alpha = .24) + geom_line(data = d, aes(y = median), colour = colour, linewidth = .65)

run_v2_2_full_stress_diagnostics <- function() {
  fit_path <- Sys.getenv("RENEWAL_V2_2_FULL_STRESS_FIT", here::here("02_Script/stan/renewal_ceara_v2_2_full_stress_fit.rds"))
  if (!file.exists(fit_path)) stop("Full-stress fit bundle not found: ", fit_path)
  b <- readRDS(fit_path)
  if (!isTRUE(b$config$stress_test) || !inherits(b$fit, "stanfit")) stop("Expected v2.2 full-stress fit bundle")
  table_dir <- here::here("03_Output/tables/renewal_v2_2/full_period_stress_test")
  figure_dir <- here::here("03_Output/figures/renewal_v2_2/full_period_stress_test")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dates <- as.Date(b$weekly_data$week_start)
  years <- as.integer(format(dates, "%Y"))
  fit <- b$fit
  # HMC convergence applies to declared parameters explored by the sampler,
  # not transformed or generated quantities such as C_pred.
  sampler_pars <- c(
    "sigma_R", "log_infection_scale", "log_hazard_relative", "p_symp",
    "rho_sym_brazil", "logit_rho_sym_ceara_mid", "reporting_trend",
    "sigma_geo", "z_site", "phi_obs"
  )
  raw <- rstan::extract(fit, pars = sampler_pars, permuted = FALSE, inc_warmup = FALSE)
  ds <- posterior::as_draws_array(raw)
  cs <- as.data.frame(posterior::summarise_draws(ds, posterior::rhat, posterior::ess_bulk, posterior::ess_tail))
  names(cs) <- sub("^posterior::", "", names(cs))
  if (".variable" %in% names(cs)) names(cs)[names(cs) == ".variable"] <- "parameter"
  bfmi <- rstan::get_bfmi(fit)
  hmc <- tibble(metric = c("divergences", "max_treedepth_hits", "maximum_Rhat", "n_Rhat_gt_1.01", "minimum_bulk_ESS", "minimum_tail_ESS", "elapsed_seconds", paste0("BFMI_chain_", seq_along(bfmi))), value = c(rstan::get_num_divergent(fit), rstan::get_num_max_treedepth(fit), max(cs$rhat, na.rm = TRUE), sum(cs$rhat > 1.01, na.rm = TRUE), min(cs$ess_bulk, na.rm = TRUE), min(cs$ess_tail, na.rm = TRUE), b$config$elapsed_seconds, bfmi))
  gate <- hmc$value[hmc$metric == "divergences"] == 0 && hmc$value[hmc$metric == "max_treedepth_hits"] == 0 && hmc$value[hmc$metric == "maximum_Rhat"] <= 1.01 && hmc$value[hmc$metric == "minimum_bulk_ESS"] >= 100 && hmc$value[hmc$metric == "minimum_tail_ESS"] >= 100 && all(bfmi >= .3)
  write.csv(hmc, file.path(table_dir, "hmc_diagnostics.csv"), row.names = FALSE)
  write.csv(head(arrange(cs, desc(rhat)), 30), file.path(table_dir, "worst_30_rhat.csv"), row.names = FALSE)
  write.csv(head(arrange(cs, ess_bulk), 30), file.path(table_dir, "worst_30_bulk_ess.csv"), row.names = FALSE)
  write.csv(head(arrange(cs, ess_tail), 30), file.path(table_dir, "worst_30_tail_ess.csv"), row.names = FALSE)
  pars <- c("C_pred", "X", "S", "U", "immune_prop", "R0_t", "R_eff_t", "rho_sym_t", "overall_detection", "expected_reported_cases", "p_symp", "p_site", "state_attack_window")
  z <- rstan::extract(fit, pars = pars)
  s_prop <- z$S / (z$S + z$U)
  sm <- bind_rows(summ(z$C_pred, dates, "posterior_predictive_cases"), summ(z$X, dates, "latent_infections"), summ(s_prop, dates, "susceptible_proportion"), summ(z$immune_prop, dates, "immune_proportion"), summ(z$R0_t, dates, "R0"), summ(z$R_eff_t, dates, "Reff"), summ(z$rho_sym_t, dates, "symptomatic_reporting"), summ(z$overall_detection, dates, "overall_detection"), summ(z$expected_reported_cases, dates, "expected_reported_cases"))
  get <- function(x) filter(sm, variable == x)
  case_df <- tibble(week_start = dates, observed = b$weekly_data$cases)
  annual <- bind_rows(lapply(sort(unique(years)), function(y) {
    idx <- which(years == y)
    xs <- rowSums(z$X[, idx, drop = FALSE])
    det <- rowSums(z$X[, idx, drop = FALSE] * z$overall_detection[, idx, drop = FALSE]) / xs
    attack <- xs / (z$S[, idx[1]] + z$U[, idx[1]])
    data.frame(year = y, observed_cases = sum(b$weekly_data$cases[idx]), posterior_infections_q025 = quantile(xs, .025), posterior_infections_median = median(xs), posterior_infections_q975 = quantile(xs, .975), overall_detection_q025 = quantile(det, .025), overall_detection_median = median(det), overall_detection_q975 = quantile(det, .975), S_start_q025 = quantile(s_prop[, idx[1]], .025), S_start_median = median(s_prop[, idx[1]]), S_start_q975 = quantile(s_prop[, idx[1]], .975), S_end_q025 = quantile(s_prop[, idx[length(idx)]], .025), S_end_median = median(s_prop[, idx[length(idx)]]), S_end_q975 = quantile(s_prop[, idx[length(idx)]], .975), attack_increment_q025 = quantile(attack, .025), attack_increment_median = median(attack), attack_increment_q975 = quantile(attack, .975))
  }))
  write.csv(annual, file.path(table_dir, "annual_attack_decomposition.csv"), row.names = FALSE)
  cumulative <- t(apply(z$X, 1, cumsum))
  write.csv(bind_rows(summ(cumulative, dates, "cumulative_infections"), sm), file.path(table_dir, "trajectory_summary.csv"), row.names = FALSE)
  low_r0 <- get("R0") |>
    mutate(observed_cases = b$weekly_data$cases, r0_width = q975 - q025) |>
    filter(observed_cases <= 5, median > 10 | r0_width > 10)
  write.csv(low_r0, file.path(table_dir, "low_incidence_extreme_or_weak_R0_weeks.csv"), row.names = FALSE)
  xscale <- scale_x_date(breaks = seq(as.Date("2015-01-01"), as.Date("2026-01-01"), by = "2 years"), date_labels = "%Y")
  p_cases <- add_band(ggplot(get("posterior_predictive_cases"), aes(week_start)), get("posterior_predictive_cases"), "#56B4E9") + geom_point(data = case_df, aes(week_start, observed), size = .35) + xscale + labs(title = "Reported cases and posterior prediction", y = "Weekly cases", x = NULL) + theme_v2()
  p_x <- add_band(ggplot(get("latent_infections"), aes(week_start)), get("latent_infections"), "#D55E00") + xscale + labs(title = "Latent infection incidence", y = "Infections", x = NULL) + theme_v2()
  states <- bind_rows(mutate(get("susceptible_proportion"), state = "Susceptible"), mutate(get("immune_proportion"), state = "Infection-derived immune"))
  p_state <- ggplot(states, aes(week_start, median, colour = state, fill = state)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .13, colour = NA) +
    geom_line() +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    xscale +
    labs(title = "Population susceptibility and immunity", y = "Proportion", x = NULL) +
    theme_v2()
  rr <- bind_rows(mutate(get("R0"), number = "R0"), mutate(get("Reff"), number = "Reff")) |> filter(match(week_start, dates) > b$stan_data$seed_weeks)
  p_r <- ggplot(rr, aes(week_start, median, colour = number, fill = number)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .12, colour = NA) +
    geom_line() +
    xscale +
    labs(title = "Time-varying reproduction numbers", y = "Reproduction number", x = NULL) +
    theme_v2()
  reporting <- bind_rows(mutate(get("symptomatic_reporting"), component = "Symptomatic reporting"), mutate(get("overall_detection"), component = "Overall detection"))
  p_q <- ggplot(reporting, aes(week_start, median, colour = component, fill = component)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .12, colour = NA) +
    geom_line() +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    xscale +
    labs(title = "Case ascertainment", y = "Probability", x = NULL) +
    theme_v2()
  sero_obs <- tibble(site = c("Juazeiro do Norte", "Quixada"), observed = c(103 / 404, 289 / 409), lower = c(binom.test(103, 404)$conf.int[1], binom.test(289, 409)$conf.int[1]), upper = c(binom.test(103, 404)$conf.int[2], binom.test(289, 409)$conf.int[2]))
  sero_mod <- bind_rows(lapply(1:2, function(j) tibble(site = sero_obs$site[j], q025 = quantile(z$p_site[, j], .025), median = median(z$p_site[, j]), q975 = quantile(z$p_site[, j], .975))))
  p_sero <- ggplot(left_join(sero_obs, sero_mod, by = "site"), aes(site, median)) +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = .15) +
    geom_point(colour = "#76558F") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = .15, position = position_nudge(x = .12)) +
    geom_point(aes(y = observed), shape = 18, position = position_nudge(x = .12), colour = "#E69F00") +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    labs(title = "Serology posterior fit", y = "Cumulative proportion", x = NULL) +
    theme_v2()
  subtitle <- if (gate) "Unchanged v2.2 stress test, 2015<U+2013>2025" else "UNCONVERGED STRESS TEST <U+2014> DO NOT INTERPRET POSTERIOR TRAJECTORIES"
  fig <- (p_cases | p_x) / (p_state | p_r) / (p_q | p_sero) + plot_annotation(title = "Chikungunya transmission with explicit demographic accounting, Cear<U+00E1>", subtitle = subtitle, tag_levels = "A", theme = theme(plot.subtitle = element_text(colour = if (gate) "grey30" else "#D55E00", face = if (gate) "plain" else "bold")))
  ggsave(file.path(figure_dir, paste0("renewal_v2_2_trajectories_ceara_full_stress", if (gate) "" else "_UNCONVERGED", ".pdf")), fig, width = 183, height = 225, units = "mm", device = cairo_pdf)
  pp <- ggplot(tibble(value = c(rbeta(50000, 30, 28), z$p_symp), source = rep(c("Prior Beta(30,28)", "Posterior"), c(50000, nrow(z$p_symp)))), aes(value, colour = source)) +
    geom_density() +
    labs(title = "p_symp prior versus posterior", x = "p_symp", y = "Density") +
    theme_v2()
  ggsave(file.path(figure_dir, "p_symp_prior_vs_posterior.pdf"), pp, width = 160, height = 100, units = "mm", device = cairo_pdf)
  message("[gate] ", if (gate) "PASS" else "FAIL <U+2014> figures labelled non-converged; no scientific interpretation")
}
if (sys.nframe() == 0) run_v2_2_full_stress_diagnostics()
