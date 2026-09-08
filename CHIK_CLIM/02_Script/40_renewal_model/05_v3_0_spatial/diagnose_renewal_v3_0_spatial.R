# =============================================================================
# diagnose_renewal_v3_0_spatial.R
#
# Diagnostics for a v3.0 fit bundle. It accepts either the small debugging
# pilot produced by fit_renewal_v3_0_spatial.R or a future all-Ceara fit
# bundle with the same list structure. It does not start sampling.
#
# Optional environment override:
#   RENEWAL_V3_0_FIT_PATH=path/to/fit_bundle.rds
# =============================================================================

required_packages <- c("here", "rstan")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(rstan)
})

posterior_interval <- function(x) {
  c(
    median = stats::median(x),
    q2.5 = stats::quantile(x, 0.025, names = FALSE),
    q97.5 = stats::quantile(x, 0.975, names = FALSE)
  )
}

summarise_time_draws <- function(draws, dates, variable, observed = NULL) {
  if (is.null(dim(draws))) draws <- matrix(draws, ncol = length(dates))
  summaries <- t(apply(draws, 2, posterior_interval))
  data.frame(
    week_start = dates,
    variable = variable,
    observed = if (is.null(observed)) NA_real_ else observed,
    median = summaries[, "median"],
    q2.5 = summaries[, "q2.5"],
    q97.5 = summaries[, "q97.5"],
    row.names = NULL
  )
}

summarise_v3_0_fit <- function(bundle) {
  fit <- bundle$fit
  prepared <- bundle$prepared
  dates <- prepared$week_dates
  draws <- rstan::extract(
    fit,
    pars = c(
      "C_ceara", "C_ceara_pred", "X_ceara_total",
      "S_ceara_prop", "immune_ceara_prop",
      "p_symp", "rho_sym_global", "overall_detection",
      "p_sero", "sero_positive_pred"
    ),
    permuted = TRUE
  )

  time_summary <- rbind(
    summarise_time_draws(draws$C_ceara_pred, dates, "Ceara posterior-predicted cases",
                         observed = draws$C_ceara[1, ]),
    summarise_time_draws(draws$X_ceara_total, dates, "Ceara latent infections"),
    summarise_time_draws(draws$S_ceara_prop, dates, "Ceara susceptible proportion"),
    summarise_time_draws(draws$immune_ceara_prop, dates, "Ceara immune proportion")
  )

  reporting <- do.call(rbind, lapply(
    c("p_symp", "rho_sym_global", "overall_detection"),
    function(name) {
      interval <- posterior_interval(draws[[name]])
      data.frame(parameter = name, t(interval), row.names = NULL)
    }
  ))

  serology <- prepared$serology
  serology_summary <- do.call(rbind, lapply(seq_len(nrow(serology)), function(j) {
    p_interval <- posterior_interval(draws$p_sero[, j])
    count_interval <- posterior_interval(draws$sero_positive_pred[, j])
    data.frame(
      site = serology$site[j],
      muni6 = serology$muni6[j],
      observed_positive = serology$positive[j],
      tested = serology$n[j],
      observed_prevalence = serology$positive[j] / serology$n[j],
      predicted_prevalence_median = p_interval["median"],
      predicted_prevalence_q2.5 = p_interval["q2.5"],
      predicted_prevalence_q97.5 = p_interval["q97.5"],
      predicted_positive_median = count_interval["median"],
      predicted_positive_q2.5 = count_interval["q2.5"],
      predicted_positive_q97.5 = count_interval["q97.5"],
      row.names = NULL
    )
  }))

  summary_matrix <- summary(fit)$summary
  rhat <- summary_matrix[, "Rhat"]
  ess <- summary_matrix[, "n_eff"]
  finite_rhat <- rhat[is.finite(rhat)]
  finite_ess <- ess[is.finite(ess)]
  hmc <- data.frame(
    metric = c(
      "divergences", "maximum_treedepth_hits", "minimum_bfmi",
      "maximum_rhat", "parameters_rhat_gt_1.01", "minimum_ess"
    ),
    value = c(
      rstan::get_num_divergent(fit),
      rstan::get_num_max_treedepth(fit),
      min(rstan::get_bfmi(fit)),
      if (length(finite_rhat)) max(finite_rhat) else NA_real_,
      sum(rhat > 1.01, na.rm = TRUE),
      if (length(finite_ess)) min(finite_ess) else NA_real_
    )
  )

  scale_draws <- rstan::extract(
    fit,
    pars = c("sigma_R", "sigma_muni_R", "p_symp", "rho_sym_global", "phi_obs"),
    permuted = TRUE
  )
  scale_matrix <- do.call(cbind, scale_draws)
  colnames(scale_matrix) <- c("sigma_R", "sigma_muni_R", "p_symp", "rho_sym_global", "phi_obs")

  list(
    time_summary = time_summary,
    reporting = reporting,
    serology = serology_summary,
    hmc = hmc,
    scale_correlation = stats::cor(scale_matrix),
    fortaleza_note = "No Fortaleza cohort count, denominator, or collection window is stored in v3.0. It is excluded from the likelihood; add an external cohort specification before calculating a numerical PPC."
  )
}

plot_v3_0_diagnostics <- function(diagnostics, output_path) {
  time <- diagnostics$time_summary
  grDevices::pdf(output_path, width = 11, height = 8.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  case <- time[time$variable == "Ceara posterior-predicted cases", ]
  graphics::plot(case$week_start, case$observed, type = "h", col = "grey50",
                 xlab = "", ylab = "weekly cases",
                 main = "Ceara observed and posterior-predicted cases")
  graphics::lines(case$week_start, case$median, col = "#0072B2", lwd = 2)
  graphics::lines(case$week_start, case$q2.5, col = "#0072B2", lty = 2)
  graphics::lines(case$week_start, case$q97.5, col = "#0072B2", lty = 2)

  latent <- time[time$variable == "Ceara latent infections", ]
  graphics::plot(latent$week_start, latent$median, type = "l", col = "#D55E00", lwd = 2,
                 ylim = range(c(latent$q2.5, latent$q97.5)),
                 xlab = "", ylab = "weekly latent infections",
                 main = "Ceara latent infections")
  graphics::lines(latent$week_start, latent$q2.5, col = "#D55E00", lty = 2)
  graphics::lines(latent$week_start, latent$q97.5, col = "#D55E00", lty = 2)

  susceptible <- time[time$variable == "Ceara susceptible proportion", ]
  immune <- time[time$variable == "Ceara immune proportion", ]
  graphics::plot(susceptible$week_start, susceptible$median, type = "l", col = "#009E73",
                 ylim = c(0, 1), lwd = 2, xlab = "", ylab = "proportion",
                 main = "Population-weighted Ceará susceptibility and immunity")
  graphics::lines(susceptible$week_start, susceptible$q2.5, col = "#009E73", lty = 2)
  graphics::lines(susceptible$week_start, susceptible$q97.5, col = "#009E73", lty = 2)
  graphics::lines(immune$week_start, immune$median, col = "#CC79A7", lwd = 2)
  graphics::legend("right", c("susceptible", "immune"), col = c("#009E73", "#CC79A7"),
                   lty = 1, bty = "n")

  reporting <- diagnostics$reporting
  graphics::barplot(reporting$median, names.arg = reporting$parameter, las = 2,
                    ylim = c(0, max(reporting$q97.5, 0.1)),
                    ylab = "probability", main = "Pooled reporting parameters")
  graphics::arrows(
    x0 = seq_len(nrow(reporting)), y0 = reporting$q2.5,
    x1 = seq_len(nrow(reporting)), y1 = reporting$q97.5,
    angle = 90, code = 3, length = 0.05
  )
}

run_v3_0_diagnostics <- function() {
  default_path <- here::here(
    "02_Script/stan/renewal_ceara_v3_0_spatial_debug_pilot.rds"
  )
  fit_path <- Sys.getenv("RENEWAL_V3_0_FIT_PATH", default_path)
  if (!file.exists(fit_path)) {
    stop(
      "No v3.0 fit bundle found at ", fit_path,
      ". Run the explicitly opt-in debug pilot or provide a completed fit bundle; this diagnostic script does not sample."
    )
  }

  bundle <- readRDS(fit_path)
  if (!all(c("fit", "prepared") %in% names(bundle))) {
    stop("Fit bundle must contain both fit and prepared objects")
  }
  diagnostics <- summarise_v3_0_fit(bundle)

  table_dir <- here::here("03_Output/tables/renewal_v3_0")
  figure_dir <- here::here("03_Output/figures/renewal_v3_0")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  utils::write.csv(
    diagnostics$time_summary,
    file.path(table_dir, "renewal_v3_0_spatial_time_diagnostics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    diagnostics$reporting,
    file.path(table_dir, "renewal_v3_0_spatial_reporting_diagnostics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    diagnostics$serology,
    file.path(table_dir, "renewal_v3_0_spatial_serology_ppc.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    diagnostics$hmc,
    file.path(table_dir, "renewal_v3_0_spatial_hmc_diagnostics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    diagnostics$scale_correlation,
    file.path(table_dir, "renewal_v3_0_spatial_scale_parameter_correlation.csv")
  )
  plot_v3_0_diagnostics(
    diagnostics,
    file.path(figure_dir, "renewal_v3_0_spatial_diagnostics.pdf")
  )
  message("[Fortaleza] ", diagnostics$fortaleza_note)
  invisible(diagnostics)
}

if (sys.nframe() == 0) {
  run_v3_0_diagnostics()
}
