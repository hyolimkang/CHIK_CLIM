# =============================================================================
# diagnose_renewal_v3_1_minimal_spatial.R
#
# Creates six-unit posterior and sampler diagnostics from the short pilot or a
# compatible future fit bundle. It never starts sampling.
# =============================================================================

required_packages <- c("here", "rstan", "posterior")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(posterior)
})

posterior_interval <- function(x) {
  c(
    median = stats::median(x),
    q2.5 = stats::quantile(x, 0.025, names = FALSE),
    q97.5 = stats::quantile(x, 0.975, names = FALSE)
  )
}

summarise_time <- function(x, dates, variable, observed = NULL) {
  if (is.null(dim(x))) x <- matrix(x, ncol = length(dates))
  q <- t(apply(x, 2, posterior_interval))
  data.frame(
    week_start = dates, variable = variable,
    observed = if (is.null(observed)) NA_real_ else observed,
    median = q[, "median"], q2.5 = q[, "q2.5"], q97.5 = q[, "q97.5"],
    row.names = NULL
  )
}

fit_hmc_diagnostics_v3_1 <- function(fit) {
  draws_array <- posterior::as_draws_array(
    rstan::extract(fit, permuted = FALSE, inc_warmup = FALSE)
  )
  convergence <- as.data.frame(posterior::summarise_draws(
    draws_array, posterior::rhat, posterior::ess_bulk, posterior::ess_tail
  ))
  names(convergence) <- sub("^posterior::", "", names(convergence))
  if (".variable" %in% names(convergence)) {
    names(convergence)[names(convergence) == ".variable"] <- "variable"
  }
  finite_rhat <- convergence$rhat[is.finite(convergence$rhat)]
  list(
    divergences = rstan::get_num_divergent(fit),
    max_treedepth_hits = rstan::get_num_max_treedepth(fit),
    bfmi = rstan::get_bfmi(fit),
    max_rhat = if (length(finite_rhat)) max(finite_rhat) else NA_real_,
    min_bulk_ess = min(convergence$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(convergence$ess_tail, na.rm = TRUE),
    convergence = convergence
  )
}

summarise_v3_1_fit <- function(bundle) {
  fit <- bundle$fit
  prepared <- bundle$prepared
  dates <- prepared$week_dates
  draws <- rstan::extract(
    fit,
    pars = c(
      "C_ceara", "C_ceara_pred", "X_ceara_total",
      "S_ceara_prop", "immune_ceara_prop", "R0_shared", "R_eff_ceara",
      "S_prop", "X", "p_symp", "rho_sym", "overall_detection",
      "p_sero", "sero_positive_pred"
    ),
    permuted = TRUE
  )

  time <- rbind(
    summarise_time(draws$C_ceara_pred, dates, "Ceara posterior-predicted cases",
      observed = draws$C_ceara[1, ]
    ),
    summarise_time(draws$X_ceara_total, dates, "Ceara latent infections"),
    summarise_time(draws$S_ceara_prop, dates, "Ceara susceptible proportion"),
    summarise_time(draws$immune_ceara_prop, dates, "Ceara immune proportion"),
    summarise_time(draws$R0_shared, dates, "Shared R0"),
    summarise_time(draws$R_eff_ceara, dates, "Population-weighted Ceara Reff")
  )

  unit_susceptibility <- do.call(rbind, lapply(seq_len(prepared$stan_data$K), function(k) {
    out <- summarise_time(
      draws$S_prop[, k, ], dates,
      paste0("Susceptible proportion: ", prepared$unit_labels[k])
    )
    out$unit_index <- k
    out$unit <- prepared$unit_labels[k]
    out
  }))

  endpoint_dates <- as.Date(c("2017-12-31", "2019-12-29"))
  endpoint_index <- match(endpoint_dates, dates)
  endpoint <- do.call(rbind, lapply(seq_len(prepared$stan_data$K), function(k) {
    do.call(rbind, lapply(seq_along(endpoint_index), function(i) {
      q <- posterior_interval(draws$S_prop[, k, endpoint_index[i]])
      data.frame(
        unit_index = k, unit = prepared$unit_labels[k],
        date = endpoint_dates[i],
        susceptible_median = q["median"],
        susceptible_q2.5 = q["q2.5"],
        susceptible_q97.5 = q["q97.5"],
        credible_interval_width = q["q97.5"] - q["q2.5"],
        row.names = NULL
      )
    }))
  }))

  cumulative_infections <- apply(draws$X, c(1, 2), sum)
  infection_intervals <- t(apply(cumulative_infections, 2, posterior_interval))
  contribution <- prepared$unit_summary
  contribution$inferred_infections_median <- infection_intervals[, "median"]
  contribution$inferred_infections_q2.5 <- infection_intervals[, "q2.5"]
  contribution$inferred_infections_q97.5 <- infection_intervals[, "q97.5"]
  contribution$inferred_infection_contribution <-
    contribution$inferred_infections_median / sum(contribution$inferred_infections_median)

  reporting <- do.call(rbind, lapply(
    c("p_symp", "rho_sym", "overall_detection"),
    function(name) {
      data.frame(
        parameter = name, t(posterior_interval(draws[[name]])),
        row.names = NULL
      )
    }
  ))

  serology <- do.call(rbind, lapply(seq_len(nrow(prepared$serology)), function(j) {
    p <- posterior_interval(draws$p_sero[, j])
    n <- posterior_interval(draws$sero_positive_pred[, j])
    data.frame(
      site = prepared$serology$site[j],
      observed_positive = prepared$serology$positive[j],
      tested = prepared$serology$n[j],
      observed_prevalence = prepared$serology$positive[j] / prepared$serology$n[j],
      predicted_prevalence_median = p["median"],
      predicted_prevalence_q2.5 = p["q2.5"],
      predicted_prevalence_q97.5 = p["q97.5"],
      predicted_positive_median = n["median"],
      predicted_positive_q2.5 = n["q2.5"],
      predicted_positive_q97.5 = n["q97.5"],
      row.names = NULL
    )
  }))

  hmc <- if (!is.null(bundle$pilot_diagnostics) &&
    all(c("convergence", "max_rhat") %in% names(bundle$pilot_diagnostics))) {
    bundle$pilot_diagnostics
  } else {
    fit_hmc_diagnostics_v3_1(fit)
  }
  worst <- hmc$convergence[order(hmc$convergence$rhat,
    decreasing = TRUE,
    na.last = NA
  ), , drop = FALSE]
  list(
    time = time,
    unit_susceptibility = unit_susceptibility,
    endpoints = endpoint,
    contribution = contribution,
    reporting = reporting,
    serology = serology,
    hmc_summary = data.frame(
      metric = c(
        "divergences", "max_treedepth_hits", "min_bfmi",
        "max_rhat", "min_bulk_ess", "min_tail_ess"
      ),
      value = c(
        hmc$divergences, hmc$max_treedepth_hits, min(hmc$bfmi),
        hmc$max_rhat, hmc$min_bulk_ess, hmc$min_tail_ess
      )
    ),
    hmc_convergence = hmc$convergence,
    worst_mixing = utils::head(worst, 10),
    fortaleza_note = "Fortaleza is represented as a transmission unit but has no serology likelihood or cohort values for a numerical external PPC."
  )
}

plot_v3_1_diagnostics <- function(d, output_path) {
  grDevices::pdf(output_path, width = 11, height = 8.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  cases <- d$time[d$time$variable == "Ceara posterior-predicted cases", ]
  graphics::plot(cases$week_start, cases$observed,
    type = "h", col = "grey50",
    xlab = "", ylab = "weekly cases",
    main = "Ceara observed and posterior-predicted cases"
  )
  graphics::lines(cases$week_start, cases$median, col = "#0072B2", lwd = 2)
  graphics::lines(cases$week_start, cases$q2.5, col = "#0072B2", lty = 2)
  graphics::lines(cases$week_start, cases$q97.5, col = "#0072B2", lty = 2)

  latent <- d$time[d$time$variable == "Ceara latent infections", ]
  graphics::plot(latent$week_start, latent$median,
    type = "l", col = "#D55E00", lwd = 2,
    ylim = range(c(latent$q2.5, latent$q97.5)),
    xlab = "", ylab = "weekly latent infections",
    main = "Total latent infections"
  )
  graphics::lines(latent$week_start, latent$q2.5, col = "#D55E00", lty = 2)
  graphics::lines(latent$week_start, latent$q97.5, col = "#D55E00", lty = 2)

  s <- d$time[d$time$variable == "Ceara susceptible proportion", ]
  u <- d$time[d$time$variable == "Ceara immune proportion", ]
  graphics::plot(s$week_start, s$median,
    type = "l", col = "#009E73", ylim = c(0, 1),
    xlab = "", ylab = "proportion",
    main = "Population-weighted Ceara susceptibility and immunity"
  )
  graphics::lines(s$week_start, s$q2.5, col = "#009E73", lty = 2)
  graphics::lines(s$week_start, s$q97.5, col = "#009E73", lty = 2)
  graphics::lines(u$week_start, u$median, col = "#CC79A7", lwd = 2)
  graphics::legend("right", c("susceptible", "immune"),
    col = c("#009E73", "#CC79A7"), lty = 1, bty = "n"
  )

  r0 <- d$time[d$time$variable == "Shared R0", ]
  reff <- d$time[d$time$variable == "Population-weighted Ceara Reff", ]
  graphics::plot(r0$week_start, r0$median,
    type = "l", col = "#0072B2", lwd = 2,
    ylim = range(c(r0$q2.5, r0$q97.5, reff$q2.5, reff$q97.5)),
    xlab = "", ylab = "reproduction number",
    main = "Shared R0 and weighted effective reproduction number"
  )
  graphics::lines(r0$week_start, r0$q2.5, col = "#0072B2", lty = 2)
  graphics::lines(r0$week_start, r0$q97.5, col = "#0072B2", lty = 2)
  graphics::lines(reff$week_start, reff$median, col = "#E69F00", lwd = 2)
  graphics::abline(h = 1, lty = 3)
  graphics::legend("topright", c("R0", "Reff"),
    col = c("#0072B2", "#E69F00"),
    lty = 1, bty = "n"
  )

  graphics::par(mfrow = c(2, 3))
  for (k in 1:6) {
    x <- d$unit_susceptibility[d$unit_susceptibility$unit_index == k, ]
    graphics::plot(x$week_start, x$median,
      type = "l", ylim = c(0, 1),
      xlab = "", ylab = "S/N", main = x$unit[1]
    )
    graphics::lines(x$week_start, x$q2.5, lty = 2)
    graphics::lines(x$week_start, x$q97.5, lty = 2)
  }
  graphics::par(mfrow = c(1, 1))

  graphics::barplot(d$reporting$median,
    names.arg = d$reporting$parameter, las = 2,
    ylim = c(0, max(d$reporting$q97.5, 0.1)),
    ylab = "probability", main = "Pooled detection parameters"
  )
}

run_v3_1_diagnostics <- function() {
  default_path <- here::here(
    "02_Script/stan/renewal_ceara_v3_1_minimal_spatial_pilot.rds"
  )
  fit_path <- Sys.getenv("RENEWAL_V3_1_FIT_PATH", default_path)
  if (!file.exists(fit_path)) {
    stop("No v3.1 fit bundle found at ", fit_path, ". This script does not sample.")
  }

  diagnostics <- summarise_v3_1_fit(readRDS(fit_path))
  table_dir <- here::here("03_Output/tables/renewal_v3_1")
  figure_dir <- here::here("03_Output/figures/renewal_v3_1")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  for (name in c(
    "time", "unit_susceptibility", "endpoints", "contribution",
    "reporting", "serology", "hmc_summary", "hmc_convergence",
    "worst_mixing"
  )) {
    utils::write.csv(diagnostics[[name]], file.path(
      table_dir, paste0("renewal_v3_1_", name, ".csv")
    ), row.names = FALSE)
  }
  plot_v3_1_diagnostics(
    diagnostics,
    file.path(figure_dir, "renewal_v3_1_minimal_spatial_diagnostics.pdf")
  )
  message("[Fortaleza] ", diagnostics$fortaleza_note)
  invisible(diagnostics)
}

if (sys.nframe() == 0) {
  run_v3_1_diagnostics()
}
