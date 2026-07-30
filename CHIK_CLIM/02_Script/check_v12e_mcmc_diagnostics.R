# ---------------------------------------------------------------------------
# check_v12e_mcmc_diagnostics.R
#
# MCMC diagnostics for fit_chik_ceara_stan_weekly.R.
#
# Run after fit_chik_ceara_stan_weekly.R has completed in the same R session.
# Required object:
#   - fit_v12e
#
# This script prints diagnostics and plots in the active R session only.
# ---------------------------------------------------------------------------

required_pkgs <- c("dplyr", "ggplot2", "posterior", "bayesplot")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(posterior)
  library(bayesplot)
})

if (!exists("fit_v12e")) {
  stop("Missing fit_v12e. Run source('02_Script/fit_chik_ceara_stan_weekly.R') first.")
}

key_pars <- c(
  "rho", "mu_log_beta", "sigma_season", "S0_frac", "I0", "phi_cases"
)

draws_array <- as.array(fit_v12e)
draws_df <- as_draws_df(draws_array)

# ---- Numeric summary -------------------------------------------------------
fit_summary <- if ("stanfit" %in% class(fit_v12e)) {
  as.data.frame(rstan::summary(fit_v12e)$summary)
} else {
  posterior::summarise_draws(draws_df) |>
    as.data.frame()
}

if (!"parameter" %in% names(fit_summary)) {
  fit_summary$parameter <- rownames(fit_summary)
}
rownames(fit_summary) <- NULL

# ---- HMC sampler diagnostics ----------------------------------------------
if ("stanfit" %in% class(fit_v12e)) {
  sampler_params <- rstan::get_sampler_params(fit_v12e, inc_warmup = FALSE)
  n_divergent <- sum(vapply(
    sampler_params,
    function(x) sum(x[, "divergent__"]),
    numeric(1)
  ))
  n_treedepth <- sum(vapply(
    sampler_params,
    function(x) sum(x[, "treedepth__"] >= 13),
    numeric(1)
  ))
  ebfmi <- vapply(
    sampler_params,
    function(x) {
      energy <- x[, "energy__"]
      mean(diff(energy)^2) / stats::var(energy)
    },
    numeric(1)
  )

  sampler_diag <- tibble::tibble(
    diagnostic = c("divergent_transitions", "max_treedepth_hits", paste0("ebfmi_chain_", seq_along(ebfmi))),
    value = c(n_divergent, n_treedepth, ebfmi)
  )

  message("Divergences: ", n_divergent)
  message("Max treedepth hits: ", n_treedepth)
  message("E-BFMI by chain: ", paste(round(ebfmi, 3), collapse = ", "))
  print(sampler_diag)
}

# ---- Trace and density plots ----------------------------------------------
available_key_pars <- intersect(key_pars, dimnames(draws_array)$parameters)

if (length(available_key_pars) > 0) {
  p_trace <- bayesplot::mcmc_trace(draws_array, pars = available_key_pars) +
    ggplot2::ggtitle("MCMC trace plots: key scalar parameters")
  print(p_trace)

  p_dens <- bayesplot::mcmc_dens_overlay(draws_array, pars = available_key_pars) +
    ggplot2::ggtitle("MCMC posterior density overlay by chain")
  print(p_dens)
}

# ---- Rhat and ESS overview -------------------------------------------------
if ("Rhat" %in% names(fit_summary)) {
  rhat_df <- fit_summary |>
    dplyr::filter(!is.na(.data$Rhat), is.finite(.data$Rhat))

  p_rhat <- ggplot(rhat_df, aes(x = .data$Rhat)) +
    geom_histogram(bins = 50, fill = "#3182bd", colour = "white") +
    geom_vline(xintercept = 1.01, linetype = "dashed", colour = "red") +
    labs(
      x = "Rhat",
      y = "Number of parameters",
      title = "Rhat distribution",
      subtitle = "Most parameters should be <= 1.01; investigate anything > 1.05"
    ) +
    theme_bw(base_size = 12)
  print(p_rhat)
}

ess_col <- intersect(c("n_eff", "ess_bulk"), names(fit_summary))
if (length(ess_col) > 0) {
  ess_col <- ess_col[[1]]
  ess_df <- fit_summary |>
    dplyr::filter(!is.na(.data[[ess_col]]), is.finite(.data[[ess_col]]))

  p_ess <- ggplot(ess_df, aes(x = .data[[ess_col]])) +
    geom_histogram(bins = 50, fill = "#31a354", colour = "white") +
    labs(
      x = ess_col,
      y = "Number of parameters",
      title = "Effective sample size distribution"
    ) +
    theme_bw(base_size = 12)
  print(p_ess)
}

# ---- Print worst parameters ------------------------------------------------
if ("Rhat" %in% names(fit_summary)) {
  message("\nWorst Rhat parameters:")
  print(
    fit_summary |>
      dplyr::filter(!is.na(.data$Rhat)) |>
      dplyr::arrange(dplyr::desc(.data$Rhat)) |>
      dplyr::select(parameter, mean, sd, Rhat, dplyr::any_of(c("n_eff", "ess_bulk", "ess_tail"))) |>
      utils::head(20)
  )
}

if (length(ess_col) > 0) {
  message("\nLowest ESS parameters:")
  print(
    fit_summary |>
      dplyr::filter(!is.na(.data[[ess_col]])) |>
      dplyr::arrange(.data[[ess_col]]) |>
      dplyr::select(parameter, mean, sd, Rhat, dplyr::any_of(c("n_eff", "ess_bulk", "ess_tail"))) |>
      utils::head(20)
  )
}


draws_df <- as.data.frame(fit_v12e)

get_med <- function(var) {
  cols <- grep(paste0("^", var, "\\["), names(draws_df), value = TRUE)
  apply(draws_df[, cols, drop = FALSE], 2, median)
}

local_inf <- get_med("local_infectious_frac_out")
import_inf <- get_med("import_frac_out")
pressure <- get_med("infectious_pressure_out")

plot(local_inf, type = "l", col = "blue", ylab = "infectious pressure")
lines(import_inf, col = "red")
lines(pressure, col = "black")
legend("topright", c("local I/N", "import", "total"), col = c("blue", "red", "black"), lty = 1)

reff <- get_med("Reff_t")
new_inf <- get_med("new_inf")

plot(reff, type = "l", col = "darkgreen", ylab = "Reff")
abline(h = 1, lty = 2)

