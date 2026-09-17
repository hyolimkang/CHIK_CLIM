# =============================================================================
# Fit pre-specified Ceará annual FOI models.
# M0: static annual FOI. M1: stationary AR(1) annual FOI departures.
# Q1-Q3: pre-specified constant overall-detection prior scenarios.
# =============================================================================

required_packages <- c("here", "rstan", "dplyr", "ggplot2", "patchwork", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(ggplot2); library(patchwork); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

script_project_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM R-project root.")
}
source(file.path(script_project_root(), "02_Script", "40_susceptibility", "prepare_dynamic_annual_foi_data.R"))
env_integer <- function(name, default) {
  value <- Sys.getenv(name, "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) stop(name, " must be a positive integer.")
  parsed
}

simulate_annual_path <- function(annual, log_mu, q, rho, sigma, dynamic, phi) {
  y <- nrow(annual); delta <- numeric(y)
  if (dynamic) {
    delta[1] <- rnorm(1, 0, sigma / sqrt(1 - rho^2))
    for (i in 2:y) delta[i] <- rho * delta[i - 1] + rnorm(1, 0, sigma)
  }
  lambda <- exp(log_mu + delta); attack <- -expm1(-lambda)
  S_start <- S_end <- U_start <- U_end <- latent <- expected <- numeric(y)
  S_start[1] <- annual$N_start[1]
  for (i in seq_len(y)) {
    latent[i] <- S_start[i] * attack[i]
    s_after <- S_start[i] - latent[i]; u_after <- U_start[i] + latent[i]
    frac_s <- s_after / annual$N_start[i]
    S_end[i] <- s_after + annual$births[i] - annual$deaths[i] * frac_s + annual$reconciliation[i] * frac_s
    U_end[i] <- u_after - annual$deaths[i] * (1 - frac_s) + annual$reconciliation[i] * (1 - frac_s)
    expected[i] <- q * latent[i]
    if (i < y) { S_start[i + 1] <- S_end[i]; U_start[i + 1] <- U_end[i] }
  }
  tibble(year = annual$year, lambda, attack_prob = attack, S_end_prop = S_end / annual$N_end,
         expected_cases = expected, predicted_cases = rnbinom(y, mu = expected + 1e-9, size = phi))
}

save_prior_predictive <- function(prepared, figure_path, n_draws = 1000L) {
  q1 <- filter(prepared$q_scenarios, scenario == "Q1_primary_product_anchor")
  set.seed(20260911)
  sims <- bind_rows(lapply(seq_len(n_draws), function(i) {
    rho <- max(-.94, min(.94, rnorm(1, 0, .45)))
    simulate_annual_path(
      prepared$annual,
      rnorm(1, prepared$stan_base_data$log_foi_prior_mean, prepared$stan_base_data$log_foi_prior_sd),
      plogis(rnorm(1, q1$q_logit_prior_mean, q1$q_logit_prior_sd)),
      rho, abs(rnorm(1, 0, .75)), TRUE, rgamma(1, 2, .1)
    ) |> mutate(draw = i)
  }))
  interval <- sims |> group_by(year) |> summarise(
    lambda_q025 = quantile(lambda, .025), lambda_mid = median(lambda), lambda_q975 = quantile(lambda, .975),
    cases_q025 = quantile(predicted_cases, .025), cases_mid = median(predicted_cases), cases_q975 = quantile(predicted_cases, .975),
    s_q025 = quantile(S_end_prop, .025), s_mid = median(S_end_prop), s_q975 = quantile(S_end_prop, .975), .groups = "drop")
  p1 <- ggplot(interval, aes(year, lambda_mid)) + geom_ribbon(aes(ymin = lambda_q025, ymax = lambda_q975), fill = "#56B4E9", alpha = .35) + geom_line() + scale_y_continuous(trans = "log10") + labs(title = "Prior predictive annual FOI", y = "Annual FOI (log scale)") + theme_classic()
  p2 <- ggplot(interval, aes(year, cases_mid)) + geom_ribbon(aes(ymin = cases_q025, ymax = cases_q975), fill = "#E69F00", alpha = .35) + geom_line() + labs(title = "Prior predictive reported cases", y = "Cases") + theme_classic()
  p3 <- ggplot(interval, aes(year, s_mid)) + geom_ribbon(aes(ymin = s_q025, ymax = s_q975), fill = "#009E73", alpha = .35) + geom_line() + labs(title = "Prior predictive susceptible proportion", y = "S/N") + theme_classic()
  ggsave(figure_path, p1 / p2 / p3, width = 190, height = 240, units = "mm", device = cairo_pdf)
}

fit_dynamic_annual_foi <- function() {
  prepared <- build_dynamic_annual_foi_data(write_outputs = TRUE)
  figure_dir <- project_path("03_Output", "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  save_prior_predictive(prepared, file.path(figure_dir, "dynamic_annual_foi_ceara_prior_predictive.pdf"))

  n_iter <- env_integer("DYNAMIC_ANNUAL_FOI_ITER", 2000L)
  n_warmup <- env_integer("DYNAMIC_ANNUAL_FOI_WARMUP", 1000L)
  n_chains <- env_integer("DYNAMIC_ANNUAL_FOI_CHAINS", 4L)
  if (n_warmup >= n_iter) stop("DYNAMIC_ANNUAL_FOI_WARMUP must be less than iterations.")
  selected <- strsplit(Sys.getenv("DYNAMIC_ANNUAL_FOI_MODELS", "M0_Q1,M1_Q1,M1_Q2,M1_Q3,M1_Q1_foi_sd1_5,M1_Q1_foi_sd2"), ",", fixed = TRUE)[[1]]
  configs <- tibble(
    fit_id = c("M0_Q1", "M1_Q1", "M1_Q2", "M1_Q3", "M1_Q1_foi_sd1_5", "M1_Q1_foi_sd2"),
    model = c("M0_static", "M1_dynamic_AR1", "M1_dynamic_AR1", "M1_dynamic_AR1", "M1_dynamic_AR1", "M1_dynamic_AR1"),
    use_dynamic = c(0L, 1L, 1L, 1L, 1L, 1L),
    scenario = c("Q1_primary_product_anchor", "Q1_primary_product_anchor", "Q2_broader_product_anchor", "Q3_lower_detection", "Q1_primary_product_anchor", "Q1_primary_product_anchor"),
    foi_sd_multiplier = c(1, 1, 1, 1, 1.5, 2)
  ) |> filter(fit_id %in% selected)
  if (!nrow(configs)) stop("No recognised fit IDs selected.")

  stan_file <- project_path("02_Script", "stan", "chik_dynamic_annual_foi.stan")
  model <- rstan::stan_model(stan_file)
  fit_dir <- project_path("02_Script", "stan")
  fit_manifest <- vector("list", nrow(configs))
  for (i in seq_len(nrow(configs))) {
    q_prior <- filter(prepared$q_scenarios, scenario == configs$scenario[i])
    stan_base <- prepared$stan_base_data
    stan_base$log_foi_prior_sd <- stan_base$log_foi_prior_sd * configs$foi_sd_multiplier[i]
    stan_data <- c(stan_base, list(
      q_logit_prior_mean = q_prior$q_logit_prior_mean,
      q_logit_prior_sd = q_prior$q_logit_prior_sd,
      use_dynamic = configs$use_dynamic[i]
    ))
    message("[fit] ", configs$fit_id[i], " | ", configs$scenario[i])
    started <- Sys.time()
    fit <- rstan::sampling(model, data = stan_data, chains = n_chains, iter = n_iter,
      warmup = n_warmup, seed = 20260912L + i, control = list(adapt_delta = .999, max_treedepth = 12), refresh = max(1L, floor(n_iter / 10)))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    bundle <- list(fit = fit, prepared = prepared, stan_data = stan_data,
      config = as.list(configs[i, ]), elapsed_seconds = elapsed)
    path <- file.path(fit_dir, paste0("dynamic_annual_foi_ceara_", configs$fit_id[i], ".rds"))
    saveRDS(bundle, path)
    fit_manifest[[i]] <- tibble(fit_id = configs$fit_id[i], model = configs$model[i], scenario = configs$scenario[i], elapsed_seconds = elapsed, path = path)
  }
  manifest <- bind_rows(fit_manifest)
  write_csv(manifest, project_path("03_Output", "tables", "dynamic_annual_foi", "dynamic_annual_foi_fit_manifest.csv"))
  invisible(manifest)
}

if (sys.nframe() == 0L) fit_dynamic_annual_foi()
