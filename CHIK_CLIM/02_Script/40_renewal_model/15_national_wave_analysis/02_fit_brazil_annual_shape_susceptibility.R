# State-specific batch extension of the validated annual SHAPE FOI model.
#
# The scientific model is the existing chik_dynamic_annual_foi_shape_v2.stan:
# one long-term FOI scale, independently varying normalised annual FOI shape,
# conditional annual case allocation, and explicit S/U demographic accounting.
# The only extension is rebuilding the same input contract for each UF.

required_shape_packages <- c("here", "rstan", "posterior")
missing_shape_packages <- required_shape_packages[
  !vapply(required_shape_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_shape_packages)) stop("Missing package(s): ", paste(missing_shape_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(posterior) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

shape_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}
source(file.path(shape_root(), "02_Script", "40_renewal_model", "15_national_wave_analysis", "00_national_wave_analysis_helpers.R"))

shape_settings <- list(
  years = 2015:2025,
  sigma_year_prior_scale = .75,
  iter = as.integer(Sys.getenv("NATIONAL_ANNUAL_SHAPE_ITER", "2000")),
  warmup = as.integer(Sys.getenv("NATIONAL_ANNUAL_SHAPE_WARMUP", "1000")),
  chains = as.integer(Sys.getenv("NATIONAL_ANNUAL_SHAPE_CHAINS", "4")),
  adapt_delta = .99,
  max_treedepth = 12L,
  # Portability is evaluated first with contrasting large, northern, and
  # northeast UFs.  No state-specific tuning follows.
  portability_states = c("CE", "AM", "SP")
)
if (shape_settings$warmup >= shape_settings$iter) stop("Warmup must be smaller than iterations.")

make_shape_input <- function(state, state_week, demography, foi_ensemble, settings) {
  annual_demography <- demography |>
    filter(state == !!state, year %in% c(settings$years, max(settings$years) + 1L)) |>
    arrange(year)
  if (nrow(annual_demography) != length(settings$years) + 1L) stop("Incomplete demographic sequence for ", state)
  boundary_dates <- as.Date(sprintf("%d-01-01", c(settings$years, max(settings$years) + 1L)))
  boundary_population <- interpolate_uf_population(annual_demography, boundary_dates)
  annual_cases <- state_week |>
    filter(state == !!state) |>
    mutate(year = year(week_start)) |>
    filter(year %in% settings$years) |>
    group_by(year) |>
    summarise(observed_cases = sum(reported_cases), .groups = "drop")
  annual <- tibble(year = settings$years, N_start = boundary_population[-length(boundary_population)],
                   N_end = boundary_population[-1]) |>
    left_join(annual_demography |> filter(year %in% settings$years) |>
                select(year, births, deaths), by = "year") |>
    left_join(annual_cases, by = "year") |>
    mutate(observed_cases = coalesce(as.integer(observed_cases), 0L),
           reconciliation = N_end - N_start - births + deaths)
  anchor <- foi_ensemble |>
    filter(uf == !!state) |>
    pull(foi_equiv)
  if (length(anchor) != 100L || any(!is.finite(anchor) | anchor <= 0)) stop("Invalid 100-member FOI ensemble for ", state)
  if (anyNA(annual) || any(annual$N_start <= 0 | annual$N_end <= 0 | annual$births < 0 | annual$deaths < 0)) {
    stop("Invalid annual state input for ", state)
  }
  list(
    annual = annual,
    stan_data = list(
      Y = nrow(annual), C = annual$observed_cases, C_total = sum(annual$observed_cases),
      N_start = annual$N_start, N_end = annual$N_end, births = annual$births,
      deaths = annual$deaths, reconciliation = annual$reconciliation,
      log_foi_prior_mean = mean(log(anchor)), log_foi_prior_sd = sd(log(anchor)),
      sigma_year_prior_scale = settings$sigma_year_prior_scale
    ),
    anchor = tibble(state = state, foi_equiv_log_mean = mean(log(anchor)), foi_equiv_log_sd = sd(log(anchor)),
                    foi_equiv_median = median(anchor), foi_equiv_q025 = quantile(anchor, .025), foi_equiv_q975 = quantile(anchor, .975))
  )
}

fit_one_state_shape <- function(model, state_input, state, settings) {
  started <- Sys.time()
  fit <- tryCatch(
    sampling(
      model, data = state_input$stan_data, chains = settings$chains,
      iter = settings$iter, warmup = settings$warmup,
      seed = 20260930L + match(state, sort(unique(read_national_demography()$state))),
      control = list(adapt_delta = settings$adapt_delta, max_treedepth = settings$max_treedepth), refresh = 0
    ),
    error = function(e) e
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (inherits(fit, "error")) return(list(status = "error", reason = conditionMessage(fit), elapsed_seconds = elapsed))
  summary <- as.data.frame(rstan::summary(fit)$summary) |> tibble::rownames_to_column("parameter")
  hmc_summary <- summary |> filter(grepl("^(log_lambda_bar|log_sigma_year|z_year|log_kappa)", parameter))
  hmc <- hmc_gate(hmc_summary, rstan::get_sampler_params(fit, inc_warmup = FALSE), settings$max_treedepth, fit,
                  pars = c("log_lambda_bar", "log_sigma_year", "z_year", "log_kappa"))
  list(status = "fitted", fit = fit, hmc = hmc, elapsed_seconds = elapsed)
}

summarise_shape_fit <- function(state, state_input, fitted) {
  draws <- rstan::extract(fitted$fit, pars = c("lambda", "attack_prob", "X", "S_start", "S_end", "q_implied"), permuted = TRUE)
  annual_summary <- bind_rows(lapply(seq_len(nrow(state_input$annual)), function(y) {
    bind_cols(
      tibble(state = state, year = state_input$annual$year[y], observed_cases = state_input$annual$observed_cases[y]),
      summarise_draws(draws$lambda[, y], "lambda"),
      summarise_draws(draws$attack_prob[, y], "attack_prob"),
      summarise_draws(draws$X[, y], "infections"),
      summarise_draws(draws$S_start[, y] / state_input$annual$N_start[y], "S_start_prop"),
      summarise_draws(draws$S_end[, y] / state_input$annual$N_end[y], "S_end_prop")
    )
  }))
  list(
    annual_summary = annual_summary,
    q_summary = bind_cols(tibble(state = state), summarise_draws(draws$q_implied, "q_implied"),
                           tibble(q_implied_gt_one_draw_fraction = mean(draws$q_implied > 1))),
    posterior_draws = draws
  )
}

run_national_annual_shape <- function() {
  paths <- ensure_national_output_dirs()
  state_week <- read_national_state_week()
  demography <- read_national_demography()
  foi_ensemble <- read_csv(national_path("01_Data", "brazil_uf_foi_ensemble.csv"), show_col_types = FALSE)
  states <- sort(unique(demography$state))
  requested <- Sys.getenv("NATIONAL_ANNUAL_SHAPE_STATES", "")
  if (nzchar(requested)) states <- strsplit(requested, ",", fixed = TRUE)[[1]]
  model <- stan_model(national_path("02_Script", "stan", "chik_dynamic_annual_foi_shape_v2.stan"))
  fits <- list(); audit <- list(); annual_out <- list(); q_out <- list(); anchors <- list()

  # Fixed portability gate before the remaining state batch.
  run_order <- unique(c(intersect(shape_settings$portability_states, states), states))
  for (state in run_order) {
    message("[annual-SHAPE] ", state)
    input <- make_shape_input(state, state_week, demography, foi_ensemble, shape_settings)
    result <- fit_one_state_shape(model, input, state, shape_settings)
    audit[[state]] <- if (identical(result$status, "fitted")) {
      bind_cols(tibble(state = state, fit_status = "fitted", elapsed_seconds = result$elapsed_seconds), result$hmc)
    } else {
      tibble(state = state, fit_status = "error", elapsed_seconds = result$elapsed_seconds, fit_reason = result$reason,
             divergences = NA_integer_, treedepth_hits = NA_integer_, max_rhat = NA_real_, min_bulk_ess = NA_real_,
             min_tail_ess = NA_real_, min_bfmi = NA_real_, bfmi_by_chain = NA_character_, hmc_pass = FALSE)
    }
    anchors[[state]] <- input$anchor
    if (identical(result$status, "fitted")) {
      fits[[state]] <- list(fit = result$fit, input = input, hmc = result$hmc, elapsed_seconds = result$elapsed_seconds)
      extracted <- summarise_shape_fit(state, input, result)
      annual_out[[state]] <- extracted$annual_summary
      q_out[[state]] <- extracted$q_summary
    }
    # Do not silently continue to all states if the pre-specified contrasting
    # portability set fails the same HMC gate.
    completed_portability <- intersect(shape_settings$portability_states, names(audit))
    if (length(completed_portability) == length(intersect(shape_settings$portability_states, states)) &&
        any(!vapply(audit[completed_portability], function(x) x$hmc_pass, logical(1)))) {
      break
    }
  }
  hmc_table <- bind_rows(audit)
  write_csv(hmc_table, file.path(paths$table, "brazil_chik_annual_shape_hmc_gate.csv"))
  write_csv(bind_rows(anchors), file.path(paths$table, "brazil_chik_annual_shape_foi_anchor.csv"))
  if (length(annual_out)) write_csv(bind_rows(annual_out), file.path(paths$table, "brazil_chik_annual_susceptibility_summary.csv"))
  if (length(q_out)) write_csv(bind_rows(q_out), file.path(paths$table, "brazil_chik_annual_shape_q_implied.csv"))
  saveRDS(list(settings = shape_settings, fits = fits, hmc = hmc_table,
               demographic_source = "IBGE 2024 population projection: NASC_T and OBT_T; annual reconciliation is exact population closure.",
               initial_condition = "S_start[2015] = N_start[2015]; U_start[2015] = 0"),
          file.path(paths$fit, "brazil_chik_annual_shape_fits.rds"))
  invisible(hmc_table)
}

if (sys.nframe() == 0L) run_national_annual_shape()
