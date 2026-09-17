# Bottom-up 5-km susceptibility reconstruction for Brazil.
#
# This is post-processing only.  It does not refit Stan, alter the frozen wave
# census, use q, or use early-Re.  The absolute spatial FOI scale comes only
# from the coherent 100-member allfoi_s1 grid ensemble.  The annual SHAPE
# posterior contributes only the within-UF temporal multiplier lambda/lambda_bar.

required_grid_packages <- c("here", "dplyr", "readr", "tibble", "rstan", "matrixStats", "sf")
missing_grid_packages <- required_grid_packages[!vapply(required_grid_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_grid_packages)) stop("Missing package(s): ", paste(missing_grid_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(rstan)
  library(matrixStats)
})

grid_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

grid_settings <- list(
  years = 2015:2025,
  selected_years = c(2015L, 2017L, 2020L, 2022L, 2024L, 2025L),
  # Each RF surface is used once.  A common realization index is paired with
  # one posterior draw from every independently fitted UF, preserving the
  # 11-year covariance within each UF trajectory.
  n_realizations = as.integer(Sys.getenv("GRID_SUSC_N_REALIZATIONS", "100")),
  quantiles = c(.025, .5, .975),
  population_tolerance_relative = 1e-8,
  aggregation_tolerance = 1e-8,
  seed = 20261002L
)

grid_paths <- function(root = grid_root()) {
  output_tag <- Sys.getenv("GRID_SUSC_OUTPUT_TAG", "")
  table_root <- file.path(root, "03_Output", "tables", "national_grid_susceptibility")
  figure_root <- file.path(root, "03_Output", "figures", "national_grid_susceptibility")
  if (nzchar(output_tag)) {
    table_root <- file.path(table_root, output_tag)
    figure_root <- file.path(figure_root, output_tag)
  }
  list(
    root = root,
    table = table_root,
    figure = figure_root,
    grid = file.path(table_root, "grid_selected_years"),
    fit = file.path(root, "02_Script", "stan", "national_wave_analysis", "brazil_chik_annual_shape_fits.rds"),
    uf_ensemble = file.path(root, "01_Data", "brazil_uf_foi_ensemble.csv")
  )
}

ensure_grid_paths <- function(paths) {
  invisible(lapply(paths[c("table", "figure", "grid")], dir.create, recursive = TRUE, showWarnings = FALSE))
  paths
}

# The established aggregation script is sourced for the validated Brazil-grid
# reader, official UF assignment, and UF lookup.  Its main block is guarded,
# so sourcing does not regenerate or overwrite the historical aggregation.
load_grid_helpers <- function(root) {
  source(file.path(root, "02_Script", "90_development_archive", "09_misc_archive", "15_build_brazil_uf_foi_ensemble.R"))
  invisible(TRUE)
}

summary_draws <- function(x, prefix) {
  tibble(
    !!paste0(prefix, "_q025") := as.numeric(quantile(x, .025)),
    !!paste0(prefix, "_median") := median(x),
    !!paste0(prefix, "_q975") := as.numeric(quantile(x, .975))
  )
}

weighted_distribution <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (!length(x)) return(rep(NA_real_, 9L))
  order_x <- order(x)
  x <- x[order_x]
  w <- w[order_x]
  cumulative_weight <- cumsum(w) / sum(w)
  weighted_quantile <- function(p) x[which(cumulative_weight >= p)[1L]]
  mean_x <- sum(w * x) / sum(w)
  setNames(
    unname(c(
      weighted_quantile(.10), weighted_quantile(.25),
      weighted_quantile(.50), weighted_quantile(.75),
      weighted_quantile(.90),
      sqrt(sum(w * (x - mean_x)^2) / sum(w)),
      weighted_quantile(.75) - weighted_quantile(.25),
      sum(w[x < .25]) / sum(w),
      sum(w[x < .50]) / sum(w)
    )),
    c("q10", "q25", "q50", "q75", "q90", "weighted_sd", "weighted_iqr", "prop_lt_025", "prop_lt_050")
  )
}

# Keep the required <0.75 and >=0.90 fractions separate for readable output.
weighted_heterogeneity <- function(s_prop, population) {
  base <- weighted_distribution(s_prop, population)
  keep <- is.finite(s_prop) & is.finite(population) & population > 0
  population <- population[keep]
  s_prop <- s_prop[keep]
  c(
    within_state_S_q10 = base[["q10"]], within_state_S_q25 = base[["q25"]],
    within_state_S_q50 = base[["q50"]], within_state_S_q75 = base[["q75"]],
    within_state_S_q90 = base[["q90"]], within_state_S_weighted_sd = base[["weighted_sd"]],
    within_state_S_weighted_iqr = base[["weighted_iqr"]],
    prop_population_S_lt_025 = base[["prop_lt_025"]],
    prop_population_S_lt_050 = base[["prop_lt_050"]],
    prop_population_S_lt_075 = sum(population[s_prop < .75]) / sum(population),
    prop_population_S_ge_090 = sum(population[s_prop >= .90]) / sum(population)
  )
}

summarise_heterogeneity_mc <- function(values) {
  result <- list()
  for (metric in colnames(values)) {
    result[[metric]] <- median(values[, metric])
    result[[paste0(metric, "_q025")]] <- as.numeric(quantile(values[, metric], .025))
    result[[paste0(metric, "_q975")]] <- as.numeric(quantile(values[, metric], .975))
  }
  as_tibble(result)
}

append_grid_quantiles <- function(grid_summary, value, prefix) {
  q <- matrixStats::rowQuantiles(value, probs = grid_settings$quantiles, na.rm = TRUE)
  grid_summary[[paste0(prefix, "_q025")]] <- q[, 1L]
  grid_summary[[paste0(prefix, "_median")]] <- q[, 2L]
  grid_summary[[paste0(prefix, "_q975")]] <- q[, 3L]
  grid_summary
}

draw_pairing <- function(fit_bundle, settings) {
  valid_draw_counts <- vapply(fit_bundle$fits, function(entry) {
    length(rstan::extract(entry$fit, pars = "lambda_bar", permuted = TRUE)$lambda_bar)
  }, integer(1))
  n_available <- min(valid_draw_counts)
  if (settings$n_realizations > 100L) stop("Only 100 coherent RF surfaces are available.")
  if (n_available < settings$n_realizations) stop("Too few SHAPE posterior draws for requested pairing.")
  set.seed(settings$seed)
  tibble(
    realization = seq_len(settings$n_realizations),
    ensemble_id = seq_len(settings$n_realizations),
    posterior_draw = sample.int(n_available, settings$n_realizations, replace = FALSE)
  )
}

validate_temporal_weights <- function(lambda, lambda_bar, state) {
  weight <- sweep(lambda, 1L, lambda_bar, "/")
  deviation <- max(abs(rowMeans(weight) - 1))
  if (!is.finite(deviation) || deviation > 1e-10) {
    stop("Annual temporal weights do not have arithmetic mean one for ", state, ".")
  }
  list(weight = weight, max_mean_deviation = deviation)
}

state_foi_qc <- function(state_grid, member_columns, state, existing_ensemble) {
  foi <- as.matrix(state_grid[, member_columns, drop = FALSE])
  population <- state_grid$tot
  denominator <- sum(population)
  pwmean <- colSums(sweep(foi, 1L, population, "*")) / denominator
  p_mean <- colSums(sweep(-expm1(-foi), 1L, population, "*")) / denominator
  equivalent <- -log1p(-p_mean)
  calculated <- tibble(
    state = state, ensemble_id = seq_along(member_columns),
    foi_pwmean_grid = pwmean, foi_equiv_grid = equivalent
  )
  expected <- existing_ensemble |>
    dplyr::filter(uf == state) |>
    transmute(state = uf, ensemble_id, foi_pwmean_existing = foi_pwmean, foi_equiv_existing = foi_equiv)
  calculated |>
    left_join(expected, by = c("state", "ensemble_id")) |>
    mutate(
      pwmean_difference = foi_pwmean_grid - foi_pwmean_existing,
      equivalent_difference = foi_equiv_grid - foi_equiv_existing,
      pwmean_abs_difference = abs(pwmean_difference),
      equivalent_abs_difference = abs(equivalent_difference)
    )
}

summarise_state_year <- function(state, year, annual, S, U, X, foi, N_end, old_S_end, heterogeneity) {
  population_grid <- N_end
  s_pop <- colSums(S) / sum(population_grid)
  denominator_foi <- colSums(sweep(foi, 1L, population_grid, "*"))
  s_foi <- colSums(S * foi) / denominator_foi
  infection_total <- colSums(X)
  old_s <- as.numeric(old_S_end[, 1L]) / sum(population_grid)
  bind_cols(
    tibble(state = state, year = year, observed_cases = annual$observed_cases),
    summary_draws(s_pop, "S_pop"),
    summary_draws(s_foi, "S_foiweighted"),
    summary_draws(infection_total, "latent_infections_total"),
    summary_draws(old_s, "S_old_homogeneous"),
    summary_draws(s_pop - old_s, "S_pop_minus_old_homogeneous"),
    summarise_heterogeneity_mc(heterogeneity)
  )
}

reconstruct_state_grid <- function(state, grid, member_columns, assignment, fit_bundle, pairing, paths, settings, existing_ensemble) {
  entry <- fit_bundle$fits[[state]]
  if (is.null(entry)) stop("Missing annual SHAPE fit for ", state)
  if (!isTRUE(entry$hmc$hmc_pass)) {
    return(list(
      state = state, status = "annual_shape_hmc_failed", n_cells = 0L,
      annual = tibble(state = state, year = settings$years, reconstruction_status = "annual_shape_hmc_failed"),
      qc = tibble(state = state, aggregation_status = "not_run_annual_hmc_failed")
    ))
  }
  index <- which(assignment$uf_abbr == state)
  state_grid <- grid[index, , drop = FALSE]
  state_grid <- state_grid[is.finite(state_grid$tot) & state_grid$tot > 0, , drop = FALSE]
  if (!nrow(state_grid)) stop("No populated valid FOI cells for ", state)
  annual <- entry$input$annual
  if (!identical(as.integer(annual$year), settings$years)) stop("Unexpected annual sequence for ", state)
  draws <- rstan::extract(entry$fit, pars = c("lambda", "lambda_bar", "S_end"), permuted = TRUE)
  temporal <- validate_temporal_weights(draws$lambda, draws$lambda_bar, state)
  selected_draws <- pairing$posterior_draw
  temporal_weight <- temporal$weight[selected_draws, , drop = FALSE]
  old_S_end <- draws$S_end[selected_draws, , drop = FALSE]
  if (any(!is.finite(temporal_weight) | temporal_weight <= 0)) stop("Invalid temporal weight for ", state)

  population_share <- state_grid$tot / sum(state_grid$tot)
  foi <- as.matrix(state_grid[, member_columns, drop = FALSE])
  if (ncol(foi) != nrow(pairing)) stop("One RF ensemble surface is required per realization.")
  qc <- state_foi_qc(state_grid, member_columns, state, existing_ensemble) |>
    mutate(max_temporal_weight_mean_deviation = temporal$max_mean_deviation)
  if (max(qc$pwmean_abs_difference, qc$equivalent_abs_difference, na.rm = TRUE) > settings$aggregation_tolerance) {
    stop("5-km FOI aggregation no longer reproduces the validated UF ensemble for ", state)
  }

  n_cell <- nrow(state_grid)
  n_mc <- nrow(pairing)
  S <- tcrossprod(population_share * annual$N_start[1L], rep(1, n_mc))
  U <- matrix(0, nrow = n_cell, ncol = n_mc)
  cumulative_X <- matrix(0, nrow = n_cell, ncol = n_mc)
  grid_summary <- tibble(
    cell_id = state_grid$cell_id,
    x = state_grid$x,
    y = state_grid$y,
    state = state,
    pop_share = population_share
  )
  grid_summary <- append_grid_quantiles(grid_summary, foi, "longterm_foi")
  annual_rows <- vector("list", length(settings$years))

  for (year_index in seq_along(settings$years)) {
    year <- settings$years[year_index]
    lambda <- sweep(foi, 2L, temporal_weight[, year_index], "*")
    attack_probability <- -expm1(-lambda)
    X <- S * attack_probability
    if (any(X < -1e-8 | X - S > 1e-8)) stop("Invalid grid infections in ", state, " ", year)
    cumulative_X <- cumulative_X + X
    S_after <- S - X
    U_after <- U + X
    composition_S <- S_after / (S_after + U_after)
    birth_grid <- population_share * annual$births[year_index]
    death_grid <- population_share * annual$deaths[year_index]
    reconciliation_grid <- population_share * annual$reconciliation[year_index]
    S <- sweep(S_after, 1L, birth_grid, "+") -
      sweep(composition_S, 1L, death_grid, "*") +
      sweep(composition_S, 1L, reconciliation_grid, "*")
    U <- U_after - sweep(1 - composition_S, 1L, death_grid, "*") +
      sweep(1 - composition_S, 1L, reconciliation_grid, "*")
    N_end <- population_share * annual$N_end[year_index]
    max_population_error <- max(abs(sweep(S + U, 1L, N_end, "-")))
    if (max_population_error > settings$population_tolerance_relative * max(annual$N_end[year_index], 1)) {
      stop("Grid demographic accounting error in ", state, " ", year, ": ", signif(max_population_error, 5))
    }
    if (any(S < -1e-7 | U < -1e-7)) stop("Negative grid stock in ", state, " ", year)

    heterogeneity <- t(vapply(seq_len(n_mc), function(m) {
      weighted_heterogeneity(S[, m] / N_end, N_end)
    }, numeric(11L)))
    annual_rows[[year_index]] <- summarise_state_year(
      state, year, annual[year_index, , drop = FALSE], S, U, X, foi, N_end,
      old_S_end[, year_index, drop = FALSE], heterogeneity
    ) |>
      mutate(
        reconstruction_status = "ok",
        population_accounting_max_abs_error = max_population_error,
        annual_weight_mean_deviation = temporal$max_mean_deviation
      )

    if (year %in% settings$selected_years) {
      suffix <- as.character(year)
      grid_summary <- append_grid_quantiles(grid_summary, lambda, paste0("annual_lambda_", suffix))
      grid_summary <- append_grid_quantiles(grid_summary, attack_probability, paste0("attack_probability_", suffix))
      grid_summary <- append_grid_quantiles(grid_summary, X, paste0("latent_infections_", suffix))
      grid_summary <- append_grid_quantiles(grid_summary, S, paste0("susceptible_count_", suffix))
      grid_summary <- append_grid_quantiles(grid_summary, sweep(S, 1L, N_end, "/"), paste0("susceptible_fraction_", suffix))
      grid_summary <- append_grid_quantiles(grid_summary, sweep(U, 1L, N_end, "/"), paste0("immune_fraction_", suffix))
    }
    if (year == max(settings$years)) {
      baseline_population <- population_share * annual$N_start[1L]
      grid_summary <- append_grid_quantiles(
        grid_summary, sweep(cumulative_X, 1L, baseline_population, "/"),
        "cumulative_infected_fraction_2025"
      )
    }
  }
  grid_file <- file.path(paths$grid, paste0("brazil_chik_grid_susceptibility_", state, "_selected_years.rds"))
  saveRDS(grid_summary, grid_file, compress = "xz")
  list(
    state = state, status = "ok", n_cells = n_cell, grid_file = grid_file,
    annual = bind_rows(annual_rows), qc = qc
  )
}

run_brazil_grid_susceptibility <- function() {
  root <- grid_root()
  paths <- ensure_grid_paths(grid_paths(root))
  if (!file.exists(paths$fit) || !file.exists(paths$uf_ensemble)) {
    stop("Required annual SHAPE fit bundle or UF FOI ensemble is missing.")
  }
  load_grid_helpers(root)
  fit_bundle <- readRDS(paths$fit)
  pairing <- draw_pairing(fit_bundle, grid_settings)
  existing_ensemble <- read_csv(paths$uf_ensemble, show_col_types = FALSE)
  allfoi_input <- read_allfoi_brazil()
  boundary_input <- make_uf_boundaries()
  assignment_input <- assign_uf(allfoi_input$brazil, boundary_input$boundaries)
  member_columns <- allfoi_input$member_columns
  valid_member <- rowSums(is.finite(as.matrix(allfoi_input$brazil[, member_columns, drop = FALSE]))) == length(member_columns)
  valid_population <- is.finite(allfoi_input$brazil$tot) & allfoi_input$brazil$tot > 0
  keep <- valid_member & valid_population
  grid <- allfoi_input$brazil[keep, c("x", "y", "tot", member_columns), drop = FALSE]
  assignment <- assignment_input$assignment[keep, , drop = FALSE]
  grid$cell_id <- assignment$cell_id
  states <- sort(unique(assignment$uf_abbr))
  # A state filter supports inexpensive engineering smoke tests.  The default
  # is all UFs and is the only setting used for the scientific output.
  requested_states <- trimws(strsplit(Sys.getenv("GRID_SUSC_STATES", ""), ",", fixed = TRUE)[[1L]])
  requested_states <- requested_states[nzchar(requested_states)]
  if (length(requested_states)) {
    unknown_states <- setdiff(requested_states, states)
    if (length(unknown_states)) stop("Unknown requested state(s): ", paste(unknown_states, collapse = ", "))
    states <- intersect(states, requested_states)
  }
  results <- vector("list", length(states))
  names(results) <- states
  for (state in states) {
    message("[grid susceptibility] ", state)
    results[[state]] <- reconstruct_state_grid(
      state, grid, member_columns, assignment, fit_bundle, pairing, paths,
      grid_settings, existing_ensemble
    )
    gc(verbose = FALSE)
  }
  annual <- bind_rows(lapply(results, `[[`, "annual"))
  qc <- bind_rows(lapply(results, `[[`, "qc"))
  manifest <- bind_rows(lapply(results, function(x) {
    tibble(state = x$state, reconstruction_status = x$status, n_cells = x$n_cells,
           grid_file = if (is.null(x$grid_file)) NA_character_ else x$grid_file)
  }))
  qc_summary <- qc |>
    summarise(
      max_pwmean_abs_difference = max(pwmean_abs_difference, na.rm = TRUE),
      max_equivalent_abs_difference = max(equivalent_abs_difference, na.rm = TRUE),
      max_temporal_weight_mean_deviation = max(max_temporal_weight_mean_deviation, na.rm = TRUE)
    )
  if (any(!is.finite(unlist(qc_summary)))) stop("FOI QC summary is invalid.")
  write_csv(pairing, file.path(paths$table, "brazil_chik_grid_susceptibility_mc_pairing.csv"))
  write_csv(qc, file.path(paths$table, "brazil_chik_grid_foi_aggregation_qc.csv"))
  write_csv(qc_summary, file.path(paths$table, "brazil_chik_grid_foi_aggregation_qc_summary.csv"))
  write_csv(manifest, file.path(paths$table, "brazil_chik_grid_susceptibility_selected_years_manifest.csv"))
  write_csv(annual, file.path(paths$table, "brazil_chik_uf_bottomup_susceptibility.csv"))
  # Keep the requested heterogeneity output separate while retaining all
  # uncertainty columns in the principal UF table.
  heterogeneity_columns <- unique(c("state", "year", grep("^(within_state_S_|prop_population_S_)", names(annual), value = TRUE)))
  write_csv(annual |> select(all_of(heterogeneity_columns)), file.path(paths$table, "brazil_chik_uf_susceptibility_heterogeneity.csv"))
  write_csv(
    annual |> select(state, year, observed_cases, starts_with("latent_infections_total"), reconstruction_status),
    file.path(paths$table, "brazil_chik_bottomup_annual_infection_timing.csv")
  )
  failed_states <- manifest |>
    dplyr::filter(reconstruction_status != "ok") |>
    pull(state)
  hmc_gate_note <- if (length(failed_states)) {
    paste0("Excluded from numerical reconstruction because their annual SHAPE HMC gates failed: ",
           paste(failed_states, collapse = ", "))
  } else {
    "All reconstructed annual SHAPE fits passed their HMC gates."
  }
  write_csv(tibble(
    initial_condition = "S_grid_start_2015 = N_grid_start_2015; U_grid_start_2015 = 0",
    spatial_scale = "allfoi_s1 cell-specific foi1...foi100; no UF lambda_bar is used as a grid scale",
    temporal_component = "Annual SHAPE posterior lambda/lambda_bar, with its within-UF annual covariance retained",
    population = "No annual 5-km population input was used; baseline tot shares allocate official UF annual population and demographic flows",
    ensemble_pairing = "100 coherent national RF surfaces, one per realization, paired with one common-index posterior draw per independently fitted UF",
    q_use = "diagnostic only; not used here",
    early_re_use = "none",
    hmc_gate = hmc_gate_note
  ), file.path(paths$table, "brazil_chik_grid_susceptibility_metadata.csv"))
  invisible(list(annual = annual, qc = qc, manifest = manifest))
}

if (sys.nframe() == 0L) run_brazil_grid_susceptibility()
