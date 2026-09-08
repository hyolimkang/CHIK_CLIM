# Prepare validated official-macroregion 8-unit inputs for v3.2.
required_packages <- c("here", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr) })
source(here::here("02_Script/40_renewal_model/05_v3_0_spatial/prepare_renewal_v3_0_spatial_data.R"))

v3_2_units <- c("Fortaleza", "Fortaleza_remainder", "Norte_Sobral",
                "Cariri_remainder", "Juazeiro_do_Norte",
                "Sertao_Central_remainder", "Quixada",
                "Litoral_Leste_Jaguaribe")

aggregate_v3_2 <- function(x, unit_index) {
  out <- rowsum(as.matrix(x), group = unit_index, reorder = FALSE)
  out[as.character(seq_along(v3_2_units)), , drop = FALSE]
}

derive_q_prior <- function(n_draws = 1000000L, seed = 32052018L) {
  set.seed(seed)
  q_draws <- stats::rbeta(n_draws, 30, 28) * stats::rbeta(n_draws, 20, 60)
  mu <- mean(q_draws); va <- stats::var(q_draws)
  concentration <- mu * (1 - mu) / va - 1
  list(a = mu * concentration, b = (1 - mu) * concentration,
       mean = mu, sd = sqrt(va), n_draws = n_draws, seed = seed)
}

prepare_renewal_v3_2_coarse_regional_data <- function() {
  base <- prepare_renewal_v3_0_spatial_data(debug_subset = FALSE, min_total_cases = NULL)
  map_path <- here::here("01_Data/reference/ceara_municipality_final_model_unit_2018.csv")
  mapping <- utils::read.csv(map_path, stringsAsFactors = FALSE, colClasses = "character",
                             check.names = FALSE)
  needed <- c("muni6", "municipality_name", "health_region_22",
              "official_macroregion_5", "final_model_unit", "source_reference")
  if (!identical(names(mapping), needed) || nrow(mapping) != 184L ||
      anyDuplicated(mapping$muni6) || anyNA(mapping[needed])) {
    stop("Official regional mapping is incomplete, duplicated, or has an unexpected schema")
  }
  if (!setequal(mapping$muni6, base$municipality_ids) ||
      !setequal(unique(mapping$final_model_unit), v3_2_units)) {
    stop("Official mapping does not give a one-to-one complete assignment to the eight v3.2 units")
  }
  mapping$unit_index <- match(mapping$final_model_unit, v3_2_units)
  mapping <- mapping[match(base$municipality_ids, mapping$muni6), ]
  if (anyNA(mapping$unit_index)) stop("A municipality lacks final model-unit assignment")

  C <- aggregate_v3_2(base$stan_data$C, mapping$unit_index)
  N_start <- aggregate_v3_2(base$stan_data$N_start, mapping$unit_index)
  N_end <- aggregate_v3_2(base$stan_data$N_end, mapping$unit_index)
  births <- aggregate_v3_2(base$stan_data$births, mapping$unit_index)
  deaths <- aggregate_v3_2(base$stan_data$deaths, mapping$unit_index)
  if (any(C < 0) || any(abs(C - round(C)) > 1e-8) ||
      any(births < 0) || any(deaths < 0)) stop("Invalid aggregated cases or demographic flows")
  if (max(abs(N_end[, -ncol(N_end)] - N_start[, -1])) > 1e-6)
    stop("Regional N_end/N_start sequence is discontinuous")
  accounting_error <- check_demographic_accounting(N_start, N_end, births, deaths)
  if (!identical(colSums(C), colSums(base$stan_data$C)) ||
      max(abs(colSums(N_start) - colSums(base$stan_data$N_start))) > 1e-6)
    stop("Regional aggregation does not reproduce statewide cases/population")

  serology <- base$serology
  sero_unit <- mapping$unit_index[match(serology$muni6, mapping$muni6)]
  if (!identical(as.integer(sero_unit), c(5L, 7L)))
    stop("Juazeiro and Quixada serology failed direct final-unit mapping")
  q_prior <- derive_q_prior()

  table_dir <- here::here("03_Output/tables/renewal_v3_2_coarse_regional")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  assignment_path <- file.path(table_dir, "renewal_v3_2_coarse_regional_assignment.csv")
  utils::write.csv(mapping, assignment_path, row.names = FALSE)

  unit_summary <- data.frame(
    unit_index = seq_along(v3_2_units), unit = v3_2_units,
    municipalities = as.integer(table(factor(mapping$unit_index, levels = 1:8))),
    mean_population = rowMeans(N_start),
    population_share = rowMeans(N_start) / sum(rowMeans(N_start)),
    cumulative_cases = rowSums(C),
    case_share = rowSums(C) / sum(C)
  )
  stan_data <- list(
    K = 8L, N = base$stan_data$N, G = base$stan_data$G,
    seed_weeks = base$stan_data$seed_weeks,
    C = unname(array(as.integer(C), dim = dim(C))), w = base$stan_data$w,
    seed_profile = rep(1 / base$stan_data$seed_weeks, base$stan_data$seed_weeks),
    N_start = unname(N_start), N_end = unname(N_end),
    births = unname(births), deaths = unname(deaths),
    q_prior_a = q_prior$a, q_prior_b = q_prior$b,
    J = nrow(serology), sero_unit = as.integer(sero_unit),
    sero_window_start = base$stan_data$sero_window_start,
    sero_window_end = base$stan_data$sero_window_end,
    sero_positive = base$stan_data$sero_positive, sero_n = base$stan_data$sero_n
  )
  list(stan_data = stan_data, week_dates = base$week_dates, mapping = mapping,
       unit_labels = v3_2_units, unit_summary = unit_summary, serology = serology,
       assignment_path = assignment_path, q_prior = q_prior,
       diagnostics = list(municipalities = nrow(mapping), units = 8L,
         weeks = base$stan_data$N, accounting_error = accounting_error,
         population_gap = base$diagnostics$maximum_population_relative_gap,
         births_source = base$diagnostics$births_source,
         deaths_source = base$diagnostics$deaths_source))
}
if (sys.nframe() == 0) {
  x <- prepare_renewal_v3_2_coarse_regional_data()
  out <- here::here("02_Script/stan/renewal_ceara_v3_2_coarse_regional_data.rds")
  saveRDS(x, out)
  message(sprintf("[save] %s | %d units, %d weeks, S+U error %.3g", out,
                  x$diagnostics$units, x$diagnostics$weeks, x$diagnostics$accounting_error))
}
