# =============================================================================
# prepare_renewal_v3_1_minimal_spatial_data.R
#
# Builds the six-unit data contract for v3.1. Burden thresholds are fixed
# before fitting and stored with the prepared data; changing them for a later
# sensitivity analysis changes only this R-layer grouping, never the Stan
# model. The underlying municipality-week demographic construction is reused
# unchanged from the validated v3.0 preparation layer.
# =============================================================================

required_packages <- c("here", "dplyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
})

source(here::here(
  "02_Script/40_renewal_model/05_v3_0_spatial/prepare_renewal_v3_0_spatial_data.R"
))

unit_labels_v3_1 <- c(
  "Fortaleza",
  "Juazeiro do Norte",
  "Quixada",
  "High-burden pooled municipalities",
  "Moderate-burden pooled municipalities",
  "Low-burden pooled municipalities"
)

make_v3_1_assignment <- function(municipality_ids, cases, high_burden_min,
                                 moderate_burden_min) {
  if (!is.finite(high_burden_min) || !is.finite(moderate_burden_min) ||
    high_burden_min <= moderate_burden_min || moderate_burden_min < 0) {
    stop("Require fixed thresholds high_burden_min > moderate_burden_min >= 0")
  }

  total_cases <- rowSums(cases)
  assignment <- data.frame(
    muni6 = as.character(municipality_ids),
    total_reported_cases_2015_2019 = as.numeric(total_cases),
    unit_index = NA_integer_,
    unit = NA_character_,
    stringsAsFactors = FALSE
  )

  # Official six-digit IBGE municipality codes. The three named units are
  # removed before burden classification, irrespective of their case totals.
  named_codes <- c("230440", "230730", "231130")
  named_indices <- match(named_codes, assignment$muni6)
  if (anyNA(named_indices)) {
    stop("Fortaleza, Juazeiro do Norte, or Quixada is absent from the Cear<U+00E1> panel")
  }
  assignment$unit_index[named_indices] <- 1:3

  pooled <- is.na(assignment$unit_index)
  assignment$unit_index[pooled &
    assignment$total_reported_cases_2015_2019 >= high_burden_min] <- 4L
  assignment$unit_index[pooled &
    assignment$total_reported_cases_2015_2019 >= moderate_burden_min &
    assignment$total_reported_cases_2015_2019 < high_burden_min] <- 5L
  assignment$unit_index[pooled & is.na(assignment$unit_index)] <- 6L

  assignment$unit <- unit_labels_v3_1[assignment$unit_index]
  if (!identical(sort(unique(assignment$unit_index)), 1:6)) {
    stop("At least one of the six fixed epidemiological units is empty")
  }
  assignment
}

aggregate_unit_matrix <- function(x, assignment) {
  out <- rowsum(as.matrix(x), group = assignment$unit_index, reorder = FALSE)
  out[as.character(1:6), , drop = FALSE]
}

write_v3_1_assignment <- function(assignment, high_burden_min,
                                  moderate_burden_min) {
  table_dir <- here::here("03_Output/tables/renewal_v3_1")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  assignment_path <- file.path(
    table_dir,
    sprintf(
      "renewal_v3_1_minimal_spatial_unit_assignment_high%s_moderate%s.csv",
      format(high_burden_min, trim = TRUE, scientific = FALSE),
      format(moderate_burden_min, trim = TRUE, scientific = FALSE)
    )
  )
  utils::write.csv(assignment, assignment_path, row.names = FALSE)
  assignment_path
}

prepare_renewal_v3_1_minimal_spatial_data <- function(
  high_burden_min = 250,
  moderate_burden_min = 25
) {
  # v3.0 supplies explicitly constructed municipality N_start/N_end, births,
  # deaths, and its complete-panel checks. No municipality data are invented
  # here; this function only sums them within pre-specified strata.
  base <- prepare_renewal_v3_0_spatial_data(
    debug_subset = FALSE,
    min_total_cases = NULL
  )

  assignment <- make_v3_1_assignment(
    base$municipality_ids,
    base$stan_data$C,
    high_burden_min = high_burden_min,
    moderate_burden_min = moderate_burden_min
  )
  assignment_path <- write_v3_1_assignment(
    assignment, high_burden_min, moderate_burden_min
  )

  C <- aggregate_unit_matrix(base$stan_data$C, assignment)
  N_start <- aggregate_unit_matrix(base$stan_data$N_start, assignment)
  N_end <- aggregate_unit_matrix(base$stan_data$N_end, assignment)
  births <- aggregate_unit_matrix(base$stan_data$births, assignment)
  deaths <- aggregate_unit_matrix(base$stan_data$deaths, assignment)

  if (any(C < 0) || any(abs(C - round(C)) > 1e-8)) {
    stop("Aggregated cases must be non-negative integers")
  }
  if (any(births < 0) || any(deaths < 0)) {
    stop("Aggregated births and deaths must be non-negative")
  }
  if (max(abs(N_end[, -ncol(N_end), drop = FALSE] -
    N_start[, -1, drop = FALSE])) > 1e-6) {
    stop("Aggregated demographic stocks are discontinuous between weeks")
  }
  accounting_error <- check_demographic_accounting(N_start, N_end, births, deaths)

  if (!identical(colSums(C), colSums(base$stan_data$C))) {
    stop("Six-unit aggregation does not reproduce statewide weekly cases")
  }
  if (max(abs(colSums(N_start) - colSums(base$stan_data$N_start))) > 1e-6) {
    stop("Six-unit aggregation does not reproduce statewide population")
  }

  serology <- base$serology
  sero_unit <- match(serology$muni6, c("230440", "230730", "231130"))
  # Juazeiro and Quixada only. Fortaleza is deliberately not a likelihood row.
  if (!identical(as.integer(sero_unit), c(2L, 3L))) {
    stop("Primary serology codes must map directly to the Juazeiro and Quixada units")
  }

  population_contribution <- rowMeans(N_start) / sum(rowMeans(N_start))
  case_contribution <- rowSums(C) / sum(C)
  unit_summary <- data.frame(
    unit_index = 1:6,
    unit = unit_labels_v3_1,
    municipalities = as.integer(table(factor(assignment$unit_index, levels = 1:6))),
    mean_population = rowMeans(N_start),
    population_contribution = population_contribution,
    reported_cases = rowSums(C),
    reported_case_contribution = case_contribution,
    stringsAsFactors = FALSE
  )

  stan_data <- list(
    K = 6L,
    N = base$stan_data$N,
    G = base$stan_data$G,
    seed_weeks = base$stan_data$seed_weeks,
    C = unname(array(as.integer(C), dim = dim(C))),
    w = base$stan_data$w,
    N_start = unname(N_start),
    N_end = unname(N_end),
    births = unname(births),
    deaths = unname(deaths),
    log_seed_prior_mean = base$stan_data$log_seed_prior_mean,
    seed_prior_sd = base$stan_data$seed_prior_sd,
    J = nrow(serology),
    sero_unit = as.integer(sero_unit),
    sero_window_start = base$stan_data$sero_window_start,
    sero_window_end = base$stan_data$sero_window_end,
    sero_positive = base$stan_data$sero_positive,
    sero_n = base$stan_data$sero_n
  )

  list(
    stan_data = stan_data,
    unit_labels = unit_labels_v3_1,
    unit_assignment = assignment,
    unit_summary = unit_summary,
    assignment_path = assignment_path,
    week_dates = base$week_dates,
    serology = serology,
    fortaleza_ppc = base$fortaleza_ppc,
    diagnostics = list(
      municipalities_all_ceara = base$diagnostics$municipalities_all_ceara,
      units = 6L,
      weeks = base$diagnostics$weeks,
      maximum_accounting_error = accounting_error,
      maximum_population_relative_gap =
        base$diagnostics$maximum_population_relative_gap,
      births_source = base$diagnostics$births_source,
      deaths_source = base$diagnostics$deaths_source,
      deaths_are_approximation = base$diagnostics$deaths_are_approximation
    ),
    config = list(
      date_start = base$config$date_start,
      date_end = base$config$date_end,
      high_burden_min = high_burden_min,
      moderate_burden_min = moderate_burden_min,
      threshold_rule = "High: total cases >= high_burden_min; Moderate: moderate_burden_min <= total cases < high_burden_min; Low: total cases < moderate_burden_min. Fortaleza, Juazeiro do Norte, and Quixada are excluded before pooling.",
      fortaleza_use = "External comparison only; excluded from serology likelihood"
    )
  )
}

if (sys.nframe() == 0) {
  high_burden_min <- as.numeric(Sys.getenv("RENEWAL_V3_1_HIGH_BURDEN_MIN", "250"))
  moderate_burden_min <- as.numeric(Sys.getenv("RENEWAL_V3_1_MODERATE_BURDEN_MIN", "25"))
  prepared <- prepare_renewal_v3_1_minimal_spatial_data(
    high_burden_min = high_burden_min,
    moderate_burden_min = moderate_burden_min
  )
  out_path <- here::here(
    "02_Script/stan/renewal_ceara_v3_1_minimal_spatial_data.rds"
  )
  saveRDS(prepared, out_path)
  message(sprintf(
    "[save] %s\n[data] %d units, %d weeks; max |S+U-N|=%.3g; assignment=%s",
    out_path, prepared$diagnostics$units, prepared$diagnostics$weeks,
    prepared$diagnostics$maximum_accounting_error, prepared$assignment_path
  ))
}
