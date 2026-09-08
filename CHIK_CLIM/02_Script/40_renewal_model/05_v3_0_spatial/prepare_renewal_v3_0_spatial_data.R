# =============================================================================
# prepare_renewal_v3_0_spatial_data.R
#
# Builds the municipality x week Stan contract for renewal v3.0. This script
# deliberately does not fit a model. It keeps all demographic construction in
# R, makes every approximation explicit, and fails before sampling on an
# incomplete or incoherent panel.
# =============================================================================

required_packages <- c("here", "dplyr", "data.table", "lubridate")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(data.table)
  library(lubridate)
})

discretize_gamma_generation_interval <- function(mean_weeks, sd_weeks, G) {
  shape <- (mean_weeks / sd_weeks)^2
  rate <- mean_weeks / sd_weeks^2
  interval_mass <- diff(stats::pgamma(0:G, shape = shape, rate = rate))
  interval_mass / sum(interval_mass)
}

matrix_from_muni_week <- function(data, muni_ids, week_dates, value) {
  index <- match(
    paste(rep(muni_ids, each = length(week_dates)),
          rep(as.character(week_dates), times = length(muni_ids))),
    paste(data$muni6, as.character(data$week_start))
  )
  values <- data[[value]][index]
  if (anyNA(values)) {
    stop("Missing ", value, " after municipality-week matrix construction")
  }
  matrix(values, nrow = length(muni_ids), ncol = length(week_dates), byrow = TRUE)
}

read_muni_births_from_sinasc_cache <- function(muni_ids, date_start, date_end) {
  # Births are aggregated by maternal residence municipality from the same
  # cached SINASC microdata used by v2.2. No values are fabricated in Stan.
  cache_dir <- here::here("01_Data/sinasc_cache")
  needed_years <- seq(lubridate::year(date_start), lubridate::year(date_end))
  pieces <- vector("list", length(needed_years))

  for (i in seq_along(needed_years)) {
    year_i <- needed_years[i]
    zip_path <- file.path(cache_dir, sprintf("SINASC_%d_csv.zip", year_i))
    csv_name <- sprintf("SINASC_%d.csv", year_i)
    if (!file.exists(zip_path)) {
      stop(
        "Required cached SINASC file is absent: ", zip_path,
        ". Supply an explicit municipality-week birth table rather than ",
        "substituting zero births."
      )
    }

    extract_cmd <- sprintf('unzip -p "%s" "%s"', zip_path, csv_name)
    header_names <- names(data.table::fread(cmd = extract_cmd, nrows = 0))
    find_col <- function(name) {
      found <- header_names[toupper(header_names) == name]
      if (!length(found)) stop("Column ", name, " not found in ", csv_name)
      found[1]
    }
    col_date <- find_col("DTNASC")
    col_muni <- find_col("CODMUNRES")

    births_i <- data.table::fread(
      cmd = extract_cmd,
      select = c(col_date, col_muni),
      colClasses = "character",
      sep = ";", quote = "\"", encoding = "Latin-1"
    )
    data.table::setnames(births_i, c(col_date, col_muni), c("date_raw", "muni_raw"))
    births_i[, date := as.Date(date_raw, format = "%d%m%Y")]
    births_i[, muni6 := substr(muni_raw, 1, 6)]
    births_i <- births_i[
      !is.na(date) & date >= date_start & date <= date_end & muni6 %in% muni_ids
    ]
    births_i[, week_start := lubridate::floor_date(date, unit = "week", week_start = 7)]
    pieces[[i]] <- births_i[, .(births = .N), by = .(muni6, week_start)]
  }

  births <- data.table::rbindlist(pieces)[
    , .(births = sum(births)), by = .(muni6, week_start)
  ]
  births[, muni6 := as.character(muni6)]
  as.data.frame(births)
}

interpolate_muni_population <- function(population, muni_ids, boundary_dates) {
  # The available municipality series is annual. Annual values are treated
  # as 1-January anchors and linearly interpolated only between adjacent
  # observed annual stocks; no population series is created in Stan.
  population <- population |>
    dplyr::transmute(muni6 = as.character(muni6), year = as.integer(year),
                     population = as.numeric(population)) |>
    dplyr::filter(muni6 %in% muni_ids)

  if (anyDuplicated(population[c("muni6", "year")])) {
    stop("Municipality annual population table has duplicate municipality-year rows")
  }

  required_years <- sort(unique(lubridate::year(boundary_dates)))
  coverage_ok <- vapply(muni_ids, function(m) {
    all(required_years %in% population$year[population$muni6 == m])
  }, logical(1))
  if (any(!coverage_ok)) {
    stop("Annual population coverage is incomplete for ", sum(!coverage_ok),
         " municipalities; cannot construct N_start/N_end")
  }

  anchor_dates <- as.Date(sprintf("%d-01-01", population$year))
  out <- vapply(muni_ids, function(m) {
    x <- population |>
      dplyr::filter(muni6 == m) |>
      dplyr::arrange(year)
    stats::approx(
      x = as.numeric(as.Date(sprintf("%d-01-01", x$year))),
      y = x$population,
      xout = as.numeric(boundary_dates),
      method = "linear",
      rule = 1
    )$y
  }, numeric(length(boundary_dates)))
  if (anyNA(out) || any(out <= 0)) {
    stop("Interpolated municipality population has missing or non-positive values")
  }
  t(out)
}

weekly_state_deaths <- function(uf_projection, week_dates) {
  annual <- uf_projection |>
    dplyr::filter(uf_code == 23L) |>
    dplyr::transmute(year = as.integer(year), deaths_total = as.numeric(deaths_total))
  daily_rate <- setNames(
    annual$deaths_total / ifelse(lubridate::leap_year(annual$year), 366, 365),
    annual$year
  )
  vapply(week_dates, function(ws) {
    dates <- seq(as.Date(ws), as.Date(ws) + 6, by = "day")
    sum(daily_rate[as.character(lubridate::year(dates))])
  }, numeric(1))
}

check_demographic_accounting <- function(N_start, N_end, births, deaths) {
  # A deterministic non-zero-infection recursion verifies the algebra used
  # in Stan. Infection transfers S to U only, so S + U must still equal the
  # externally supplied population stock at every week.
  M <- nrow(N_start)
  N <- ncol(N_start)
  S <- N_start[, 1]
  U <- rep(0, M)
  max_error <- 0

  for (t in seq_len(N)) {
    max_error <- max(max_error, max(abs(S + U - N_start[, t])))
    if (t < N) {
      X <- 0.001 * S
      s_after <- S - X
      u_after <- U + X
      frac_s <- s_after / (s_after + u_after)
      reconciliation <- N_end[, t] - N_start[, t] - births[, t] + deaths[, t]
      S <- s_after + births[, t] - deaths[, t] * frac_s + reconciliation * frac_s
      U <- u_after - deaths[, t] * (1 - frac_s) + reconciliation * (1 - frac_s)
    }
  }
  if (max_error > 1e-5) {
    stop("Municipality S + U accounting failed; maximum absolute error = ", max_error)
  }
  max_error
}

prepare_renewal_v3_0_spatial_data <- function(debug_subset = FALSE, min_total_cases = NULL) {
  UF_CODE <- "23"
  DATE_START <- as.Date("2015-01-04")
  DATE_END <- as.Date("2019-12-29")
  G <- 8L
  SEED_WEEKS <- 8L

  serology <- data.frame(
    site = c("Juazeiro do Norte", "Quixada"),
    muni6 = c("230730", "231130"),
    window_start = as.Date(c("2018-06-03", "2018-06-03")),
    window_end = as.Date(c("2018-12-30", "2019-12-29")),
    positive = c(103L, 289L),
    n = c(404L, 409L),
    stringsAsFactors = FALSE
  )
  fortaleza_ppc <- data.frame(
    site = "Fortaleza cohort (external PPC only)",
    muni6 = "230440",
    stringsAsFactors = FALSE
  )

  panel <- readRDS(here::here("01_Data/chik_dlnm_panel_muni_week_2015_2025.rds")) |>
    dplyr::transmute(
      muni6 = as.character(muni6),
      week_start = as.Date(week_start),
      cases = as.numeric(cases_confirmed)
    ) |>
    dplyr::filter(
      substr(muni6, 1, 2) == UF_CODE,
      week_start >= DATE_START,
      week_start <= DATE_END
    ) |>
    dplyr::arrange(muni6, week_start)

  week_dates <- seq(DATE_START, DATE_END, by = "week")
  muni_ids_all <- sort(unique(panel$muni6))
  expected_rows <- length(muni_ids_all) * length(week_dates)
  if (nrow(panel) != expected_rows) {
    stop("Ceara panel is not a complete municipality x week rectangle: got ",
         nrow(panel), ", expected ", expected_rows)
  }
  if (anyDuplicated(panel[c("muni6", "week_start")])) {
    stop("Duplicated municipality-week rows in the Ceara case panel")
  }
  if (!identical(sort(unique(panel$week_start)), week_dates)) {
    stop("Week dates are incomplete or not in the intended Sunday-start sequence")
  }
  if (anyNA(panel$cases) || any(panel$cases < 0) ||
      any(abs(panel$cases - round(panel$cases)) > 1e-8)) {
    stop("Reported municipality-week cases must be complete non-negative integers")
  }

  # This construction explicitly repeats each municipality's ordered weeks;
  # the subsequent check prevents a renewal lag from crossing municipality IDs.
  cases_all <- matrix_from_muni_week(panel, muni_ids_all, week_dates, "cases")
  if (any(vapply(seq_len(nrow(cases_all)), function(i) {
    !identical(
      panel$week_start[panel$muni6 == muni_ids_all[i]],
      week_dates
    )
  }, logical(1)))) {
    stop("Municipality-specific week ordering failed; lag histories could cross IDs")
  }

  population <- readRDS(here::here("01_Data/ibge_pop_muni_year_2015_2025.rds"))
  boundary_dates <- c(week_dates, max(week_dates) + 7)
  population_boundary <- interpolate_muni_population(
    population, muni_ids_all, boundary_dates
  )
  N_start_all <- population_boundary[, seq_along(week_dates), drop = FALSE]
  N_end_all <- population_boundary[, seq_along(week_dates) + 1L, drop = FALSE]

  births_long <- read_muni_births_from_sinasc_cache(
    muni_ids_all, DATE_START, DATE_END + 6
  )
  birth_index <- match(
    paste(rep(muni_ids_all, each = length(week_dates)),
          rep(as.character(week_dates), times = length(muni_ids_all))),
    paste(births_long$muni6, as.character(births_long$week_start))
  )
  birth_values <- births_long$births[birth_index]
  birth_values[is.na(birth_values)] <- 0
  births_all <- matrix(
    birth_values, nrow = length(muni_ids_all), ncol = length(week_dates),
    byrow = TRUE
  )

  uf_projection <- readRDS(
    here::here("01_Data/ibge_population_projection_uf_2024revision.rds")
  )
  state_deaths <- weekly_state_deaths(uf_projection, week_dates)

  # Municipal all-cause deaths are not available in the current project
  # inputs. The R layer therefore allocates the documented IBGE state annual
  # all-cause total to municipalities by start-of-week population share. This
  # is explicit, preserves the state weekly total exactly, and is a stated
  # limitation rather than an unrecorded Stan-side assumption.
  deaths_all <- sweep(N_start_all, 2, colSums(N_start_all), "/") *
    matrix(rep(state_deaths, each = nrow(N_start_all)),
           nrow = nrow(N_start_all), ncol = length(state_deaths))

  if (any(births_all < 0) || any(deaths_all < 0) ||
      any(!is.finite(births_all)) || any(!is.finite(deaths_all))) {
    stop("Births and deaths must be finite and non-negative")
  }
  if (max(abs(N_end_all[, -ncol(N_end_all), drop = FALSE] -
              N_start_all[, -1, drop = FALSE])) > 1e-6) {
    stop("N_end[m,t] does not equal N_start[m,t+1]; demographic stocks are discontinuous")
  }
  accounting_error <- check_demographic_accounting(
    N_start_all, N_end_all, births_all, deaths_all
  )

  state_demography <- readRDS(here::here("01_Data/ceara_weekly_demography.rds")) |>
    dplyr::filter(week_start %in% week_dates) |>
    dplyr::arrange(week_start)
  if (!identical(state_demography$week_start, week_dates)) {
    stop("State demographic contract lacks the requested 2015-2019 weeks")
  }
  population_relative_gap <- (colSums(N_start_all) - state_demography$N_start) /
    state_demography$N_start
  if (max(abs(population_relative_gap)) > 0.10) {
    stop("Municipality population sum differs from the state stock by more than 10%")
  }

  if (!all(serology$muni6 %in% muni_ids_all) ||
      !fortaleza_ppc$muni6 %in% muni_ids_all) {
    stop("One or more hard-coded IBGE municipality codes are absent from the Ceara panel")
  }
  sero_start_all <- match(serology$window_start, week_dates)
  sero_end_all <- match(serology$window_end, week_dates)
  if (anyNA(sero_start_all) || anyNA(sero_end_all) || any(sero_start_all > sero_end_all)) {
    stop("Serology collection windows do not map to complete intended weekly intervals")
  }

  if (debug_subset && !is.null(min_total_cases)) {
    stop("debug_subset and min_total_cases are mutually exclusive selection modes")
  }
  total_cases <- rowSums(cases_all)
  required <- c(serology$muni6, fortaleza_ppc$muni6)

  if (debug_subset) {
    remaining <- setdiff(muni_ids_all, required)
    high <- remaining[order(total_cases[match(remaining, muni_ids_all)],
                            decreasing = TRUE)][seq_len(min(5L, length(remaining)))]
    remaining <- setdiff(remaining, high)
    low <- remaining[order(total_cases[match(remaining, muni_ids_all)],
                           decreasing = FALSE)][seq_len(min(5L, length(remaining)))]
    muni_ids <- unique(c(required, high, low))
  } else if (!is.null(min_total_cases)) {
    # Case counts are extremely concentrated in Ceara: 56 of 184 municipalities
    # already account for 94% of all reported cases over 2015-2019, and most
    # of the remainder have single- or double-digit cumulative case counts,
    # far too sparse to identify a municipality-specific trajectory. This
    # keeps every municipality clearing the threshold, plus the mandatory
    # serology/PPC municipalities regardless of their own case count (e.g.
    # Juazeiro do Norte has only 196 cases, below a 250-case threshold, but
    # is required by the serology likelihood).
    above_threshold <- muni_ids_all[total_cases > min_total_cases]
    muni_ids <- unique(c(required, above_threshold))
  } else {
    muni_ids <- muni_ids_all
  }
  muni_idx <- match(muni_ids, muni_ids_all)

  cases <- cases_all[muni_idx, , drop = FALSE]
  N_start <- N_start_all[muni_idx, , drop = FALSE]
  N_end <- N_end_all[muni_idx, , drop = FALSE]
  births <- births_all[muni_idx, , drop = FALSE]
  deaths <- deaths_all[muni_idx, , drop = FALSE]
  sero_muni <- match(serology$muni6, muni_ids)
  if (anyNA(sero_muni)) stop("Municipality subset omitted a primary serology municipality")

  stan_data <- list(
    M = length(muni_ids),
    N = length(week_dates),
    G = G,
    seed_weeks = SEED_WEEKS,
    C = unname(array(as.integer(cases), dim = dim(cases))),
    w = as.vector(discretize_gamma_generation_interval(2, 1, G)),
    N_start = unname(N_start),
    N_end = unname(N_end),
    births = unname(births),
    deaths = unname(deaths),
    log_seed_prior_mean = rep(log(10), SEED_WEEKS),
    seed_prior_sd = 1.5,
    J = nrow(serology),
    sero_muni = as.integer(sero_muni),
    sero_window_start = as.integer(sero_start_all),
    sero_window_end = as.integer(sero_end_all),
    sero_positive = as.integer(serology$positive),
    sero_n = as.integer(serology$n)
  )

  list(
    stan_data = stan_data,
    municipality_ids = muni_ids,
    week_dates = week_dates,
    serology = serology,
    fortaleza_ppc = fortaleza_ppc,
    diagnostics = list(
      panel_rows = nrow(panel),
      municipalities_all_ceara = length(muni_ids_all),
      municipalities_in_stan_data = length(muni_ids),
      weeks = length(week_dates),
      duplicate_muni_week_rows = 0L,
      maximum_accounting_error = accounting_error,
      maximum_population_relative_gap = max(abs(population_relative_gap)),
      population_relative_gap_range = range(population_relative_gap),
      births_source = "SINASC individual records by maternal residence, weekly aggregation",
      deaths_source = "IBGE Ceara annual all-cause deaths, day-weighted then allocated by municipality start-week population share",
      deaths_are_approximation = TRUE
    ),
    config = list(
      date_start = DATE_START,
      date_end = DATE_END,
      debug_subset = debug_subset,
      min_total_cases = min_total_cases,
      generation_interval_mean_weeks = 2,
      generation_interval_sd_weeks = 1,
      primary_serology_sites = serology$site,
      fortaleza_use = "External posterior predictive comparison only; excluded from likelihood"
    )
  )
}

if (sys.nframe() == 0) {
  debug_subset <- tolower(Sys.getenv("RENEWAL_V3_0_DEBUG_SUBSET", "false")) %in%
    c("1", "true", "yes")
  min_total_cases_env <- Sys.getenv("RENEWAL_V3_0_MIN_CASES", "")
  min_total_cases <- if (nzchar(min_total_cases_env)) as.numeric(min_total_cases_env) else NULL
  prepared <- prepare_renewal_v3_0_spatial_data(
    debug_subset = debug_subset, min_total_cases = min_total_cases
  )
  suffix <- if (debug_subset) {
    "_debug"
  } else if (!is.null(min_total_cases)) {
    sprintf("_min%d_cases", min_total_cases)
  } else {
    "_all_ceara"
  }
  out_path <- here::here(
    "02_Script/stan", paste0("renewal_ceara_v3_0_spatial_data", suffix, ".rds")
  )
  saveRDS(prepared, out_path)
  message(sprintf(
    "[save] %s\n[data] %d municipalities in Stan data (%d all-Ceara), %d weeks; max |S+U-N|=%.3g; max state-population gap=%.3f%%",
    out_path,
    prepared$diagnostics$municipalities_in_stan_data,
    prepared$diagnostics$municipalities_all_ceara,
    prepared$diagnostics$weeks,
    prepared$diagnostics$maximum_accounting_error,
    100 * prepared$diagnostics$maximum_population_relative_gap
  ))
}
