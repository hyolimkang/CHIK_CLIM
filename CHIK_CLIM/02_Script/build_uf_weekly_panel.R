# build_uf_weekly_panel.R
# Aggregate muni-week SINAN panel to state x week, join population and births.
# Set uf_code / uf_name in the config block when running interactively.

resolve_data_path <- function(candidates) {
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop("None of these paths exist:\n", paste("  ", candidates, collapse = "\n"))
  }
  found[[1]]
}

smooth_population_for_transmission_model <- function(pop_year) {
  pop_year <- pop_year |>
    dplyr::arrange(year) |>
    dplyr::mutate(population_raw = population)

  # IBGE 2022 is a Census benchmark, while surrounding years are annual
  # estimates. For weekly SIR dynamics, avoid treating this benchmark revision
  # as a sudden demographic shock by interpolating 2022/2023 between the
  # consistent estimate years 2021 and 2024.
  if (all(c(2021L, 2024L) %in% pop_year$year)) {
    pop_2021 <- pop_year$population[pop_year$year == 2021L][1]
    pop_2024 <- pop_year$population[pop_year$year == 2024L][1]

    for (yr in c(2022L, 2023L)) {
      if (yr %in% pop_year$year) {
        pop_year$population[pop_year$year == yr] <- round(
          pop_2021 + (pop_2024 - pop_2021) * ((yr - 2021) / (2024 - 2021))
        )
      }
    }
  }

  pop_year
}

load_weekly_panel <- function(panel_path = NULL) {
  if (is.null(panel_path)) {
    panel_path <- resolve_data_path(c(
      here::here("01_Data/chik_brazil_muni_week_2015_2024.rds"),
      "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_MORBID/CHIK_MORBID/01_Data/chik_brazil_muni_week_2015_2024.rds",
      "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_MORBID/CHIK_MORBID/01_Data/chik_brazil_muni_week_2015_2024.csv"
    ))
  }

  if (grepl("\\.rds$", panel_path, ignore.case = TRUE)) {
    panel <- readRDS(panel_path)
  } else {
    panel <- readr::read_csv(panel_path, show_col_types = FALSE)
  }

  panel |>
    dplyr::mutate(muni6 = sprintf("%06d", as.integer(muni6)))
}

build_uf_weekly_panel <- function(panel_week,
                                    uf_code,
                                    uf_name,
                                    year_start,
                                    year_end,
                                    case_var = "cases_confirmed",
                                    deaths_var = "deaths_from_chik",
                                    pop_path = NULL,
                                    births_path = NULL) {
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate")
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
  if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
  if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")

  uf_weekly <- panel_week |>
    dplyr::mutate(week_start = as.Date(week_start)) |>
    dplyr::filter(substr(muni6, 1, 2) == uf_code) |>
    dplyr::group_by(week_start) |>
    dplyr::summarise(
      cases  = sum(.data[[case_var]], na.rm = TRUE),
      deaths = sum(.data[[deaths_var]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::filter(
      week_start >= as.Date(sprintf("%d-01-01", year_start)),
      week_start <= as.Date(sprintf("%d-12-31", year_end))
    ) |>
    dplyr::arrange(week_start) |>
    dplyr::mutate(
      year         = lubridate::year(week_start),
      week_of_year = lubridate::isoweek(week_start),
      t            = dplyr::row_number()
    )

  if (is.null(pop_path)) {
    pop_path <- resolve_data_path(c(
      here::here(sprintf("01_Data/ibge_pop_muni_year_%d_%d.csv", year_start, year_end)),
      here::here("01_Data/ibge_pop_muni_year_2015_2024.csv"),
      "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_MORBID/CHIK_MORBID/01_Data/ibge_pop_muni_year_2015_2024.csv"
    ))
  }

  pop <- readr::read_csv(pop_path, show_col_types = FALSE)
  pop_year <- pop |>
    dplyr::mutate(muni6 = sprintf("%06d", as.integer(muni6))) |>
    dplyr::filter(substr(muni6, 1, 2) == uf_code) |>
    dplyr::group_by(year) |>
    dplyr::summarise(population = sum(population), .groups = "drop") |>
    smooth_population_for_transmission_model()

  uf_weekly <- uf_weekly |>
    dplyr::left_join(pop_year, by = "year")

  if (is.null(births_path)) {
    births_path <- here::here(sprintf("01_Data/%s_births_annual_ibge.csv", uf_name))
  }
  if (!file.exists(births_path)) {
    stop("Births file not found: ", births_path,
         "\nRun fetch_ce_births_ibge.R for this UF first.")
  }

  births_annual <- readr::read_csv(births_path, show_col_types = FALSE)
  if (year_end > max(births_annual$year)) {
    births_annual <- dplyr::bind_rows(
      births_annual,
      tibble::tibble(
        year   = year_end,
        births = births_annual$births[births_annual$year == max(births_annual$year)],
        source = sprintf(
          "IBGE SIDRA t2609 UF%s (%d carried forward)",
          uf_code, max(births_annual$year)
        )
      )
    )
  }

  uf_weekly <- uf_weekly |>
    dplyr::left_join(births_annual |> dplyr::select(year, births), by = "year") |>
    dplyr::mutate(births_weekly = births / 52)

  if (any(is.na(uf_weekly$births_weekly))) {
    stop("Missing births for some years — check births join for UF", uf_code, ".")
  }

  uf_weekly
}

# Ensure year / week_of_year exist (handles old RDS that used `week` only).
standardize_uf_weekly <- function(df) {
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate")

  if (!"week_start" %in% names(df)) {
    stop("`week_start` column missing from uf weekly panel.")
  }

  df <- df |>
    dplyr::mutate(week_start = as.Date(week_start))

  if (!"year" %in% names(df)) {
    df <- df |> dplyr::mutate(year = lubridate::year(week_start))
  }

  if (!"week_of_year" %in% names(df)) {
    if ("week" %in% names(df)) {
      df <- df |> dplyr::mutate(week_of_year = as.integer(week))
    } else {
      df <- df |> dplyr::mutate(week_of_year = lubridate::isoweek(week_start))
    }
  }

  df
}

# ISO week index for Stan (length must match nrow(ce_fit)).
make_week_id <- function(df, W = 52L) {
  df <- standardize_uf_weekly(df)
  as.integer(pmin(df$week_of_year, W))
}

check_stan_panel <- function(stan_data, label = "stan_data") {
  N <- stan_data$N
  req_len <- c(
    week_id = N,
    cases = N,
    year_id = N,
    pop = N,
    births_weekly = N
  )
  bad <- vapply(names(req_len), function(nm) {
    if (!nm %in% names(stan_data)) return(TRUE)
    length(stan_data[[nm]]) != req_len[[nm]]
  }, logical(1))
  if (any(bad)) {
    msgs <- vapply(names(bad)[bad], function(nm) {
      have <- if (nm %in% names(stan_data)) length(stan_data[[nm]]) else NA_integer_
      sprintf("  %s: expected length %d, got %s", nm, N, have)
    }, character(1))
    stop(
      sprintf("%s: Stan data length mismatch (N = %d):\n%s",
              label, N, paste(msgs, collapse = "\n")),
      call. = FALSE
    )
  }
  invisible(stan_data)
}

if (sys.nframe() == 0) {
  # ============================================================
  # USER CONFIG
  # ============================================================
  UF_CODE    <- "23"    # CE=23, MG=31, SP=35, RJ=33, BA=29, ...
  UF_NAME    <- "ce"    # lowercase prefix for 01_Data file names
  YEAR_START <- 2014
  YEAR_END   <- 2025
  CASE_VAR   <- "cases_confirmed"
  # ============================================================

  panel_week <- load_weekly_panel()
  uf_weekly  <- build_uf_weekly_panel(
    panel_week = panel_week,
    uf_code    = UF_CODE,
    uf_name    = UF_NAME,
    year_start = YEAR_START,
    year_end   = YEAR_END,
    case_var   = CASE_VAR
  )

  out_path <- here::here(sprintf("01_Data/%s_weekly_%d_%d.rds", UF_NAME, YEAR_START, YEAR_END))
  saveRDS(uf_weekly, out_path)
  message("Saved: ", out_path, " (", nrow(uf_weekly), " weeks)")
}
