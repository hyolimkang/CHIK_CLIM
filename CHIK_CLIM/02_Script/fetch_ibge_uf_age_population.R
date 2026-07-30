# ---------------------------------------------------------------------------
# fetch_ibge_uf_age_population.R
#
# Fetch IBGE population projections by UF, year, and single age from SIDRA
# table 7358, then aggregate to the 20 age groups used in downstream burden
# calculations.
#
# Data source:
#   IBGE SIDRA table 7358, variable 606
#   "Populacao, por sexo e idade" (Projecao da Populacao, revision 2018)
#
# Outputs:
#   01_Data/ibge_pop_uf_age_group_2015_2024.rds / .csv
#     uf_code, uf_name, year, age_group, age_mid, age_lower, age_upper,
#     population
#
#   01_Data/ibge_pop_uf_age_group_wide_2015_2024.rds / .csv
#     one row per UF-year, one column per age group population
#
#   01_Data/ibge_pop_uf_age_group_split12_2015_2024.rds / .csv
#   01_Data/ibge_pop_uf_age_group_split12_wide_2015_2024.rds / .csv
#     same outputs, but with age 12 as its own group and 13-17 separated
#
#   01_Data/ibge_pop_uf_single_age_expanded_2015_2024.rds / .csv
#     uf_code, uf_name, year, age, population
#
# The grouped output uses:
#   <1, 1-4, 5-9, 10-11, 12-17, 18-19, ..., 80-84, 85+
# ---------------------------------------------------------------------------

for (p in c("dplyr", "readr", "tidyr", "httr2", "jsonlite", "here", "stringr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(httr2)
  library(jsonlite)
  library(here)
  library(stringr)
})

# ---- Config ----------------------------------------------------------------
YEAR_START <- 2015L
YEAR_END <- 2024L
YEARS <- YEAR_START:YEAR_END

# Set UF_CODES before sourcing this script if you only need selected states,
# e.g. UF_CODES <- c("23", "31").
if (!exists("UF_CODES")) {
  UF_CODES <- c(
    "11", "12", "13", "14", "15", "16", "17",
    "21", "22", "23", "24", "25", "26", "27", "28", "29",
    "31", "32", "33", "35",
    "41", "42", "43",
    "50", "51", "52", "53"
  )
}

age_gr_levels <- c(
  "<1", "1-4", "5-9", "10-11", "12-17", "18-19", "20-24", "25-29",
  "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
  "60-64", "65-69", "70-74", "75-79", "80-84", "85+"
)

age_group_map <- tibble::tibble(
  age_group = factor(age_gr_levels, levels = age_gr_levels),
  age_mid = c(
    mean(0:1), mean(1:4), mean(5:9), mean(10:11), mean(12:17),
    mean(18:19), mean(20:24), mean(25:29), mean(30:34), mean(35:39),
    mean(40:44), mean(45:49), mean(50:54), mean(55:59), mean(60:64),
    mean(65:69), mean(70:74), mean(75:79), mean(80:84), mean(85:89)
  ),
  age_lower = c(0L, 1L, 5L, 10L, 12L, 18L, 20L, 25L, 30L, 35L,
                40L, 45L, 50L, 55L, 60L, 65L, 70L, 75L, 80L, 85L),
  age_upper = c(0L, 4L, 9L, 11L, 17L, 19L, 24L, 29L, 34L, 39L,
                44L, 49L, 54L, 59L, 64L, 69L, 74L, 79L, 84L, Inf)
)

age_groups <- age_group_map$age_mid

age_gr_levels_split12 <- c(
  "<1", "1-4", "5-9", "10-11", "12", "13-17", "18-19", "20-24",
  "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
  "60-64", "65-69", "70-74", "75-79", "80-84", "85+"
)

age_group_map_split12 <- tibble::tibble(
  age_group = factor(age_gr_levels_split12, levels = age_gr_levels_split12),
  age_mid = c(
    mean(0:1), mean(1:4), mean(5:9), mean(10:11), 12,
    mean(13:17), mean(18:19), mean(20:24), mean(25:29), mean(30:34),
    mean(35:39), mean(40:44), mean(45:49), mean(50:54), mean(55:59),
    mean(60:64), mean(65:69), mean(70:74), mean(75:79), mean(80:84),
    mean(85:89)
  ),
  age_lower = c(0L, 1L, 5L, 10L, 12L, 13L, 18L, 20L, 25L, 30L,
                35L, 40L, 45L, 50L, 55L, 60L, 65L, 70L, 75L, 80L, 85L),
  age_upper = c(0L, 4L, 9L, 11L, 12L, 17L, 19L, 24L, 29L, 34L,
                39L, 44L, 49L, 54L, 59L, 64L, 69L, 74L, 79L, 84L, Inf)
)

age_groups_split12 <- age_group_map_split12$age_mid

sidra_year_codes <- c(
  `2015` = 49031, `2016` = 49032, `2017` = 49033, `2018` = 49034,
  `2019` = 49035, `2020` = 49036, `2021` = 49037, `2022` = 49038,
  `2023` = 49039, `2024` = 49040
)

sidra_age_codes <- tibble::tibble(
  age = c(0:89, 90L),
  age_code = c(
    49111, 6558:6566, 6567:6576, 6577:6581, 6582, 6656:6659,
    6583:6587, 6588:6592, 6593:6597, 6598:6602, 6603:6607,
    6608:6612, 6613:6617, 6618:6622, 6623:6627, 6628:6632,
    6633:6642, 49110
  )
)

if (!all(as.character(YEARS) %in% names(sidra_year_codes))) {
  stop("YEARS contains values not mapped in sidra_year_codes.", call. = FALSE)
}

DATA_DIR <- here::here("01_Data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_GROUP_RDS <- file.path(
  DATA_DIR,
  sprintf("ibge_pop_uf_age_group_%d_%d.rds", YEAR_START, YEAR_END)
)
OUT_GROUP_CSV <- sub("\\.rds$", ".csv", OUT_GROUP_RDS)

OUT_GROUP_WIDE_RDS <- file.path(
  DATA_DIR,
  sprintf("ibge_pop_uf_age_group_wide_%d_%d.rds", YEAR_START, YEAR_END)
)
OUT_GROUP_WIDE_CSV <- sub("\\.rds$", ".csv", OUT_GROUP_WIDE_RDS)

OUT_GROUP_SPLIT12_RDS <- file.path(
  DATA_DIR,
  sprintf("ibge_pop_uf_age_group_split12_%d_%d.rds", YEAR_START, YEAR_END)
)
OUT_GROUP_SPLIT12_CSV <- sub("\\.rds$", ".csv", OUT_GROUP_SPLIT12_RDS)

OUT_GROUP_SPLIT12_WIDE_RDS <- file.path(
  DATA_DIR,
  sprintf("ibge_pop_uf_age_group_split12_wide_%d_%d.rds", YEAR_START, YEAR_END)
)
OUT_GROUP_SPLIT12_WIDE_CSV <- sub("\\.rds$", ".csv", OUT_GROUP_SPLIT12_WIDE_RDS)

OUT_SINGLE_RDS <- file.path(
  DATA_DIR,
  sprintf("ibge_pop_uf_single_age_expanded_%d_%d.rds", YEAR_START, YEAR_END)
)
OUT_SINGLE_CSV <- sub("\\.rds$", ".csv", OUT_SINGLE_RDS)

# ---- Helpers ---------------------------------------------------------------
find_col <- function(x, pattern, last = FALSE) {
  hit <- grep(pattern, names(x), ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0) stop("No column matching pattern: ", pattern, call. = FALSE)
  if (last) hit[[length(hit)]] else hit[[1]]
}

parse_sidra_value <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", x, fixed = TRUE)))
}

parse_age_label <- function(age_label) {
  age_label <- stringr::str_squish(age_label)
  dplyr::case_when(
    stringr::str_detect(age_label, "^\\d+\\s+anos?$") ~
      as.integer(stringr::str_extract(age_label, "^\\d+")),
    stringr::str_detect(age_label, "^90\\s+anos\\s+ou\\s+mais$") ~ 90L,
    TRUE ~ NA_integer_
  )
}

age_to_group <- function(age, group_map) {
  out <- rep(NA_character_, length(age))
  for (i in seq_len(nrow(group_map))) {
    lo <- group_map$age_lower[i]
    hi <- group_map$age_upper[i]
    valid_age <- !is.na(age)
    in_group <- if (is.infinite(hi)) {
      valid_age & age >= lo
    } else {
      valid_age & age >= lo & age <= hi
    }
    out[in_group] <- as.character(group_map$age_group[i])
  }
  factor(out, levels = levels(group_map$age_group))
}

fetch_one_uf <- function(uf_code, years) {
  year_arg <- paste(unname(sidra_year_codes[as.character(years)]), collapse = ",")
  age_arg <- paste(sidra_age_codes$age_code, collapse = ",")
  api_url <- sprintf(
    "https://apisidra.ibge.gov.br/values/t/7358/n3/%s/v/606/p/2018/c2/6794/c287/%s/c1933/%s/f/u?formato=json",
    uf_code, age_arg, year_arg
  )

  message(sprintf("[sidra] UF %s, years %s-%s ...", uf_code, min(years), max(years)))

  req <- httr2::request(api_url) |>
    httr2::req_timeout(90)

  raw_txt <- httr2::req_perform(req) |>
    httr2::resp_body_string()

  raw <- jsonlite::fromJSON(raw_txt, simplifyDataFrame = TRUE)
  if (nrow(raw) <= 1) stop("SIDRA returned no data for UF ", uf_code, call. = FALSE)

  # First row is the SIDRA header. Table 7358 returns two "Ano" columns:
  # period year (= 2018) and projected population year (= 2015, 2016, ...).
  # make.unique() avoids duplicate-name failures; we then take the last Ano.
  names(raw) <- make.unique(as.character(unlist(raw[1, ], use.names = FALSE)))
  raw <- raw[-1, , drop = FALSE]

  uf_code_col <- find_col(raw, "Unidade da Federa..o.*C.digo|UF.*C.digo|NC")
  uf_name_col <- find_col(raw, "Unidade da Federa..o$|UF$|Nível Territorial")
  age_col <- find_col(raw, "^Idade$")
  year_col <- find_col(raw, "^Ano(\\.\\d+)?$", last = TRUE)
  value_col <- find_col(raw, "^Valor$")

  raw |>
    transmute(
      uf_code = sprintf("%02s", .data[[uf_code_col]]),
      uf_name = .data[[uf_name_col]],
      year = as.integer(.data[[year_col]]),
      age_label = .data[[age_col]],
      age = parse_age_label(.data[[age_col]]),
      population = round(parse_sidra_value(.data[[value_col]]))
    ) |>
    filter(!is.na(.data$age), !is.na(.data$population), .data$population > 0)
}

expand_single_age_for_simulator <- function(single_age, max_age = 100L) {
  single_age |>
    mutate(
      age_upper_expand = ifelse(.data$age >= 90L, max_age, .data$age),
      age_lower_expand = .data$age,
      n_age = .data$age_upper_expand - .data$age_lower_expand + 1L
    ) |>
    rowwise() |>
    mutate(
      age_expanded = list(seq.int(.data$age_lower_expand, .data$age_upper_expand)),
      population_each = .data$population / .data$n_age
    ) |>
    ungroup() |>
    select(.data$uf_code, .data$uf_name, .data$year, .data$age_expanded, .data$population_each) |>
    tidyr::unnest(.data$age_expanded) |>
    transmute(
      uf_code = .data$uf_code,
      uf_name = .data$uf_name,
      year = .data$year,
      age = as.integer(.data$age_expanded),
      population = .data$population_each
    )
}

aggregate_age_groups <- function(single_age, group_map) {
  single_age |>
    filter(!is.na(.data$age)) |>
    mutate(age_group = age_to_group(.data$age, group_map)) |>
    filter(!is.na(.data$age_group)) |>
    group_by(.data$uf_code, .data$uf_name, .data$year, .data$age_group) |>
    summarise(population = sum(.data$population, na.rm = TRUE), .groups = "drop") |>
    left_join(group_map, by = "age_group") |>
    select(
      .data$uf_code, .data$uf_name, .data$year, .data$age_group,
      .data$age_mid, .data$age_lower, .data$age_upper, .data$population
    ) |>
    arrange(.data$uf_code, .data$year, .data$age_lower)
}

make_age_group_wide <- function(grouped_pop) {
  grouped_pop |>
    select(.data$uf_code, .data$uf_name, .data$year, .data$age_group, .data$population) |>
    tidyr::pivot_wider(
      names_from = .data$age_group,
      values_from = .data$population
    ) |>
    arrange(.data$uf_code, .data$year)
}

# ---- Fetch and aggregate ---------------------------------------------------
raw_single_age <- dplyr::bind_rows(lapply(UF_CODES, fetch_one_uf, years = YEARS))

pop_age_group <- aggregate_age_groups(raw_single_age, age_group_map)

pop_age_group_wide <- make_age_group_wide(pop_age_group)

pop_age_group_split12 <- aggregate_age_groups(raw_single_age, age_group_map_split12)

pop_age_group_split12_wide <- make_age_group_wide(pop_age_group_split12)

pop_single_expanded <- expand_single_age_for_simulator(raw_single_age) |>
  arrange(.data$uf_code, .data$year, .data$age)

# ---- Save ------------------------------------------------------------------
saveRDS(pop_age_group, OUT_GROUP_RDS)
readr::write_csv(pop_age_group, OUT_GROUP_CSV)

saveRDS(pop_age_group_wide, OUT_GROUP_WIDE_RDS)
readr::write_csv(pop_age_group_wide, OUT_GROUP_WIDE_CSV)

saveRDS(pop_age_group_split12, OUT_GROUP_SPLIT12_RDS)
readr::write_csv(pop_age_group_split12, OUT_GROUP_SPLIT12_CSV)

saveRDS(pop_age_group_split12_wide, OUT_GROUP_SPLIT12_WIDE_RDS)
readr::write_csv(pop_age_group_split12_wide, OUT_GROUP_SPLIT12_WIDE_CSV)

saveRDS(pop_single_expanded, OUT_SINGLE_RDS)
readr::write_csv(pop_single_expanded, OUT_SINGLE_CSV)

message(sprintf(
  "[save] UF-year-age-group population\n  %s\n  %s rows | %d UFs | %d years | %d age groups",
  OUT_GROUP_RDS,
  format(nrow(pop_age_group), big.mark = ","),
  dplyr::n_distinct(pop_age_group$uf_code),
  dplyr::n_distinct(pop_age_group$year),
  dplyr::n_distinct(pop_age_group$age_group)
))

message(sprintf(
  "[save] Expanded single-age population for simulator initialisation\n  %s\n  %s rows",
  OUT_SINGLE_RDS,
  format(nrow(pop_single_expanded), big.mark = ",")
))

message(sprintf(
  "[save] Wide UF-year N matrix by age group\n  %s\n  %s rows",
  OUT_GROUP_WIDE_RDS,
  format(nrow(pop_age_group_wide), big.mark = ",")
))

message(sprintf(
  "[save] Split-12 UF-year-age-group population\n  %s\n  %s rows | %d age groups",
  OUT_GROUP_SPLIT12_RDS,
  format(nrow(pop_age_group_split12), big.mark = ","),
  dplyr::n_distinct(pop_age_group_split12$age_group)
))

message(sprintf(
  "[save] Split-12 wide UF-year N matrix by age group\n  %s\n  %s rows",
  OUT_GROUP_SPLIT12_WIDE_RDS,
  format(nrow(pop_age_group_split12_wide), big.mark = ",")
))
