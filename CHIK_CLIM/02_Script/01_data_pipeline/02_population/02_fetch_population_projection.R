# ===========================================================================
# 07_fetch_ibge_population_projection_uf.R
#
# Purpose
# -------
# UF-level annual population for the Ceara renewal-model track (state-level
# N_pop_t), sourced from a single, internally consistent IBGE revision
# instead of the muni-level annual estimates used elsewhere in this repo
# (01_data_pipeline/02_population/01_fetch_population.R), which mix several vintages:
# SIDRA table 6579 annual estimates for most years, the 2022 Demographic
# Census (table 4714) for 2022, and linear interpolation for 2023. IBGE
# itself recommends against mixing revisions in one series, and mixing them
# here produces a visible artifact: Ceara's muni-summed population drops
# from 9,240,580 (2021) to 8,794,957 (2022, Census year) before rebounding
# to 9,014,301 (2023, interpolated) -- an artificial ~450k (~5%) dip-and-
# recover that has nothing to do with real demographic change. Since the
# renewal model (40_renewal_model/) drives weekly susceptible-pool changes
# directly off population deltas, that artifact would be read by the model
# as a mass in/out-migration event.
#
# Source
# ------
# IBGE "Projecoes da Populacao: Brasil e Unidades da Federacao", revision
# 2024 (3rd edition, released 2025-11-07), covering 2000-2070 by UF. Table 4
# ("indicadores") carries one consistent annual series per UF for total
# population plus the compensating-equation components IBGE publishes
# directly: births (NASC) and deaths (OBT). Net migration is not published
# as its own column in this table; it is implicit in how POP_T evolves
# (P[t+1] = P[t] + births - deaths +/- net migration) rather than break-outable
# from these columns alone. Age-structured UF population elsewhere in this
# repo (01_data_pipeline/02_population/03_fetch_age_population.R) still uses the older
# SIDRA table 7358 / revision 2018 -- a separate, pre-existing inconsistency
# in the SIR and vaccine-impact tracks that this script does not touch.
#   https://ftp.ibge.gov.br/Projecao_da_Populacao/Projecao_da_Populacao_2024/projecoes_2024_tab4_indicadores.xlsx
#
# Output
# ------
#   01_Data/ibge_population_projection_uf_2024revision.rds (+ .csv)
#     columns: uf_code, uf_sigla, uf_name, year, pop_total, pop_male,
#              pop_female, births_total, deaths_total, natural_growth_rate
#     coverage: all 27 UF, 2000-2070
# ===========================================================================

for (p in c("here", "dplyr", "readr", "readxl")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr); library(readxl)
})

SOURCE_URL <- paste0(
  "https://ftp.ibge.gov.br/Projecao_da_Populacao/",
  "Projecao_da_Populacao_2024/projecoes_2024_tab4_indicadores.xlsx"
)

fetch_ibge_population_projection_uf <- function(cache_path) {
  if (!file.exists(cache_path)) {
    message("[download] ", SOURCE_URL)
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    utils::download.file(SOURCE_URL, cache_path, mode = "wb", quiet = FALSE)
  } else {
    message("[cache] using existing download: ", cache_path)
  }

  # Row 7 of the "4) INDICADORES" sheet holds the real header (ANO, CÓD.,
  # SIGLA, LOCAL, POP_T, ...); rows 1-6 are a title block, data starts row 8.
  # UF rows use the standard 2-digit IBGE state codes (11-53); code 0 is
  # Brasil and 1-5 are the five macro-regions, both excluded here.
  raw <- readxl::read_excel(
    cache_path, sheet = "4) INDICADORES", col_names = FALSE, progress = FALSE
  )
  header <- as.character(raw[7, ])
  data <- raw[8:nrow(raw), ]
  names(data) <- header

  data |>
    dplyr::transmute(
      uf_code = as.integer(.data[["CÓD."]]),
      uf_sigla = .data[["SIGLA"]],
      uf_name = .data[["LOCAL"]],
      year = as.integer(.data[["ANO"]]),
      pop_total = as.numeric(.data[["POP_T"]]),
      pop_male = as.numeric(.data[["POP_H"]]),
      pop_female = as.numeric(.data[["POP_M"]]),
      births_total = as.numeric(.data[["NASC_T"]]),
      deaths_total = as.numeric(.data[["OBT_T"]]),
      natural_growth_rate = as.numeric(.data[["TCV"]])
    ) |>
    dplyr::filter(.data$uf_code >= 11, .data$uf_code <= 53) |>
    dplyr::arrange(.data$uf_code, .data$year)
}

if (sys.nframe() == 0) {
  cache_path <- here::here("01_Data/ibge_cache/projecoes_2024_tab4_indicadores.xlsx")
  population <- fetch_ibge_population_projection_uf(cache_path)

  out_path <- here::here("01_Data/ibge_population_projection_uf_2024revision.rds")
  saveRDS(population, out_path)
  readr::write_csv(population, sub("\\.rds$", ".csv", out_path))

  message(sprintf(
    "[save] %s\n  %s rows | %d UFs | %d-%d",
    out_path, format(nrow(population), big.mark = ","),
    dplyr::n_distinct(population$uf_code),
    min(population$year), max(population$year)
  ))
}
