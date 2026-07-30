# fetch_ce_births_ibge.R
# Annual live births from IBGE SIDRA Registro Civil, table 2609.
# Variable: "Nascidos vivos registrados no ano" (total, all maternal ages).
# Generic: works for any Brazilian UF — set uf_code and uf_name in the
# config block at the bottom before running interactively.

fetch_uf_births_annual <- function(uf_code = "23", uf_name = "ce",
                                    year_start = 2014, year_end = 2024) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    install.packages("jsonlite")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  if (!requireNamespace("here", quietly = TRUE)) {
    install.packages("here")
  }
  if (!requireNamespace("readr", quietly = TRUE)) {
    install.packages("readr")
  }

  period <- sprintf("%d-%d", year_start, year_end)
  url <- sprintf(
    "https://apisidra.ibge.gov.br/values/t/2609/n3/%s/p/%s/v/allxp/f/u?formato=json",
    uf_code, period
  )

  raw <- jsonlite::fromJSON(url)
  if (nrow(raw) < 2) {
    stop(sprintf("SIDRA returned no rows for UF%s births.", uf_code))
  }

  out <- raw[-1, ] |>
    dplyr::transmute(
      year   = as.integer(D2N),
      births = as.integer(V),
      source = sprintf("IBGE SIDRA t2609 UF%s", uf_code)
    ) |>
    dplyr::arrange(year)

  out
}

if (sys.nframe() == 0) {
  # ============================================================
  # USER CONFIG 
  # ============================================================
  UF_CODE    <- "31"    # CE=23, MG=31, SP=35, RJ=33, BA=29, ...
  UF_NAME    <- "MG"    
  YEAR_START <- 2014
  YEAR_END   <- 2024
  # ============================================================

  births  <- fetch_uf_births_annual(UF_CODE, UF_NAME, YEAR_START, YEAR_END)
  out_csv <- here::here(sprintf("01_Data/%s_births_annual_ibge.csv", UF_NAME))
  readr::write_csv(births, out_csv)
  saveRDS(births, sub("\\.csv$", ".rds", out_csv))
  message("Saved: ", out_csv)
  print(births)
}
