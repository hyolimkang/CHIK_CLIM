# National (27-UF) weekly climate table, built by looping the existing
# single-state build_state_weekly() helper over every Brazilian UF.
#
# This is pure aggregation of already-validated inputs (the national DLNM
# muni-week panel + the IBGE 2024-revision UF population projection) -- no
# new data source, no modification of build_state_weekly() itself.

required_climate_packages <- c("here", "dplyr", "readr")
missing_climate_packages <- required_climate_packages[
  !vapply(required_climate_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_climate_packages)) stop("Missing package(s): ", paste(missing_climate_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr) })

climate_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}
source(file.path(climate_root(), "02_Script", "40_renewal_model", "00_shared", "01_build_ceara_state_weekly.R"))
source(file.path(climate_root(), "02_Script", "40_renewal_model", "15_national_wave_analysis", "00_national_wave_analysis_helpers.R"))

run_build_brazil_uf_weekly_climate <- function() {
  paths <- ensure_national_output_dirs()
  panel <- readRDS(file.path(climate_root(), "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  uf_lookup <- readRDS(file.path(climate_root(), "01_Data", "ibge_population_projection_uf_2024revision.rds")) |>
    transmute(uf_code = sprintf("%02d", as.integer(uf_code)), state = uf_sigla) |>
    distinct() |>
    arrange(state)
  stopifnot(nrow(uf_lookup) == 27L)

  message("[uf-weekly-climate] building for ", nrow(uf_lookup), " UFs ...")
  climate <- bind_rows(lapply(seq_len(nrow(uf_lookup)), function(i) {
    state <- uf_lookup$state[i]
    uf_code <- uf_lookup$uf_code[i]
    message("  ", state, " (", uf_code, ")")
    build_state_weekly(panel, uf_code) |>
      transmute(state = state, week_start, year, Tmean, PRCP)
  }))

  if (n_distinct(climate$state) != 27L) stop("Expected 27 states, got ", n_distinct(climate$state))
  if (anyNA(climate[c("Tmean", "PRCP")])) {
    n_na <- sum(is.na(climate$Tmean) | is.na(climate$PRCP))
    warning(n_na, " state-weeks have NA Tmean/PRCP -- likely municipalities with no climate coverage; will propagate as NA downstream.")
  }

  out_path <- file.path(paths$table, "brazil_chik_uf_weekly_climate.csv")
  write_csv(climate, out_path)
  message("[uf-weekly-climate] saved: ", out_path, " (", nrow(climate), " rows, ", n_distinct(climate$state), " states)")
  invisible(climate)
}

if (sys.nframe() == 0L) run_build_brazil_uf_weekly_climate()
