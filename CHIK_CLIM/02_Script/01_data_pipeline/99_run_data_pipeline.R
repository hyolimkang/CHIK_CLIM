# ============================================================
# DATA PIPELINE -- TOP-LEVEL ORCHESTRATOR
# ============================================================
#
# Sources each stage's fetch/build scripts in order. Every fetch script
# guards against re-fetching when its destination file already exists
# (e.g. 01_surveillance/01_fetch_sinan.R checks file.exists(dest) before
# hitting the network), so re-running this is safe and will not
# unnecessarily repeat expensive downloads.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this pipeline.")
  source(setup_path, local = .GlobalEnv)
}

STAGE_DIR <- project_path("02_Script/01_data_pipeline")

message("=== 01/05: surveillance (SINAN fetch + clean) ===")
source(file.path(STAGE_DIR, "01_surveillance/01_fetch_sinan.R"), local = .GlobalEnv)
source(file.path(STAGE_DIR, "01_surveillance/02_clean_sinan.R"), local = .GlobalEnv)

message("\n=== 02/05: population (IBGE population, projection, age population) ===")
source(file.path(STAGE_DIR, "02_population/01_fetch_population.R"), local = .GlobalEnv)
source(file.path(STAGE_DIR, "02_population/02_fetch_population_projection.R"), local = .GlobalEnv)
source(file.path(STAGE_DIR, "02_population/03_fetch_age_population.R"), local = .GlobalEnv)

message("\n=== 03/05: demography (per-state SINASC births -> weekly demography) ===")
for (f in list.files(file.path(STAGE_DIR, "03_demography"), pattern = "^01_.*\\.R$", full.names = TRUE)) source(f, local = .GlobalEnv)
for (f in list.files(file.path(STAGE_DIR, "03_demography"), pattern = "^03_.*\\.R$", full.names = TRUE)) source(f, local = .GlobalEnv)

message("\n=== 04/05: geography (municipality centroids) ===")
source(file.path(STAGE_DIR, "04_geography/01_fetch_municipality_centroids.R"), local = .GlobalEnv)

message("\n=== 05/05: climate (ERA5-Land fetch -> national weekly climate -> DLNM panel) ===")
source(file.path(STAGE_DIR, "05_climate/01_fetch_era5_temperature.R"), local = .GlobalEnv)
source(file.path(STAGE_DIR, "05_climate/02_fetch_era5_precipitation.R"), local = .GlobalEnv)
source(file.path(STAGE_DIR, "05_climate/03_build_national_weekly_climate.R"), local = .GlobalEnv)
source(file.path(STAGE_DIR, "05_climate/04_build_climate_panel.R"), local = .GlobalEnv)

message("\n[data_pipeline] Complete.")
