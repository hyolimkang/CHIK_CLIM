# ===========================================================================
# 06_fetch_muni_centroids.R
#
# Purpose
# -------
# Download and save Brazil's full municipality (município) boundaries:
#   - full polygons (municipality boundary) — used for climate zonal
#     extraction (30_climate_covariates_dlnm/03_build_national_climate_weekly.R
#     uses exactextractr::exact_extract() to compute area-weighted means)
#   - centroids (lon/lat) — kept as a reference/visualization convenience
#
# Why polygon area-weighted mean instead of centroid
# ----------------------------------------------------
# Sampling a single centroid point per municipality can fail to represent
# the actual population/terrain, especially for large, irregularly-shaped
# municipalities like those in the Amazon. So the main analysis instead
# uses a zonal approach that weights ERA5-Land grid cells by how much area
# they share with the municipality boundary:
#   T_{m,t} = sum(A_mg * T_g,t) / sum(A_mg).
#
# Data source
# -----------
# Downloads the official IBGE malha municipal shapefile directly
# (geoftp.ibge.gov.br). The `geobr` package was tried first, but:
#   - geobr >= 2.0.0 internally uses `duckspatial::ddbs_open_dataset()`,
#     which references base R's `%||%` (added only in R 4.4.0), so it
#     breaks on R 4.3.x with "could not find function %||%".
#   - geobr 1.9.x (older, gpkg-based) points at a metadata server
#     (www.ipea.gov.br/geobr/metadata/metadata_gpkg.csv) that now 404s —
#     no longer served.
# Both paths are dead ends, so fetching the IBGE source directly is more
# reliable.
#
# Output
# ------
#   01_Data/ibge_muni_polygons.rds
#     sf object, columns: muni7, muni6, geometry (EPSG:4326, same CRS as
#     the ERA5-Land grid — exact_extract can use it as-is)
#   01_Data/ibge_muni_centroids.rds (+ .csv)
#     columns : muni7, muni6, lon, lat (centroid coordinates, reference only)
# ===========================================================================

for (p in c("here", "sf", "dplyr", "readr", "curl")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(sf); library(dplyr); library(readr)
})

fetch_muni_boundaries <- function(
    year = 2020,
    shp_url = sprintf(
      "https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_%d/Brasil/BR/BR_Municipios_%d.zip",
      year, year
    )
) {
  zip_path <- tempfile(fileext = ".zip")
  message("[ibge] downloading municipal shapefile (year = ", year, ") ...")
  h <- curl::new_handle(ssl_verifypeer = FALSE)
  curl::curl_download(shp_url, zip_path, handle = h, quiet = TRUE)

  unzip_dir <- tempfile("br_muni_")
  dir.create(unzip_dir)
  utils::unzip(zip_path, exdir = unzip_dir)

  shp_file <- list.files(unzip_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)[1]
  if (is.na(shp_file)) stop("No .shp found inside downloaded zip: ", shp_url)

  message("[sf] reading shapefile: ", basename(shp_file))
  muni_sf <- sf::st_read(shp_file, quiet = TRUE)

  code_col <- intersect(c("CD_MUN", "CD_GEOCMU", "code_muni"), names(muni_sf))[1]
  if (is.na(code_col)) {
    stop("Could not find a municipality-code column among: ", paste(names(muni_sf), collapse = ", "))
  }

  muni_sf <- muni_sf %>%
    transmute(
      muni7 = sprintf("%07s", as.character(.data[[code_col]])),
      muni6 = substr(muni7, 1, 6)
    ) %>%
    arrange(muni6)

  # ERA5-Land NetCDF grids are lon/lat EPSG:4326 — keep polygons in the
  # same CRS so exact_extract() doesn't need an on-the-fly reprojection.
  if (is.na(sf::st_crs(muni_sf)) || sf::st_crs(muni_sf) != sf::st_crs(4326)) {
    muni_sf <- sf::st_transform(muni_sf, 4326)
  }

  message("[sf] computing centroids ...")
  centroid_coords <- sf::st_coordinates(suppressWarnings(sf::st_centroid(muni_sf)))

  centroids <- muni_sf %>%
    sf::st_drop_geometry() %>%
    mutate(lon = centroid_coords[, "X"], lat = centroid_coords[, "Y"])

  list(polygons = muni_sf, centroids = centroids)
}

if (sys.nframe() == 0) {
  boundaries <- fetch_muni_boundaries()

  poly_path <- here::here("01_Data/ibge_muni_polygons.rds")
  saveRDS(boundaries$polygons, poly_path)

  centroid_rds <- here::here("01_Data/ibge_muni_centroids.rds")
  centroid_csv <- sub("\\.rds$", ".csv", centroid_rds)
  saveRDS(boundaries$centroids, centroid_rds)
  readr::write_csv(boundaries$centroids, centroid_csv)

  message(sprintf(
    "[save] municipality boundaries\n  %s (polygons)\n  %s\n  %s\n  %d municipalities",
    poly_path, centroid_rds, centroid_csv, nrow(boundaries$centroids)
  ))
}
