# ------------------------------------------------------------------------------
# Pipeline-wide libraries
# Order matters: load `raster` LAST so its functions take precedence over
# `terra` for this script (which uses raster::raster/projectRaster/resample/stack).
# Other scripts that prefer terra should load only `terra` (no `raster`).
# ------------------------------------------------------------------------------
# Data wrangling
library(dplyr)
library(readr)
library(data.table)
library(tibble)

# Spatial (sf, terra optional; raster goes last)
library(sf)
library(geodata)
library(terra)
library(raster)
library(FNN)

# FOI modelling pipeline
library(randomForest)

# Parallel execution (foreach + %dopar%)
library(foreach)
library(doParallel)
library(parallel)

# Plotting
library(ggplot2)

options(scipen = 999)

# ------------------------------------------------------------------------------
# 0) Paths and helpers
# ------------------------------------------------------------------------------
# Auto-detect project root even when run from 01_Data.
cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
project_root <- if (basename(cwd) == "01_Data") dirname(cwd) else cwd

resolve_path <- function(candidates) {
  candidates <- unlist(candidates, recursive = TRUE, use.names = FALSE)
  if (is.null(candidates) || length(candidates) == 0) {
    stop("Path candidates are NULL or empty.")
  }
  if (!is.character(candidates)) {
    stop("Path candidates must be a character vector.")
  }
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0) {
    stop(
      "No file found among candidates:\n",
      paste0(" - ", candidates, collapse = "\n")
    )
  }
  normalizePath(hits[[1]], winslash = "/", mustWork = TRUE)
}

raster_scale <- function(ras, log_transform = TRUE) {
  vec <- as.vector(ras)
  if (log_transform) {
    zinf <- min(vec[vec > 0], na.rm = TRUE)
    vec <- log(vec + 0.5 * zinf)
  }
  vec <- scale(vec)
  values(ras) <- vec
  ras
}

raster_scale_neg <- function(ras, log_transform = TRUE) {
  vec <- as.vector(ras)
  if (log_transform) {
    min_val <- min(vec, na.rm = TRUE)
    offset <- abs(min_val) + 1
    vec <- log(vec + offset)
  }
  vec <- scale(vec)
  values(ras) <- vec
  ras
}

apply_transform <- function(ras, transform_key) {
  if (transform_key == "scale_log") return(raster_scale(ras, log_transform = TRUE))
  if (transform_key == "scale_no_log") return(raster_scale(ras, log_transform = FALSE))
  if (transform_key == "scale_neg_log") return(raster_scale_neg(ras, log_transform = TRUE))
  if (transform_key == "none") return(ras)
  stop("Unknown transform key: ", transform_key)
}

source("02_Script/Functions/fixNAs_adj.R")
source("02_Script/Functions/thinning_chik.R")

# ------------------------------------------------------------------------------
# 1) Covariate catalog (current baseline)
# ------------------------------------------------------------------------------
# Keep this table as the single source of truth for current baseline covariates.
covariate_catalog <- tibble::tribble(
  ~covariate_name, ~raster_name, ~type, ~transform, ~resample_method, ~candidate_paths,
  "Tsuit", "tsuit", "climate", "scale_log", "template", list(c(
    file.path(project_root, "01_Data", "final", "Dengue_temperature_suitaiblity_masked_.tif")
  )),
  "Tmin", "tmin", "climate", "scale_no_log", "project_bilinear", list(c(
    file.path(project_root, "01_Data", "final", "Tmin_TerraClim_2010_2020_005dg_masked_.tif")
  )),
  "PRCP", "precip", "climate", "scale_log", "project_bilinear", list(c(
    file.path(project_root, "01_Data", "final", "PRCP_TerraClim_2010_2020_005dg_masked_.tif")
  )),
  "Pop_dens", "pop_dens", "demography", "scale_log", "project_bilinear", list(c(
    file.path(project_root, "01_Data", "final", "landscan_global_2022_masked.tif")
  )),
  "GDP", "gdp_nat", "socioeconomic", "scale_log", "project_bilinear", list(c(
    file.path(project_root, "01_Data", "final", "GDP_2009_2019_National_masked.tif"),
    file.path(project_root, "01_Data", "final", "GDP_2009_2019_National_masked_.tif")
  )),
  "Albo", "albo", "vector", "scale_no_log", "project_bilinear", list(c(
    file.path(project_root, "01_Data", "final", "Albopictus_mean_2020_rcp60_spreadXsuit_masked_.tif")
  )),
  "Aegyp", "aegyp", "vector", "scale_no_log", "project_bilinear", list(c(
    file.path(project_root, "01_Data", "final", "Aegypti_mean_2020_rcp60_spreadXsuit_masked_.tif")
  )),
  "CHIKRisk", "chikrisk", "risk_proxy", "scale_log", "resample_bilinear", list(c(
    file.path(project_root, "01_Data", "final", "CHIK_riskmap_wmean_masked.tif"),
    file.path(project_root, "01_Data", "final", "CHIK_rangemap_100pred_2024_02_06.tif")
  )),
  "GDP_cap", "gdp_per", "socioeconomic", "scale_log", "resample_bilinear", list(c(
    file.path(project_root, "01_Data", "final", "GDP_2009_2019_1km_masked_.tif")
  )),
  "pop_dens", "pop_dens_dup", "demography", "scale_log", "project_bilinear", list(c(
    file.path(project_root, "01_Data", "final", "landscan_global_2022_masked.tif")
  )),
  "NDVI", "ndvi", "environment", "scale_no_log", "project_bilinear", list(c(
    file.path(project_root, "01_Data", "final", "NDVI_2010_2020_005dg_masked_.tif")
  )),
  "DHI", "dhi", "health_system", "none", "resample_bilinear", list(c(
    file.path(project_root, "01_Data", "final", "DHI_global_clusters_14c_sv20_masked_.tif")
  ))
)

# Export a readable inventory for planning / future matching.
write_csv(
  covariate_catalog %>% dplyr::select(covariate_name, type, transform, resample_method),
  file.path(project_root, "01_Data", "current_covariate_catalog.csv")
)

# ------------------------------------------------------------------------------
# 2) Load + transform covariates
# ------------------------------------------------------------------------------
covariates <- list()

for (i in seq_len(nrow(covariate_catalog))) {
  row_i <- covariate_catalog[i, ]
  raster_path <- resolve_path(row_i$candidate_paths[[1]])
  message("Loading covariate: ", row_i$covariate_name, " | ", raster_path)
  ras <- raster::raster(raster_path)
  ras <- apply_transform(ras, row_i$transform)
  covariates[[row_i$covariate_name]] <- ras
}

# Set template from Tsuit
template <- covariates[["Tsuit"]]

# ------------------------------------------------------------------------------
# 3) Align rasters to template
# ------------------------------------------------------------------------------
align_to_template <- function(ras, method_key, template_raster) {
  if (method_key == "template") return(ras)
  if (method_key == "project_bilinear") return(projectRaster(ras, template_raster, method = "bilinear"))
  if (method_key == "resample_bilinear") return(resample(ras, template_raster, method = "bilinear"))
  stop("Unknown resample method: ", method_key)
}

for (i in seq_len(nrow(covariate_catalog))) {
  nm <- covariate_catalog$covariate_name[i]
  method_i <- covariate_catalog$resample_method[i]
  covariates[[nm]] <- align_to_template(covariates[[nm]], method_i, template)
}

# ------------------------------------------------------------------------------
# 4) Create final stack (matches old covlist_upd order)
# ------------------------------------------------------------------------------
final_order <- c("Tsuit", "Tmin", "PRCP", "GDP", "Albo", "Aegyp", "CHIKRisk", "GDP_cap", "pop_dens", "NDVI", "DHI")
covlist_upd <- stack(covariates[final_order])
names(covlist_upd) <- final_order
plot(covlist_upd)

save(
  covlist_upd,
  file = file.path(project_root, "01_Data", "current_covariates_upd.RData")
)

message("Saved:")
message(" - ", file.path(project_root, "01_Data", "current_covariate_catalog.csv"))
message(" - ", file.path(project_root, "01_Data", "current_covariates_upd.RData"))



# ------------------------------------------------------------------------------
# 5) FOI data 
# ------------------------------------------------------------------------------
# extract covariate values
load("01_Data/current_covariates_upd.RData")

# extract FOI and occurrence values
chik_foi<- read.csv("01_Data/chik_foi.csv")
load("01_Data/chik_occ.RData")

# create only occurrence + foi
chik_foi <- chik_foi[,c(5:9, 15)]
new_order <- c("long", "lat", "mfoi")
chik_foi <- chik_foi[,new_order]
chik_occ <- chik_occ[,c(1:2)]
names(chik_occ) <- c("long", "lat")
chik_foi_merge <- chik_foi[,c(1:2)]
chik_occ <- rbind(chik_foi_merge, chik_occ)

template = tsuit

p_covs <- extract(covlist_upd, data.frame(chik_foi$long, chik_foi$lat), df= T, na.rm=TRUE)
p_covs$FOI <- chik_foi$mfoi
p_covs$FOI_lo <- chik_foi$mfoi_lo
p_covs$FOI_hi <- chik_foi$mfoi_hi
p_covs$Longitude <- chik_foi$long
p_covs$Latitude <- chik_foi$lat
p_covs$NumTested <- chik_foi$NumTested
p_covs$region <- chik_foi$region
p_covs$study_no <- chik_foi$study_no
p_covs$country <- chik_foi$country
p_covs <- fixNAs(p_covs, covlist_upd)
(point_idx <- which(apply(p_covs[, names(covlist_upd)], 1, function(row) any(is.na(row)))))


# ------------------------------------------------------------------------------
# 6) kNN imputation for any FOI point that still has NA covariates
#
# For each row with NA in any covariate column, find the k nearest neighbours
# (in lon/lat space) and fill missing covariates with the neighbour mean.
# This was the imputation step in old_covariates.R, generalised here so it is
# robust to NAs appearing in any row (not only row 50).
# ------------------------------------------------------------------------------
p_covs_sf <- st_as_sf(p_covs, coords = c("Longitude", "Latitude"), remove = FALSE)

cov_cols <- names(covlist_upd)
na_rows  <- which(apply(p_covs[, cov_cols], 1, function(row) any(is.na(row))))

if (length(na_rows) > 0) {
  message(sprintf("kNN imputing %d row(s) with NA covariates: %s",
                  length(na_rows), paste(na_rows, collapse = ", ")))

  coords_matrix <- st_coordinates(p_covs_sf)
  k_neighbours  <- 3

  for (na_row in na_rows) {
    # Use only rows that are themselves complete as the reference pool.
    complete_idx <- setdiff(seq_len(nrow(p_covs)), na_rows)
    nn_set <- get.knnx(
      data  = coords_matrix[complete_idx, , drop = FALSE],
      query = coords_matrix[na_row, , drop = FALSE],
      k     = k_neighbours
    )
    nn_global_idx <- complete_idx[nn_set$nn.index[1, ]]

    for (cv in cov_cols) {
      if (is.na(p_covs[[cv]][na_row])) {
        p_covs[[cv]][na_row] <- mean(p_covs[[cv]][nn_global_idx], na.rm = TRUE)
      }
    }
  }

  remaining <- which(apply(p_covs[, cov_cols], 1, function(row) any(is.na(row))))
  if (length(remaining) == 0) {
    message("All covariate NAs successfully imputed.")
  } else {
    warning(sprintf("Remaining NA rows after kNN imputation: %s",
                    paste(remaining, collapse = ", ")))
  }
} else {
  message("No NA covariates in p_covs; imputation step skipped.")
}

# Refresh sf object with imputed values
p_covs_sf <- st_as_sf(p_covs, coords = c("Longitude", "Latitude"), remove = FALSE)

ggplot(data = p_covs_sf) +
  geom_sf() +
  geom_sf_text(aes(label = ID), check_overlap = TRUE, nudge_y = 0.01) +
  theme_minimal()

