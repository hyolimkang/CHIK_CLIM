library(geodata)
library(terra)

# ------------------------------------------------------------------------------
# WorldClim downloads for CHIK FOI future projection 
#
# Variable rationale (matched to current_covariates_baseline.R):
#   - Tsuit (Dengue temperature suitability) -> need monthly tavg, apply same
#                                               suitability function as current
#   - Tmin  (TerraClim 2010-2020)            -> monthly tmin, then annual mean
#   - PRCP  (TerraClim 2010-2020)            -> monthly prec, then annual sum
#
# We also pull `bioc` for sensitivity / extrapolation checks (bio1, bio6, bio12).
#
# Current TerraClim baseline = 2010-2020, but WorldClim historical = 1970-2000.
# We download WorldClim historical so the future projection can be done as a
# *change-factor (delta)* applied on top of the TerraClim baseline:
#
#   Tmin_future_proj  = Tmin_TerraClim + (Tmin_WC_future - Tmin_WC_hist)
#   PRCP_future_proj  = PRCP_TerraClim * (PRCP_WC_future / PRCP_WC_hist)
#   Tsuit_future_proj = suitability_fn(tavg_WC_future)   # recomputed
# ------------------------------------------------------------------------------

out_root <- "01_Data/worldclim_cmip6"
hist_root <- "01_Data/worldclim_historical"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(hist_root, recursive = TRUE, showWarnings = FALSE)

# ---- Config ------------------------------------------------------------------
gcms     <- c("MIROC6", "MPI-ESM1-2-HR", "GFDL-ESM4")  # MVP: 3-GCM ensemble
ssps     <- c("245", "585")                            # main scenarios
periods  <- c("2041-2060")                             # 2050s window
res_deg  <- 2.5                                        # arc-minutes (~5 km eq.)

# Variables we actually need:
#   "tavg" -> Tsuit recompute
#   "tmin" -> Tmin matching
#   "prec" -> PRCP matching
#   "bioc" -> sensitivity (bio1/bio6/bio12 etc.)
vars_future <- c("tavg", "tmin", "prec", "bioc")

# ---- Future downloads --------------------------------------------------------
download_future <- function(model, ssp, time, var, res, path) {
  message(sprintf("Downloading future: %s | SSP%s | %s | %s", model, ssp, time, var))
  cmip6_world(
    model = model,
    ssp   = ssp,
    time  = time,
    var   = var,
    res   = res,
    path  = path
  )
}

future_stacks <- list()
for (gcm in gcms) {
  for (ssp in ssps) {
    for (per in periods) {
      for (v in vars_future) {
        key <- paste(gcm, ssp, per, v, sep = "|")
        future_stacks[[key]] <- tryCatch(
          download_future(gcm, ssp, per, v, res_deg, out_root),
          error = function(e) {
            message("  FAILED: ", conditionMessage(e))
            NULL
          }
        )
      }
    }
  }
}

# ---- Historical baseline (1970-2000) for delta calculation -------------------
# These are the WorldClim "current" climatologies that match the GCM historical.
vars_hist <- c("tavg", "tmin", "prec", "bio")  # `bio` for historical bioclim
hist_stacks <- list()
for (v in vars_hist) {
  message(sprintf("Downloading historical: %s", v))
  hist_stacks[[v]] <- tryCatch(
    worldclim_global(var = v, res = res_deg, path = hist_root),
    error = function(e) {
      message("  FAILED: ", conditionMessage(e))
      NULL
    }
  )
}

message("\nDone. Future files in: ", normalizePath(out_root))
message("Historical files in: ",      normalizePath(hist_root))
