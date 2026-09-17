library(geodata)
library(terra)

# ------------------------------------------------------------------------------
# Pilot: build future PRCP raster aligned to TerraClim 2010-2020 baseline
#
# Goal of this pilot
#   Verify the *alignment + delta downscaling pipeline* end-to-end on ONE
#   variable (PRCP) before scaling to the full ensemble.
#
# Pipeline
#   1. Load current PRCP (TerraClim 2010-2020) and current template (Tsuit).
#   2. Download WorldClim historical monthly precipitation (1970-2000).
#   3. Download WorldClim CMIP6 future monthly precipitation
#      (MIROC6, SSP245 and SSP585, 2041-2060).
#   4. Aggregate monthly -> annual total (Σ over 12 months) for both periods.
#   5. Compute multiplicative delta = future_annual / hist_annual on WC grid.
#   6. Resample delta to current template.
#   7. Future PRCP = TerraClim PRCP × delta.   (still in raw mm/year)
#   8. Save future PRCP rasters for SSP245 and SSP585.
#
# Notes
#   - Multiplicative delta is the standard for precipitation
#     (avoids negative values, preserves spatial dryness pattern).
#   - This pilot keeps PRCP in raw units (mm/year). Scaling/log transform for
#     the RF model is applied later, using the TRAINING mean/sd, not future's.
#   - Diagnostics are printed/plotted at the end so we can sanity check the
#     pipeline before extending to other variables.
# ------------------------------------------------------------------------------

# ---- Paths -------------------------------------------------------------------
cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
project_root <- if (basename(cwd) == "01_Data") dirname(cwd) else cwd

terraclim_prcp_path <- file.path(
  project_root, "01_Data", "final",
  "PRCP_TerraClim_2010_2020_005dg_masked_.tif"
)
template_path <- file.path(
  project_root, "01_Data", "final",
  "Dengue_temperature_suitaiblity_masked_.tif"
)

wc_hist_dir   <- file.path(project_root, "01_Data", "worldclim_historical")
wc_future_dir <- file.path(project_root, "01_Data", "worldclim_cmip6")
out_dir       <- file.path(project_root, "01_Data", "future_pilot")

dir.create(wc_hist_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(wc_future_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir,       recursive = TRUE, showWarnings = FALSE)

# ---- Pilot config ------------------------------------------------------------
gcm     <- "MIROC6"
ssps    <- c("245", "585")
period  <- "2041-2060"
res_deg <- 2.5  # arc-minutes (~5 km)

# ---- 1) Current PRCP + template ---------------------------------------------
message("[1/8] Loading current PRCP and template...")
prcp_current <- rast(terraclim_prcp_path)
template     <- rast(template_path)

# ---- 2) Historical monthly prec (1970-2000) ---------------------------------
message("[2/8] Downloading WorldClim historical monthly prec (1970-2000)...")
wc_hist_monthly <- worldclim_global(
  var  = "prec",
  res  = res_deg,
  path = wc_hist_dir
)  # 12-layer SpatRaster: prec_01..prec_12

# ---- 3) Future monthly prec (one GCM × two SSPs × 2041-2060) ----------------
download_future_prec <- function(model, ssp, time, res, path) {
  message(sprintf("[3/8] Downloading future prec: %s | SSP%s | %s ...",
                  model, ssp, time))
  cmip6_world(
    model = model, ssp = ssp, time = time,
    var = "prec", res = res, path = path
  )
}
wc_fut_monthly <- list()
for (ssp in ssps) {
  wc_fut_monthly[[ssp]] <- download_future_prec(
    gcm, ssp, period, res_deg, wc_future_dir
  )
}

# ---- 4) Aggregate monthly -> annual total -----------------------------------
message("[4/8] Aggregating monthly -> annual total precipitation...")
hist_annual <- sum(wc_hist_monthly)            # mm/year (climatology)
fut_annual  <- lapply(wc_fut_monthly, function(x) sum(x))

# ---- 5) Multiplicative delta on WC grid -------------------------------------
message("[5/8] Computing multiplicative delta = future_annual / hist_annual...")
# Avoid divide-by-zero: clamp historical at a tiny floor.
hist_floor <- 1  # mm/year
hist_annual_safe <- ifel(hist_annual < hist_floor, hist_floor, hist_annual)

delta_list <- lapply(fut_annual, function(fa) fa / hist_annual_safe)
names(delta_list) <- ssps

# ---- 6) Resample delta to current template ----------------------------------
message("[6/8] Resampling delta to current template...")
delta_on_template <- lapply(delta_list, function(d) {
  resample(d, template, method = "bilinear")
})

# ---- 7) Apply delta to TerraClim baseline -----------------------------------
message("[7/8] Computing future PRCP = TerraClim × delta ...")
prcp_current_on_template <- resample(prcp_current, template, method = "bilinear")
prcp_future_list <- lapply(delta_on_template, function(d) prcp_current_on_template * d)
names(prcp_future_list) <- ssps

# ---- 8) Save outputs --------------------------------------------------------
message("[8/8] Saving future PRCP rasters...")
for (ssp in ssps) {
  out_path <- file.path(
    out_dir,
    sprintf("PRCP_future_%s_SSP%s_%s.tif", gcm, ssp, period)
  )
  writeRaster(prcp_future_list[[ssp]], out_path, overwrite = TRUE)
  message("  saved: ", out_path)
}

# ---- Diagnostics (base R, namespace-safe) -----------------------------------
summary_vec <- function(r, label) {
  v <- terra::values(r)
  v <- v[!is.na(v)]
  cat(sprintf("  %-22s min=%.3f  median=%.3f  mean=%.3f  max=%.3f\n",
              label, min(v), stats::median(v), mean(v), max(v)))
}

message("\n=== Diagnostics ===")
for (ssp in ssps) {
  cat(sprintf("\nSSP%s:\n", ssp))
  summary_vec(delta_on_template[[ssp]],     "delta (future/hist)")
  summary_vec(prcp_current_on_template,     "PRCP current (mm/yr)")
  summary_vec(prcp_future_list[[ssp]],      "PRCP future  (mm/yr)")
}

# ---- Delta visualization (saved to PNG) -------------------------------------
message("\n[viz] Saving delta map and histogram PNGs...")
viz_dir <- file.path(out_dir, "viz")
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

for (ssp in ssps) {
  d  <- delta_on_template[[ssp]]
  pc <- prcp_current_on_template
  pf <- prcp_future_list[[ssp]]

  # 1) Delta map + histogram
  png(file.path(viz_dir, sprintf("delta_SSP%s.png", ssp)),
      width = 1400, height = 900, res = 130)
  par(mfrow = c(2, 2), mar = c(3, 3, 3, 4))
  plot(d, main = sprintf("Multiplicative delta (SSP%s)", ssp),
       col = hcl.colors(50, "RdBu"))
  hist(terra::values(d), breaks = 80,
       main = sprintf("Delta distribution (SSP%s)", ssp),
       xlab = "future / historical")
  abline(v = 1, col = "red", lwd = 2, lty = 2)
  plot(pc, main = "PRCP current (TerraClim 2010-2020)",
       col = hcl.colors(50, "Blues", rev = TRUE))
  plot(pf, main = sprintf("PRCP future (SSP%s, 2041-2060)", ssp),
       col = hcl.colors(50, "Blues", rev = TRUE))
  dev.off()
}

message("Saved viz to: ", viz_dir)
