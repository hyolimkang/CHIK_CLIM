library(terra)
library(randomForest)

# ------------------------------------------------------------------------------
# Pilot prediction: current vs future PRCP, configurable ensemble size & SSPs.
#
# Quick mode (default)  : N_MODELS_PILOT = 1   SSPS_PILOT = "245"
# Full  mode            : N_MODELS_PILOT = 100 SSPS_PILOT = c("245", "585")
#
# Inputs assumed in environment OR loadable:
#   - rf_mod_hyper : list of trained RF models (built in 02_Script/precip_pilot.R)
#   - covlist_upd  : current covariate stack (saved by 01_Data/current_covariates_baseline.R)
#   - future PRCP rasters from 01_Data/build_future_PRCP_pilot.R
#   - raw TerraClim PRCP raster (for recovering training scaling parameters)
# ------------------------------------------------------------------------------

# ============================== CONFIG ========================================
N_MODELS_PILOT <- 1            # 1 for quick test, 100 for full ensemble
SSPS_PILOT     <- c("245")     # c("245") for quick, c("245","585") for full
GCM            <- "MIROC6"     # only one GCM in this pilot
PERIOD         <- "2041-2060"
# ==============================================================================

# ---- Paths -------------------------------------------------------------------
cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
project_root <- if (basename(cwd) %in% c("01_Data", "02_Script")) dirname(cwd) else cwd

current_stack_path <- file.path(project_root, "01_Data", "current_covariates_upd.RData")
prcp_raw_path      <- file.path(project_root, "01_Data", "final",
                                "PRCP_TerraClim_2010_2020_005dg_masked_.tif")
fut_dir            <- file.path(project_root, "01_Data", "future_pilot")
out_dir            <- file.path(project_root, "01_Data", "future_pilot", "predict")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tag <- sprintf("n%d_%s", N_MODELS_PILOT, paste0("SSP", SSPS_PILOT, collapse = "_"))

# ---- 1) Make sure rf_mod_hyper exists, then subset --------------------------
if (!exists("rf_mod_hyper")) {
  stop(
    "`rf_mod_hyper` not found in environment.\n",
    "  Run 02_Script/precip_pilot.R first to train the RF models,\n",
    "  or load a saved RData that contains it."
  )
}
if (!is.list(rf_mod_hyper) || length(rf_mod_hyper) < 1) {
  stop("rf_mod_hyper must be a non-empty list of trained RF models.")
}
n_avail <- length(rf_mod_hyper)
if (N_MODELS_PILOT > n_avail) {
  stop(sprintf("N_MODELS_PILOT=%d but only %d models in rf_mod_hyper.",
               N_MODELS_PILOT, n_avail))
}
models_pilot <- rf_mod_hyper[seq_len(N_MODELS_PILOT)]
message(sprintf("[1/7] Using %d / %d trained RF model(s).", N_MODELS_PILOT, n_avail))

# ---- 2) Load current covariate stack ----------------------------------------
message("[2/7] Loading current covariate stack...")
load(current_stack_path)                     # loads `covlist_upd`
covlist_current <- terra::rast(covlist_upd)  # convert to SpatRaster

req_vars <- c("Tsuit", "PRCP", "GDP", "Albo", "Aegyp", "CHIKRisk")
missing <- setdiff(req_vars, names(covlist_current))
if (length(missing) > 0) stop("Missing covariate layers: ", paste(missing, collapse = ", "))
covlist_current <- covlist_current[[req_vars]]   # keep only what RF expects

# ---- 3) Recover PRCP scaling parameters from raw TerraClim PRCP -------------
# IMPORTANT: must use TRAINING parameters, NOT future's own mean/sd.
message("[3/7] Recovering PRCP scaling parameters from raw TerraClim raster...")
prcp_raw <- terra::rast(prcp_raw_path)
v_raw    <- terra::values(prcp_raw); v_raw <- v_raw[!is.na(v_raw)]
zinf     <- min(v_raw[v_raw > 0])
v_log    <- log(v_raw + 0.5 * zinf)
prcp_mean <- mean(v_log)
prcp_sd   <- stats::sd(v_log)
cat(sprintf("   zinf=%.4f  log_mean=%.4f  log_sd=%.4f\n", zinf, prcp_mean, prcp_sd))

scale_like_training <- function(prcp_raster_raw_units) {
  log_r <- terra::app(prcp_raster_raw_units,
                      fun = function(x) log(x + 0.5 * zinf))
  (log_r - prcp_mean) / prcp_sd
}

# ---- 4) Build future stacks (PRCP swapped, scaled identically) --------------
build_future_stack <- function(ssp_label) {
  fut_path <- file.path(
    fut_dir,
    sprintf("PRCP_future_%s_SSP%s_%s.tif", GCM, ssp_label, PERIOD)
  )
  if (!file.exists(fut_path)) {
    stop("Future PRCP file not found: ", fut_path,
         "\n  Run 01_Data/build_future_PRCP_pilot.R first.")
  }
  prcp_fut_raw     <- terra::rast(fut_path)
  prcp_fut_scaled  <- scale_like_training(prcp_fut_raw)
  prcp_fut_aligned <- terra::resample(prcp_fut_scaled,
                                      covlist_current[["PRCP"]],
                                      method = "bilinear")
  out <- covlist_current
  out[["PRCP"]] <- prcp_fut_aligned
  out
}

message(sprintf("[4/7] Building future stacks for: %s",
                paste0("SSP", SSPS_PILOT, collapse = ", ")))
future_stacks <- lapply(SSPS_PILOT, build_future_stack)
names(future_stacks) <- SSPS_PILOT

# ---- 5) Ensemble prediction (chunked, summary stats) ------------------------
# For each grid cell, run all `models` predictions and summarise into:
#   median, lo (2.5%), hi (97.5%), sd
# This mirrors the old workflow that built foi_mid / foi_lo / foi_hi / foi_sd
# from a 100-column matrix of per-model predictions.
# Works for any number of models (with N=1, lo=hi=median and sd=NA).
predict_ensemble_summary <- function(stk, models, chunk_size = 50000) {
  df <- terra::as.data.frame(stk, xy = FALSE, na.rm = FALSE)
  n_cells <- nrow(df)
  med <- lo <- hi <- sd_v <- rep(NA_real_, n_cells)

  ok_idx <- which(stats::complete.cases(df))
  if (length(ok_idx) == 0) return(list(median = med, lo = lo, hi = hi, sd = sd_v))

  chunks <- split(ok_idx, ceiling(seq_along(ok_idx) / chunk_size))
  pb <- txtProgressBar(min = 0, max = length(chunks), style = 3)
  for (k in seq_along(chunks)) {
    rows  <- chunks[[k]]
    chunk <- df[rows, , drop = FALSE]
    # 100-column (or N-column) matrix of per-model predictions.
    preds_mat <- do.call(cbind, lapply(models, function(m) {
      predict(m, newdata = chunk, type = "response")
    }))

    med[rows]  <- apply(preds_mat, 1, stats::median, na.rm = TRUE)
    if (ncol(preds_mat) >= 2) {
      lo[rows]   <- apply(preds_mat, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
      hi[rows]   <- apply(preds_mat, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
      sd_v[rows] <- apply(preds_mat, 1, stats::sd,        na.rm = TRUE)
    } else {
      lo[rows]   <- med[rows]
      hi[rows]   <- med[rows]
      sd_v[rows] <- NA_real_
    }
    setTxtProgressBar(pb, k)
  }
  close(pb)
  list(median = med, lo = lo, hi = hi, sd = sd_v)
}

raster_from_vec <- function(template_stk, vec, name) {
  r <- terra::rast(template_stk[[1]])
  terra::values(r) <- vec
  names(r) <- name
  r
}

summary_to_stack <- function(template_stk, summary_list, prefix) {
  # Build a 4-layer SpatRaster: prefix_median / prefix_lo / prefix_hi / prefix_sd
  rs <- list(
    raster_from_vec(template_stk, summary_list$median, paste0(prefix, "_median")),
    raster_from_vec(template_stk, summary_list$lo,     paste0(prefix, "_lo")),
    raster_from_vec(template_stk, summary_list$hi,     paste0(prefix, "_hi")),
    raster_from_vec(template_stk, summary_list$sd,     paste0(prefix, "_sd"))
  )
  do.call(c, rs)
}

message("[5/7] Predicting current FOI summary (median / lo / hi / sd) ...")
sum_curr   <- predict_ensemble_summary(covlist_current, models_pilot)
foi_current_stack <- summary_to_stack(covlist_current, sum_curr, "FOI_current")
foi_current       <- foi_current_stack[["FOI_current_median"]]   # primary FOI map

foi_future_stacks <- list()
foi_future_list   <- list()                                       # backwards-compat for viz
for (ssp in SSPS_PILOT) {
  message(sprintf("[5/7] Predicting future FOI SSP%s summary ...", ssp))
  sum_fut <- predict_ensemble_summary(future_stacks[[ssp]], models_pilot)
  stk_fut <- summary_to_stack(covlist_current, sum_fut, sprintf("FOI_SSP%s", ssp))
  foi_future_stacks[[ssp]] <- stk_fut
  foi_future_list[[ssp]]   <- stk_fut[[sprintf("FOI_SSP%s_median", ssp)]]
}

# ---- 6) Save outputs --------------------------------------------------------
# Defensive: if user runs this section interactively without sourcing from top,
# rebuild `tag` from CONFIG so writeRaster always has a valid filename suffix.
if (!exists("tag")) {
  if (!exists("N_MODELS_PILOT")) N_MODELS_PILOT <- 1
  if (!exists("SSPS_PILOT"))     SSPS_PILOT     <- c("245")
  tag <- sprintf("n%d_%s", N_MODELS_PILOT,
                 paste0("SSP", SSPS_PILOT, collapse = "_"))
  message("Recreated tag (was missing): ", tag)
}

message("[6/7] Saving rasters ...")

# Save the full 4-layer summary stack (median/lo/hi/sd) AND a single-layer
# median raster for downstream use (mirrors old foi_layer1 / foi_raster).
terra::writeRaster(
  foi_current_stack,
  file.path(out_dir, sprintf("FOI_current_summary_%s.tif", tag)),
  overwrite = TRUE
)
terra::writeRaster(
  foi_current,
  file.path(out_dir, sprintf("FOI_current_median_%s.tif", tag)),
  overwrite = TRUE
)

for (ssp in SSPS_PILOT) {
  terra::writeRaster(
    foi_future_stacks[[ssp]],
    file.path(out_dir, sprintf("FOI_SSP%s_summary_%s.tif", ssp, tag)),
    overwrite = TRUE
  )
  terra::writeRaster(
    foi_future_list[[ssp]],
    file.path(out_dir, sprintf("FOI_SSP%s_median_%s.tif", ssp, tag)),
    overwrite = TRUE
  )
}

# ---- 7) Visualization -------------------------------------------------------
message("[7/7] Saving comparison map ...")
n_ssp   <- length(SSPS_PILOT)
n_panel <- 1 + 2 * n_ssp + 1   # current + (future + diff per SSP) + hist

png(file.path(out_dir, sprintf("pilot_FOI_comparison_%s.png", tag)),
    width = 1500, height = 350 * ceiling(n_panel / 2), res = 130)
par(mfrow = c(ceiling(n_panel / 2), 2), mar = c(3, 3, 3, 4))

plot(foi_current,
     main = sprintf("FOI current (n=%d RF)", N_MODELS_PILOT),
     col = hcl.colors(50, "viridis"))

last_diff <- NULL
for (ssp in SSPS_PILOT) {
  fut <- foi_future_list[[ssp]]
  dif <- fut - foi_current
  last_diff <- dif

  plot(fut,
       main = sprintf("FOI SSP%s (PRCP swapped)", ssp),
       col = hcl.colors(50, "viridis"))
  plot(dif,
       main = sprintf("ΔFOI SSP%s - current", ssp),
       col = hcl.colors(50, "RdBu"))
}

if (!is.null(last_diff)) {
  v <- terra::values(last_diff); v <- v[!is.na(v)]
  hist(v, breaks = 80,
       main = sprintf("ΔFOI SSP%s distribution",
                      utils::tail(SSPS_PILOT, 1)),
       xlab = "ΔFOI")
  abline(v = 0, col = "red", lwd = 2, lty = 2)
}

dev.off()

message("Saved: ",
        file.path(out_dir, sprintf("pilot_FOI_comparison_%s.png", tag)))
