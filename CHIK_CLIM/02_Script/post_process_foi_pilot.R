# =============================================================================
# Post-process pilot FOI predictions: masking + population linkage
#
# Reproduces the relevant logic from the user's old code:
#   1) Apply two CHIK range masks
#        chik_binary : cell-level CHIK binmap (CHIK_binmap_2024_04_24.tif)
#        chik_occ    : country-level CHIK binmap (CHIK_binmap_by_country2024_04_24.tif)
#      A cell is "in range" iff both masks are > 0.
#   2) Optionally zero-out USA and China at the country level (rnaturalearth).
#   3) Merge cell-wise FOI with age-stratified female + male population
#      (f_pop / m_pop) to build foi_comb_all_df, suitable for downstream
#      burden / exposure analyses.
#
# Inputs
#   01_Data/CHIK_binmap_2024_04_24.tif
#   01_Data/CHIK_binmap_by_country2024_04_24.tif
#   01_Data/f_pop.RData                # data.frame: x, y, "f 0" ... "f 18"
#   01_Data/m_pop.RData                # data.frame: x, y, "m 0" ... "m 18"
#   01_Data/future_pilot/predict/FOI_current_median_<tag>.tif
#   01_Data/future_pilot/predict/FOI_SSP<ssp>_median_<tag>.tif
#
# Outputs (under 01_Data/future_pilot/predict/)
#   CHIK_cmask_combined.tif                          # 0/1 combined range mask
#   FOI_<scenario>_median_masked_<tag>.tif           # masked FOI raster
#   foi_comb_all_df_<scenario>_<tag>.RData           # FOI x population table
#
# Run AFTER 02_Script/predict_future_PRCP_pilot.R has produced the median
# rasters for the active CONFIG (N_MODELS_PILOT, SSPS_PILOT).
# =============================================================================

# ---- Libraries --------------------------------------------------------------
suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
})

# ---- Bump R's vector memory cap (macOS default is 16 Gb) -------------------
# f_pop and m_pop expand to ~10 Gb each in memory; loading both at once
# overruns the default cap.  We try to lift the cap to 64 Gb (virtual; if the
# machine has less RAM, the OS will swap and the script just runs slower).
#
# IMPORTANT: if R was started with a hard limit (e.g. via ~/.Renviron's
# R_MAX_VSIZE or `--max-vsize=...` on the command line), mid-session calls to
# mem.maxVSize() will be silently ignored and the cap stays at the startup
# value.  In that case the only fix is to set R_MAX_VSIZE in ~/.Renviron and
# restart R.  The diagnostic below prints the cap before/after the bump so
# you can tell whether it actually took effect.
if (.Platform$OS.type != "windows") {
  before_mb <- tryCatch(utils::mem.maxVSize(), error = function(e) NA_real_)
  tryCatch(
    utils::mem.maxVSize(vsize = 64000),   # 64 Gb (units: Mb)
    error = function(e) {
      Sys.setenv(R_MAX_VSIZE = "64Gb")
      message("[mem] mem.maxVSize() errored; set R_MAX_VSIZE=64Gb instead.")
    }
  )
  after_mb <- tryCatch(utils::mem.maxVSize(), error = function(e) NA_real_)
  message(sprintf("[mem] vsize cap: before = %.0f Mb, after = %.0f Mb",
                  before_mb, after_mb))
  if (!is.na(after_mb) && after_mb < 32000) {
    message("[mem] WARNING: cap is still < 32 Gb after the bump.\n",
            "      The session likely has a hard startup limit that ",
            "mid-session calls cannot raise.\n",
            "      Either set DO_POP_MERGE <- FALSE below to skip the ",
            "memory-heavy step,\n",
            "      or add the line `R_MAX_VSIZE=64Gb` to ~/.Renviron and ",
            "restart R.")
  }
}

# ---- Paths ------------------------------------------------------------------
cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
project_root <- if (basename(cwd) %in% c("01_Data", "02_Script")) dirname(cwd) else cwd

data_dir    <- file.path(project_root, "01_Data")
predict_dir <- file.path(data_dir, "future_pilot", "predict")
dir.create(predict_dir, recursive = TRUE, showWarnings = FALSE)

# ---- CONFIG -----------------------------------------------------------------
# Match the predict script's tag.  If those globals are not set (e.g. running
# this script standalone) we default to the same n=1 / SSP245 pilot.
if (!exists("N_MODELS_PILOT")) N_MODELS_PILOT <- 1
if (!exists("SSPS_PILOT"))     SSPS_PILOT     <- c("245")
if (!exists("tag")) {
  tag <- sprintf("n%d_%s", N_MODELS_PILOT,
                 paste0("SSP", SSPS_PILOT, collapse = "_"))
}

# Toggle expensive steps (population merge needs ~1.5 GB RAM).
APPLY_COUNTRY_EXCLUSION <- TRUE                          # USA, China -> 0
EXCLUDED_COUNTRIES      <- c("United States of America", "China")
DO_POP_MERGE            <- FALSE                         # merge f_pop + m_pop
# NOTE: f_pop + m_pop together can need >16 Gb in R's vector heap.  If you hit
# `Error: vector memory limit ... reached`, leave DO_POP_MERGE = FALSE for now.
# The masking and FOI raster outputs still come out fine; the population merge
# (used downstream for burden) can be re-enabled once R's vsize cap is raised
# (see the [mem] diagnostic above), or in a separate, smaller workspace.

message(sprintf("[post-process] tag = %s", tag))

# ---- 1. Load FOI rasters for the active tag --------------------------------
foi_paths <- list(
  current = file.path(predict_dir, sprintf("FOI_current_median_%s.tif", tag))
)
for (ssp in SSPS_PILOT) {
  foi_paths[[paste0("SSP", ssp)]] <-
    file.path(predict_dir, sprintf("FOI_SSP%s_median_%s.tif", ssp, tag))
}
missing <- foi_paths[!vapply(foi_paths, file.exists, logical(1))]
if (length(missing) > 0) {
  stop("Missing FOI median rasters:\n  ",
       paste(unlist(missing), collapse = "\n  "),
       "\nRun 02_Script/predict_future_PRCP_pilot.R first.")
}
foi_rasters <- lapply(foi_paths, terra::rast)

# Reference grid = current FOI raster (covariate grid).
ref_grid <- foi_rasters$current

# ---- 2. Load and align mask rasters ----------------------------------------
chik_bin_path <- file.path(data_dir, "CHIK_binmap_2024_04_24.tif")
chik_occ_path <- file.path(data_dir, "CHIK_binmap_by_country2024_04_24.tif")
stopifnot(file.exists(chik_bin_path), file.exists(chik_occ_path))

chik_bin <- terra::rast(chik_bin_path)
chik_occ <- terra::rast(chik_occ_path)

# Resample to the FOI grid using nearest neighbour (binary masks).
chik_bin_aligned <- terra::resample(chik_bin, ref_grid, method = "near")
chik_occ_aligned <- terra::resample(chik_occ, ref_grid, method = "near")

# Combined mask: 1 where both binmaps say CHIK is plausibly present.
cmask_combined <- (chik_bin_aligned > 0) & (chik_occ_aligned > 0)
cmask_combined <- terra::ifel(is.na(cmask_combined), 0, cmask_combined)
names(cmask_combined) <- "cmask"

# ---- 2b. Optional country-level exclusion (USA, China) ---------------------
if (APPLY_COUNTRY_EXCLUSION) {
  countries <- ne_countries(scale = "medium", returnclass = "sf")
  excl_sf <- countries[countries$name_long %in% EXCLUDED_COUNTRIES |
                         countries$sovereignt %in% EXCLUDED_COUNTRIES, ]
  if (nrow(excl_sf) > 0) {
    excl_vec <- terra::vect(excl_sf)
    excl_rast <- terra::rasterize(excl_vec, ref_grid, field = 1, background = 0)
    cmask_combined <- cmask_combined * (excl_rast == 0)
    message(sprintf("[mask] zeroed-out countries: %s",
                    paste(EXCLUDED_COUNTRIES, collapse = ", ")))
  }
}

terra::writeRaster(cmask_combined,
                   file.path(predict_dir, "CHIK_cmask_combined.tif"),
                   overwrite = TRUE)
message("[mask] saved CHIK_cmask_combined.tif")

# ---- 3. Apply mask to each FOI raster --------------------------------------
masked_paths <- character(0)
for (scen in names(foi_rasters)) {
  foi_r <- foi_rasters[[scen]]
  foi_masked <- foi_r * cmask_combined
  out_path <- file.path(predict_dir,
                        sprintf("FOI_%s_median_masked_%s.tif", scen, tag))
  terra::writeRaster(foi_masked, out_path, overwrite = TRUE)
  masked_paths[scen] <- out_path
  rng <- terra::global(foi_masked, c("min", "max"), na.rm = TRUE)
  message(sprintf("[mask] %s -> [%g, %g]   %s",
                  scen, rng[1, 1], rng[1, 2], basename(out_path)))
}

# ---- 4. Population merge (foi_comb_all_df) ---------------------------------
# Mirrors the old code:
#   all_pop <- f_pop[,3:21] + m_pop[,3:21]; cbind(f_pop[,1:2], all_pop)
#   colnames(all_pop)[3:20] <- 1:18
#   foi_comb_all <- merge(foi_mat, all_pop, by = c("x","y"))
#   foi_comb_all$mask_val_bin <- extract(chik_binary, ...)
#   foi_comb_all$mask_val_occ <- extract(chik_occ, ...)
#   foi_comb_all <- foi_comb_all[!is.na(foi_comb_all$"1"), ]
if (DO_POP_MERGE) {
  # Memory-frugal pattern: never hold both raw f_pop and m_pop at once.
  # Step 1: load f_pop -> extract xy + age matrix -> drop f_pop
  # Step 2: load m_pop -> add age columns into f's matrix in place -> drop m_pop
  # Step 3: drop NA-pop rows early (the cap on size for everything downstream)

  load_one_object <- function(path) {
    e <- new.env()
    load(path, envir = e)
    obj <- get(ls(e)[1], envir = e)
    rm(e); gc()
    obj
  }

  message("[pop] loading f_pop ...")
  f_pop <- load_one_object(file.path(data_dir, "f_pop.RData"))
  stopifnot(all(c("x", "y") %in% names(f_pop)))
  pop_xy_df <- f_pop[, c("x", "y")]
  age_idx_f <- 3:ncol(f_pop)
  age_mat <- as.matrix(f_pop[, age_idx_f])
  storage.mode(age_mat) <- "double"
  rm(f_pop); gc()

  message("[pop] loading m_pop and summing in place ...")
  m_pop <- load_one_object(file.path(data_dir, "m_pop.RData"))
  stopifnot(all(c("x", "y") %in% names(m_pop)))
  age_idx_m <- 3:ncol(m_pop)
  if (length(age_idx_m) != length(age_idx_f)) {
    stop("f_pop and m_pop have different number of age columns.")
  }
  for (j in seq_along(age_idx_m)) {
    age_mat[, j] <- age_mat[, j] +
      as.numeric(m_pop[[age_idx_m[j]]])
  }
  rm(m_pop); gc()

  # Match old column naming: "1","2",...  (use character so $"1" still works).
  colnames(age_mat) <- as.character(seq_len(ncol(age_mat)))

  # Filter early: drop rows whose first age column is NA (no population).
  keep <- !is.na(age_mat[, 1])
  age_mat   <- age_mat[keep, , drop = FALSE]
  pop_xy_df <- pop_xy_df[keep, , drop = FALSE]
  rm(keep); gc()
  message(sprintf("[pop] kept %d population cells (NA-filtered).",
                  nrow(age_mat)))

  pop_xy <- as.matrix(pop_xy_df)

  # Mask values from the *raw* binmaps at pop xy (matches old code semantics).
  mask_val_bin <- terra::extract(chik_bin, pop_xy)[, 1]
  mask_val_occ <- terra::extract(chik_occ, pop_xy)[, 1]

  for (scen in names(foi_rasters)) {
    foi_masked_r <- terra::rast(masked_paths[[scen]])
    foi_at_pop <- terra::extract(foi_masked_r, pop_xy)[, 1]
    rm(foi_masked_r); gc()

    foi_comb_all_df <- data.frame(
      x = pop_xy_df$x,
      y = pop_xy_df$y,
      age_mat,
      foi_mid      = foi_at_pop,
      mask_val_bin = mask_val_bin,
      mask_val_occ = mask_val_occ,
      check.names  = FALSE
    )

    out_rdata <- file.path(predict_dir,
                           sprintf("foi_comb_all_df_%s_%s.RData", scen, tag))
    save(foi_comb_all_df, file = out_rdata)
    message(sprintf("[pop] %s -> %d rows  %s",
                    scen, nrow(foi_comb_all_df), basename(out_rdata)))
    rm(foi_comb_all_df, foi_at_pop); gc()
  }

  rm(age_mat, pop_xy, pop_xy_df, mask_val_bin, mask_val_occ); gc()
} else {
  message("[pop] DO_POP_MERGE = FALSE; skipped population merge.")
}

message("[post-process] done.")
