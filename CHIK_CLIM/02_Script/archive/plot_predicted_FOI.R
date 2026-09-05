library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(scales)
library(dplyr)
library(viridisLite)

# ------------------------------------------------------------------------------
# Plot predicted FOI rasters (current / future) using the old-style log10 scale.
#
# Tidyterra-free: rasters are converted to data.frames and drawn via geom_raster.
# Mirrors the old workflow:
#   df <- as.data.frame(foi_layer, xy = TRUE)
#   geom_tile(filter(df, val == 0), fill = "beige") +
#   geom_raster(filter(df, val >  0), aes(fill = val))
#
# Inputs (auto-loaded from disk if not in environment):
#   01_Data/future_pilot/predict/FOI_current_median_<tag>.tif
#   01_Data/future_pilot/predict/FOI_SSP<ssp>_median_<tag>.tif
# ------------------------------------------------------------------------------

# ---- Paths -------------------------------------------------------------------
cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
project_root <- if (basename(cwd) %in% c("01_Data", "02_Script")) dirname(cwd) else cwd

predict_dir <- file.path(project_root, "01_Data", "future_pilot", "predict")
plot_dir    <- file.path(project_root, "01_Data", "future_pilot", "predict", "maps")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Pull tag from CONFIG (or use the same default as predict script) -------
if (!exists("tag")) {
  if (!exists("N_MODELS_PILOT")) N_MODELS_PILOT <- 1
  if (!exists("SSPS_PILOT"))     SSPS_PILOT     <- c("245")
  tag <- sprintf("n%d_%s", N_MODELS_PILOT,
                 paste0("SSP", SSPS_PILOT, collapse = "_"))
}

# ---- World boundaries (country-level) ---------------------------------------
world_boundaries <- ne_countries(scale = "medium", returnclass = "sf")

# ---- Memory-safe raster -> plotting data.frame -----------------------------
# A native 0.05deg global raster has ~26M cells; keeping the full grid in
# a data.frame plus dplyr filter copies overruns macOS R's 16Gb cap.
# We drop NA cells at conversion time and (if still too large) aggregate
# by an integer factor before plotting.  This is purely for display; the
# rasters on disk remain at native resolution.
MAX_PLOT_CELLS <- 4e6

raster_to_plot_df <- function(r, fun = "max") {
  if (terra::ncell(r) > MAX_PLOT_CELLS) {
    fact <- ceiling(sqrt(terra::ncell(r) / MAX_PLOT_CELLS))
    r <- terra::aggregate(r, fact = fact, fun = fun, na.rm = TRUE)
    message(sprintf("[plot] aggregated raster by factor %d -> %d cells",
                    fact, terra::ncell(r)))
  }
  df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  val_col <- setdiff(names(df), c("x", "y"))[1]
  names(df)[names(df) == val_col] <- "val"
  df
}

# ---- Default range mask path (c_mask) used to grey-out non-range countries --
# Priority order:
#   1) The combined (cell-level x country-level) mask produced by
#      02_Script/archive/post_process_foi_pilot.R, if it exists.
#   2) Otherwise the cell-level CHIK binmap on disk, which matches the old
#      chik_binary used in the original masking step.
DEFAULT_CMASK_PATH <- {
  combined <- file.path(project_root, "01_Data", "future_pilot",
                        "predict", "CHIK_cmask_combined.tif")
  fallback <- file.path(project_root, "01_Data", "CHIK_binmap_2024_04_24.tif")
  if (file.exists(combined)) combined else fallback
}

# ---- FOI plotting helper (log10 scale, 0 -> beige, outside-range -> grey) ---
plot_foi_map <- function(foi_raster,
                         title       = "FOI",
                         c_mask_path = DEFAULT_CMASK_PATH,
                         xlim        = c(-180, 180),
                         ylim        = c(-50, 75),
                         palette     = rev(viridisLite::rocket(5))) {
  # max-aggregate preserves FOI hot spots when downsampling for display.
  df_all  <- raster_to_plot_df(foi_raster, fun = "max")
  df_zero <- df_all[df_all$val == 0, , drop = FALSE]
  df_pos  <- df_all[df_all$val >  0, , drop = FALSE]
  rm(df_all); gc()

  # Optional range mask: cells where mask <= 0 are drawn in light grey
  # to mark "outside chikungunya geographic range".  Aggregating with
  # fun = "max" on a 0/1 mask means a block is greyed only when *every*
  # original cell inside it was outside the CHIK range.
  df_mask_zero <- NULL
  if (!is.null(c_mask_path) && file.exists(c_mask_path)) {
    cmask_r      <- terra::rast(c_mask_path)
    cmask_r      <- terra::resample(cmask_r, foi_raster, method = "near")
    df_mask      <- raster_to_plot_df(cmask_r, fun = "max")
    df_mask_zero <- df_mask[df_mask$val <= 0, , drop = FALSE]
    rm(df_mask, cmask_r); gc()
  }

  log_min    <- log10(0.001)
  log_max    <- log10(0.1)
  log_breaks <- seq(log_min, log_max, length.out = 6)
  breaks     <- 10 ^ log_breaks
  ext_labels <- scales::label_number(accuracy = 0.001)(breaks)

  p <- ggplot()

  # Layer order matches the old code:
  #   1) FOI = 0 -> beige
  #   2) c_mask = 0 -> lightgrey (overlays beige in outside-range cells)
  #   3) FOI > 0 -> log10 colour scale
  #   4) world boundaries on top
  p <- p + geom_tile(data = df_zero,
                     aes(x = x, y = y), fill = "beige", alpha = 1)
  if (!is.null(df_mask_zero)) {
    p <- p + geom_tile(data = df_mask_zero,
                       aes(x = x, y = y), fill = "lightgrey", alpha = 1)
  }

  p <- p +
    geom_raster(data = df_pos,
                aes(x = x, y = y, fill = val)) +
    geom_sf(data = world_boundaries,
            fill = NA, color = "black", linewidth = 0.3) +
    scale_fill_gradientn(
      name     = "FOI",
      colors   = palette,
      trans    = "log10",
      breaks   = breaks,
      labels   = ext_labels,
      limits   = c(min(breaks), max(breaks)),
      na.value = "transparent",
      oob      = scales::squish
    ) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      legend.position  = "right"
    )

  p
}

# ---- ΔFOI plotting helper (diverging linear, quantile-clipped) -------------
# Uses |delta|'s 95th percentile as symmetric limit so the central 90% of
# variation is amplified visually. Larger absolute changes are squished to
# the colour endpoints (oob = squish).
plot_diff_map <- function(diff_raster,
                          title       = "ΔFOI",
                          c_mask_path = DEFAULT_CMASK_PATH,
                          xlim        = c(-180, 180),
                          ylim        = c(-50, 75),
                          quant_clip  = 0.95) {
  # mean-aggregate is appropriate for a continuous (signed) signal.
  df <- raster_to_plot_df(diff_raster, fun = "mean")

  # Optional range mask
  df_mask_zero <- NULL
  if (!is.null(c_mask_path) && file.exists(c_mask_path)) {
    cmask_r      <- terra::rast(c_mask_path)
    cmask_r      <- terra::resample(cmask_r, diff_raster, method = "near")
    df_mask      <- raster_to_plot_df(cmask_r, fun = "max")
    df_mask_zero <- df_mask[df_mask$val <= 0, , drop = FALSE]
    rm(df_mask, cmask_r); gc()
  }

  q <- stats::quantile(abs(df$val), probs = quant_clip, na.rm = TRUE)
  if (!is.finite(q) || q == 0) q <- max(abs(df$val), na.rm = TRUE)
  lims <- c(-q, q)

  p <- ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = val))
  if (!is.null(df_mask_zero)) {
    p <- p + geom_tile(data = df_mask_zero,
                       aes(x = x, y = y), fill = "lightgrey")
  }

  p <- p +
    geom_sf(data = world_boundaries, fill = NA, color = "black",
            linewidth = 0.3) +
    scale_fill_gradient2(
      name     = "ΔFOI",
      low      = "#2166AC",
      mid      = "white",
      high     = "#B2182B",
      midpoint = 0,
      limits   = lims,
      oob      = scales::squish,
      na.value = "transparent"
    ) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title    = title,
         subtitle = sprintf("Color clipped at ±%.1e (|ΔFOI| %.0fth pct)",
                            q, quant_clip * 100),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank()
    )

  p
}

# ---- Load + plot rasters ----------------------------------------------------
load_or_skip <- function(path) {
  if (!file.exists(path)) {
    message("Missing (skipped): ", path); return(NULL)
  }
  terra::rast(path)
}

# Prefer the masked raster produced by post_process_foi_pilot.R if present;
# otherwise fall back to the raw predicted raster.  The masked version has
# FOI = 0 outside the CHIK range mask, which then renders as beige and is
# overlaid by the c_mask grey -> matches the old foi_layer1 / foi_layer2 look.
pick_raster_path <- function(scen_name) {
  masked_p <- file.path(predict_dir,
                        sprintf("FOI_%s_median_masked_%s.tif", scen_name, tag))
  raw_p    <- file.path(predict_dir,
                        sprintf("FOI_%s_median_%s.tif", scen_name, tag))
  if (file.exists(masked_p)) masked_p else raw_p
}

foi_current_path <- pick_raster_path("current")
foi_current      <- load_or_skip(foi_current_path)
if (!is.null(foi_current)) {
  message("Using current raster: ", foi_current_path)
}

if (!is.null(foi_current)) {
  p_curr <- plot_foi_map(foi_current,
                         title = sprintf("FOI current (median, %s)", tag))
  ggsave(file.path(plot_dir, sprintf("map_FOI_current_%s.png", tag)),
         p_curr, width = 10, height = 5, dpi = 300)
  message("Saved: ", file.path(plot_dir, sprintf("map_FOI_current_%s.png", tag)))
}

ssps_to_plot <- if (exists("SSPS_PILOT")) SSPS_PILOT else c("245")

for (ssp in ssps_to_plot) {
  fut_path <- pick_raster_path(sprintf("SSP%s", ssp))
  fut_r <- load_or_skip(fut_path)
  if (is.null(fut_r)) next
  message("Using SSP", ssp, " raster: ", fut_path)

  p_fut <- plot_foi_map(
    fut_r,
    title = sprintf("FOI future SSP%s (median, %s)", ssp, tag)
  )
  ggsave(file.path(plot_dir,
                   sprintf("map_FOI_SSP%s_%s.png", ssp, tag)),
         p_fut, width = 10, height = 5, dpi = 300)
  message("Saved: ",
          file.path(plot_dir, sprintf("map_FOI_SSP%s_%s.png", ssp, tag)))

  if (!is.null(foi_current)) {
    diff_r <- fut_r - foi_current
    p_diff <- plot_diff_map(diff_r,
                            title = sprintf("ΔFOI (SSP%s - current)", ssp))
    ggsave(file.path(plot_dir,
                     sprintf("map_dFOI_SSP%s_%s.png", ssp, tag)),
           p_diff, width = 10, height = 5, dpi = 300)
    message("Saved: ",
            file.path(plot_dir, sprintf("map_dFOI_SSP%s_%s.png", ssp, tag)))
  }
}
