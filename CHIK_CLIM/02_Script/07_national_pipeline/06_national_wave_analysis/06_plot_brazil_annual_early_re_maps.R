# Annual maps of statewide early Re from the frozen major-wave census.
#
# A state can have multiple major waves beginning in one year.  For a single
# state-year fill, the map displays the maximum usable Re6 posterior median;
# if no usable Re6 exists it displays the earliest major wave as unavailable.
# An n>1 label identifies every state-year with multiple major waves.  This is
# a display rule only: all waves remain separate rows in the master table.

required_map_packages <- c("here", "sf", "dplyr", "readr", "ggplot2", "scales", "patchwork", "cowplot", "png")
missing_map_packages <- required_map_packages[
  !vapply(required_map_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_map_packages)) stop("Missing package(s): ", paste(missing_map_packages, collapse = ", "))
suppressPackageStartupMessages({ library(sf); library(dplyr); library(readr); library(ggplot2); library(scales) })

map_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}

build_state_geometry <- function(root) {
  municipalities <- readRDS(file.path(root, "01_Data", "ibge_muni_polygons.rds"))
  lookup <- readRDS(file.path(root, "01_Data", "ibge_muni_name_lookup.rds")) |>
    transmute(muni6 = as.character(muni6), state = uf)
  municipal_geometry <- municipalities |>
    mutate(muni6 = as.character(muni6)) |>
    inner_join(lookup, by = "muni6")
  # This tolerance is solely a cartographic simplification for a national map.
  # It does not alter state membership or any epidemic quantities.
  municipal_geometry <- st_simplify(municipal_geometry, dTolerance = 0.02, preserveTopology = TRUE)
  # Do not union high-resolution municipal coastlines: that operation is very
  # slow and is unnecessary for a state-filled display.  Municipal polygons
  # are drawn without strokes, so their state-level fill reads as one region.
  centre_coordinates <- st_coordinates(st_centroid(municipal_geometry))
  state_points <- municipal_geometry |>
    st_drop_geometry() |>
    mutate(x = centre_coordinates[, 1], y = centre_coordinates[, 2]) |>
    group_by(state) |>
    summarise(x = mean(x), y = mean(y), .groups = "drop") |>
    st_as_sf(coords = c("x", "y"), crs = st_crs(municipal_geometry))
  list(municipal = municipal_geometry, state_points = state_points)
}

select_state_year_event <- function(master, years) {
  events <- master |>
    mutate(calendar_year = as.integer(calendar_year), usable = re6_analysis_usable & !is.na(Re_early_6_median)) |>
    dplyr::filter(calendar_year %in% years)
  count_table <- events |> count(calendar_year, state, name = "n_major_waves")
  representative <- bind_rows(lapply(split(events, interaction(events$calendar_year, events$state, drop = TRUE)), function(x) {
    usable <- dplyr::filter(x, usable)
    if (nrow(usable)) slice_max(usable, Re_early_6_median, n = 1, with_ties = FALSE) else slice_min(x, onset_week, n = 1, with_ties = FALSE)
  })) |>
    left_join(count_table, by = c("calendar_year", "state")) |>
    mutate(event_status = if_else(usable, if_else(is_recurrence, "recurrent_usable", "first_usable"), "major_re_unavailable"))
  all_states <- sort(unique(master$state))
  tidyr::expand_grid(calendar_year = years, state = all_states) |>
    left_join(representative, by = c("calendar_year", "state")) |>
    mutate(
      n_major_waves = coalesce(n_major_waves, 0L),
      event_status = coalesce(event_status, "no_major_wave"),
      map_Re = if_else(event_status %in% c("first_usable", "recurrent_usable"), Re_early_6_median, NA_real_)
    )
}

annual_map <- function(year_value, geometry, state_year, fill_limits, show_legend = FALSE) {
  event_data <- dplyr::filter(state_year, calendar_year == year_value)
  data <- geometry$municipal |> left_join(event_data, by = "state")
  state_points <- geometry$state_points |> left_join(event_data, by = "state")
  unavailable <- dplyr::filter(data, event_status == "major_re_unavailable")
  first <- dplyr::filter(state_points, event_status == "first_usable")
  recurrent <- dplyr::filter(state_points, event_status == "recurrent_usable")
  multiple <- dplyr::filter(state_points, n_major_waves > 1L)

  ggplot() +
    geom_sf(data = data, aes(fill = map_Re), colour = NA) +
    geom_sf(data = unavailable, fill = "white", colour = "grey25", linewidth = .04) +
    geom_sf(data = recurrent, shape = 1, colour = "black", size = 1.8, stroke = .55) +
    geom_sf(data = first, shape = 21, fill = "white", colour = "black", size = 2.8, stroke = .9) +
    geom_sf_text(data = multiple, aes(label = paste0("n=", n_major_waves)), size = 2.1, fontface = "bold") +
    scale_fill_viridis_c(option = "C", limits = fill_limits, trans = "sqrt", oob = squish,
                         name = expression(R[e]^early), na.value = "grey92") +
    coord_sf(datum = NA) +
    labs(title = as.character(year_value)) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", hjust = .5, size = 11),
          legend.position = if (show_legend) "bottom" else "none")
}

run_annual_early_re_maps <- function() {
  root <- map_root()
  table_path <- file.path(root, "03_Output", "tables", "national_wave_analysis", "brazil_chik_wave_analysis_master.csv")
  figure_dir <- file.path(root, "03_Output", "figures", "national_wave_analysis", "annual_early_re_maps")
  if (!file.exists(table_path)) stop("Master wave table is missing; run 04_build_brazil_wave_analysis_master.R first.")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  years <- 2015:2025
  master <- read_csv(table_path, show_col_types = FALSE)
  state_year <- select_state_year_event(master, years)
  geometry <- build_state_geometry(root)
  range_re <- range(state_year$map_Re, na.rm = TRUE)
  maps <- lapply(seq_along(years), function(i) annual_map(
    years[i], geometry = geometry, state_year = state_year, fill_limits = range_re,
    show_legend = i == 1L
  ))
  combined <- patchwork::wrap_plots(maps, ncol = 4) +
    patchwork::plot_annotation(
      title = "Brazilian chikungunya early Re by major-wave onset year",
      subtitle = "Fill: maximum usable six-week Re among waves beginning in the UF-year; light grey: no major wave; dashed outline: major wave but Re6 unavailable; thick outline: first major wave; n>1: multiple major waves began that year.",
      caption = "All panels use the same square-root Re colour scale. This is a descriptive map, not a climate or recurrence model."
    )
  output_path <- file.path(figure_dir, "brazil_chik_annual_early_re_maps_2015_2025.png")
  # A raster device keeps the 11-panel map responsive while retaining every
  # municipal polygon used to create the state fills.
  png(output_path, width = 3600, height = 2700, res = 300, type = "cairo")
  print(combined)
  grDevices::dev.off()
  # Preserve the responsive raster map in a standards-compliant one-page PDF
  # as well. This avoids fragile enormous vector PDFs from municipal geometry.
  pdf_path <- file.path(figure_dir, "brazil_chik_annual_early_re_maps_2015_2025.pdf")
  image <- png::readPNG(output_path)
  grDevices::pdf(pdf_path, width = 12, height = 9, onefile = TRUE)
  grid::grid.newpage(); grid::grid.raster(image, interpolate = TRUE)
  grDevices::dev.off()
  write_csv(st_drop_geometry(geometry$state_points) |> left_join(state_year, by = "state"),
            file.path(root, "03_Output", "tables", "national_wave_analysis", "brazil_chik_annual_early_re_map_display_data.csv"))
  message("[save] ", output_path)
  invisible(output_path)
}

if (sys.nframe() == 0L) run_annual_early_re_maps()
