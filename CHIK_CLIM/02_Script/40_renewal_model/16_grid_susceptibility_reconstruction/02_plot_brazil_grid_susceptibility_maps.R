# Fine-scale maps from the compact outputs of 01_reconstruct_brazil_grid_susceptibility.R.
#
# This is a visualisation-only post-processing step.  It maps posterior/ensemble
# medians and does not alter FOI, susceptibility, demographic accounting, or
# the annual SHAPE posterior draws.

required_map_packages <- c("here", "dplyr", "readr", "tibble", "ggplot2", "sf", "scales")
missing_map_packages <- required_map_packages[
  !vapply(required_map_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_map_packages)) stop("Missing package(s): ", paste(missing_map_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(sf)
  library(scales)
})

map_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

map_paths <- function(root = map_root()) {
  output_tag <- Sys.getenv("GRID_SUSC_OUTPUT_TAG", "")
  table_root <- file.path(root, "03_Output", "tables", "national_grid_susceptibility")
  figure_root <- file.path(root, "03_Output", "figures", "national_grid_susceptibility")
  if (nzchar(output_tag)) {
    table_root <- file.path(table_root, output_tag)
    figure_root <- file.path(figure_root, output_tag)
  }
  list(
    table = table_root,
    grid = file.path(table_root, "grid_selected_years"),
    figure = figure_root,
    boundary = file.path(root, "01_Data", "ibge_muni_polygons.rds")
  )
}

read_grid_map_data <- function(grid_directory) {
  grid_files <- sort(list.files(
    grid_directory,
    pattern = "^brazil_chik_grid_susceptibility_[A-Z]{2}_selected_years\\.rds$",
    full.names = TRUE
  ))
  if (!length(grid_files)) stop("No completed state grid summaries found in: ", grid_directory)
  bind_rows(lapply(grid_files, readRDS))
}

get_grid_spacing <- function(x) {
  spacing <- median(diff(sort(unique(x))))
  if (!is.finite(spacing) || spacing <= 0) stop("Could not determine regular 5-km grid spacing.")
  spacing
}

read_uf_outline <- function(path) {
  if (!file.exists(path)) stop("IBGE municipality boundary layer is unavailable: ", path)
  municipalities <- readRDS(path)
  if (is.na(st_crs(municipalities))) stop("IBGE municipality boundary layer has no CRS.")
  st_transform(st_union(st_make_valid(municipalities)), 4326)
}

make_timevarying_susceptibility_map <- function(grid, outline, paths, years) {
  values <- bind_rows(lapply(years, function(year) {
    column <- paste0("susceptible_fraction_", year, "_median")
    if (!column %in% names(grid)) stop("Missing susceptibility field: ", column)
    transmute(grid, x, y, year = factor(year, levels = years), susceptible_fraction = .data[[column]])
  }))
  dx <- get_grid_spacing(values$x)
  dy <- get_grid_spacing(values$y)
  plot <- ggplot(values, aes(x, y, fill = susceptible_fraction)) +
    geom_tile(width = dx, height = dy) +
    geom_sf(data = outline, inherit.aes = FALSE, fill = NA, colour = "grey20", linewidth = .08) +
    coord_sf(expand = FALSE) +
    facet_wrap(~year, ncol = 3) +
    scale_fill_viridis_c(
      option = "C", limits = c(0, 1), oob = scales::squish,
      labels = scales::label_percent(accuracy = 1), name = "Susceptible\nfraction (S/N)"
    ) +
    labs(
      title = "Brazil chikungunya susceptibility at 5-km resolution",
      subtitle = "Posterior/ensemble median; all panels share a 0-100% scale",
      x = NULL, y = NULL
    ) +
    theme_void(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  output <- file.path(paths$figure, "brazil_chik_grid_susceptibility_timevarying_selected_years.png")
  ggsave(output, plot, width = 330, height = 225, units = "mm", dpi = 300, bg = "white")
  output
}

make_cumulative_infection_map <- function(grid, outline, paths) {
  column <- "cumulative_infected_fraction_2025_median"
  if (!column %in% names(grid)) stop("Missing cumulative-infection field: ", column)
  dx <- get_grid_spacing(grid$x)
  dy <- get_grid_spacing(grid$y)
  upper <- max(grid[[column]], na.rm = TRUE)
  if (!is.finite(upper) || upper <= 0) stop("Invalid cumulative infected fraction.")
  plot <- ggplot(grid, aes(x, y, fill = .data[[column]])) +
    geom_tile(width = dx, height = dy) +
    geom_sf(data = outline, inherit.aes = FALSE, fill = NA, colour = "grey20", linewidth = .08) +
    coord_sf(expand = FALSE) +
    scale_fill_viridis_c(
      option = "B", limits = c(0, upper), oob = scales::squish,
      labels = scales::label_percent(accuracy = 1), name = "Cumulative\ninfected fraction"
    ) +
    labs(
      title = "Cumulative chikungunya infection-derived immunity by 2025",
      subtitle = "5-km posterior/ensemble median; cumulative infections divided by the 2015 baseline cell population",
      x = NULL, y = NULL
    ) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
  output <- file.path(paths$figure, "brazil_chik_grid_cumulative_infected_fraction_2025.png")
  ggsave(output, plot, width = 190, height = 225, units = "mm", dpi = 300, bg = "white")
  output
}

make_susceptibility_change_map <- function(grid, outline, paths) {
  start_column <- "susceptible_fraction_2015_median"
  end_column <- "susceptible_fraction_2025_median"
  if (!all(c(start_column, end_column) %in% names(grid))) stop("Missing 2015 or 2025 susceptibility fields.")
  grid <- mutate(grid, susceptibility_change = .data[[end_column]] - .data[[start_column]])
  limit <- max(abs(grid$susceptibility_change), na.rm = TRUE)
  if (!is.finite(limit) || limit <= 0) stop("Invalid susceptibility change.")
  dx <- get_grid_spacing(grid$x)
  dy <- get_grid_spacing(grid$y)
  plot <- ggplot(grid, aes(x, y, fill = susceptibility_change)) +
    geom_tile(width = dx, height = dy) +
    geom_sf(data = outline, inherit.aes = FALSE, fill = NA, colour = "grey20", linewidth = .08) +
    coord_sf(expand = FALSE) +
    scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
      limits = c(-limit, limit), labels = scales::label_percent(accuracy = 1),
      name = expression(Delta * " susceptible fraction")
    ) +
    labs(
      title = "Change in 5-km chikungunya susceptibility, 2015–2025",
      subtitle = "Posterior/ensemble median; red denotes greater depletion",
      x = NULL, y = NULL
    ) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
  output <- file.path(paths$figure, "brazil_chik_grid_susceptibility_change_2015_2025.png")
  ggsave(output, plot, width = 190, height = 225, units = "mm", dpi = 300, bg = "white")
  output
}

run_grid_susceptibility_maps <- function() {
  paths <- map_paths()
  dir.create(paths$figure, recursive = TRUE, showWarnings = FALSE)
  grid <- read_grid_map_data(paths$grid)
  outline <- read_uf_outline(paths$boundary)
  selected_years <- c(2015L, 2017L, 2020L, 2022L, 2024L, 2025L)
  output <- c(
    timevarying_susceptibility = make_timevarying_susceptibility_map(grid, outline, paths, selected_years),
    cumulative_infected_2025 = make_cumulative_infection_map(grid, outline, paths),
    susceptibility_change_2015_2025 = make_susceptibility_change_map(grid, outline, paths)
  )
  write_csv(
    tibble(figure = names(output), path = unname(output), summary = "posterior/ensemble median"),
    file.path(paths$table, "brazil_chik_grid_susceptibility_map_manifest.csv")
  )
  invisible(output)
}

if (sys.nframe() == 0L) run_grid_susceptibility_maps()
