# Build Brazilian UF-level long-term FOI summaries from the 100-member
# allfoi ensemble. This script is deliberately an aggregation/audit workflow;
# it does not fit a dynamic FOI, susceptibility, reporting, or climate model.

required_packages <- c("here", "sf", "dplyr", "tibble", "readr", "ggplot2", "scales", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

project_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM R-project root.")
}

project_path <- function(...) file.path(project_root(), ...)

uf_lookup <- tibble(
  uf = c(
    "11", "12", "13", "14", "15", "16", "17", "21", "22", "23", "24", "25", "26", "27",
    "28", "29", "31", "32", "33", "35", "41", "42", "43", "50", "51", "52", "53"
  ),
  uf_abbr = c(
    "RO", "AC", "AM", "RR", "PA", "AP", "TO", "MA", "PI", "CE", "RN", "PB", "PE", "AL",
    "SE", "BA", "MG", "ES", "RJ", "SP", "PR", "SC", "RS", "MS", "MT", "GO", "DF"
  ),
  state_name = c(
    "Rondonia", "Acre", "Amazonas", "Roraima", "Para", "Amapa", "Tocantins", "Maranhao", "Piaui", "Ceara", "Rio Grande do Norte", "Paraiba", "Pernambuco", "Alagoas",
    "Sergipe", "Bahia", "Minas Gerais", "Espirito Santo", "Rio de Janeiro", "Sao Paulo", "Parana", "Santa Catarina", "Rio Grande do Sul", "Mato Grosso do Sul", "Mato Grosso", "Goias", "Distrito Federal"
  )
)

make_uf_boundaries <- function() {
  polygon_path <- project_path("01_Data", "ibge_muni_polygons.rds")
  if (!file.exists(polygon_path)) {
    stop("Official IBGE municipal boundary layer is unavailable: ", polygon_path)
  }
  municipalities <- readRDS(polygon_path)
  if (!inherits(municipalities, "sf") || !all(c("muni6", "geometry") %in% names(municipalities))) {
    stop("IBGE municipal boundary file must be an sf object with muni6 and geometry.")
  }
  if (is.na(st_crs(municipalities))) stop("IBGE municipal boundaries have no CRS.")
  municipalities <- st_transform(municipalities, 4326) |>
    mutate(uf = substr(muni6, 1, 2)) |>
    inner_join(uf_lookup, by = "uf")
  boundaries <- municipalities |>
    group_by(uf, uf_abbr, state_name) |>
    summarise(geometry = st_union(geometry), .groups = "drop") |>
    st_make_valid()
  if (nrow(boundaries) != 27L || anyDuplicated(boundaries$uf)) {
    stop("Dissolving official IBGE municipalities did not create exactly 27 UFs.")
  }
  list(boundaries = boundaries, polygon_path = polygon_path)
}

read_allfoi_brazil <- function() {
  allfoi_path <- project_path("01_Data", "allfoi_s1.RData")
  if (!file.exists(allfoi_path)) stop("allfoi input does not exist: ", allfoi_path)
  load(allfoi_path) # Loads object `allfoi` in this function's environment.
  required <- c("x", "y", "country", "iso3", "tot", "foi_mid", "foi_lo", "foi_hi")
  missing <- setdiff(required, names(allfoi))
  if (length(missing) > 0L) stop("allfoi is missing: ", paste(missing, collapse = ", "))
  member_columns <- paste0("foi", seq_len(100L))
  missing_members <- setdiff(member_columns, names(allfoi))
  if (length(missing_members) > 0L) {
    stop("allfoi does not contain all 100 expected ensemble members: ", paste(missing_members, collapse = ", "))
  }
  if (!any(allfoi$iso3 == "BRA", na.rm = TRUE)) {
    stop("No iso3 == 'BRA' rows found; do not substitute an imprecise Brazil filter.")
  }
  keep_columns <- c("x", "y", "country", "iso3", "tot", member_columns, "foi_mid", "foi_lo", "foi_hi")
  brazil <- allfoi[allfoi$iso3 == "BRA", keep_columns, drop = FALSE]
  rm(allfoi)
  gc(verbose = FALSE)

  if (any(!is.finite(brazil$x) | !is.finite(brazil$y))) stop("Brazil grid has non-finite coordinates.")
  if (any(brazil$x < -180 | brazil$x > 180 | brazil$y < -90 | brazil$y > 90)) {
    stop("allfoi x/y are not longitude/latitude degrees; CRS needs explicit revision before joining.")
  }
  if (any(brazil$tot < 0, na.rm = TRUE)) stop("Grid population `tot` contains negative values.")

  list(
    brazil = brazil,
    member_columns = member_columns,
    allfoi_path = allfoi_path,
    n_missing_population = sum(!is.finite(brazil$tot))
  )
}

assign_uf <- function(brazil, boundaries) {
  # `allfoi` coordinates are documented/validated here as WGS84 lon/lat.
  points <- st_as_sf(brazil[, c("x", "y")], coords = c("x", "y"), crs = 4326, remove = FALSE)
  membership <- st_within(points, boundaries, sparse = TRUE)
  n_membership <- lengths(membership)
  direct_index <- vapply(membership, function(x) if (length(x) == 1L) x else NA_integer_, integer(1))
  unmatched <- which(n_membership == 0L)
  ambiguous <- which(n_membership > 1L)
  if (length(ambiguous) > 0L) {
    stop("Some grid centroids join to more than one UF; inspect before aggregation.")
  }

  assignment_method <- rep("within", nrow(brazil))
  nearest_distance_m <- rep(0, nrow(brazil))
  if (length(unmatched) > 0L) {
    # Preserve, quantify, and explicitly label this fallback instead of
    # silently dropping populated coastal/border points. Compute exact
    # fallback distances in small projected chunks: a single all-cell
    # distance call can allocate several GB for coastal grid cells.
    nearest <- st_nearest_feature(points[unmatched, ], boundaries)
    boundaries_m <- st_transform(boundaries, 5880)
    chunk_size <- 5000L
    chunk_starts <- seq.int(1L, length(unmatched), by = chunk_size)
    for (chunk_start in chunk_starts) {
      chunk_position <- seq.int(
        chunk_start,
        min(chunk_start + chunk_size - 1L, length(unmatched))
      )
      chunk_points <- st_transform(points[unmatched[chunk_position], ], 5880)
      nearest_distance_m[unmatched[chunk_position]] <- as.numeric(st_distance(
        chunk_points, boundaries_m[nearest[chunk_position], ], by_element = TRUE
      ))
    }
    direct_index[unmatched] <- nearest
    assignment_method[unmatched] <- "nearest_uf_fallback"
  }
  if (anyNA(direct_index)) stop("Some cells remain without a UF assignment.")

  # Flag a cell as near a UF boundary when one of its eight immediately
  # adjacent cells on the original 5-km grid belongs to another UF. This
  # avoids a costly point-to-complex-boundary distance calculation, is tied
  # directly to the grid resolution, and is strictly an audit flag: it never
  # changes an assignment or drops a cell.
  x_index <- match(brazil$x, sort(unique(brazil$x)))
  y_index <- match(brazil$y, sort(unique(brazil$y)))
  grid_key <- paste(x_index, y_index, sep = "_")
  if (anyDuplicated(grid_key)) stop("Brazil input contains duplicate grid coordinates.")
  near_uf_boundary <- rep(FALSE, nrow(brazil))
  neighbour_offsets <- expand.grid(dx = -1:1, dy = -1:1) |>
    filter(dx != 0 | dy != 0)
  for (offset_index in seq_len(nrow(neighbour_offsets))) {
    neighbour_key <- paste(
      x_index + neighbour_offsets$dx[offset_index],
      y_index + neighbour_offsets$dy[offset_index],
      sep = "_"
    )
    neighbour_index <- match(neighbour_key, grid_key)
    valid_neighbour <- which(!is.na(neighbour_index))
    different_uf <- valid_neighbour[
      direct_index[valid_neighbour] != direct_index[neighbour_index[valid_neighbour]]
    ]
    near_uf_boundary[different_uf] <- TRUE
  }
  assignment <- tibble(
    cell_id = seq_len(nrow(brazil)),
    uf = boundaries$uf[direct_index],
    uf_abbr = boundaries$uf_abbr[direct_index],
    state_name = boundaries$state_name[direct_index],
    assignment_method = assignment_method,
    nearest_distance_m = nearest_distance_m,
    near_uf_boundary = near_uf_boundary
  )
  list(
    assignment = assignment,
    audit = list(
      n_direct = sum(assignment_method == "within"),
      n_raw_unmatched = length(unmatched),
      n_raw_ambiguous = length(ambiguous),
      n_nearest_fallback = sum(assignment_method == "nearest_uf_fallback"),
      n_near_uf_boundary = sum(assignment$near_uf_boundary, na.rm = TRUE),
      maximum_nearest_distance_m = if (length(unmatched)) max(nearest_distance_m[unmatched]) else 0
    )
  )
}

weighted_state_mean <- function(value, weight, group, levels) {
  valid <- is.finite(value) & is.finite(weight) & weight > 0
  numerator <- rowsum(ifelse(valid, weight * value, 0), group = factor(group, levels = levels), reorder = TRUE)
  denominator <- rowsum(ifelse(valid, weight, 0), group = factor(group, levels = levels), reorder = TRUE)
  result <- as.numeric(numerator[, 1] / denominator[, 1])
  result[denominator[, 1] == 0] <- NA_real_
  list(mean = result, valid_population = as.numeric(denominator[, 1]), valid = valid)
}

aggregate_ensemble <- function(brazil, assignment, member_columns) {
  state_levels <- uf_lookup$uf
  group <- assignment$uf
  weight <- brazil$tot
  valid_all_members <- rowSums(
    is.finite(as.matrix(brazil[, member_columns, drop = FALSE]))
  ) == length(member_columns)
  valid_population_all_members <- is.finite(weight) & weight > 0 & valid_all_members
  total_population <- as.numeric(rowsum(
    ifelse(is.finite(weight) & weight > 0, weight, 0),
    group = factor(group, levels = state_levels), reorder = TRUE
  )[, 1])
  valid_population_all <- as.numeric(rowsum(
    ifelse(valid_population_all_members, weight, 0),
    group = factor(group, levels = state_levels), reorder = TRUE
  )[, 1])
  n_cells <- as.integer(table(factor(group, levels = state_levels)))
  # Count prediction availability separately from positive-population cells.
  # The coverage quantities below are the corresponding population-weighted checks.
  n_valid_cells_all <- as.integer(rowsum(
    as.integer(valid_all_members), group = factor(group, levels = state_levels), reorder = TRUE
  )[, 1])
  if (any(n_valid_cells_all > n_cells)) {
    stop("UF valid-cell counts exceed total cells; check UF aggregation order.")
  }

  ensemble_rows <- vector("list", length(member_columns))
  member_coverage <- vector("list", length(member_columns))
  for (member_index in seq_along(member_columns)) {
    values <- brazil[[member_columns[member_index]]]
    arithmetic <- weighted_state_mean(values, weight, group, state_levels)
    probability <- -expm1(-values)
    equivalent_probability <- weighted_state_mean(probability, weight, group, state_levels)
    if (any(equivalent_probability$mean >= 1, na.rm = TRUE)) {
      stop("A UF-level annual infection probability is >= 1; equivalent FOI is undefined.")
    }
    equivalent <- -log1p(-equivalent_probability$mean)
    ensemble_rows[[member_index]] <- tibble(
      uf = state_levels,
      ensemble_id = member_index,
      foi_pwmean = arithmetic$mean,
      foi_equiv = equivalent
    )
    member_coverage[[member_index]] <- tibble(
      uf = state_levels,
      ensemble_id = member_index,
      population_valid_member = arithmetic$valid_population
    )
  }
  ensemble <- bind_rows(ensemble_rows) |>
    left_join(uf_lookup, by = "uf") |>
    select(uf = uf_abbr, state_name, ensemble_id, foi_pwmean, foi_equiv)
  coverage_by_member <- bind_rows(member_coverage)

  state_population <- uf_lookup |>
    mutate(
      n_cells = n_cells,
      n_valid_cells = n_valid_cells_all,
      population_total = total_population,
      population_valid_foi = valid_population_all,
      population_coverage = population_valid_foi / population_total
    ) |>
    left_join(
      coverage_by_member |>
        group_by(uf) |>
        summarise(
          population_coverage_min_member = min(population_valid_member / total_population[match(uf, state_levels)]),
          population_coverage_max_member = max(population_valid_member / total_population[match(uf, state_levels)]),
          .groups = "drop"
        ),
      by = "uf"
    )
  list(ensemble = ensemble, state_population = state_population, coverage_by_member = coverage_by_member)
}

summarise_ensemble <- function(ensemble, state_population) {
  distribution_summary <- function(x) {
    c(
      mid = median(x), mean = mean(x), sd = sd(x), lo = quantile(x, .025),
      q25 = quantile(x, .25), q75 = quantile(x, .75), hi = quantile(x, .975)
    )
  }
  summary <- ensemble |>
    group_by(uf, state_name) |>
    summarise(
      across(
        c(foi_pwmean, foi_equiv),
        list(
          mid = median, mean = mean, sd = sd, lo = ~ quantile(.x, .025), q25 = ~ quantile(.x, .25),
          q75 = ~ quantile(.x, .75), hi = ~ quantile(.x, .975)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    ) |>
    left_join(state_population |> select(uf = uf_abbr, state_name, n_cells, n_valid_cells, population_total, population_valid_foi, population_coverage, population_coverage_min_member, population_coverage_max_member), by = c("uf", "state_name")) |>
    mutate(
      arithmetic_minus_equivalent = foi_pwmean_mid - foi_equiv_mid,
      arithmetic_vs_equivalent_relative_difference = arithmetic_minus_equivalent / foi_equiv_mid
    ) |>
    select(
      uf, state_name, n_cells, n_valid_cells, population_total, population_valid_foi, population_coverage,
      population_coverage_min_member, population_coverage_max_member,
      starts_with("foi_pwmean"), starts_with("foi_equiv"), arithmetic_minus_equivalent,
      arithmetic_vs_equivalent_relative_difference
    )
  if (nrow(summary) != 27L || any(!is.finite(as.matrix(summary |> select(starts_with("foi_")))))) {
    stop("Invalid state-level ensemble summary.")
  }
  summary
}

old_cell_summary_aggregation <- function(brazil, assignment) {
  state_levels <- uf_lookup$uf
  group <- assignment$uf
  weight <- brazil$tot
  old_mid <- weighted_state_mean(brazil$foi_mid, weight, group, state_levels)$mean
  old_lo <- weighted_state_mean(brazil$foi_lo, weight, group, state_levels)$mean
  old_hi <- weighted_state_mean(brazil$foi_hi, weight, group, state_levels)$mean
  tibble(
    uf = uf_lookup$uf_abbr,
    old_foi_mid_wmean = old_mid,
    old_foi_lo_wmean = old_lo,
    old_foi_hi_wmean = old_hi,
    old_interval_width = old_hi - old_lo
  )
}

make_figures <- function(summary, ensemble, boundaries, figure_path) {
  map_data <- boundaries |>
    left_join(summary, by = c("uf_abbr" = "uf", "state_name"))
  map_theme <- theme_void(base_size = 8) + theme(plot.title = element_text(face = "bold"))
  p1 <- ggplot(map_data) +
    geom_sf(aes(fill = foi_pwmean_mid), colour = "white", linewidth = .08) +
    scale_fill_viridis_c(option = "C", labels = scales::label_number(accuracy = .001)) +
    labs(title = "A. UF arithmetic population-weighted FOI", fill = "FOI") + map_theme
  p2 <- ggplot(map_data) +
    geom_sf(aes(fill = foi_equiv_mid), colour = "white", linewidth = .08) +
    scale_fill_viridis_c(option = "C", labels = scales::label_number(accuracy = .001)) +
    labs(title = "B. UF equivalent FOI", fill = "FOI") + map_theme
  p3 <- ggplot(summary, aes(foi_pwmean_mid, foi_equiv_mid, label = uf)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
    geom_point(colour = "#0072B2", size = 1.8) +
    geom_text(size = 2.1, nudge_y = .00015, check_overlap = TRUE) +
    labs(title = "C. Arithmetic versus equivalent UF FOI", x = "Arithmetic population-weighted FOI", y = "Equivalent FOI") +
    theme_classic(base_size = 8.5)
  selected <- unique(c(
    "CE",
    summary$uf[which.max(summary$foi_equiv_mid)],
    summary$uf[which.min(summary$foi_equiv_mid)],
    summary$uf[which.max(summary$population_total)]
  ))
  p4 <- ensemble |>
    filter(uf %in% selected) |>
    ggplot(aes(reorder(uf, foi_equiv), foi_equiv)) +
    geom_boxplot(fill = "#56B4E9", width = .55, outlier.size = .5) +
    labs(title = "D. Selected UF equivalent-FOI ensemble distributions", x = "UF", y = "Equivalent FOI") +
    theme_classic(base_size = 8.5)
  p5 <- ggplot(summary, aes(reorder(uf, population_coverage), population_coverage)) +
    geom_col(fill = "#009E73") +
    geom_hline(yintercept = c(.95, .99), linetype = c(2, 3), colour = c("#D55E00", "#CC79A7")) +
    coord_flip() +
    scale_y_continuous(labels = scales::label_percent(accuracy = .1), limits = c(0, 1)) +
    labs(title = "E. All-member valid-population coverage", x = "UF", y = "Population coverage") +
    theme_classic(base_size = 8.5)
  ggsave(
    figure_path, (p1 | p2) / (p3 | p4) / p5,
    width = 190, height = 255, units = "mm", device = cairo_pdf
  )
}

run_brazil_uf_foi_aggregation <- function() {
  output_ensemble <- project_path("01_Data", "brazil_uf_foi_ensemble.csv")
  output_summary <- project_path("01_Data", "brazil_uf_foi_longterm.csv")
  table_directory <- project_path("03_Output", "tables", "brazil_uf_foi_aggregation")
  figure_directory <- project_path("03_Output", "figures", "brazil_uf_foi_aggregation")
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

  allfoi_input <- read_allfoi_brazil()
  boundary_input <- make_uf_boundaries()
  assignment_input <- assign_uf(allfoi_input$brazil, boundary_input$boundaries)
  aggregation <- aggregate_ensemble(allfoi_input$brazil, assignment_input$assignment, allfoi_input$member_columns)
  summary <- summarise_ensemble(aggregation$ensemble, aggregation$state_population)
  old <- old_cell_summary_aggregation(allfoi_input$brazil, assignment_input$assignment)
  comparison <- summary |>
    left_join(old, by = "uf") |>
    mutate(
      new_foi_pwmean_ensemble_width = foi_pwmean_hi - foi_pwmean_lo,
      old_minus_new_midpoint = old_foi_mid_wmean - foi_pwmean_mid,
      old_minus_new_lower = old_foi_lo_wmean - foi_pwmean_lo,
      old_minus_new_upper = old_foi_hi_wmean - foi_pwmean_hi,
      old_minus_new_interval_width = old_interval_width - new_foi_pwmean_ensemble_width
    )
  audit <- tibble(
    metric = c(
      "brazil_grid_cells", "brazil_population_represented", "cells_missing_population", "cells_with_all_100_member_predictions",
      "cells_missing_any_member_prediction", "raw_unmatched_cells", "raw_ambiguous_cells", "nearest_fallback_cells",
      "near_uf_boundary_cells", "maximum_nearest_fallback_distance_m", "minimum_state_coverage", "median_state_coverage",
      "n_states_below_99pct_coverage", "n_states_below_95pct_coverage"
    ),
    value = c(
      nrow(allfoi_input$brazil), sum(allfoi_input$brazil$tot, na.rm = TRUE), allfoi_input$n_missing_population,
      sum(rowSums(is.finite(as.matrix(allfoi_input$brazil[, allfoi_input$member_columns, drop = FALSE]))) == 100L),
      sum(rowSums(is.finite(as.matrix(allfoi_input$brazil[, allfoi_input$member_columns, drop = FALSE]))) < 100L),
      assignment_input$audit$n_raw_unmatched, assignment_input$audit$n_raw_ambiguous, assignment_input$audit$n_nearest_fallback,
      assignment_input$audit$n_near_uf_boundary, assignment_input$audit$maximum_nearest_distance_m,
      min(summary$population_coverage), median(summary$population_coverage),
      sum(summary$population_coverage < .99), sum(summary$population_coverage < .95)
    )
  )
  state_assignment_audit <- assignment_input$assignment |>
    group_by(uf_abbr, state_name, assignment_method) |>
    summarise(
      n_cells = n(),
      n_near_uf_boundary = sum(near_uf_boundary, na.rm = TRUE),
      .groups = "drop"
    )

  if (nrow(aggregation$ensemble) != 2700L || anyDuplicated(aggregation$ensemble[c("uf", "ensemble_id")])) {
    stop("Expected 27 x 100 unique UF-ensemble rows.")
  }
  if (any(summary$population_coverage <= 0 | summary$population_coverage > 1)) stop("Invalid population coverage.")
  if (any(aggregation$ensemble$foi_pwmean < 0 | aggregation$ensemble$foi_equiv < 0)) stop("Negative FOI found after aggregation.")

  readr::write_csv(aggregation$ensemble, output_ensemble)
  readr::write_csv(summary, output_summary)
  readr::write_csv(comparison, file.path(table_directory, "old_vs_ensemble_first_aggregation.csv"))
  readr::write_csv(aggregation$coverage_by_member |> left_join(uf_lookup |> select(uf, uf_abbr), by = "uf") |> select(uf = uf_abbr, ensemble_id, population_valid_member), file.path(table_directory, "population_coverage_by_member.csv"))
  readr::write_csv(audit, file.path(table_directory, "aggregation_audit.csv"))
  readr::write_csv(state_assignment_audit, file.path(table_directory, "uf_cell_assignment_audit.csv"))
  readr::write_csv(summary |> filter(population_coverage < .99), file.path(table_directory, "ufs_below_99pct_population_coverage.csv"))
  make_figures(summary, aggregation$ensemble, boundary_input$boundaries, file.path(figure_directory, "brazil_uf_foi_aggregation_diagnostics.pdf"))

  ce_audit <- comparison |> filter(uf == "CE")
  write_csv(ce_audit, file.path(table_directory, "ceara_foi_aggregation_audit.csv"))
  message("[save] ", output_ensemble)
  message("[save] ", output_summary)
  message("[Ceara] old midpoint = ", signif(ce_audit$old_foi_mid_wmean, 6),
          "; ensemble arithmetic median = ", signif(ce_audit$foi_pwmean_mid, 6),
          "; equivalent median = ", signif(ce_audit$foi_equiv_mid, 6))
  invisible(list(summary = summary, comparison = comparison, audit = audit, ce_audit = ce_audit))
}


if (sys.nframe() == 0L) {
  run_brazil_uf_foi_aggregation()
}
