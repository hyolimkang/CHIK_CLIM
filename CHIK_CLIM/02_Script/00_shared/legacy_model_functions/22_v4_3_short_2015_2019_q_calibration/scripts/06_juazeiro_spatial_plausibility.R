# Mode audit Part 5: spatial plausibility of the high-q mode's implied
# Juazeiro/Ceara contrast (~17x), using data INDEPENDENT of the Juazeiro
# serology likelihood -- confirmed CHIKV incidence (same case definition as
# C[t]: cases_confirmed, CLASSI_FIN==13) from the muni-week SINAN panel,
# NOT the serosurvey itself. Pre-survey period: 2015-01-04 through the
# Juazeiro survey midpoint (2018-09-15).

required_packages <- c("here", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tibble) })

spatial_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

run_juazeiro_spatial_plausibility <- function() {
  root <- spatial_root()
  out_dir <- file.path(root, "03_Output", "model_fits", "ceara", "v4_3_short_2015_2019_q_calibration")

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  panel <- panel |> mutate(muni6 = as.character(muni6), week_start = as.Date(week_start))
  survey_midpoint <- as.Date("2018-09-15")
  ceara_munis <- panel |> filter(substr(muni6, 1, 2) == "23")

  pre_survey <- ceara_munis |> filter(week_start <= survey_midpoint)

  # Cumulative confirmed cases and average population per municipality over
  # the pre-survey period (same case definition as C[t]: cases_confirmed).
  by_muni <- pre_survey |> group_by(muni6) |> summarise(
    cumulative_confirmed = sum(cases_confirmed, na.rm = TRUE),
    mean_population = mean(population, na.rm = TRUE), .groups = "drop"
  ) |> mutate(incidence_per_capita = cumulative_confirmed / mean_population)

  juazeiro_muni6 <- "230730"
  juazeiro_row <- by_muni |> filter(muni6 == juazeiro_muni6)
  if (!nrow(juazeiro_row)) stop("Juazeiro do Norte (muni6=230730) not found in panel.")

  ceara_state_incidence <- sum(by_muni$cumulative_confirmed) / sum(by_muni$mean_population)
  juazeiro_incidence <- juazeiro_row$incidence_per_capita
  ratio <- juazeiro_incidence / ceara_state_incidence

  by_muni <- by_muni |> mutate(ratio_to_state = incidence_per_capita / ceara_state_incidence) |>
    arrange(desc(incidence_per_capita)) |>
    mutate(percentile_rank = 100 * (1 - (row_number() - 1) / (n() - 1)))
  juazeiro_percentile <- by_muni |> filter(muni6 == juazeiro_muni6) |> pull(percentile_rank)
  p90_ratio <- quantile(by_muni$ratio_to_state, .90, names = FALSE)
  p95_ratio <- quantile(by_muni$ratio_to_state, .95, names = FALSE)
  max_ratio <- max(by_muni$ratio_to_state)

  result <- tibble(
    metric = c("juazeiro_incidence_per_capita", "ceara_state_incidence_per_capita", "juazeiro_over_ceara_ratio",
               "juazeiro_percentile_rank_among_ceara_munis", "ceara_munis_p90_ratio_to_state",
               "ceara_munis_p95_ratio_to_state", "ceara_munis_max_ratio_to_state", "n_ceara_munis",
               "high_q_mode_required_ratio"),
    value = c(juazeiro_incidence, ceara_state_incidence, ratio, juazeiro_percentile, p90_ratio, p95_ratio, max_ratio, nrow(by_muni), 25.5 / 1.5)
  )
  write_csv(result, file.path(out_dir, "juazeiro_spatial_plausibility.csv"))
  write_csv(by_muni, file.path(out_dir, "ceara_muni_incidence_ranking_presurvey.csv"))

  message("[spatial plausibility] pre-survey period: 2015-01-04 to ", survey_midpoint)
  message(sprintf("[spatial plausibility] Juazeiro incidence/capita = %.5f | Ceara state incidence/capita = %.5f | ratio = %.2fx",
                   juazeiro_incidence, ceara_state_incidence, ratio))
  message(sprintf("[spatial plausibility] Juazeiro percentile rank among %d Ceara munis (by pre-survey confirmed-case incidence) = %.1f%%",
                   nrow(by_muni), juazeiro_percentile))
  message(sprintf("[spatial plausibility] Ceara munis' 90th/95th/max percentile ratio-to-state = %.2fx / %.2fx / %.2fx", p90_ratio, p95_ratio, max_ratio))
  message(sprintf("[spatial plausibility] high-q mode REQUIRES Juazeiro to be ~%.1fx the Ceara state average", 25.5 / 1.5))
  message(sprintf("[spatial plausibility] Juazeiro's OWN case-based ratio-to-state is %.2fx -- %s the high-q mode's required %.1fx",
                   ratio, if (ratio < 25.5/1.5 * 0.5) "FAR BELOW" else if (ratio < 25.5/1.5) "below" else "consistent with or above", 25.5/1.5))
  invisible(result)
}

if (sys.nframe() == 0L) run_juazeiro_spatial_plausibility()
