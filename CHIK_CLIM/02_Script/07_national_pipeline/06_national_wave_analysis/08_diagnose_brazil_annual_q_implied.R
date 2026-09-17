# Posterior-draw diagnostic of implied confirmed-case detection by UF and year.
# q is never fed back into the conditional annual-SHAPE likelihood.

required_q_packages <- c("here", "rstan", "dplyr", "readr", "ggplot2", "scales")
missing_q_packages <- required_q_packages[!vapply(required_q_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_q_packages)) stop("Missing package(s): ", paste(missing_q_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr); library(ggplot2); library(scales) })

q_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}

summarise_q_draws <- function(x, prefix = "q_implied") {
  tibble(
    !!paste0(prefix, "_q025") := quantile(x, .025),
    !!paste0(prefix, "_median") := median(x),
    !!paste0(prefix, "_q975") := quantile(x, .975),
    !!paste0(prefix, "_prob_gt_one") := mean(x > 1)
  )
}

run_annual_q_implied_diagnostic <- function() {
  root <- q_root()
  table_dir <- file.path(root, "03_Output", "07_national_pipeline", "tables", "national_wave_analysis", "q_implied")
  figure_dir <- file.path(root, "03_Output", "07_national_pipeline", "figures", "national_wave_analysis", "q_implied")
  fit_path <- file.path(root, "03_Output", "07_national_pipeline", "model_fits", "national_wave_analysis", "brazil_chik_annual_shape_fits.rds")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(fit_path)) stop("Annual SHAPE fit bundle is missing.")
  bundle <- readRDS(fit_path)
  hmc <- read_csv(file.path(root, "03_Output", "07_national_pipeline", "tables", "national_wave_analysis", "brazil_chik_annual_shape_hmc_gate.csv"), show_col_types = FALSE) |>
    select(state, annual_hmc_pass = hmc_pass)

  # Existing externally informed reference: p_symp ~ Beta(30,28) and
  # symptomatic reporting ~ Beta(20,60).  Its central 95% interval defines
  # transparent high/low plausibility flags; it is not used in fitting.
  set.seed(20261001L)
  q_reference <- rbeta(500000L, 30, 28) * rbeta(500000L, 20, 60)
  reference_low <- quantile(q_reference, .025)
  reference_high <- quantile(q_reference, .975)

  annual_rows <- list(); total_rows <- list()
  draw_rows <- list(); total_draw_rows <- list()
  for (state in names(bundle$fits)) {
    entry <- bundle$fits[[state]]
    draws_x <- rstan::extract(entry$fit, pars = "X", permuted = TRUE)$X
    annual <- entry$input$annual
    # sweep(x, "/") would calculate X/C.  The requested implied detection
    # probability is C/X for every posterior draw and calendar year.
    q_year <- sweep(draws_x, 2L, annual$observed_cases, FUN = function(x, c) c / x)
    q_total <- sum(annual$observed_cases) / rowSums(draws_x)
    annual_rows[[state]] <- bind_rows(lapply(seq_len(ncol(q_year)), function(y) {
      summary <- summarise_q_draws(q_year[, y])
      bind_cols(tibble(state = state, year = annual$year[y], observed_cases = annual$observed_cases[y],
                       reference_q025 = reference_low, reference_q975 = reference_high,
                       prob_below_reference = mean(q_year[, y] < reference_low),
                       prob_above_reference = mean(q_year[, y] > reference_high)), summary)
    }))
    total_rows[[state]] <- bind_cols(
      tibble(state = state, observed_cases_total = sum(annual$observed_cases),
             reference_q025 = reference_low, reference_q975 = reference_high,
             prob_below_reference = mean(q_total < reference_low), prob_above_reference = mean(q_total > reference_high)),
      summarise_q_draws(q_total, "q_total_implied")
    )
    total_draw_rows[[state]] <- tibble(
      state = state,
      draw = seq_along(q_total),
      q_total_implied = q_total
    )
    draw_rows[[state]] <- as_tibble(q_year, .name_repair = "minimal") |>
      setNames(paste0("year_", annual$year)) |>
      mutate(draw = row_number(), state = state) |>
      tidyr::pivot_longer(starts_with("year_"), names_to = "year", values_to = "q_implied") |>
      mutate(year = as.integer(sub("year_", "", year)))
  }
  annual_summary <- bind_rows(annual_rows) |>
    left_join(hmc, by = "state") |>
    mutate(
      flag_q_gt_one = q_implied_prob_gt_one > .025,
      flag_implausibly_high = prob_above_reference >= .95,
      flag_implausibly_low = prob_below_reference >= .95,
      flag_label = case_when(
        !annual_hmc_pass ~ "annual HMC fail",
        flag_q_gt_one ~ "Pr(q > 1) > 2.5%",
        flag_implausibly_high ~ "high vs external 95% reference",
        flag_implausibly_low ~ "low vs external 95% reference",
        TRUE ~ "none"
      )
    )
  total_summary <- bind_rows(total_rows) |>
    left_join(hmc, by = "state") |>
    mutate(
      flag_q_gt_one = q_total_implied_prob_gt_one > .025,
      flag_implausibly_high = prob_above_reference >= .95,
      flag_implausibly_low = prob_below_reference >= .95
    )
  write_csv(annual_summary, file.path(table_dir, "brazil_chik_annual_q_implied_summary.csv"))
  write_csv(total_summary, file.path(table_dir, "brazil_chik_total_q_implied_summary.csv"))
  saveRDS(bind_rows(draw_rows), file.path(table_dir, "brazil_chik_annual_q_implied_posterior_draws.rds"), compress = "xz")
  saveRDS(bind_rows(total_draw_rows), file.path(table_dir, "brazil_chik_total_q_implied_posterior_draws.rds"), compress = "xz")

  heatmap_data <- annual_summary |> mutate(fill_value = if_else(annual_hmc_pass, q_implied_median, NA_real_),
    text = case_when(flag_q_gt_one ~ "!", flag_implausibly_high ~ "H", flag_implausibly_low ~ "L", TRUE ~ ""))
  p <- ggplot(heatmap_data, aes(year, reorder(state, q_implied_median), fill = fill_value)) +
    geom_tile(colour = "white", linewidth = .25) +
    geom_text(aes(label = text), fontface = "bold", size = 3) +
    # Square-root scaling retains zero-case years (q = 0) while still making
    # low detection values legible; a log scale would drop those valid cells.
    scale_fill_viridis_c(trans = "sqrt", na.value = "grey75", name = expression(median(q[implied]))) +
    scale_x_continuous(breaks = 2015:2025) +
    labs(title = "Posterior implied detection by UF and year", subtitle = "!: Pr(q > 1) > 2.5%; H/L: posterior probability >=95% above/below external 95% detection reference; grey: annual HMC gate failed.", x = NULL, y = NULL) +
    theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(figure_dir, "brazil_chik_annual_q_implied_heatmap.pdf"), p, width = 250, height = 180, units = "mm", device = cairo_pdf)
  invisible(list(annual = annual_summary, total = total_summary))
}

if (sys.nframe() == 0L) run_annual_q_implied_diagnostic()
