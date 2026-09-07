# =============================================================================
# 06_extract_renewal_v2_1_full_period_attack_decomposition.R
#
# Extract annual attack-rate decomposition from the v2.1 2015-2025 stress-test
# fit. This is a read-only post-processing script: it neither recompiles nor
# samples the Stan model.
#
# Definitions
# -----------
# * posterior infections: annual sum of latent X[t].
# * implied case detection: annual observed cases / annual latent infections.
# * S/N start: susceptible proportion at the first model week of a calendar year.
# * S/N end: the next year's first-week S/N (after the preceding year's final
#   infection and demographic update); for 2025, post-infection S/N in its last
#   model week.
# * attack increment: annual sum of X[t] divided by population in the first
#   model week of that calendar year. It is not S/N loss when population changes.
#
# Environment overrides:
#   RENEWAL_V2_1_FULL_DECOMP_FIT, RENEWAL_V2_1_FULL_DECOMP_TAG
# =============================================================================

required_packages <- c("here", "rstan", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(rstan)
})

summarise_draws <- function(draws) {
  quantiles <- stats::quantile(draws, probs = c(0.025, 0.5, 0.975))
  c(
    q025 = unname(quantiles[1]),
    median = unname(quantiles[2]),
    q975 = unname(quantiles[3])
  )
}

format_count_interval <- function(summary) {
  formatC(summary["median"], format = "f", digits = 0, big.mark = ",") |>
    paste0(
      " [", formatC(summary["q025"], format = "f", digits = 0, big.mark = ","),
      ", ", formatC(summary["q975"], format = "f", digits = 0, big.mark = ","), "]"
    )
}

format_percent_interval <- function(summary) {
  paste0(
    formatC(100 * summary["median"], format = "f", digits = 1), "% [",
    formatC(100 * summary["q025"], format = "f", digits = 1), "%, ",
    formatC(100 * summary["q975"], format = "f", digits = 1), "%]"
  )
}

if (sys.nframe() == 0) {
  default_fit <- here::here("02_Script/stan/renewal_ceara_v2_1_full_period_fit.rds")
  fit_path <- Sys.getenv("RENEWAL_V2_1_FULL_DECOMP_FIT", unset = default_fit)
  if (!file.exists(fit_path)) stop("Stress-test fit bundle not found: ", fit_path)
  bundle <- readRDS(fit_path)
  if (!isTRUE(bundle$config$stress_test) || !inherits(bundle$fit, "stanfit")) {
    stop("Expected a v2.1 full-period stress-test fit bundle")
  }

  weekly <- bundle$weekly_data
  dates <- weekly$week_start
  population <- bundle$stan_data$N_pop_t
  if (length(dates) != length(population) || any(population <= 0)) {
    stop("Weekly dates and population vector are incompatible")
  }
  draws <- rstan::extract(bundle$fit, pars = c("X", "S"))
  if (ncol(draws$X) != length(dates) || ncol(draws$S) != length(dates)) {
    stop("Latent trajectories and weekly dates are incompatible")
  }

  years <- sort(unique(as.integer(format(dates, "%Y"))))
  annual <- lapply(years, function(year) {
    index <- which(as.integer(format(dates, "%Y")) == year)
    start <- min(index)
    last <- max(index)
    infections <- rowSums(draws$X[, index, drop = FALSE])
    detection <- sum(weekly$cases[index]) / infections
    s_start <- draws$S[, start] / population[start]
    s_end <- if (last < length(dates)) {
      draws$S[, last + 1L] / population[last + 1L]
    } else {
      (draws$S[, last] - draws$X[, last]) / population[last]
    }
    attack_increment <- infections / population[start]

    c(
      year = year,
      observed_cases = sum(weekly$cases[index]),
      posterior_infections = summarise_draws(infections),
      implied_case_detection = summarise_draws(detection),
      susceptible_start = summarise_draws(s_start),
      susceptible_end = summarise_draws(s_end),
      attack_increment = summarise_draws(attack_increment)
    )
  })

  numeric_table <- do.call(rbind, annual) |> as.data.frame()
  numeric_table$year <- as.integer(numeric_table$year)
  names(numeric_table) <- sub("\\.(q025|median|q975)$", "_\\1", names(numeric_table))

  markdown_table <- data.frame(
    year = numeric_table$year,
    observed_cases = formatC(numeric_table$observed_cases, format = "f", digits = 0, big.mark = ","),
    posterior_infections = vapply(seq_len(nrow(numeric_table)), function(i) {
      format_count_interval(c(
        q025 = numeric_table$posterior_infections_q025[i],
        median = numeric_table$posterior_infections_median[i],
        q975 = numeric_table$posterior_infections_q975[i]
      ))
    }, character(1)),
    implied_case_detection = vapply(seq_len(nrow(numeric_table)), function(i) {
      format_percent_interval(c(
        q025 = numeric_table$implied_case_detection_q025[i],
        median = numeric_table$implied_case_detection_median[i],
        q975 = numeric_table$implied_case_detection_q975[i]
      ))
    }, character(1)),
    `S/N start` = vapply(seq_len(nrow(numeric_table)), function(i) {
      format_percent_interval(c(
        q025 = numeric_table$susceptible_start_q025[i],
        median = numeric_table$susceptible_start_median[i],
        q975 = numeric_table$susceptible_start_q975[i]
      ))
    }, character(1)),
    `S/N end` = vapply(seq_len(nrow(numeric_table)), function(i) {
      format_percent_interval(c(
        q025 = numeric_table$susceptible_end_q025[i],
        median = numeric_table$susceptible_end_median[i],
        q975 = numeric_table$susceptible_end_q975[i]
      ))
    }, character(1)),
    `attack increment` = vapply(seq_len(nrow(numeric_table)), function(i) {
      format_percent_interval(c(
        q025 = numeric_table$attack_increment_q025[i],
        median = numeric_table$attack_increment_median[i],
        q975 = numeric_table$attack_increment_q975[i]
      ))
    }, character(1)),
    check.names = FALSE
  )

  output_dir <- here::here("03_Output/tables")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tag <- Sys.getenv("RENEWAL_V2_1_FULL_DECOMP_TAG", unset = "")
  suffix <- if (nzchar(tag)) paste0("_", tag) else ""
  numeric_path <- file.path(output_dir, paste0("renewal_v2_1_full_period_attack_decomposition", suffix, ".csv"))
  markdown_path <- file.path(output_dir, paste0("renewal_v2_1_full_period_attack_decomposition", suffix, ".md"))
  utils::write.csv(numeric_table, numeric_path, row.names = FALSE)

  markdown_lines <- c(
    "# Annual attack-rate decomposition: v2.1 full-period stress test",
    "",
    "Posterior values are median [2.5%, 97.5%]. See script header for definitions and the 2025 end-of-window convention.",
    "",
    paste0("| ", paste(names(markdown_table), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(markdown_table)), collapse = "|"), "|"),
    apply(markdown_table, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  )
  writeLines(markdown_lines, markdown_path)
  message("[save] ", numeric_path)
  message("[save] ", markdown_path)
}
