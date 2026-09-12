# Posterior predictive checks per q node and for the combined modular
# ensemble (task Section 14). Explicit rows for 2017/2021/2022/2023, plus
# expected/observed ratio for 2017 and 2022. Does NOT tune q nodes to fix
# any visible mismatch -- if all plausible q nodes underpredict 2022, that
# is reported as a STRUCTURAL_MODEL_MISMATCH, not something to patch here.

required_packages <- c("here", "rstan", "dplyr", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(readr) })

ppc_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

HIGHLIGHT_YEARS <- c(2017L, 2021L, 2022L, 2023L)

annual_ppc_for_node <- function(base_dir, node_id) {
  fit_path <- file.path(base_dir, "outputs", node_id, paste0("renewal_ceara_v4_1_modular_q_fit_", node_id, ".rds"))
  bundle <- readRDS(fit_path)
  years <- as.integer(format(as.Date(bundle$weekly_data$week_start), "%Y"))
  draws <- rstan::extract(bundle$fit, pars = "expected_reported_cases")$expected_reported_cases
  bind_rows(lapply(sort(unique(years)), function(year) {
    idx <- which(years == year)
    expected <- rowSums(draws[, idx, drop = FALSE])
    observed <- sum(bundle$weekly_data$cases[idx])
    tibble(
      node_id = node_id, q_value = bundle$config$q_external, year = year, observed_cases = observed,
      expected_cases_median = median(expected), expected_cases_q025 = quantile(expected, .025), expected_cases_q975 = quantile(expected, .975),
      ratio_expected_over_observed = median(expected) / observed
    )
  }))
}

run_make_ppc <- function(node_ids = c("q10", "q30", "q50", "q70", "q90")) {
  root <- ppc_root()
  base_dir <- file.path(root, "02_Script", "40_renewal_model", "19_v4_1_modular_q")

  available <- node_ids[file.exists(file.path(base_dir, "outputs", node_ids, paste0("renewal_ceara_v4_1_modular_q_fit_", node_ids, ".rds")))]
  if (!length(available)) stop("No completed node fits found.")
  message("[ppc] computing annual PPC for nodes: ", paste(available, collapse = ", "))

  per_node <- bind_rows(lapply(available, function(n) annual_ppc_for_node(base_dir, n)))
  write_csv(per_node, file.path(base_dir, "annual_ppc_per_node.csv"))

  highlight <- filter(per_node, year %in% HIGHLIGHT_YEARS) |>
    select(node_id, q_value, year, observed_cases, expected_cases_median, expected_cases_q025, expected_cases_q975, ratio_expected_over_observed)
  write_csv(highlight, file.path(base_dir, "annual_ppc_highlight_years_per_node.csv"))
  message("[ppc] highlight years (2017/2021/2022/2023) by node:")
  print(as.data.frame(highlight))

  # Combined-ensemble annual PPC, if the combined scalar/weekly ensemble
  # exists (requires 05_combine_modular_posterior.R to have been run).
  combined_weekly_path <- file.path(base_dir, "combined_modular_ensemble_weekly.csv")
  if (file.exists(combined_weekly_path)) {
    weekly <- read_csv(combined_weekly_path, show_col_types = FALSE) |> mutate(year = as.integer(format(week_start, "%Y")))
    # X (latent infections), not expected_reported_cases, is stored in the
    # weekly ensemble; approximate combined expected reported cases using
    # each node's own q, so re-derive from per-node expected cases instead:
    # the per-node table already has this; combined PPC is the q-weighted
    # (0.20 each) pool of per-node expected annual cases.
    combined_annual <- per_node |> group_by(year) |> summarise(
      observed_cases = first(observed_cases),
      expected_cases_median = median(expected_cases_median), # pooled across nodes, equal weight
      expected_cases_q025 = min(expected_cases_q025), expected_cases_q975 = max(expected_cases_q975),
      .groups = "drop"
    ) |> mutate(ratio_expected_over_observed = expected_cases_median / observed_cases)
    write_csv(combined_annual, file.path(base_dir, "annual_ppc_combined_ensemble.csv"))
    combined_highlight <- filter(combined_annual, year %in% HIGHLIGHT_YEARS)
    message("[ppc] combined-ensemble highlight years:")
    print(as.data.frame(combined_highlight))
  } else {
    message("[ppc] combined ensemble not built yet -- run 05_combine_modular_posterior.R for the pooled PPC.")
  }

  underpredicts_2022 <- all(filter(per_node, year == 2022L)$ratio_expected_over_observed < 0.9)
  if (underpredicts_2022 && length(available) == 5) {
    message("[ppc] *** STRUCTURAL_MODEL_MISMATCH: ALL FIVE q nodes underpredict 2022 (expected/observed < 0.9). This is a remaining structural limitation -- do NOT respond with time-varying q. ***")
  }
  invisible(list(per_node = per_node, highlight = highlight))
}

if (sys.nframe() == 0L) run_make_ppc()
