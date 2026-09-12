# Combines the five HMC-passing conditional q-node fits into the final
# modular ensemble (task Section 12). Run only after ALL required nodes
# have passed the strict HMC gate (04_check_hmc_gate.R).
#
# This does NOT average point estimates. It samples an equal number of
# posterior draws from each node's fit, tags each retained draw with its
# node_id and q value, and pools them -- since each of the 5 strata
# represents 20% prior probability, this pooling gives each node equal
# (0.20) weight in the final ensemble by construction. The resulting
# ensemble mixes (a) transmission-posterior uncertainty conditional on q
# and (b) external q uncertainty. It is a discrete numerical-quadrature
# approximation to a continuous external q distribution -- NOT a fully
# joint Bayesian posterior in q (the case likelihood never updates q; see
# DEVELOPMENT_LOG.md).

required_packages <- c("here", "rstan", "dplyr", "tibble", "readr", "purrr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(readr); library(purrr) })

combine_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

DRAWS_PER_NODE <- as.integer(Sys.getenv("MODULAR_Q_DRAWS_PER_NODE", "1000"))

load_node_draws <- function(base_dir, node_id, n_draws) {
  out_dir <- file.path(base_dir, "outputs", node_id)
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_1_modular_q_fit_", node_id, ".rds"))
  hmc_path <- file.path(out_dir, paste0("hmc_diagnostics_", node_id, ".csv"))
  if (!file.exists(fit_path)) stop("Missing fit for node ", node_id, ": ", fit_path)
  if (!file.exists(hmc_path)) stop("Run 04_check_hmc_gate.R for node ", node_id, " first.")
  hmc <- read_csv(hmc_path, show_col_types = FALSE)
  hmc_pass <- hmc$value[hmc$metric == "divergences"] == 0 && hmc$value[hmc$metric == "max_treedepth_hits"] == 0 &&
    hmc$value[hmc$metric == "maximum_Rhat"] <= 1.01 && hmc$value[hmc$metric == "minimum_bulk_ESS"] >= 100
  if (!hmc_pass) stop("Node ", node_id, " does NOT pass the strict HMC gate -- per Section 19, stop and report; do not combine a failed node.")

  bundle <- readRDS(fit_path)
  years <- as.integer(format(as.Date(bundle$weekly_data$week_start), "%Y"))
  draws <- rstan::extract(bundle$fit, pars = c("S_prop", "immune_prop", "X", "R0_t", "R_eff_t"))
  n_available <- nrow(draws$X)
  n_take <- min(n_draws, n_available)
  set.seed(20260911L + match(node_id, c("q10", "q30", "q50", "q70", "q90")))
  idx <- sample.int(n_available, n_take)

  idx_2017 <- which(years == 2017L)[1]; idx_2022 <- which(years == 2022L)[1]; idx_2025 <- which(years == 2025L)[1]
  tibble(
    node_id = node_id, q_value = bundle$config$q_external, draw = idx,
    S_2017 = draws$S_prop[idx, idx_2017], S_2022 = draws$S_prop[idx, idx_2022], S_2025 = draws$S_prop[idx, idx_2025],
    R0_2017 = draws$R0_t[idx, idx_2017], R0_2022 = draws$R0_t[idx, idx_2022],
    cumulative_infections = rowSums(draws$X[idx, , drop = FALSE])
  ) |> mutate(
    weekly = list(list(S_prop = draws$S_prop[idx, ], immune_prop = draws$immune_prop[idx, ], X = draws$X[idx, ], R0_t = draws$R0_t[idx, ], R_eff_t = draws$R_eff_t[idx, ], week_start = bundle$weekly_data$week_start))
  )
}

run_combine_modular_posterior <- function(node_ids = c("q10", "q30", "q50", "q70", "q90")) {
  root <- combine_root()
  base_dir <- file.path(root, "02_Script", "40_renewal_model", "19_v4_1_modular_q")
  q_nodes <- read_csv(file.path(base_dir, "q_nodes.csv"), show_col_types = FALSE)

  node_data <- map(node_ids, ~ load_node_draws(base_dir, .x, DRAWS_PER_NODE))
  names(node_data) <- node_ids

  scalar_ensemble <- bind_rows(lapply(node_data, function(d) select(d, -weekly))) |>
    left_join(select(q_nodes, node_id, prior_weight), by = "node_id")
  write_csv(scalar_ensemble, file.path(base_dir, "combined_modular_ensemble_scalars.csv"))

  # Weekly trajectory ensemble: stack all nodes' sampled draws (equal count
  # per node = equal 0.20 weight per node in the pooled ensemble).
  week_start <- node_data[[1]]$weekly[[1]]$week_start
  weekly_stack <- function(field) do.call(rbind, lapply(node_data, function(d) d$weekly[[1]][[field]]))
  S_all <- weekly_stack("S_prop"); U_all <- weekly_stack("immune_prop"); X_all <- weekly_stack("X")
  R0_all <- weekly_stack("R0_t"); Reff_all <- weekly_stack("R_eff_t")

  weekly_summary <- tibble(
    week_start = week_start,
    S_prop_median = apply(S_all, 2, median), S_prop_q025 = apply(S_all, 2, quantile, .025), S_prop_q975 = apply(S_all, 2, quantile, .975),
    U_prop_median = apply(U_all, 2, median), U_prop_q025 = apply(U_all, 2, quantile, .025), U_prop_q975 = apply(U_all, 2, quantile, .975),
    X_median = apply(X_all, 2, median), X_q025 = apply(X_all, 2, quantile, .025), X_q975 = apply(X_all, 2, quantile, .975),
    R0_median = apply(R0_all, 2, median), R0_q025 = apply(R0_all, 2, quantile, .025), R0_q975 = apply(R0_all, 2, quantile, .975),
    Reff_median = apply(Reff_all, 2, median), Reff_q025 = apply(Reff_all, 2, quantile, .025), Reff_q975 = apply(Reff_all, 2, quantile, .975)
  )
  write_csv(weekly_summary, file.path(base_dir, "combined_modular_ensemble_weekly.csv"))

  scalar_summary <- scalar_ensemble |> summarise(
    q_mean = weighted.mean(q_value, prior_weight), q_sd = sqrt(weighted.mean((q_value - weighted.mean(q_value, prior_weight))^2, prior_weight)),
    S_2017_median = median(S_2017), S_2017_q025 = quantile(S_2017, .025), S_2017_q975 = quantile(S_2017, .975),
    S_2022_median = median(S_2022), S_2022_q025 = quantile(S_2022, .025), S_2022_q975 = quantile(S_2022, .975),
    S_2025_median = median(S_2025), S_2025_q025 = quantile(S_2025, .025), S_2025_q975 = quantile(S_2025, .975),
    R0_2017_median = median(R0_2017), R0_2017_q025 = quantile(R0_2017, .025), R0_2017_q975 = quantile(R0_2017, .975),
    R0_2022_median = median(R0_2022), R0_2022_q025 = quantile(R0_2022, .025), R0_2022_q975 = quantile(R0_2022, .975),
    cumulative_infections_median = median(cumulative_infections), cumulative_infections_q025 = quantile(cumulative_infections, .025), cumulative_infections_q975 = quantile(cumulative_infections, .975)
  )
  write_csv(scalar_summary, file.path(base_dir, "combined_modular_ensemble_scalar_summary.csv"))

  message("[combine] modular ensemble built from ", length(node_ids), " nodes x ", DRAWS_PER_NODE, " draws each = ", nrow(scalar_ensemble), " total draws.")
  message("[combine] IMPORTANT: this is modular uncertainty propagation (external q distribution x conditional transmission posteriors), NOT a fully joint Bayesian posterior in q.")
  invisible(list(scalar_ensemble = scalar_ensemble, weekly_summary = weekly_summary, scalar_summary = scalar_summary))
}

if (sys.nframe() == 0L) run_combine_modular_posterior()
