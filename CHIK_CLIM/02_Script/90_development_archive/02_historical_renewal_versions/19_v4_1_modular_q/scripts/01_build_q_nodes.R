# Builds the five stratified q nodes for modular reporting-fraction
# uncertainty propagation (task spec Sections 2-3).
#
# External epidemiological information: p_symptomatic ~ Beta(30,28),
# p_detect_given_symptomatic ~ Beta(20,60), q = p_symptomatic * p_detect.
# >=1e6 Monte Carlo draws (here: 2e6, reproducible seed) of the INDUCED
# product distribution are generated; five equal-probability strata
# (0-20/20-40/40-60/60-80/80-100%) are represented by their 10th/30th/50th/
# 70th/90th percentiles, taken DIRECTLY from the induced draws -- NOT from
# the previously-checked Beta(16.2,108.7) approximation, per explicit
# instruction not to hard-code approximate values if the exact distribution
# differs. The Beta approximation is still reported alongside for
# comparison/documentation only.

required_packages <- c("here", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(tibble); library(readr) })

q_nodes_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

run_build_q_nodes <- function() {
  root <- q_nodes_root()
  out_dir <- file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "19_v4_1_modular_q")
  n_mc <- 2e6
  set.seed(20260911L)
  p_symptomatic <- rbeta(n_mc, 30, 28)
  p_detect <- rbeta(n_mc, 20, 60)
  q_induced <- p_symptomatic * p_detect

  full_probs <- c(.025, .05, .10, .25, .50, .75, .90, .95, .975)
  induced_summary <- tibble(
    statistic = c("mean", "sd", paste0("q", full_probs * 100)),
    value = c(mean(q_induced), sd(q_induced), quantile(q_induced, full_probs, names = FALSE))
  )
  write_csv(induced_summary, file.path(out_dir, "q_induced_distribution_summary.csv"))
  message("[q-nodes] induced product distribution (n=", n_mc, "):")
  print(induced_summary, n = Inf)

  set.seed(20260911L)
  beta_approx_draws <- rbeta(n_mc, 16.2, 108.7)
  beta_summary <- tibble(
    statistic = c("mean", "sd", paste0("q", full_probs * 100)),
    value = c(mean(beta_approx_draws), sd(beta_approx_draws), quantile(beta_approx_draws, full_probs, names = FALSE))
  )
  comparison <- merge(induced_summary, beta_summary, by = "statistic", suffixes = c("_induced", "_Beta_16.2_108.7"))
  write_csv(comparison, file.path(out_dir, "q_induced_vs_Beta_16.2_108.7_comparison.csv"))
  message("[q-nodes] comparison against previously-checked Beta(16.2, 108.7) written.")

  # Five equal-probability strata (0-20/20-40/40-60/60-80/80-100%),
  # represented by the 10th/30th/50th/70th/90th percentile of the INDUCED
  # distribution -- the median of each stratum under a uniform-within-
  # stratum approximation, and exactly what Section 3 specifies.
  node_quantiles <- c(0.10, 0.30, 0.50, 0.70, 0.90)
  node_values <- quantile(q_induced, node_quantiles, names = FALSE)
  q_nodes <- tibble(
    node_id = paste0("q", node_quantiles * 100),
    prior_quantile = node_quantiles,
    q_value = node_values,
    prior_weight = 0.20
  )
  write_csv(q_nodes, file.path(out_dir, "q_nodes.csv"))
  message("[q-nodes] five stratified nodes (from the exact induced distribution, not the Beta approximation):")
  print(q_nodes)
  invisible(q_nodes)
}

if (sys.nframe() == 0L) run_build_q_nodes()
