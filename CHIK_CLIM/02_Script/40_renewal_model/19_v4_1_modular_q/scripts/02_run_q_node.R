# Fits ONE conditional q-node of the v4.1 modular-q model.
#
# q_external is DATA for this fit (looked up from q_nodes.csv by node id) --
# never a sampled parameter. Everything else is byte-for-byte identical to
# the HMC-passing v4.0/td14 backbone (verified by 00_audit_v4_0_reference.R).
# The compiled Stan model is built ONCE per R session and reused across
# nodes if this script is sourced repeatedly in the same session (e.g. from
# a driver loop) rather than re-invoked via Rscript each time; when invoked
# via Rscript per node, rstan's on-disk auto_write cache (same .stan file
# hash) avoids recompilation on the 2nd-5th node regardless.
#
# Per-chain CmdStan-style CSVs are written to a PERMANENT path from the
# start of sampling (outputs/<node_id>/chain_<id>.csv), so
# 03_monitor_v4_1_fit.R can inspect live progress. save_warmup is on by
# default (confirmed via smoke test: "# save_warmup=1" appears in the
# written CSV header even without explicitly setting it), and refresh=25
# gives frequent (if not live-streamed on Windows PSOCK) console reporting.
#
# Usage: Q_NODE_ID=q50 Rscript 02_run_q_node.R

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

modular_q_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

# Reuse v4.0's exact helper functions (generation_weights, compute_hmc_gate,
# v4_initial_values, build_state_weekly) -- sourcing, not reimplementing,
# per the "do not reconstruct from memory" instruction. v4_1_root()/
# v4_1_paths()/make_v4_1_data() from the ABORTED joint-q attempt
# (18_v4_1_minimal_no_vaccine) are NOT sourced here -- this script is
# independent of that folder.
root <- modular_q_root()
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "40_renewal_model", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))

modular_q_dir <- function(root = modular_q_root()) file.path(root, "02_Script", "40_renewal_model", "19_v4_1_modular_q")

# Mirrors v4.0's make_v4_data() exactly, except the data field is named
# q_external (matching the Stan file) instead of q_fixed, and q is not
# passed at all here -- it is added by the caller from the node table.
make_v4_1_modular_data <- function(weekly, imports_per_week, year_effect_prior_sd) {
  dates <- as.Date(weekly$week_start)
  years <- sort(unique(as.integer(format(dates, "%Y"))))
  expected_dates <- seq(min(dates), max(dates), by = "week")
  if (!identical(dates, expected_dates)) stop("Weekly data must be complete and consecutive.")
  accounting_error <- max(abs(
    weekly$N_end - weekly$N_start - weekly$births + weekly$all_cause_deaths - weekly$net_population_reconciliation
  ))
  if (accounting_error > 1e-6) stop("Demographic accounting failed before Stan fitting.")
  list(
    stan_data = list(
      N = nrow(weekly), G = 8L, C = as.integer(weekly$cases), w = generation_weights(),
      N_start = as.numeric(weekly$N_start), N_end = as.numeric(weekly$N_end),
      births = as.numeric(weekly$births), deaths = as.numeric(weekly$all_cause_deaths),
      Y = length(years), year_id = match(as.integer(format(dates, "%Y")), years),
      seasonal_sin = sin(2 * pi * seq_len(nrow(weekly)) / 52.1775),
      seasonal_cos = cos(2 * pi * seq_len(nrow(weekly)) / 52.1775),
      imports_per_week = imports_per_week, year_effect_prior_sd = year_effect_prior_sd
    ),
    years = years,
    demographic_accounting_error = accounting_error
  )
}

run_q_node <- function(node_id = Sys.getenv("Q_NODE_ID", "q50")) {
  base_dir <- modular_q_dir(root)
  q_nodes_path <- file.path(base_dir, "q_nodes.csv")
  if (!file.exists(q_nodes_path)) stop("Run 01_build_q_nodes.R first: ", q_nodes_path, " not found.")
  q_nodes <- read_csv(q_nodes_path, show_col_types = FALSE)
  node_row <- filter(q_nodes, node_id == !!node_id)
  if (!nrow(node_row)) stop("Unknown node_id '", node_id, "'. Known nodes: ", paste(q_nodes$node_id, collapse = ", "))
  q_value <- node_row$q_value[1]

  out_dir <- file.path(base_dir, "outputs", node_id)
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  settings <- list(
    date_start = as.Date("2015-01-04"), date_end = as.Date("2025-12-21"),
    imports_per_week = 1, year_effect_prior_sd = 0.40,
    iter = 2000L, warmup = 1000L, chains = 4L,
    adapt_delta = .95, max_treedepth = 14L
  )

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |>
    filter(week_start >= settings$date_start, week_start <= settings$date_end)
  demography_path <- file.path(root, "01_Data", "ceara_weekly_demography_2015_2025.rds")
  weekly_demography <- readRDS(demography_path) |>
    filter(week_start >= settings$date_start, week_start <= settings$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) {
    stop("Case and demographic panels do not have the same full weekly sequence.")
  }

  prepared <- make_v4_1_modular_data(weekly, settings$imports_per_week, settings$year_effect_prior_sd)
  prepared$stan_data$q_external <- q_value

  node_seed <- 20260911L + match(node_id, q_nodes$node_id) * 1000L
  sample_file_base <- file.path(progress_dir, "chain")
  message(sprintf(
    "[%s] q_external=%.5f | %d weeks | seed=%d | sample_file base=%s",
    node_id, q_value, prepared$stan_data$N, node_seed, sample_file_base
  ))

  stan_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_1_modular_q.stan")
  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data,
    chains = settings$chains, iter = settings$iter, warmup = settings$warmup,
    seed = node_seed,
    init = v4_initial_values(prepared$stan_data),
    control = list(adapt_delta = settings$adapt_delta, max_treedepth = settings$max_treedepth),
    refresh = 25,
    sample_file = sample_file_base
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs"), settings$max_treedepth)

  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly,
    config = c(settings, list(
      model_version = "v4_1_modular_q", node_id = node_id, q_external = q_value, years = prepared$years,
      demographic_accounting_error = prepared$demographic_accounting_error,
      demography_source = normalizePath(demography_path, winslash = "/", mustWork = TRUE),
      elapsed_seconds = elapsed_seconds, stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE),
      seed = node_seed,
      q_role = "MODULAR: fixed data within this fit; varied ACROSS fits per q_nodes.csv to propagate external reporting-fraction uncertainty; not estimated jointly with transmission parameters (see DEVELOPMENT_LOG.md)."
    )),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_1_modular_q_fit_", node_id, ".rds"))
  saveRDS(bundle, fit_path)

  v4_0_td14_bundle <- readRDS(file.path(root, "02_Script", "stan", "renewal_ceara_v4_0_minimal_no_vaccine_fit_td14.rds"))
  benchmark_seconds <- v4_0_td14_bundle$config$elapsed_seconds
  runtime_ratio <- elapsed_seconds / benchmark_seconds
  runtime_flag <- if (runtime_ratio > 3) "COMPUTATIONALLY_EXPENSIVE" else "normal"
  runtime_record <- tibble(
    node_id = node_id, wall_clock_seconds = elapsed_seconds, v4_0_td14_benchmark_seconds = benchmark_seconds,
    ratio_vs_benchmark = runtime_ratio, seconds_per_iteration = elapsed_seconds / (settings$iter * settings$chains),
    flag = runtime_flag
  )
  write_csv(runtime_record, file.path(out_dir, paste0("runtime_", node_id, ".csv")))

  message(sprintf(
    "[%s] DONE | HMC gate: %s | elapsed=%.0fs (%.2fx td14 benchmark of %.0fs) | flag=%s",
    node_id, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, runtime_ratio, benchmark_seconds, runtime_flag
  ))
  if (runtime_flag == "COMPUTATIONALLY_EXPENSIVE") {
    message("[", node_id, "] WARNING: >3x the v4.0/td14 benchmark wall-clock time -- flagged COMPUTATIONALLY_EXPENSIVE even if HMC passes.")
  }
  invisible(bundle)
}

if (sys.nframe() == 0L) run_q_node()
