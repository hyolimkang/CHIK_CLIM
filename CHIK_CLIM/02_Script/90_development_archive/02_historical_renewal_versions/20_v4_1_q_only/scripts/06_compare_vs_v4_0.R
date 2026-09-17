# Direct side-by-side comparison of a v4.1 q-only stage fit against the
# HMC-passing v4.0/td14 reference (task Section 9). Quantifies exactly how
# much computational difficulty q alone introduces.
#
# Usage: STAGE=B Rscript 06_compare_vs_v4_0.R

required_packages <- c("here", "rstan", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(tibble); library(readr) })

compare_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

sampler_stats <- function(fit) {
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  treedepth <- unlist(lapply(sp, function(x) x[, "treedepth__"]))
  leapfrog <- unlist(lapply(sp, function(x) x[, "n_leapfrog__"]))
  list(
    treedepth_median = median(treedepth), treedepth_p95 = quantile(treedepth, .95, names = FALSE),
    td14_hits = sum(treedepth == 14), n_leapfrog_median = median(leapfrog), n_leapfrog_p95 = quantile(leapfrog, .95, names = FALSE),
    divergences = rstan::get_num_divergent(fit), bfmi_min = min(rstan::get_bfmi(fit)), bfmi_max = max(rstan::get_bfmi(fit))
  )
}

run_compare_vs_v4_0 <- function(stage = Sys.getenv("STAGE", "B")) {
  root <- compare_root()
  out_dir <- file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "20_v4_1_q_only", "outputs", paste0("stage", stage))
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_1_q_only_fit_stage", stage, ".rds"))
  if (!file.exists(fit_path)) stop("No fit found for stage '", stage, "': ", fit_path)
  q_only_bundle <- readRDS(fit_path)
  v4_0_bundle <- readRDS(file.path(root, "02_Script", "stan", "renewal_ceara_v4_0_minimal_no_vaccine_fit_td14.rds"))

  v4_0_stats <- sampler_stats(v4_0_bundle$fit)
  q_only_stats <- sampler_stats(q_only_bundle$fit)
  v4_0_hmc <- read_csv(file.path(root, "03_Output", "tables", "renewal_v4_0_minimal_no_vaccine_td14", "hmc_diagnostics.csv"), show_col_types = FALSE)
  q_only_hmc_path <- file.path(out_dir, paste0("hmc_diagnostics_stage", stage, ".csv"))
  q_only_hmc <- if (file.exists(q_only_hmc_path)) read_csv(q_only_hmc_path, show_col_types = FALSE) else NULL

  comparison <- tibble(
    metric = c("runtime_seconds", "median_treedepth", "p95_treedepth", "td14_hits", "median_n_leapfrog", "p95_n_leapfrog",
               "divergences", "BFMI_min", "BFMI_max", "max_Rhat", "min_bulk_ESS"),
    v4_0_td14 = c(
      v4_0_bundle$config$elapsed_seconds, v4_0_stats$treedepth_median, v4_0_stats$treedepth_p95, v4_0_stats$td14_hits,
      v4_0_stats$n_leapfrog_median, v4_0_stats$n_leapfrog_p95, v4_0_stats$divergences, v4_0_stats$bfmi_min, v4_0_stats$bfmi_max,
      v4_0_hmc$value[v4_0_hmc$metric == "maximum_Rhat"], v4_0_hmc$value[v4_0_hmc$metric == "minimum_bulk_ESS"]
    ),
    v4_1_q_only = c(
      q_only_bundle$config$elapsed_seconds, q_only_stats$treedepth_median, q_only_stats$treedepth_p95, q_only_stats$td14_hits,
      q_only_stats$n_leapfrog_median, q_only_stats$n_leapfrog_p95, q_only_stats$divergences, q_only_stats$bfmi_min, q_only_stats$bfmi_max,
      if (!is.null(q_only_hmc)) q_only_hmc$value[q_only_hmc$metric == "maximum_Rhat"] else NA_real_,
      if (!is.null(q_only_hmc)) q_only_hmc$value[q_only_hmc$metric == "minimum_bulk_ESS"] else NA_real_
    )
  )
  comparison$ratio_q_only_over_v4_0 <- comparison$v4_1_q_only / comparison$v4_0_td14

  write_csv(comparison, file.path(out_dir, paste0("comparison_vs_v4_0_stage", stage, ".csv")))
  message("[stage", stage, "] v4.0/td14 vs v4.1 q-only comparison:")
  print(as.data.frame(comparison))

  # NOTE: v4.0/td14 ran 1000+1000 iterations; this stage may have run far
  # fewer -- runtime_seconds and per-metric ratios must be read alongside
  # the iteration counts, not as a pure like-for-like unless STAGE=D.
  message(sprintf("[stage%s] iteration counts: v4.0/td14 = 1000+1000; this stage = %d+%d", stage, q_only_bundle$config$warmup, q_only_bundle$config$iter_sampling))
  invisible(comparison)
}

if (sys.nframe() == 0L) run_compare_vs_v4_0()
