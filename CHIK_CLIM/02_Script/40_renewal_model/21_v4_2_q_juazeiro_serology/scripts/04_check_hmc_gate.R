# Post-hoc strict HMC gate check for one completed v4.2 stage fit.
# Gate (unchanged from v4.0): post-warmup divergences=0, post-warmup
# max-treedepth(14) hits=0, max Rhat<=1.01, min bulk ESS>=100, BFMI>=0.30 in
# every chain.
#
# Usage: STAGE=B Rscript 04_check_hmc_gate.R

required_packages <- c("here", "rstan", "posterior", "dplyr", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(posterior); library(dplyr); library(tibble); library(readr) })

gate_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

run_check_hmc_gate <- function(stage = Sys.getenv("STAGE", "B")) {
  root <- gate_root()
  out_dir <- file.path(root, "02_Script", "40_renewal_model", "21_v4_2_q_juazeiro_serology", "outputs", paste0("stage", stage))
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_2_fit_stage", stage, ".rds"))
  if (!file.exists(fit_path)) stop("No fit found for stage '", stage, "': ", fit_path)
  bundle <- readRDS(fit_path)
  fit <- bundle$fit

  core_pars <- c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs", "q")
  raw <- rstan::extract(fit, pars = core_pars, permuted = FALSE, inc_warmup = FALSE)
  draws <- posterior::as_draws_array(raw)
  parameter_summary <- as.data.frame(posterior::summarise_draws(draws, posterior::rhat, posterior::ess_bulk, posterior::ess_tail))
  names(parameter_summary) <- sub("^posterior::", "", names(parameter_summary))
  names(parameter_summary)[names(parameter_summary) == ".variable"] <- "parameter"
  bfmi <- rstan::get_bfmi(fit)

  hmc <- tibble(
    metric = c("divergences", "max_treedepth_hits", "maximum_Rhat", "n_Rhat_gt_1.01", "minimum_bulk_ESS", "minimum_tail_ESS", "elapsed_seconds", paste0("BFMI_chain_", seq_along(bfmi))),
    value = c(
      rstan::get_num_divergent(fit), rstan::get_num_max_treedepth(fit),
      max(parameter_summary$rhat, na.rm = TRUE), sum(parameter_summary$rhat > 1.01, na.rm = TRUE),
      min(parameter_summary$ess_bulk, na.rm = TRUE), min(parameter_summary$ess_tail, na.rm = TRUE),
      bundle$config$elapsed_seconds, bfmi
    )
  )
  hmc_pass <- hmc$value[hmc$metric == "divergences"] == 0 &&
    hmc$value[hmc$metric == "max_treedepth_hits"] == 0 &&
    hmc$value[hmc$metric == "maximum_Rhat"] <= 1.01 &&
    hmc$value[hmc$metric == "minimum_bulk_ESS"] >= 100 && all(bfmi >= .3)

  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  treedepth_all <- unlist(lapply(sampler_params, function(x) x[, "treedepth__"]))
  leapfrog_all <- unlist(lapply(sampler_params, function(x) x[, "n_leapfrog__"]))
  chain_diag <- tibble(
    chain = seq_along(sampler_params),
    treedepth_median = vapply(sampler_params, function(x) median(x[, "treedepth__"]), numeric(1)),
    treedepth_hits = vapply(sampler_params, function(x) sum(x[, "treedepth__"] >= bundle$config$max_treedepth), integer(1)),
    n_leapfrog_median = vapply(sampler_params, function(x) median(x[, "n_leapfrog__"]), numeric(1)),
    n_leapfrog_p95 = vapply(sampler_params, function(x) quantile(x[, "n_leapfrog__"], .95, names = FALSE), numeric(1)),
    mean_accept_stat = vapply(sampler_params, function(x) mean(x[, "accept_stat__"]), numeric(1)),
    n_divergent = vapply(sampler_params, function(x) sum(x[, "divergent__"]), numeric(1)),
    bfmi = bfmi
  )
  treedepth_quantiles <- tibble(quantile = c("50%", "90%", "95%", "99%", "max"), treedepth = c(quantile(treedepth_all, c(.5, .9, .95, .99)), max(treedepth_all)))

  write_csv(hmc, file.path(out_dir, paste0("hmc_diagnostics_stage", stage, ".csv")))
  write_csv(arrange(parameter_summary, desc(rhat)), file.path(out_dir, paste0("parameter_diagnostics_stage", stage, ".csv")))
  write_csv(chain_diag, file.path(out_dir, paste0("chain_specific_diagnostics_stage", stage, ".csv")))
  write_csv(treedepth_quantiles, file.path(out_dir, paste0("treedepth_quantiles_stage", stage, ".csv")))

  message(sprintf("[stage%s] STRICT HMC GATE: %s", stage, if (hmc_pass) "PASS" else "FAIL"))
  message(sprintf("  divergences=%d | max_treedepth_hits=%d | max_Rhat=%.5f | min_bulk_ESS=%.1f | min_tail_ESS=%.1f | BFMI range=[%.3f, %.3f] | median n_leapfrog=%.0f",
                   hmc$value[hmc$metric == "divergences"], hmc$value[hmc$metric == "max_treedepth_hits"],
                   hmc$value[hmc$metric == "maximum_Rhat"], hmc$value[hmc$metric == "minimum_bulk_ESS"], hmc$value[hmc$metric == "minimum_tail_ESS"],
                   min(bfmi), max(bfmi), median(leapfrog_all)))
  print(as.data.frame(chain_diag))
  if (!hmc_pass) message("  *** GATE FAILED at stage ", stage, ". Per Sections 7-8: do NOT proceed automatically, do NOT loosen the gate. Report before continuing. ***")
  invisible(list(hmc = hmc, hmc_pass = hmc_pass, chain_diag = chain_diag, treedepth_quantiles = treedepth_quantiles))
}

if (sys.nframe() == 0L) run_check_hmc_gate()
