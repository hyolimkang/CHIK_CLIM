# Mode audit Part 2 (corrected): bridge sampling on the two EXPLICITLY
# truncated partition models (q<0.25, q>=0.25). Because the truncation is
# implemented via parameter bounds (not a model/prior change -- see the
# .stan file headers), log_prob() at any point equals the ORIGINAL joint
# log-density there, so bridge_sampler gives the correctly normalized
# regional integral of the ORIGINAL target directly -- no separate
# normalization correction needed.

required_packages <- c("here", "rstan", "bridgesampling", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(bridgesampling); library(tibble); library(readr) })

bridge_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

reattach_live_log_prob <- function(bundle, model) {
  throwaway <- rstan::sampling(model, data = bundle$stan_data, chains = 1, iter = 2, warmup = 1, refresh = 0)
  bundle$fit@stanmodel <- model
  bundle$fit@.MISC <- throwaway@.MISC
  bundle$fit
}

run_partition_bridge_sampling <- function() {
  root <- bridge_root()
  base_dir <- file.path(root, "03_Output", "model_fits", "ceara", "v4_3_short_2015_2019_q_calibration")

  low_bundle <- readRDS(file.path(base_dir, "outputs/partition_low/renewal_ceara_v4_3_fit_partition_low.rds"))
  high_bundle <- readRDS(file.path(base_dir, "outputs/partition_high/renewal_ceara_v4_3_fit_partition_high.rds"))

  message("[partition bridge] recompiling both truncated stanmodels...")
  low_model <- rstan::stan_model(file.path(root, "02_Script/stan/renewal_ceara_v4_3_weak_q_truncated.stan"))
  high_model <- rstan::stan_model(file.path(root, "02_Script/stan/renewal_ceara_v4_3_weak_q_truncated_high.stan"))
  low_bundle$fit <- reattach_live_log_prob(low_bundle, low_model)
  high_bundle$fit <- reattach_live_log_prob(high_bundle, high_model)

  message("[partition bridge] running bridge sampling on q<0.25 partition...")
  bridge_low <- bridgesampling::bridge_sampler(low_bundle$fit, silent = TRUE)
  message("[partition bridge] running bridge sampling on q>=0.25 partition...")
  bridge_high <- bridgesampling::bridge_sampler(high_bundle$fit, silent = TRUE)
  bf_result <- bridgesampling::bf(bridge_low, bridge_high)

  result <- tibble(
    region = c("q<0.25 (low)", "q>=0.25 (high)"),
    log_marginal_likelihood = c(bridge_low$logml, bridge_high$logml)
  )
  result$relative_posterior_mass <- exp(result$log_marginal_likelihood - max(result$log_marginal_likelihood))
  result$relative_posterior_mass <- result$relative_posterior_mass / sum(result$relative_posterior_mass)
  write_csv(result, file.path(base_dir, "mode_audit_partition_bridge_sampling_CORRECTED.csv"))

  message("[partition bridge] CORRECTED regional posterior mass (explicit q<0.25 vs q>=0.25 partition):")
  print(as.data.frame(result))
  message(sprintf("[partition bridge] Bayes factor (low region vs high region) = %.4g", as.numeric(bf_result$bf)))
  message("[partition bridge] This IS the rigorous Part-2 answer -- unlike the earlier basin-trapped-chains estimate, this correctly integrates the ORIGINAL joint target over each explicit region.")
  invisible(result)
}

if (sys.nframe() == 0L) run_partition_bridge_sampling()
