# Addendum to Step 3: divergences concentrate at LARGER tau_year
# (01_diagnostics_all_models.R finding) -- check whether the UPPER TAIL
# of the reported tau_year posterior (not just the median) is sensitive
# to excluding divergent iterations. Cheap, read-only, no new fitting.

suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble) })
source(file.path(ROOT, "02_Script/02_ceara_pipeline/14_hmc_reliability_audit/00_config_hmc_audit.R"))

bundles <- lapply(MODELS_TO_AUDIT, function(m) readRDS(m$path))
names(bundles) <- names(MODELS_TO_AUDIT)

tail_check <- bind_rows(lapply(names(bundles), function(m) {
  fit <- bundles[[m]]$fit
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  div_vec <- unlist(lapply(sp, function(x) x[, "divergent__"]))
  ext <- rstan::extract(fit, pars = "tau_year", permuted = FALSE)
  n_chain <- dim(ext)[2]
  tau_vec <- unlist(lapply(seq_len(n_chain), function(c) ext[, c, "tau_year"]))
  stopifnot(length(div_vec) == length(tau_vec))
  tau_all <- tau_vec; tau_excl <- tau_vec[div_vec == 0]
  tibble(model = m, n_divergent = sum(div_vec),
         q50_all = median(tau_all), q50_excl = median(tau_excl),
         q90_all = quantile(tau_all, .90), q90_excl = quantile(tau_excl, .90),
         q95_all = quantile(tau_all, .95), q95_excl = quantile(tau_excl, .95),
         q975_all = quantile(tau_all, .975), q975_excl = quantile(tau_excl, .975),
         pct_shift_q95 = 100 * (quantile(tau_excl, .95) - quantile(tau_all, .95)) / quantile(tau_all, .95),
         pct_shift_q975 = 100 * (quantile(tau_excl, .975) - quantile(tau_all, .975)) / quantile(tau_all, .975))
}))
write_csv(tail_check, file.path(DIR_TABLES_HA, "HMC_AUDIT_tail_bias_check.csv"))
message("\n=== Upper-tail sensitivity: tau_year quantiles, all draws vs excluding divergent iterations ===")
print(as.data.frame(tail_check), digits = 4)
