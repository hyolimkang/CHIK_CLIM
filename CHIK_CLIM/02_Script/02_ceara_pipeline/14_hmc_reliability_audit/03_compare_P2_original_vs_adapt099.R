# STEP 4 (cont'd) -- compare original P2 (adapt_delta=0.95, warmup=800,
# sampling=400) vs the adapt_delta=0.99/warmup=1500/sampling=500 rerun,
# parameter by parameter. Both fits preserved untouched; this script is
# read-only against both.

suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(posterior) })
source(file.path(ROOT, "02_Script/02_ceara_pipeline/14_hmc_reliability_audit/00_config_hmc_audit.R"))

orig <- readRDS(file.path(DIR_FITS_YV, "ce_tau_year_P2.rds"))
rerun <- readRDS(file.path(DIR_FITS_HA, "ce_tau_year_P2_adapt099.rds"))

compare_one <- function(fit, label) {
  draws_arr <- as_draws_array(fit, pars = GATE_PARS[GATE_PARS %in% names(fit@par_dims)])
  smry <- summarise_draws(draws_arr, "mean", "median", "sd", "rhat", "ess_bulk", "ess_tail", quantile2) |> mutate(model = label)
  smry
}
cmp <- bind_rows(compare_one(orig$fit, "P2_original_adapt095"), compare_one(rerun$fit, "P2_rerun_adapt099"))
write_csv(cmp, file.path(DIR_TABLES_HA, "HMC_AUDIT_P2_original_vs_adapt099_per_parameter.csv"))
message("=== Per-parameter comparison: P2 original (adapt_delta=0.95) vs rerun (adapt_delta=0.99) ===")
print(as.data.frame(cmp |> dplyr::select(model, variable, mean, median, sd, rhat, ess_bulk, ess_tail) |> dplyr::arrange(variable, model)), digits = 4)

tau_cmp <- cmp |> dplyr::filter(variable == "tau_year") |> dplyr::select(model, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail) |> as.data.frame()
message("\n=== tau_year specifically ===")
print(as.data.frame(tau_cmp), digits = 4)

sp <- rstan::get_sampler_params(rerun$fit, inc_warmup = FALSE)
per_chain_rerun <- bind_rows(lapply(seq_along(sp), function(c) tibble(chain = c, divergences = sum(sp[[c]][, "divergent__"]),
                                                                        max_treedepth_hits = sum(sp[[c]][, "treedepth__"] >= 14),
                                                                        mean_treedepth = mean(sp[[c]][, "treedepth__"]))))
write_csv(per_chain_rerun, file.path(DIR_TABLES_HA, "HMC_AUDIT_P2_adapt099_per_chain.csv"))
message("\n=== P2 rerun (adapt_delta=0.99) per-chain ==="); print(as.data.frame(per_chain_rerun), digits = 4)

message("\nOriginal P2 elapsed: ", orig$config$elapsed_seconds, "s | Rerun P2 elapsed: ", rerun$config$elapsed_seconds, "s")
message("Original P2 HMC gate pass: ", orig$hmc$hmc_pass, " | Rerun HMC gate pass: ", rerun$hmc$hmc_pass)
