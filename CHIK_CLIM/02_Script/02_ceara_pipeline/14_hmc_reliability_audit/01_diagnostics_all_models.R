# STEP 1-3 of the HMC reliability audit:
#   1. Per-model, per-chain divergences/divergence rate, Rhat (classic +
#      rank-normalised), bulk/tail ESS, E-BFMI, treedepth-saturation hits.
#   2. Diagnostic pairs plots (tau_year x alpha_R/sigma_season_year/beta_T/beta_P),
#      divergent transitions highlighted, for every model.
#   3. Whether divergent iterations concentrate at small tau_year or
#      correlate with specific other parameters; whether excluding
#      divergent draws would shift the reported tau_year summary
#      (a bias check, NOT a correction -- divergent draws are reported,
#      never silently dropped from the headline posterior).
#
# Read-only against all 5 existing fits. No new Stan fitting here.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "posterior", "bayesplot", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(posterior); library(bayesplot); library(ggplot2) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/14_hmc_reliability_audit/00_config_hmc_audit.R"))

OKABE <- c(blue = "#0072B2", orange = "#D55E00", yellow = "#E69F00", green = "#009E73", purple = "#CC79A7", grey = "#56B4E9")

# ================================================================
# Load all 5 fits once.
# ================================================================
bundles <- lapply(MODELS_TO_AUDIT, function(m) readRDS(m$path))
names(bundles) <- names(MODELS_TO_AUDIT)

# ================================================================
# 1. Per-model, per-chain diagnostic table.
# ================================================================
diagnose_one <- function(model_id, bundle) {
  fit <- bundle$fit
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  n_chains <- length(sp)
  n_iter_per_chain <- nrow(sp[[1]])
  max_td <- bundle$config$max_treedepth
  if (is.null(max_td)) max_td <- 14L

  bfmi_by_chain <- tryCatch(rstan::get_bfmi(fit), error = function(e) rep(NA_real_, n_chains))

  draws_arr <- as_draws_array(fit, pars = GATE_PARS[GATE_PARS %in% names(fit@par_dims)])
  smry <- summarise_draws(draws_arr, "rhat", "ess_bulk", "ess_tail")

  per_chain <- bind_rows(lapply(seq_len(n_chains), function(c) {
    tibble(model = model_id, chain = c,
           n_iter = n_iter_per_chain,
           divergences = sum(sp[[c]][, "divergent__"]),
           divergence_rate = mean(sp[[c]][, "divergent__"]),
           max_treedepth_hits = sum(sp[[c]][, "treedepth__"] >= max_td),
           mean_treedepth = mean(sp[[c]][, "treedepth__"]),
           max_treedepth_reached = max(sp[[c]][, "treedepth__"]),
           ebfmi = bfmi_by_chain[c])
  }))

  model_level <- tibble(model = model_id, label = MODELS_TO_AUDIT[[model_id]]$label,
                         total_post_warmup_iter = n_chains * n_iter_per_chain,
                         total_divergences = sum(per_chain$divergences),
                         overall_divergence_rate = mean(per_chain$divergences) / n_iter_per_chain,
                         total_max_treedepth_hits = sum(per_chain$max_treedepth_hits),
                         min_ebfmi = min(bfmi_by_chain, na.rm = TRUE),
                         max_rhat = max(smry$rhat, na.rm = TRUE), worst_rhat_par = smry$variable[which.max(smry$rhat)],
                         min_ess_bulk = min(smry$ess_bulk, na.rm = TRUE), worst_ess_bulk_par = smry$variable[which.min(smry$ess_bulk)],
                         min_ess_tail = min(smry$ess_tail, na.rm = TRUE), worst_ess_tail_par = smry$variable[which.min(smry$ess_tail)])

  list(per_chain = per_chain, model_level = model_level, par_summary = smry |> mutate(model = model_id))
}

all_diag <- lapply(names(bundles), function(m) diagnose_one(m, bundles[[m]]))
names(all_diag) <- names(bundles)

per_chain_all <- bind_rows(lapply(all_diag, `[[`, "per_chain"))
model_level_all <- bind_rows(lapply(all_diag, `[[`, "model_level"))
par_summary_all <- bind_rows(lapply(all_diag, `[[`, "par_summary"))

write_csv(per_chain_all, file.path(DIR_TABLES_HA, "HMC_AUDIT_per_chain.csv"))
write_csv(model_level_all, file.path(DIR_TABLES_HA, "HMC_AUDIT_model_level.csv"))
write_csv(par_summary_all, file.path(DIR_TABLES_HA, "HMC_AUDIT_per_parameter_rhat_ess.csv"))

message("\n=== HMC audit: model-level summary ===")
print(as.data.frame(model_level_all[, c("model","total_post_warmup_iter","total_divergences","overall_divergence_rate",
                                          "total_max_treedepth_hits","min_ebfmi","max_rhat","min_ess_bulk","min_ess_tail")]), digits = 4)
message("\n=== HMC audit: per-chain detail ===")
print(as.data.frame(per_chain_all), digits = 4)

# ================================================================
# 2. Diagnostic pairs plots -- divergences highlighted.
# ================================================================
for (m in names(bundles)) {
  fit <- bundles[[m]]$fit
  pars_present <- PAIRS_PARS[PAIRS_PARS %in% names(fit@par_dims)]
  p <- bayesplot::mcmc_pairs(as.array(fit), pars = pars_present, np = bayesplot::nuts_params(fit),
                               off_diag_args = list(size = 0.5, alpha = 0.3))
  ggsave(file.path(DIR_FIGURES_HA, sprintf("HMC_AUDIT_pairs_%s.png", m)), p, width = 240, height = 240, units = "mm", dpi = 300, bg = "white")
  message(sprintf("[saved] HMC_AUDIT_pairs_%s.png", m))
}

# ================================================================
# 3a. Do divergences concentrate at small tau_year, or correlate with
# other parameters? Extract (permuted=FALSE) tau_year + candidates,
# aligned with sampler_params' per-iteration divergent indicator.
# ================================================================
div_conc <- bind_rows(lapply(names(bundles), function(m) {
  fit <- bundles[[m]]$fit
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  div_vec <- unlist(lapply(sp, function(x) x[, "divergent__"]))
  ext <- rstan::extract(fit, pars = c("tau_year", "alpha_R", "sigma_season_year", "beta_T", "beta_P"), permuted = FALSE)
  # extract(permuted=FALSE) dims: [iteration, chain, parameter] -- flatten chain-major to match sp's chain-then-iteration order
  n_iter <- dim(ext)[1]; n_chain <- dim(ext)[2]
  flat <- function(par) as.vector(ext[, , par]) # column-major: iteration varies fastest within chain -- matches unlist(lapply(sp,...)) chain-major order only if we reshape consistently
  # Rebuild explicitly chain-by-chain to guarantee alignment (avoid silent transpose bugs):
  tau_vec <- unlist(lapply(seq_len(n_chain), function(c) ext[, c, "tau_year"]))
  alpha_vec <- unlist(lapply(seq_len(n_chain), function(c) ext[, c, "alpha_R"]))
  sigma_vec <- unlist(lapply(seq_len(n_chain), function(c) ext[, c, "sigma_season_year"]))
  betaT_vec <- unlist(lapply(seq_len(n_chain), function(c) ext[, c, "beta_T"]))
  betaP_vec <- unlist(lapply(seq_len(n_chain), function(c) ext[, c, "beta_P"]))
  stopifnot(length(div_vec) == length(tau_vec))
  tibble(model = m, divergent = div_vec, tau_year = tau_vec, alpha_R = alpha_vec, sigma_season_year = sigma_vec, beta_T = betaT_vec, beta_P = betaP_vec)
}))
write_csv(div_conc |> group_by(model) |> summarise(
    n_total = n(), n_divergent = sum(divergent == 1),
    median_tau_divergent = median(tau_year[divergent == 1]), median_tau_nondivergent = median(tau_year[divergent == 0]),
    p10_tau_overall = quantile(tau_year, .10), p90_tau_overall = quantile(tau_year, .90),
    median_alphaR_divergent = median(alpha_R[divergent == 1]), median_alphaR_nondivergent = median(alpha_R[divergent == 0]),
    median_sigma_divergent = median(sigma_season_year[divergent == 1]), median_sigma_nondivergent = median(sigma_season_year[divergent == 0]),
    .groups = "drop"),
  file.path(DIR_TABLES_HA, "HMC_AUDIT_divergence_concentration.csv"))

message("\n=== Where do divergences concentrate? (median parameter value, divergent vs non-divergent iterations) ===")
print(as.data.frame(read_csv(file.path(DIR_TABLES_HA, "HMC_AUDIT_divergence_concentration.csv"), show_col_types = FALSE)), digits = 4)

p_tau_div <- ggplot(div_conc, aes(x = tau_year, fill = factor(divergent))) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, alpha = 0.6, position = "identity") +
  facet_wrap(~model, scales = "free", ncol = 2) +
  scale_fill_manual(name = NULL, values = c("0" = OKABE["blue"], "1" = OKABE["orange"]), labels = c("non-divergent", "divergent")) +
  labs(title = "Where do divergent transitions fall in tau_year?", x = "tau_year", y = "Density") +
  theme_classic(base_size = 9) + theme(legend.position = "top")
ggsave(file.path(DIR_FIGURES_HA, "HMC_AUDIT_divergence_vs_tau_year.png"), p_tau_div, width = 220, height = 180, units = "mm", dpi = 300, bg = "white")

# ================================================================
# 3b. Bias check: does excluding divergent draws shift the headline
# tau_year posterior summary? (diagnostic only -- divergent draws are
# NEVER silently dropped from any reported headline posterior)
# ================================================================
bias_check <- div_conc |> group_by(model) |> summarise(
  median_tau_all_draws = median(tau_year),
  median_tau_excl_divergent = median(tau_year[divergent == 0]),
  mean_tau_all_draws = mean(tau_year),
  mean_tau_excl_divergent = mean(tau_year[divergent == 0]),
  pct_shift_median = 100 * (median(tau_year[divergent == 0]) - median(tau_year)) / median(tau_year),
  .groups = "drop")
write_csv(bias_check, file.path(DIR_TABLES_HA, "HMC_AUDIT_bias_check_exclude_divergent.csv"))
message("\n=== Bias check: tau_year summary, all draws vs excluding divergent iterations ===")
print(as.data.frame(bias_check), digits = 4)

message("\n[HMC audit Step 1-3] Complete -> ", DIR_TABLES_HA, " / ", DIR_FIGURES_HA)
