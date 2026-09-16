# Diagnostics for the Pernambuco v4.9 global-q ridge-reparameterised model.

required_packages <- c("rstan", "posterior", "dplyr", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(posterior); library(dplyr); library(tibble); library(readr) })

pe_project_root <- function() {
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(current, "CHIK_CLIM.Rproj"))) return(current)
  candidate <- file.path(current, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Run from the outer repository or inner CHIK_CLIM project root.")
}

root <- pe_project_root()
tag <- Sys.getenv("TAG", "modelC_U14")
stage <- Sys.getenv("STAGE", "pilot")
out_dir <- file.path(root, "02_Script", "40_renewal_model", "33_pernambuco_v4_9_replication", "outputs", "modelC_global_q_ridge_reparam")
table_dir <- file.path(root, "03_Output", "tables", "pernambuco_v4_9_ridge_reparam")
fit_path <- file.path(out_dir, paste0("renewal_pe_v4_9_global_q_ridge_reparam_", tag, "_", stage, ".rds"))
if (!file.exists(fit_path)) stop("Fit not found: ", fit_path)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(fit_path)
fit <- bundle$fit
raw <- rstan::extract(fit, pars = c("gamma_ridge", "alpha_R", "z_year_contrast", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2", "sigma_season_year", "z_season_year", "phi_obs", "logit_q"), permuted = FALSE, inc_warmup = FALSE)
parameter_summary <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(raw), posterior::rhat, posterior::ess_bulk, posterior::ess_tail))
names(parameter_summary)[names(parameter_summary) == ".variable"] <- "parameter"

q <- as.vector(rstan::extract(fit, "q")$q)
alpha <- as.vector(rstan::extract(fit, "alpha_R")$alpha_R)
gamma <- as.vector(rstan::extract(fit, "gamma_ridge")$gamma_ridge)
coordinate_diagnostics <- tibble(
  quantity = c("cor_logitq_alpha_R", "cor_logitq_gamma_ridge", "q_median", "q_lo95", "q_hi95", "alpha_R_median"),
  value = c(cor(qlogis(q), alpha), cor(qlogis(q), gamma), median(q), quantile(q, .025), quantile(q, .975), median(alpha))
)
bfmi <- rstan::get_bfmi(fit)
hmc <- tibble(
  metric = c("divergences", "max_treedepth_hits", "maximum_Rhat", "minimum_bulk_ESS", "minimum_tail_ESS", paste0("BFMI_chain_", seq_along(bfmi))),
  value = c(rstan::get_num_divergent(fit), rstan::get_num_max_treedepth(fit), max(parameter_summary$rhat), min(parameter_summary$ess_bulk), min(parameter_summary$ess_tail), bfmi)
)
write_csv(parameter_summary, file.path(table_dir, paste0("parameter_diagnostics_", tag, "_", stage, ".csv")))
write_csv(coordinate_diagnostics, file.path(table_dir, paste0("ridge_coordinate_diagnostics_", tag, "_", stage, ".csv")))
write_csv(hmc, file.path(table_dir, paste0("hmc_diagnostics_", tag, "_", stage, ".csv")))
message("[", tag, "] stored HMC gate: ", if (bundle$hmc$hmc_pass) "PASS" else "FAIL",
  " | cor(logit q, alpha_R)=", sprintf("%.3f", coordinate_diagnostics$value[1]),
  " | cor(logit q, gamma_ridge)=", sprintf("%.3f", coordinate_diagnostics$value[2]))
