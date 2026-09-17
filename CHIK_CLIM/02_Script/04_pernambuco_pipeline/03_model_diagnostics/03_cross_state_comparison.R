# Pernambuco external replication -- Section 11: descriptive cross-state
# table (Ceara, Bahia, Pernambuco). No pooling, no hierarchical model.
#
# IMPORTANT ASYMMETRY: Ceara's frozen v4.9 uses FIXED q (data, not
# estimated) -- there is no Ceara "global-q posterior" to report on the
# same footing as Bahia/PE. This is stated explicitly rather than silently
# treating a fixed value as if it were an estimated posterior.

required_packages <- c("rstan", "dplyr", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(readr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/tables/pernambuco_v4_9_replication")

ce <- readRDS(file.path(root, "03_Output/02_ceara_pipeline/model_fits/baseline/v4_9_hierarchical_seasonality/outputs/q0.05/renewal_ceara_v4_9_fit_q0.05.rds"))
ba <- readRDS(file.path(root, "03_Output/03_bahia_pipeline/model_fits/v4_9_replication/outputs_global_q_multisite_serology/on_full2/renewal_bahia_global_q_fit_on_full2.rds"))
pe <- readRDS(file.path(root, "03_Output/04_pernambuco_pipeline/model_fits/v4_9_replication/outputs/modelB_U14/renewal_pe_global_q_fit_modelB_U14.rds"))

immune_at <- function(bundle, year) {
  dates <- as.Date(bundle$weekly_data$week_start)
  idx <- if (year == "2025") nrow(bundle$weekly_data) else max(which(format(dates, "%Y") == year))
  ip <- rstan::extract(bundle$fit, "immune_prop")$immune_prop[, idx]
  sprintf("%.2f%% [%.2f,%.2f]", 100*median(ip), 100*quantile(ip,.025), 100*quantile(ip,.975))
}

eta_geo_summary <- function(bundle) {
  eg <- tryCatch(rstan::extract(bundle$fit, "eta_geo")$eta_geo, error = function(e) NULL)
  if (is.null(eg)) return("n/a")
  sprintf("median |eta_geo| = %.2f (range %.2f-%.2f)", median(abs(apply(eg, 2, median))), min(apply(eg,2,median)), max(apply(eg,2,median)))
}

comparison <- tibble(
  state = c("Ceara", "Bahia", "Pernambuco"),
  epidemic_regime = c("Large early epidemic (2016-17), substantial apparent depletion, delayed 2022 recurrence (required conditioned seed)",
                       "Recurrent state-level waves nearly every year, limited apparent statewide depletion, strong documented sub-regional spatial turnover",
                       "Intermediate: one large early epidemic (2015-17, ~171/100k) followed by recurrent annual/biennial waves without a multi-year gap"),
  global_q_posterior = c("q FIXED at 0.05 (reference; not estimated -- NOT directly comparable to Bahia/PE)",
                          sprintf("median=%.3f [%.3f,%.3f] (95%% CrI)", median(as.vector(rstan::extract(ba$fit,"q")$q)),
                                  quantile(as.vector(rstan::extract(ba$fit,"q")$q),.025), quantile(as.vector(rstan::extract(ba$fit,"q")$q),.975)),
                          sprintf("median=%.3f [%.3f,%.3f] (95%% CrI)", median(as.vector(rstan::extract(pe$fit,"q")$q)),
                                  quantile(as.vector(rstan::extract(pe$fit,"q")$q),.025), quantile(as.vector(rstan::extract(pe$fit,"q")$q),.975))),
  immune_fraction_2018 = c(immune_at(ce, "2018"), immune_at(ba, "2018"), immune_at(pe, "2018")),
  immune_fraction_2022 = c(immune_at(ce, "2022"), immune_at(ba, "2022"), immune_at(pe, "2022")),
  immune_fraction_2025 = c(immune_at(ce, "2025"), immune_at(ba, "2025"), immune_at(pe, "2025")),
  serology_observations_used = c("1 (Juazeiro anchor, single-site)", "6 (multisite, geo-adjusted)", "1 (U14, Recife, geo-adjusted; U07/U08 reserved)"),
  geographic_offset_magnitude = c("n/a (single-site anchor, no geographic offset)", eta_geo_summary(ba), eta_geo_summary(pe)),
  hmc_diagnostics = c(sprintf("divergences=%d, treedepth=%d, max_rhat=%.3f, PASS=%s", ce$hmc$divergences, ce$hmc$max_treedepth_hits, ce$hmc$maximum_rhat, ce$hmc$hmc_pass),
                       sprintf("divergences=%d, treedepth=%d, max_rhat=%.3f, PASS=%s", ba$hmc$divergences, ba$hmc$max_treedepth_hits, ba$hmc$maximum_rhat, ba$hmc$hmc_pass),
                       sprintf("divergences=%d, treedepth=%d, max_rhat=%.3f, PASS=%s", pe$hmc$divergences, pe$hmc$max_treedepth_hits, pe$hmc$maximum_rhat, pe$hmc$hmc_pass))
)

message("=== CE / BA / PE descriptive comparison (Section 11) ===")
print(as.data.frame(comparison), digits = 3)
write_csv(comparison, file.path(table_dir, "CE_BA_PE_comparison.csv"))
message("\n[saved] ", file.path(table_dir, "CE_BA_PE_comparison.csv"))
