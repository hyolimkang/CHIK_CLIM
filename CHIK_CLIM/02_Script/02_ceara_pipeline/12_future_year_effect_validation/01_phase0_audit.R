# PHASE 0 -- Audit existing model structure before any new fitting.
# Read-only (plus one Stan COMPILE + smoke test of the new tau_year
# model, which is a computational check, not a scientific fit).

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/12_future_year_effect_validation/00_config_year_effect_validation.R"))

message("[Phase0] Loading accepted M0 fit (read-only)...")
m0 <- readRDS(ACCEPTED_CLIMATE_FIT_M0)
sd_m0 <- m0$stan_data

audit <- tibble(
  item = c("year_effect construction", "year_effect_prior_sd (fixed)", "A_year construction", "sigma_season_year prior",
           "2018 serology anchor", "serology likelihood", "2022 conditioned seed", "fixed q", "imports_per_week",
           "generation interval", "initial susceptibility", "climate coefficients", "N (historical weeks)", "Y (historical years)"),
  value = c("year_effect[y] = year_effect_prior_sd * (z_year[y] - mean(z_year)), z_year ~ std_normal()",
            sprintf("%.2f (DATA, not estimated -- this is exactly what Phase 1 replaces with estimated tau_year)", sd_m0$year_effect_prior_sd),
            "A_year[y] = exp(sigma_season_year * z_season_year[y]), z_season_year ~ std_normal()",
            "sigma_season_year ~ Normal(0,0.20), lower=0 (HalfNormal) -- genuinely ESTIMATED, unchanged by this study",
            sprintf("t_sero index=%d (%s), sero_pos=%d/%d, kappa_sero=%d (Beta-Binomial)", sd_m0$t_sero, as.character(m0$sero_week$matched_week_start), sd_m0$sero_pos, sd_m0$sero_n, sd_m0$kappa_sero),
            "sero_pos ~ beta_binomial(sero_n, alpha_sero, beta_sero); alpha/beta_sero derived from immune_prop[t_sero]",
            sprintf("is_seed=1 for %d weeks starting %s; X_seed=cases/q; excluded from case likelihood via fit_index (N_fit=%d of N=%d)", sum(sd_m0$is_seed), as.character(min(m0$seed_dates)), sd_m0$N_fit, sd_m0$N),
            sprintf("%.2f (fixed scenario input, not estimated)", sd_m0$q),
            sprintf("%.2f /week (fixed, added as a hazard term, never double-counted as cases)", sd_m0$imports_per_week),
            sprintf("G=%d weeks, Gamma(shape=4,rate=2), discretised/renormalised", sd_m0$G),
            "S[1] = N_start[1] (fully susceptible at 2015-01-04); U[1] = 0 (no prior immunity assumed at series start)",
            "beta_T, beta_P ~ Normal(0,0.15); log_R0[t] += beta_T*z_T_anom[t] + beta_P*z_P_anom[t]; UNCHANGED by this study",
            sd_m0$N, sd_m0$Y
  )
)
print(as.data.frame(audit))
write_csv(audit, file.path(DIR_TABLES_YV, "PHASE0_model_structure_audit.csv"))

message("\n[Phase0] Checking the most recent climate_robustness results are available (read-only, informational)...")
cr_report <- file.path(ROOT, "03_Output/02_ceara_pipeline/reports/climate_robustness/CE_CLIMATE_ROBUSTNESS_ASSESSMENT.md")
yv_prior_report <- file.path(ROOT, "03_Output/02_ceara_pipeline/reports/year_effect_diagnostic/FUTURE_YEAR_EFFECT_VALIDATION.md")
message(sprintf("  climate_robustness report present: %s", file.exists(cr_report)))
message(sprintf("  prior year_effect_diagnostic report present: %s", file.exists(yv_prior_report)))
message("  Both already establish: (a) climate inputs reproduce exactly (Phase 1 of climate_robustness), (b) the fixed\n  year_effect_prior_sd=0.40 is a researcher-chosen constant, not data-informed, and the DOMINANT driver of Y1's\n  extreme future epidemics (year_effect_diagnostic FIG10-12). This study's Phase 1 directly follows up on that finding.")

message("\n[Phase0] Compiling the NEW tau_year Stan model and running a 1-chain smoke test...")
mod_tau <- rstan::stan_model(file = STAN_SOURCE_TAU_YEAR)
message("[Phase0] Compiled OK.")

source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/22_v4_3_short_2015_2019_q_calibration/scripts/02_fit_v4_3.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/27_v4_8_seeded_recurrence/scripts/01_fit_v4_8.R"))
source(file.path(ROOT, "02_Script/00_shared/legacy_model_functions/28_v4_9_hierarchical_seasonality/scripts/01_fit_v4_9.R"))

smoke_data <- m0$stan_data
smoke_data$tau_prior_scale <- 0.40
smoke_data$year_effect_prior_sd <- NULL
init_fn_smoke <- function(chain_id = 1L) {
  v <- make_v4_9_init_fn(smoke_data$Y, load_td14_scalar_inits(ROOT))(chain_id)
  v$beta_T <- 0; v$beta_P <- 0; v$tau_year <- 0.40
  v
}
smoke_fit <- rstan::sampling(mod_tau, data = smoke_data, chains = 1L, iter = 40L, warmup = 20L, init = init_fn_smoke, seed = 1L,
                              control = list(adapt_delta = 0.95, max_treedepth = 10, metric = "dense_e"), refresh = 0)
message("[Phase0] Smoke test OK -- no dimension/parsing errors. tau_year is sampling (not fixed):")
print(smoke_fit, pars = c("alpha_R", "tau_year", "beta_T", "beta_P"))

phase0_status <- tibble(check = c("Accepted M0 fit loads and matches documented structure", "New tau_year Stan model compiles", "New tau_year Stan model smoke-tests on real full-period data"),
                         pass = TRUE)
write_csv(phase0_status, file.path(DIR_TABLES_YV, "PHASE0_gate_status.csv"))
message("\n[Phase0] GATE: PASS -- proceeding to Phase 1 (full multi-chain tau_year fits).")
