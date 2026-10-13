# HMC reliability audit for the tau_year hierarchical fits (P1-P3,
# full period; H1/H2, training-only) -- requested as a PREREQUISITE
# gate before continuing the retrospective resurgence diagnostic
# (13_retrospective_resurgence_diagnostic). Preserves all existing
# fits; does not overwrite or commit anything accepted.

if (!exists("ROOT")) stop("ROOT not found -- source CHIK_CLIM/00_project_setup.R first.")

DIR_FITS_YV <- file.path(ROOT, "03_Output/02_ceara_pipeline/model_fits/future_year_effect_validation")
for (nm in c("ce_tau_year_P1.rds", "ce_tau_year_P2.rds", "ce_tau_year_P3.rds", "ce_tau_year_training_H1.rds", "ce_tau_year_training_H2.rds")) {
  f <- file.path(DIR_FITS_YV, nm)
  if (!file.exists(f)) stop("Missing required fit, cannot proceed: ", f)
}

DIR_SCRIPTS_HA <- file.path(ROOT, "02_Script/02_ceara_pipeline/14_hmc_reliability_audit")
DIR_REPORTS_HA <- file.path(ROOT, "03_Output/02_ceara_pipeline/reports/hmc_reliability_audit")
DIR_TABLES_HA  <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/hmc_reliability_audit")
DIR_FIGURES_HA <- file.path(ROOT, "03_Output/02_ceara_pipeline/figures/hmc_reliability_audit")
DIR_FITS_HA    <- file.path(ROOT, "03_Output/02_ceara_pipeline/model_fits/hmc_reliability_audit") # NEW fits only (e.g. P2 rerun) -- never overwrites DIR_FITS_YV
for (d in c(DIR_REPORTS_HA, DIR_TABLES_HA, DIR_FIGURES_HA, DIR_FITS_HA)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

MODELS_TO_AUDIT <- list(
  P1 = list(path = file.path(DIR_FITS_YV, "ce_tau_year_P1.rds"), label = "P1: full period, HalfNormal(0,0.20)"),
  P2 = list(path = file.path(DIR_FITS_YV, "ce_tau_year_P2.rds"), label = "P2: full period, HalfNormal(0,0.40) [PRIMARY]"),
  P3 = list(path = file.path(DIR_FITS_YV, "ce_tau_year_P3.rds"), label = "P3: full period, HalfNormal(0,0.80)"),
  H1 = list(path = file.path(DIR_FITS_YV, "ce_tau_year_training_H1.rds"), label = "H1: training-only 2015-2019"),
  H2 = list(path = file.path(DIR_FITS_YV, "ce_tau_year_training_H2.rds"), label = "H2: training-only 2015-2021")
)

GATE_PARS <- c("alpha_R", "tau_year", "sigma_season_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2", "beta_T", "beta_P", "phi_obs")
PAIRS_PARS <- c("tau_year", "alpha_R", "sigma_season_year", "beta_T", "beta_P")

message("[hmc_audit config] ready. Models: ", paste(names(MODELS_TO_AUDIT), collapse = ", "))
