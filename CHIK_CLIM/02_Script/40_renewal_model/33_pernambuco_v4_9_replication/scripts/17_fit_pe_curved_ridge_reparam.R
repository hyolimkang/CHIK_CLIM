# Pernambuco v4.9 global-q FREE-QUADRATIC ridge-reparameterised fit
# (Model D). Supersedes the affine-only Model C, which left 26 post-warmup
# divergences at the low-q/high-depletion tail (see
# PE_curved_ridge_divergence_map.png / PE_free_quadratic_ridge_fit.csv --
# GATE PASS: RMSE reduced 20.7%, residual correlation with d/d^2 ~ 0).
#
# c0, u0, b1, b2 are FIXED constants estimated externally (script 16) from
# non-divergent draws of the affine full run -- NOT estimated in Stan.

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

pe_project_root <- function() {
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(current, "CHIK_CLIM.Rproj"))) return(current)
  candidate <- file.path(current, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Run from the outer repository or inner CHIK_CLIM project root.")
}

root <- pe_project_root()
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "40_renewal_model", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script", "40_renewal_model", "22_v4_3_short_2015_2019_q_calibration", "scripts", "02_fit_v4_3.R"))
source(file.path(root, "02_Script", "40_renewal_model", "33_pernambuco_v4_9_replication", "scripts", "02_build_pe_serology_window.R"))

base_dir <- file.path(root, "02_Script", "40_renewal_model", "33_pernambuco_v4_9_replication")
out_dir <- file.path(base_dir, "outputs", "modelD") # short path (MAX_PATH lesson from Model C)
table_dir <- file.path(root, "03_Output", "tables", "pernambuco_v4_9_curved_reparam")

LOGIT_Q_PRIOR_MEAN <- qlogis(0.10)
LOGIT_Q_PRIOR_SD <- 1.0
SERO_GEOGRAPHIC_SD <- 1.0

curve_constants <- readRDS(file.path(table_dir, "PE_free_quadratic_ridge_constants.rds"))
if (!isTRUE(curve_constants$gate_pass)) stop("Free-quadratic ridge GATE did not pass (script 16) -- do not proceed to Stan.")
C0 <- curve_constants$c0; U0 <- curve_constants$u0; B1 <- curve_constants$b1; B2 <- curve_constants$b2
message(sprintf("Using fixed curved-ridge constants: c0=%.6f, u0=%.6f, b1=%.6f, b2=%.6f", C0, U0, B1, B2))

make_year_contrast_basis <- function(Y) {
  if (Y < 2L) stop("At least two calendar years are needed for a sum-to-zero basis.")
  qr_input <- cbind(rep(1, Y), diag(Y)[, seq_len(Y - 1L), drop = FALSE])
  basis <- qr.Q(qr(qr_input))[, 2:Y, drop = FALSE]
  if (max(abs(colSums(basis))) > 1e-10 || max(abs(crossprod(basis) - diag(Y - 1L))) > 1e-10) {
    stop("Invalid annual-effect contrast basis.")
  }
  basis
}

make_pe_curved_data <- function(weekly, sero_on = TRUE) {
  prepared <- make_v4_3_data(
    weekly, imports_per_week = 1, year_effect_prior_sd = 0.40,
    t_sero = 1L, sero_pos = 0L, sero_n = 1L, kappa_sero = 50
  )
  n <- nrow(weekly)
  prepared$stan_data$seasonal_sin2 <- sin(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$seasonal_cos2 <- cos(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$t_sero <- NULL; prepared$stan_data$sero_pos <- NULL; prepared$stan_data$sero_n <- NULL
  prepared$stan_data$is_seed <- integer(n)
  prepared$stan_data$X_seed <- numeric(n)
  prepared$stan_data$N_fit <- n
  prepared$stan_data$fit_index <- seq_len(n)
  prepared$stan_data$logit_q_prior_mean <- LOGIT_Q_PRIOR_MEAN
  prepared$stan_data$logit_q_prior_sd <- LOGIT_Q_PRIOR_SD
  prepared$stan_data$c0 <- C0; prepared$stan_data$u0 <- U0; prepared$stan_data$b1 <- B1; prepared$stan_data$b2 <- B2
  prepared$stan_data$year_contrast_basis <- make_year_contrast_basis(prepared$stan_data$Y)

  if (sero_on) {
    sero_built <- build_pe_sero_stan_fields(as.Date(weekly$week_start))
    prepared$stan_data <- c(prepared$stan_data, sero_built$stan_fields)
    prepared$sero_audit <- sero_built$audit
  } else {
    prepared$stan_data$J_sero <- 0L
    prepared$stan_data$sero_n_positive <- integer(0)
    prepared$stan_data$sero_n_tested <- integer(0)
    prepared$stan_data$sero_window_start_idx <- integer(0)
    prepared$stan_data$sero_window_n_weeks <- integer(0)
    prepared$sero_audit <- NULL
  }
  prepared$stan_data$sero_geographic_sd <- SERO_GEOGRAPHIC_SD
  prepared
}

make_pe_curved_init_fn <- function(Y, J_sero, scalar_inits) {
  q_starts <- c(0.006, 0.010, 0.020, 0.050)
  gamma_offsets <- c(-0.03, -0.01, 0.01, 0.03)
  function(chain_id = 1L) {
    idx <- (as.integer(chain_id) - 1L) %% 4L + 1L
    base <- scalar_inits[[idx]]
    list(
      gamma_curve = gamma_offsets[idx], # near 0 -- lies ON the fitted quadratic curve at this chain's logit_q start
      logit_q = qlogis(q_starts[idx]),
      z_year_contrast = rep(0, Y - 1L),
      beta_sin1 = base$beta_sin, beta_cos1 = base$beta_cos,
      beta_sin2 = 0, beta_cos2 = 0,
      sigma_season_year = 0.05, z_season_year = rep(0, Y),
      phi_obs = base$phi_obs,
      z_geo = if (J_sero > 0L) as.array(rep(0, J_sero)) else numeric(0)
    )
  }
}

run_pe_v4_9_curved_ridge_reparam <- function(
    sero_on = as.logical(Sys.getenv("SERO_ON", "TRUE")),
    stage = Sys.getenv("STAGE", "pilot")) {
  if (!stage %in% c("pilot", "full")) stop("STAGE must be pilot or full.")
  defaults <- if (stage == "pilot") list(warmup = 300L, sampling = 300L) else list(warmup = 1200L, sampling = 1000L)
  settings <- list(
    warmup = as.integer(Sys.getenv("WARMUP", defaults$warmup)),
    sampling = as.integer(Sys.getenv("ITER_SAMPLING", defaults$sampling)),
    chains = 4L, adapt_delta = as.numeric(Sys.getenv("ADAPT_DELTA", "0.95")), max_treedepth = 14L, metric = "dense_e"
  )
  if (settings$warmup < 1L || settings$sampling < 1L) stop("Warmup and sampling iterations must be positive.")
  if (settings$adapt_delta > 0.95 + 1e-9) message("NOTE: adapt_delta > 0.95 requested -- per instruction this should only be used after a clean 0.95 curved-ridge fit.")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
  weekly <- readRDS(file.path(root, "03_Output", "tables", "pernambuco_v4_9_replication", "pernambuco_weekly_input.rds"))
  prepared <- make_pe_curved_data(weekly, sero_on = sero_on)
  scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_pe_curved_init_fn(prepared$stan_data$Y, prepared$stan_data$J_sero, scalar_inits)
  stan_path <- file.path(root, "02_Script", "stan", "renewal_pernambuco_v4_9_global_q_curved_reparam.stan")
  sample_file <- file.path(progress_dir, paste0("chain_", if (sero_on) "U14" else "caseonly"))

  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data, chains = settings$chains,
    iter = settings$warmup + settings$sampling, warmup = settings$warmup,
    seed = 20260914L, init = init_fn,
    control = list(adapt_delta = settings$adapt_delta, max_treedepth = settings$max_treedepth, metric = settings$metric),
    refresh = 25, sample_file = sample_file
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit,
    c("gamma_curve", "alpha_R", "z_year_contrast", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2", "sigma_season_year", "z_season_year", "phi_obs", "logit_q"),
    settings$max_treedepth
  )
  tag <- if (sero_on) "modelD_U14" else "modelD_caseonly"
  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, sero_audit = prepared$sero_audit,
    config = c(settings, list(model_version = "pernambuco_v4_9_global_q_curved_reparam", tag = tag,
      stage = stage, sero_on = sero_on, curve_c0 = C0, curve_u0 = U0, curve_b1 = B1, curve_b2 = B2,
      logit_q_reference = U0, q_prior_mean = 0.10,
      q_prior_logit_sd = LOGIT_Q_PRIOR_SD, elapsed_seconds = elapsed_seconds,
      likelihood_note = "Same v4.9 likelihood and priors; computational curved-coordinate transform only.")),
    hmc = hmc
  )
  ad_suffix <- if (abs(settings$adapt_delta - 0.95) > 1e-9) paste0("_ad", gsub("\\.", "", sprintf("%.2f", settings$adapt_delta))) else ""
  fit_path <- file.path(out_dir, paste0("pe_curved_", if (sero_on) "U14" else "caseonly", "_", stage, ad_suffix, ".rds"))
  saveRDS(bundle, fit_path)
  saveRDS(prepared$stan_data, file.path(out_dir, paste0("stan_data_", tag, "_", stage, ad_suffix, ".rds")))
  if (!is.null(prepared$sero_audit)) {
    write_csv(prepared$sero_audit, file.path(table_dir, paste0("serology_window_", tag, ".csv")))
  }
  message(sprintf("[%s] HMC gate: %s | %.0fs | %s", tag, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_pe_v4_9_curved_ridge_reparam()
