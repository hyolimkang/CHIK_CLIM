# Pernambuco Model D -- LONG-WARMUP diagnostic (same quadratic model,
# same fixed c0/u0/b1/b2, same priors). Purpose: give the dense mass
# matrix / step-size adaptation more time before concluding the chain-4
# pattern is a structural geometry failure vs. an adaptation-quality
# issue. No model component changed.
#
# Initialisation: each chain starts from a DIFFERENT non-divergent draw of
# the existing affine posterior (Model C), transformed into the new
# quadratic coordinates (gamma_curve = alpha_R - c0 - b1*d - b2*d^2).
# All other parameters (z_year_contrast, beta_sin/cos, sigma_season_year,
# z_season_year, phi_obs, z_geo) are shared unchanged between Model C and
# Model D, so the SAME draw's values are reused directly for those.

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "40_renewal_model", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script", "40_renewal_model", "22_v4_3_short_2015_2019_q_calibration", "scripts", "02_fit_v4_3.R"))
source(file.path(root, "02_Script", "40_renewal_model", "33_pernambuco_v4_9_replication", "scripts", "02_build_pe_serology_window.R"))

base_dir <- file.path(root, "02_Script", "40_renewal_model", "33_pernambuco_v4_9_replication")
out_dir <- file.path(base_dir, "outputs", "modelD")
table_dir <- file.path(root, "03_Output", "tables", "pernambuco_v4_9_curved_reparam")

LOGIT_Q_PRIOR_MEAN <- qlogis(0.10)
LOGIT_Q_PRIOR_SD <- 1.0
SERO_GEOGRAPHIC_SD <- 1.0

curve_constants <- readRDS(file.path(table_dir, "PE_free_quadratic_ridge_constants.rds"))
C0 <- curve_constants$c0; U0 <- curve_constants$u0; B1 <- curve_constants$b1; B2 <- curve_constants$b2
message(sprintf("Using fixed curved-ridge constants: c0=%.6f, u0=%.6f, b1=%.6f, b2=%.6f", C0, U0, B1, B2))

make_year_contrast_basis <- function(Y) {
  qr_input <- cbind(rep(1, Y), diag(Y)[, seq_len(Y - 1L), drop = FALSE])
  qr.Q(qr(qr_input))[, 2:Y, drop = FALSE]
}

make_pe_curved_data <- function(weekly, sero_on = TRUE) {
  prepared <- make_v4_3_data(weekly, imports_per_week = 1, year_effect_prior_sd = 0.40,
                              t_sero = 1L, sero_pos = 0L, sero_n = 1L, kappa_sero = 50)
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
  sero_built <- build_pe_sero_stan_fields(as.Date(weekly$week_start))
  prepared$stan_data <- c(prepared$stan_data, sero_built$stan_fields)
  prepared$sero_audit <- sero_built$audit
  prepared$stan_data$sero_geographic_sd <- SERO_GEOGRAPHIC_SD
  prepared
}

# ---- Build dispersed, posterior-informed inits from the Model C (affine) full run ----
build_posterior_informed_inits <- function() {
  aff <- readRDS(file.path(out_dir, "..", "modelC", "pe_ridge_U14_full.rds"))
  fit_aff <- aff$fit
  sp <- rstan::get_sampler_params(fit_aff, inc_warmup = FALSE)
  divergent <- unlist(lapply(sp, function(x) x[, "divergent__"]))
  nondiv <- which(divergent == 0)

  # permuted=FALSE preserves chain-major [iteration, chain, param] order,
  # matching get_sampler_params()'s per-chain ordering used for `divergent`
  # -- permuted=TRUE draws are randomly shuffled and would misalign here.
  pars_shared <- c("logit_q", "alpha_R", "z_year_contrast", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                    "sigma_season_year", "z_season_year", "phi_obs", "z_geo")
  arr <- rstan::extract(fit_aff, pars = pars_shared, permuted = FALSE, inc_warmup = FALSE)
  n_iter <- dim(arr)[1]; n_chain <- dim(arr)[2]
  flat_idx <- expand.grid(iter = seq_len(n_iter), chain = seq_len(n_chain))
  flat_idx <- flat_idx[order(flat_idx$chain, flat_idx$iter), ] # chain-major, matches `divergent`
  dimnames3 <- dimnames(arr)[[3]]
  get_flat_scalar <- function(par) {
    col <- grep(paste0("^", par, "(\\[|$)"), dimnames3)
    if (length(col) == 1) mapply(function(i, c) arr[i, c, col], flat_idx$iter, flat_idx$chain)
    else NA
  }
  get_flat_vector <- function(par_prefix) {
    cols <- grep(paste0("^", par_prefix, "\\["), dimnames3)
    t(mapply(function(i, c) arr[i, c, cols], flat_idx$iter, flat_idx$chain))
  }

  logit_q_all <- get_flat_scalar("logit_q")
  alpha_R_all <- get_flat_scalar("alpha_R")
  logit_q_nd <- logit_q_all[nondiv]

  # 4 dispersed draws spanning the typical (non-divergent) posterior range.
  target_q <- quantile(logit_q_nd, probs = c(0.10, 0.35, 0.65, 0.90))
  pick_idx <- sapply(target_q, function(tq) nondiv[which.min(abs(logit_q_nd - tq))])
  message("Dispersed init draw logit_q values: ", paste(sprintf("%.4f", logit_q_all[pick_idx]), collapse = ", "))

  z_year_contrast_all <- get_flat_vector("z_year_contrast")
  z_season_year_all <- get_flat_vector("z_season_year")
  z_geo_all <- get_flat_vector("z_geo")
  beta_sin1_all <- get_flat_scalar("beta_sin1"); beta_cos1_all <- get_flat_scalar("beta_cos1")
  beta_sin2_all <- get_flat_scalar("beta_sin2"); beta_cos2_all <- get_flat_scalar("beta_cos2")
  sigma_season_year_all <- get_flat_scalar("sigma_season_year")
  phi_obs_all <- get_flat_scalar("phi_obs")

  lapply(pick_idx, function(idx) {
    d_val <- logit_q_all[idx] - U0
    gamma_curve_val <- alpha_R_all[idx] - C0 - B1 * d_val - B2 * d_val^2
    list(
      gamma_curve = gamma_curve_val,
      logit_q = logit_q_all[idx],
      z_year_contrast = as.array(z_year_contrast_all[idx, ]),
      beta_sin1 = beta_sin1_all[idx], beta_cos1 = beta_cos1_all[idx],
      beta_sin2 = beta_sin2_all[idx], beta_cos2 = beta_cos2_all[idx],
      sigma_season_year = sigma_season_year_all[idx],
      z_season_year = as.array(z_season_year_all[idx, ]),
      phi_obs = phi_obs_all[idx],
      z_geo = as.array(z_geo_all[idx, ])
    )
  })
}

run_longwarmup_diagnostic <- function() {
  settings <- list(warmup = 1200L, sampling = 500L, chains = 4L, adapt_delta = 0.95, max_treedepth = 14L, metric = "dense_e")
  progress_dir <- file.path(out_dir, "chains_longwarmup")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  weekly <- readRDS(file.path(root, "03_Output", "tables", "pernambuco_v4_9_replication", "pernambuco_weekly_input.rds"))
  prepared <- make_pe_curved_data(weekly, sero_on = TRUE)

  init_list <- tryCatch(build_posterior_informed_inits(), error = function(e) {
    message("Posterior-informed init failed (", conditionMessage(e), ") -- falling back to dispersed scalar init scheme.")
    NULL
  })
  if (is.null(init_list)) {
    scalar_inits <- load_td14_scalar_inits(root)
    q_starts <- c(0.006, 0.010, 0.020, 0.050)
    gamma_offsets <- c(-0.03, -0.01, 0.01, 0.03)
    init_list <- lapply(1:4, function(i) list(
      gamma_curve = gamma_offsets[i], logit_q = qlogis(q_starts[i]),
      z_year_contrast = rep(0, prepared$stan_data$Y - 1L),
      beta_sin1 = scalar_inits[[i]]$beta_sin, beta_cos1 = scalar_inits[[i]]$beta_cos,
      beta_sin2 = 0, beta_cos2 = 0, sigma_season_year = 0.05, z_season_year = rep(0, prepared$stan_data$Y),
      phi_obs = scalar_inits[[i]]$phi_obs, z_geo = as.array(rep(0, prepared$stan_data$J_sero))
    ))
  }
  init_fn <- function(chain_id = 1L) init_list[[(as.integer(chain_id) - 1L) %% 4L + 1L]]

  stan_path <- file.path(root, "02_Script", "stan", "renewal_pernambuco_v4_9_global_q_curved_reparam.stan")
  sample_file <- file.path(progress_dir, "chain")

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
  bundle <- list(fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, sero_audit = prepared$sero_audit,
                  config = c(settings, list(model_version = "pernambuco_v4_9_global_q_curved_reparam_longwarmup",
                                             curve_c0 = C0, curve_u0 = U0, curve_b1 = B1, curve_b2 = B2,
                                             elapsed_seconds = elapsed_seconds)),
                  hmc = hmc)
  fit_path <- file.path(out_dir, "pe_curved_U14_longwarmup.rds")
  saveRDS(bundle, fit_path)
  message(sprintf("[longwarmup] HMC gate: %s | %.0fs | %s", if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_longwarmup_diagnostic()
