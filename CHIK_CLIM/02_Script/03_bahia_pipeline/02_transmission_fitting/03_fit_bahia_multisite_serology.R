# Bahia integrated multi-site serology extension of the frozen v4.9
# structure (renewal_ceara_v4_9_bahia_multisite_serology.stan). NOT a new
# transmission model -- see the .stan file header. Fits ONLY q=0.05 and
# q=0.10 as pilots, per explicit instruction; do not run the full q sweep
# without review.
#
# Usage: QVAL=0.05 Rscript 07_fit_bahia_multisite_serology.R

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

bahia_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

root <- ROOT # from 00_project_setup.R -- bahia_root() left defined above but unused, no scientific change
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "22_v4_3_short_2015_2019_q_calibration", "scripts", "02_fit_v4_3.R"))
source(file.path(root, "02_Script", "03_bahia_pipeline", "02_transmission_fitting", "02_build_bahia_serology_windows.R"))

base_dir <- file.path(root, "03_Output", "model_fits", "bahia", "v4_9_replication")
SERO_GEOGRAPHIC_SD <- 1.0 # FIXED, not estimated (per instruction)

make_bahia_multisite_data <- function(weekly, imports_per_week, year_effect_prior_sd, q) {
  prepared <- make_v4_3_data(weekly, imports_per_week, year_effect_prior_sd,
                              t_sero = 1L, sero_pos = 0L, sero_n = 1L, kappa_sero = 50) # t_sero/sero_pos/sero_n unused by this model; kept only because make_v4_3_data requires them
  n <- nrow(weekly)
  prepared$stan_data$seasonal_sin2 <- sin(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$seasonal_cos2 <- cos(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$q <- q
  prepared$stan_data$t_sero <- NULL; prepared$stan_data$sero_pos <- NULL; prepared$stan_data$sero_n <- NULL # not in this Stan file's data block
  prepared$stan_data$is_seed <- integer(n)
  prepared$stan_data$X_seed <- numeric(n)
  prepared$stan_data$N_fit <- n
  prepared$stan_data$fit_index <- seq_len(n)

  sero_built <- build_bahia_sero_stan_fields(as.Date(weekly$week_start))
  prepared$stan_data <- c(prepared$stan_data, sero_built$stan_fields)
  prepared$stan_data$sero_geographic_sd <- SERO_GEOGRAPHIC_SD
  prepared$sero_audit <- sero_built$audit
  prepared
}

make_bahia_multisite_init_fn <- function(Y, J_sero, td14_scalar_inits) {
  offsets <- c(-0.06, -0.02, 0.02, 0.06)
  function(chain_id = 1L) {
    cid <- as.integer(chain_id)
    base <- td14_scalar_inits[[(cid - 1L) %% length(td14_scalar_inits) + 1L]]
    offset <- offsets[(cid - 1L) %% 4L + 1L]
    list(
      alpha_R = base$alpha_R, phi_obs = base$phi_obs,
      beta_sin1 = base$beta_sin, beta_cos1 = base$beta_cos,
      beta_sin2 = 0, beta_cos2 = 0,
      sigma_season_year = 0.05, z_season_year = rep(0, Y),
      z_year = rep(offset / 0.40, Y),
      z_geo = rep(0, J_sero)
    )
  }
}

run_fit_bahia_multisite <- function(qval = as.numeric(Sys.getenv("QVAL", "0.05")), stage = Sys.getenv("STAGE", "pilot")) {
  tag <- paste0("q", sprintf("%.2f", qval))
  message("[bahia-multisite-", tag, "] frozen v4.9 structure + 6-site geo-adjusted serology | q=", qval, " | stage=", stage)

  out_dir <- file.path(base_dir, "outputs_multisite_serology", tag)
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  fixed <- list(
    date_start = as.Date("2015-01-04"), date_end = as.Date("2025-12-21"),
    imports_per_week = 1, year_effect_prior_sd = 0.40,
    chains = 4L, max_treedepth = 14L, adapt_delta = .95, metric = "dense_e"
  )
  default_warmup <- if (stage == "pilot") "150" else "800"
  default_sampling <- if (stage == "pilot") "100" else "400"
  settings <- list(warmup = as.integer(Sys.getenv("WARMUP", default_warmup)),
                    iter_sampling = as.integer(Sys.getenv("ITER_SAMPLING", default_sampling)))

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "29") |> dplyr::filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  demography_path <- file.path(root, "01_Data", "bahia_weekly_demography_2015_2025.rds")
  weekly_demography <- readRDS(demography_path) |> dplyr::filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch.")

  prepared <- make_bahia_multisite_data(weekly, fixed$imports_per_week, fixed$year_effect_prior_sd, qval)

  iter_total <- settings$warmup + settings$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  td14_scalar_inits <- load_td14_scalar_inits(root)
  init_fn <- make_bahia_multisite_init_fn(prepared$stan_data$Y, prepared$stan_data$J_sero, td14_scalar_inits)

  stan_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_9_bahia_multisite_serology.stan")
  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data,
    chains = fixed$chains, iter = iter_total, warmup = settings$warmup,
    seed = 20260912L, init = init_fn,
    control = list(adapt_delta = fixed$adapt_delta, max_treedepth = fixed$max_treedepth, metric = fixed$metric),
    refresh = 25, sample_file = sample_file_base
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                                  "sigma_season_year", "z_season_year", "phi_obs", "z_geo"), fixed$max_treedepth)

  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, sero_audit = prepared$sero_audit,
    config = c(fixed, settings, list(model_version = "bahia_v4_9_multisite_serology", q = qval, tag = tag, stage = stage,
                                      sero_geographic_sd = SERO_GEOGRAPHIC_SD,
                                      years = prepared$years, elapsed_seconds = elapsed_seconds,
                                      stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE))),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_bahia_multisite_fit_", tag, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[bahia-multisite-%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", tag, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_bahia_multisite()
