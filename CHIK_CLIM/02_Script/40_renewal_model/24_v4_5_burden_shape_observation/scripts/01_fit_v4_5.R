# v4.5: FINAL 2015-2019 identifiability test. Latent weekly renewal process,
# q definition, and q's weak logit-normal calibration prior are UNCHANGED
# from v4.3/v4.4. Only the case observation model changes: (1) absolute
# yearly burden (Y=5 NB terms) ties q to total infections per year, and
# (2) weekly temporal shape within each year, conditional on the OBSERVED
# yearly total, via a manual Dirichlet-multinomial log density (one joint
# term per year -- no weekly-independence assumption at all). Juazeiro 2018
# serology is included unconditionally (same conservative treatment as
# v4.2-v4.4: beta-binomial, kappa_sero=50, fixed).
#
# Multimodal initialization: chains 1-2 near q~0.05, chains 3-4 near q~0.8,
# testing whether all four converge to one coherent posterior.
#
# Usage: Rscript 01_fit_v4_5.R

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

v4_5_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

root <- v4_5_root()
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "40_renewal_model", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script", "40_renewal_model", "22_v4_3_short_2015_2019_q_calibration", "scripts", "02_fit_v4_3.R"))

base_dir <- file.path(root, "02_Script", "40_renewal_model", "24_v4_5_burden_shape_observation")

run_fit_v4_5 <- function(fit_id = "v4_5") {
  tag <- fit_id
  message("[", tag, "] v4.5 burden+shape observation model | Juazeiro serology included unconditionally")

  out_dir <- file.path(base_dir, "outputs", tag)
  progress_dir <- file.path(out_dir, "chains")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  fixed <- list(
    date_start = as.Date("2015-01-04"), date_end = as.Date("2019-12-29"),
    imports_per_week = 1, year_effect_prior_sd = 0.40,
    chains = 4L, max_treedepth = 14L, adapt_delta = .95, metric = "dense_e",
    sero_pos = 103L, sero_n = 404L, kappa_sero = 50
  )
  settings <- list(warmup = 400L, iter_sampling = 200L)

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |> filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  demography_path <- file.path(root, "01_Data", "ceara_weekly_demography_2015_2025.rds")
  weekly_demography <- readRDS(demography_path) |> filter(week_start >= fixed$date_start, week_start <= fixed$date_end)
  weekly <- weekly_cases |> select(-population) |> inner_join(weekly_demography, by = "week_start") |> arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch.")
  message(sprintf("[%s] %s to %s | %d weeks", tag, fixed$date_start, fixed$date_end, nrow(weekly)))

  sero_week <- find_sero_week_index(as.Date(weekly$week_start))
  prepared <- make_v4_3_data(weekly, fixed$imports_per_week, fixed$year_effect_prior_sd,
                              sero_week$index, fixed$sero_pos, fixed$sero_n, fixed$kappa_sero)
  prepared$stan_data$use_serology <- NULL # v4.5's data block has no such field -- serology is unconditional

  iter_total <- settings$warmup + settings$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  td14_scalar_inits <- load_td14_scalar_inits(root)

  # Multimodal initialization (mirrors v4.4 Fit A4/B4): chains 1-2 low-q
  # region, chains 3-4 high-q region -- not to pick a winner, but to test
  # whether both regions still lead to persistent modes under the two-part
  # burden+shape observation model.
  q_init_by_chain <- c(0.04, 0.07, 0.75, 0.90)
  offsets <- c(-0.06, -0.02, 0.02, 0.06)
  init_fn <- function(chain_id = 1L) {
    cid <- as.integer(chain_id)
    base <- td14_scalar_inits[[(cid - 1L) %% length(td14_scalar_inits) + 1L]]
    offset <- offsets[(cid - 1L) %% 4L + 1L]
    c(base, list(z_year = rep(offset / fixed$year_effect_prior_sd, prepared$stan_data$Y),
                 eta_q = qlogis(q_init_by_chain[cid]),
                 phi_burden = 20, kappa_shape = 20))
  }

  stan_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_5_burden_shape.stan")
  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data,
    chains = fixed$chains, iter = iter_total, warmup = settings$warmup,
    seed = 20260912L, init = init_fn,
    control = list(adapt_delta = fixed$adapt_delta, max_treedepth = fixed$max_treedepth, metric = fixed$metric),
    refresh = 25, sample_file = sample_file_base
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_burden", "kappa_shape", "eta_q"), fixed$max_treedepth)

  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, sero_week = sero_week,
    config = c(fixed, settings, list(model_version = "v4_5_burden_shape_observation", fit_id = fit_id, tag = tag,
                                      years = prepared$years, elapsed_seconds = elapsed_seconds,
                                      stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE))),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_5_fit_", tag, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", tag, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_v4_5()
