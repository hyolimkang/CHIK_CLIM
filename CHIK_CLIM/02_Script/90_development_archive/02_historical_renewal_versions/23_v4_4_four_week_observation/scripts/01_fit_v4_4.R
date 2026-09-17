# v4.4 DIAGNOSTIC: 4-week observation-likelihood sensitivity (does weekly
# pseudo-replication create/amplify the low-q/high-q multimodality found in
# v4.3?). The latent weekly renewal/S/U/R0 process is UNCHANGED; only the
# case observation likelihood is aggregated into 78 non-overlapping 4-week
# blocks (weeks 1-312 of the 313-week 2015-2019 series). The final
# incomplete block (week 313 alone) is EXCLUDED from the likelihood in BOTH
# fits (documented rule, Section 2) via include_in_block=0 -- the latent
# process still runs through week 313 regardless (needed for S_2019).
#
# FIT A4: 4-week likelihood, weak q prior, NO serology.
# FIT B4: same + Juazeiro 2018 (103/404, beta-binomial, kappa_sero=50).
# Both use MULTIMODAL initialization: chains 1-2 in q~0.03-0.08, chains 3-4
# in q~0.70-0.90 (Section 5) -- not to pick a winner, but to test whether
# BOTH regions still lead to persistent modes under 4-week aggregation.
#
# Usage: FIT=A4 Rscript 01_fit_v4_4.R
#        FIT=B4 Rscript 01_fit_v4_4.R

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

v4_4_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

root <- v4_4_root()
setwd(root)
here::i_am("CHIK_CLIM.Rproj")
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script", "00_shared", "legacy_model_functions", "22_v4_3_short_2015_2019_q_calibration", "scripts", "02_fit_v4_3.R"))

base_dir <- file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "23_v4_4_four_week_observation")
BLOCK_SIZE <- 4L

build_block_structure <- function(weekly) {
  n_weeks <- nrow(weekly)
  n_full_blocks <- n_weeks %/% BLOCK_SIZE
  n_weeks_used <- n_full_blocks * BLOCK_SIZE
  block_id <- rep(n_full_blocks, n_weeks) # trailing leftover weeks point at last block, but are excluded via include_in_block
  block_id[seq_len(n_weeks_used)] <- rep(seq_len(n_full_blocks), each = BLOCK_SIZE)
  include_in_block <- as.integer(seq_len(n_weeks) <= n_weeks_used)
  observed_block <- vapply(seq_len(n_full_blocks), function(b) sum(weekly$cases[block_id == b & include_in_block == 1]), integer(1))
  list(n_blocks = n_full_blocks, block_id = block_id, include_in_block = include_in_block, observed_block = observed_block,
       n_weeks_excluded = n_weeks - n_weeks_used)
}

run_fit_v4_4 <- function(fit_id = Sys.getenv("FIT", "A4")) {
  use_serology <- if (fit_id == "A4") 0L else 1L
  tag <- paste0("fit", fit_id)
  message("[", tag, "] v4.4 four-week observation diagnostic | use_serology=", use_serology)

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

  sero_week <- find_sero_week_index(as.Date(weekly$week_start))
  prepared <- make_v4_3_data(weekly, fixed$imports_per_week, fixed$year_effect_prior_sd,
                              sero_week$index, fixed$sero_pos, fixed$sero_n, fixed$kappa_sero)
  prepared$stan_data$use_serology <- NULL

  blocks <- build_block_structure(weekly)
  message(sprintf("[%s] %d weeks -> %d full 4-week blocks (%d weeks used, %d trailing week(s) excluded from likelihood)",
                   tag, nrow(weekly), blocks$n_blocks, blocks$n_blocks * BLOCK_SIZE, blocks$n_weeks_excluded))
  prepared$stan_data$n_blocks <- blocks$n_blocks
  prepared$stan_data$block_id <- blocks$block_id
  prepared$stan_data$include_in_block <- blocks$include_in_block
  prepared$stan_data$observed_block <- blocks$observed_block
  prepared$stan_data$use_serology <- use_serology

  iter_total <- settings$warmup + settings$iter_sampling
  sample_file_base <- file.path(progress_dir, "chain")
  td14_scalar_inits <- load_td14_scalar_inits(root)

  # Section 5: multimodal initialization -- chains 1-2 low region, 3-4 high region.
  q_init_by_chain <- c(0.04, 0.07, 0.75, 0.90)
  offsets <- c(-0.06, -0.02, 0.02, 0.06)
  init_fn <- function(chain_id = 1L) {
    cid <- as.integer(chain_id)
    base <- td14_scalar_inits[[(cid - 1L) %% length(td14_scalar_inits) + 1L]]
    offset <- offsets[(cid - 1L) %% 4L + 1L]
    c(base, list(z_year = rep(offset / fixed$year_effect_prior_sd, prepared$stan_data$Y),
                 eta_q = qlogis(q_init_by_chain[cid]), phi_obs_block = 20))
  }

  stan_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_3_4week_diagnostic.stan")
  started <- Sys.time()
  fit <- rstan::stan(
    file = stan_path, data = prepared$stan_data,
    chains = fixed$chains, iter = iter_total, warmup = settings$warmup,
    seed = 20260912L, init = init_fn,
    control = list(adapt_delta = fixed$adapt_delta, max_treedepth = fixed$max_treedepth, metric = fixed$metric),
    refresh = 25, sample_file = sample_file_base
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs_block", "eta_q"), fixed$max_treedepth)

  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, sero_week = sero_week, blocks = blocks,
    config = c(fixed, settings, list(model_version = "v4_4_four_week_observation", fit_id = fit_id, tag = tag,
                                      use_serology = use_serology, years = prepared$years, elapsed_seconds = elapsed_seconds,
                                      stan_source = normalizePath(stan_path, winslash = "/", mustWork = TRUE))),
    hmc = hmc
  )
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_4_fit_", tag, ".rds"))
  saveRDS(bundle, fit_path)
  message(sprintf("[%s] DONE | HMC gate: %s | elapsed=%.0fs | saved: %s", tag, if (hmc$hmc_pass) "PASS" else "FAIL", elapsed_seconds, fit_path))
  invisible(bundle)
}

if (sys.nframe() == 0L) run_fit_v4_4()
