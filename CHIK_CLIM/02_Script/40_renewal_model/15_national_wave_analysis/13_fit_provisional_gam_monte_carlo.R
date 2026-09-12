# Monte Carlo uncertainty propagation for the provisional GAM checkpoint
# (Section 12). Primary climate window (-6:-2) and primary (1.0x) FOI scale
# only. For each of >=200 iterations: resample ONE Re_early_6 draw per wave
# from the real stored posterior (brazil_chik_major_wave_early_re_posterior_draws.rds)
# and ONE matching S_onset draw per wave from script 10's real raw-draw file
# (brazil_chik_wave_onset_susceptibility_draws_scale1_0.rds) -- both are
# genuine posterior draws, not approximated from quantiles. Refit the three
# parsimonious GAMs on each resampled dataset and store the distribution of
# the susceptibility effect, climate smooths, and comparison metrics.

required_mc_packages <- c("here", "dplyr", "readr", "tibble", "mgcv", "tidyr")
missing_mc_packages <- required_mc_packages[
  !vapply(required_mc_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_mc_packages)) stop("Missing package(s): ", paste(missing_mc_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tibble); library(mgcv); library(tidyr) })

mc_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}
source(file.path(mc_root(), "02_Script", "40_renewal_model", "15_national_wave_analysis", "00_national_wave_analysis_helpers.R"))

mc_settings <- list(n_iterations = as.integer(Sys.getenv("PROVISIONAL_GAM_MC_ITER", "200")), seed = 20261010L)

run_provisional_gam_monte_carlo <- function() {
  paths <- ensure_national_output_dirs()
  episode_path <- file.path(paths$table, "brazil_chik_Re_S_climate_episode_analysis_scale1_0.csv")
  re_draws_path <- file.path(paths$table, "brazil_chik_major_wave_early_re_posterior_draws.rds")
  s_draws_path <- file.path(paths$table, "brazil_chik_wave_onset_susceptibility_draws_scale1_0.rds")
  for (p in c(episode_path, re_draws_path, s_draws_path)) if (!file.exists(p)) stop("Missing required input: ", p)

  episode <- read_csv(episode_path, show_col_types = FALSE)
  re_draws <- readRDS(re_draws_path) |> filter(window_weeks == 6, wave_id %in% episode$wave_id)
  s_draws <- readRDS(s_draws_path) |> filter(wave_id %in% episode$wave_id)

  n_re_draws <- re_draws |> distinct(wave_id, draw) |> count(wave_id) |> pull(n) |> unique()
  n_s_draws <- s_draws |> distinct(wave_id, draw) |> count(wave_id) |> pull(n) |> unique()
  message("[mc] Re draws per wave: ", paste(n_re_draws, collapse = ","), "; S draws per wave: ", paste(n_s_draws, collapse = ","))

  base <- episode |> select(state, wave_id, S_prop_point = S_prop, temp = temp_pre6_2, precip = precip_pre6_2) |>
    mutate(state = factor(state), log1p_precip = log1p(precip))

  set.seed(mc_settings$seed)
  results <- vector("list", mc_settings$n_iterations)
  n_ok <- 0L
  for (iter in seq_len(mc_settings$n_iterations)) {
    re_sample <- re_draws |> group_by(wave_id) |> slice_sample(n = 1) |> ungroup() |> select(wave_id, Re_sample = Re)
    s_sample <- s_draws |> group_by(wave_id) |> slice_sample(n = 1) |> ungroup() |> select(wave_id, S_sample = S_onset)
    data_iter <- base |> inner_join(re_sample, by = "wave_id") |> inner_join(s_sample, by = "wave_id") |>
      filter(Re_sample > 0) |> mutate(logRe = log(Re_sample), S_prop = S_sample)
    if (nrow(data_iter) < nrow(base) * 0.9) next

    fit_result <- tryCatch({
      m_climate <- gam(logRe ~ s(temp, k = 4) + s(log1p_precip, k = 4) + s(state, bs = "re"), data = data_iter, method = "REML")
      m_climate_s <- gam(logRe ~ s(temp, k = 4) + s(log1p_precip, k = 4) + s(S_prop, k = 3) + s(state, bs = "re"), data = data_iter, method = "REML")
      m_gate <- gam(logRe ~ offset(log(S_prop)) + s(temp, k = 4) + s(log1p_precip, k = 4) + s(state, bs = "re"), data = data_iter, method = "REML")
      s_climate <- summary(m_climate); s_climate_s <- summary(m_climate_s); s_gate <- summary(m_gate)
      smooth_s_match <- grep("^s\\(S_prop\\)", rownames(s_climate_s$s.table), value = TRUE)
      tibble(
        iter = iter,
        AIC_climate = AIC(m_climate), AIC_climate_S = AIC(m_climate_s), AIC_gate = AIC(m_gate),
        dev_expl_climate = s_climate$dev.expl, dev_expl_climate_S = s_climate_s$dev.expl, dev_expl_gate = s_gate$dev.expl,
        edf_temp_climate = s_climate$s.table["s(temp)", "edf"], edf_precip_climate = s_climate$s.table["s(log1p_precip)", "edf"],
        edf_S = if (length(smooth_s_match)) s_climate_s$s.table[smooth_s_match[1], "edf"] else NA_real_,
        p_S = if (length(smooth_s_match)) s_climate_s$s.table[smooth_s_match[1], "p-value"] else NA_real_
      )
    }, error = function(e) NULL)
    if (!is.null(fit_result)) { results[[iter]] <- fit_result; n_ok <- n_ok + 1L }
  }

  mc_table <- bind_rows(results)
  message("[mc] successful iterations: ", n_ok, " / ", mc_settings$n_iterations)
  if (n_ok < 0.5 * mc_settings$n_iterations) stop("Fewer than half of MC iterations succeeded -- investigate before trusting this summary.")

  mc_summary <- tibble(
    metric = c("AIC_climate", "AIC_climate_S", "AIC_gate", "dev_expl_climate", "dev_expl_climate_S", "dev_expl_gate",
               "edf_temp_climate", "edf_precip_climate", "edf_S", "p_S"),
    median = vapply(mc_table[-1], median, numeric(1), na.rm = TRUE),
    q025 = vapply(mc_table[-1], quantile, numeric(1), probs = .025, na.rm = TRUE),
    q975 = vapply(mc_table[-1], quantile, numeric(1), probs = .975, na.rm = TRUE)
  )
  mc_summary$prop_p_S_below_0_05 <- NA_real_
  mc_summary$prop_p_S_below_0_05[mc_summary$metric == "p_S"] <- mean(mc_table$p_S < 0.05, na.rm = TRUE)
  mc_summary$prop_AIC_climate_S_better <- NA_real_
  mc_summary$prop_AIC_climate_S_better[mc_summary$metric == "AIC_climate_S"] <- mean(mc_table$AIC_climate_S < mc_table$AIC_climate, na.rm = TRUE)

  write_csv(mc_table, file.path(paths$table, "brazil_chik_provisional_gam_monte_carlo_draws.csv"))
  write_csv(mc_summary, file.path(paths$table, "brazil_chik_provisional_gam_monte_carlo_summary.csv"))
  message("[mc] P(AIC_climate_S < AIC_climate) across iterations = ", round(mean(mc_table$AIC_climate_S < mc_table$AIC_climate, na.rm = TRUE), 3))
  message("[mc] P(p_S < 0.05) across iterations = ", round(mean(mc_table$p_S < 0.05, na.rm = TRUE), 3))
  invisible(list(draws = mc_table, summary = mc_summary))
}

if (sys.nframe() == 0L) run_provisional_gam_monte_carlo()
