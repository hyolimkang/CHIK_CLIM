# Nationwide early-phase Re estimation for the frozen V2 major-wave census.
# The existing validated K-episode renewal Stan program is used unchanged for
# each window, retaining its shared NB2 overdispersion parameter (phi).

required_early_re_packages <- c("here", "rstan", "posterior")
missing_early_re_packages <- required_early_re_packages[!vapply(required_early_re_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_early_re_packages)) stop("Missing package(s): ", paste(missing_early_re_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(posterior) })
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

early_re_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}
source(file.path(early_re_root(), "02_Script", "40_renewal_model", "15_national_wave_analysis", "00_national_wave_analysis_helpers.R"))

early_re_settings <- list(
  gi_shape = 4, gi_rate = 2, gi_support_weeks = 8L,
  windows = c(4L, 6L, 8L), primary_window = 6L,
  iter = as.integer(Sys.getenv("NATIONAL_EARLY_RE_ITER", "2000")),
  warmup = as.integer(Sys.getenv("NATIONAL_EARLY_RE_WARMUP", "1000")),
  chains = as.integer(Sys.getenv("NATIONAL_EARLY_RE_CHAINS", "4")),
  adapt_delta = .90, max_treedepth = 11L
)
if (early_re_settings$warmup >= early_re_settings$iter) stop("Warmup must be smaller than iterations.")
early_re_settings$weights <- {
  x <- diff(pgamma(0:early_re_settings$gi_support_weeks, shape = early_re_settings$gi_shape, rate = early_re_settings$gi_rate))
  x / sum(x)
}

renewal_lambda <- function(cases, target_index, weights) {
  history <- target_index - seq_along(weights)
  if (any(history < 1L)) return(NA_real_)
  sum(weights * cases[history])
}

build_window_batch <- function(major, weekly_cases, window, settings) {
  eligible_column <- paste0("re", window, "_prefit_eligible")
  status <- major |> transmute(
    state, wave_id, eligible = .data[[eligible_column]],
    fit_status = if_else(eligible, "pending", "not_attempted"),
    fit_reason = if_else(eligible, NA_character_, "frozen census pre-fit eligibility is FALSE")
  )
  matrices <- list(); kept <- integer()
  for (i in seq_len(nrow(major))) {
    if (!status$eligible[i]) next
    wave <- major[i, ]; weekly <- filter(weekly_cases, state == wave$state)
    onset_index <- match(wave$onset_week, weekly$week_start)
    target <- onset_index + seq_len(window) - 1L
    if (is.na(onset_index) || max(target) > nrow(weekly)) {
      status$fit_status[i] <- "not_fittable"; status$fit_reason[i] <- "window is outside available weekly surveillance"; next
    }
    lambda <- vapply(target, renewal_lambda, numeric(1), cases = weekly$reported_cases, weights = settings$weights)
    if (anyNA(lambda) || any(lambda <= 0)) {
      status$fit_status[i] <- "not_fittable"; status$fit_reason[i] <- "pre-onset renewal infectiousness is unavailable or zero"; next
    }
    kept <- c(kept, i)
    matrices[[length(matrices) + 1L]] <- list(I_obs = as.integer(weekly$reported_cases[target]), Lambda = lambda)
  }
  if (!length(kept)) return(list(status = status, kept = kept, data = NULL, windows = list()))
  list(status = status, kept = kept,
       data = list(K = length(kept), W = window, I_obs = do.call(rbind, lapply(matrices, `[[`, "I_obs")), Lambda = do.call(rbind, lapply(matrices, `[[`, "Lambda"))),
       windows = matrices)
}

fit_window_batch <- function(model, batch, window, settings) {
  if (is.null(batch$data)) return(list(status = "no_fittable_waves"))
  started <- Sys.time()
  fit <- tryCatch(
    sampling(model, data = batch$data, chains = settings$chains, iter = settings$iter, warmup = settings$warmup,
             seed = 20260925L + window, control = list(adapt_delta = settings$adapt_delta, max_treedepth = settings$max_treedepth), refresh = 0),
    error = function(e) e
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (inherits(fit, "error")) return(list(status = "error", reason = conditionMessage(fit), elapsed_seconds = elapsed))
  summary <- as.data.frame(rstan::summary(fit)$summary) |> tibble::rownames_to_column("parameter")
  hmc_summary <- summary |> filter(grepl("^(Re|phi)", parameter))
  global_hmc <- hmc_gate(hmc_summary, rstan::get_sampler_params(fit, inc_warmup = FALSE), settings$max_treedepth, fit, pars = c("Re", "phi"))
  raw_re <- rstan::extract(fit, pars = "Re", permuted = FALSE, inc_warmup = FALSE)
  tail_table <- posterior::summarise_draws(posterior::as_draws_array(raw_re), posterior::ess_tail)
  list(status = "fitted", fit = fit, summary = summary, global_hmc = global_hmc,
       tail_ess = tail_table[[2]], elapsed_seconds = elapsed)
}

empty_window_columns <- function(window, status, reason = NA_character_) {
  prefix <- paste0("Re_early_", window)
  tibble(
    !!paste0(prefix, "_q025") := NA_real_, !!paste0(prefix, "_median") := NA_real_, !!paste0(prefix, "_q975") := NA_real_,
    !!paste0("re", window, "_fit_status") := status, !!paste0("re", window, "_fit_reason") := reason,
    !!paste0("re", window, "_fit_pass") := FALSE, !!paste0("re", window, "_divergences") := NA_real_,
    !!paste0("re", window, "_treedepth_hits") := NA_real_, !!paste0("re", window, "_rhat") := NA_real_,
    !!paste0("re", window, "_bulk_ess") := NA_real_, !!paste0("re", window, "_tail_ess") := NA_real_,
    !!paste0("re", window, "_min_bfmi") := NA_real_, !!paste0("re", window, "_bfmi_by_chain") := NA_character_,
    !!paste0("re", window, "_elapsed_seconds") := NA_real_,
    !!paste0("re", window, "_ppc_mean_abs_error") := NA_real_, !!paste0("re", window, "_ppc_all_weeks_95_coverage") := NA_real_
  )
}

append_window_results <- function(output, batch, fitted, window) {
  bind_cols(output, bind_rows(lapply(seq_len(nrow(output)), function(i) {
    current <- batch$status[i, ]; position <- match(i, batch$kept)
    if (is.na(position)) return(empty_window_columns(window, current$fit_status, current$fit_reason))
    if (!identical(fitted$status, "fitted")) return(empty_window_columns(window, "error", fitted$reason))
    re_row <- filter(fitted$summary, parameter == paste0("Re[", position, "]"))
    re_draws <- rstan::extract(fitted$fit, pars = "Re", permuted = TRUE)$Re[, position]
    rep_draws <- rstan::extract(fitted$fit, pars = "I_rep", permuted = TRUE)$I_rep[, position, ]
    observed <- batch$windows[[position]]$I_obs
    coverage <- mean(vapply(seq_len(ncol(rep_draws)), function(j) {
      interval <- quantile(rep_draws[, j], c(.025, .975)); observed[j] >= interval[1] && observed[j] <= interval[2]
    }, logical(1)))
    pass <- isTRUE(fitted$global_hmc$hmc_pass) && re_row$Rhat <= 1.01 && re_row$n_eff >= 400 && fitted$tail_ess[position] >= 400
    prefix <- paste0("Re_early_", window)
    bind_cols(summarise_draws(re_draws, prefix), tibble(
      !!paste0("re", window, "_fit_status") := "fitted", !!paste0("re", window, "_fit_reason") := NA_character_,
      !!paste0("re", window, "_fit_pass") := pass, !!paste0("re", window, "_divergences") := fitted$global_hmc$divergences,
      !!paste0("re", window, "_treedepth_hits") := fitted$global_hmc$treedepth_hits, !!paste0("re", window, "_rhat") := re_row$Rhat,
      !!paste0("re", window, "_bulk_ess") := re_row$n_eff, !!paste0("re", window, "_tail_ess") := fitted$tail_ess[position],
      !!paste0("re", window, "_min_bfmi") := fitted$global_hmc$min_bfmi, !!paste0("re", window, "_bfmi_by_chain") := fitted$global_hmc$bfmi_by_chain,
      !!paste0("re", window, "_elapsed_seconds") := fitted$elapsed_seconds,
      !!paste0("re", window, "_ppc_mean_abs_error") := mean(abs(colMeans(rep_draws) - observed)),
      !!paste0("re", window, "_ppc_all_weeks_95_coverage") := coverage
    ))
  })))
}

write_early_re_posterior_draws <- function(fit_bundle, paths) {
  # Compact, analysis-facing posterior draw file.  The complete stanfit
  # objects (including posterior-predictive I_rep) remain in the companion
  # fit bundle under 02_Script/stan/national_wave_analysis.
  draws <- bind_rows(lapply(names(fit_bundle$fits), function(key) {
    entry <- fit_bundle$fits[[key]]
    if (!identical(entry$fitted$status, "fitted")) return(tibble())
    re_matrix <- rstan::extract(entry$fitted$fit, pars = "Re", permuted = TRUE)$Re
    wave_lookup <- entry$batch$status[entry$batch$kept, c("state", "wave_id")]
    colnames(re_matrix) <- paste0("re_index_", seq_len(ncol(re_matrix)))
    as_tibble(re_matrix, .name_repair = "minimal") |>
      mutate(draw = row_number()) |>
      pivot_longer(-draw, names_to = "re_index", values_to = "Re") |>
      mutate(re_index = as.integer(sub("re_index_", "", re_index)),
             state = wave_lookup$state[re_index], wave_id = wave_lookup$wave_id[re_index],
             window_weeks = entry$batch$data$W) |>
      select(state, wave_id, window_weeks, draw, Re)
  }))
  saveRDS(draws, file.path(paths$table, "brazil_chik_major_wave_early_re_posterior_draws.rds"), compress = "xz")
}

run_national_early_re <- function() {
  paths <- ensure_national_output_dirs()
  major <- read_national_wave_census() |> filter(major_epidemic_primary) |> arrange(state, onset_week)
  weekly_cases <- read_national_state_week()
  model <- rstan::stan_model(national_path("02_Script", "stan", "renewal_ceara_episode_re.stan"))
  output <- major; fit_bundle <- list(config = early_re_settings, gi_weights = early_re_settings$weights, fits = list())
  for (window in early_re_settings$windows) {
    message(sprintf("[early-Re] %d-week batch", window))
    batch <- build_window_batch(major, weekly_cases, window, early_re_settings)
    fitted <- fit_window_batch(model, batch, window, early_re_settings)
    output <- append_window_results(output, batch, fitted, window)
    fit_bundle$fits[[paste0("window_", window)]] <- list(batch = batch, fitted = fitted)
  }
  output <- output |> group_by(state) |>
    mutate(weeks_since_previous_wave = as.numeric(onset_week - lag(end_week)), previous_wave_total_cases = lag(total_cases),
           weeks_since_previous_wave = if_else(wave_order == 1L, NA_real_, weeks_since_previous_wave),
           previous_wave_total_cases = if_else(wave_order == 1L, NA_real_, previous_wave_total_cases)) |>
    ungroup() |> mutate(preceding_trough_ratio = trough_ratio, re6_analysis_usable = re6_prefit_eligible & re6_fit_pass,
                         calendar_year = year(onset_week))
  write_csv(output, file.path(paths$table, "brazil_chik_major_wave_early_re.csv"))
  saveRDS(fit_bundle, file.path(paths$fit, "brazil_chik_major_wave_early_re_fits.rds"))
  write_early_re_posterior_draws(fit_bundle, paths)
  write_csv(tibble(gi_distribution = "Gamma(shape = 4, rate = 2)", gi_mean_weeks = 2, gi_sd_weeks = 1,
                   discretization = "weekly interval probabilities", truncation_weeks = 8,
                   weights = paste(signif(early_re_settings$weights, 7), collapse = ";"),
                   primary_window_weeks = 6L, sensitivity_windows_weeks = "4;8",
                   hmc_gate = "0 batch divergences; 0 batch max-treedepth hits; Re-specific Rhat <= 1.01; Re-specific bulk/tail ESS >= 400; batch BFMI >= 0.30."),
            file.path(paths$table, "brazil_chik_early_re_metadata.csv"))
  invisible(output)
}

if (sys.nframe() == 0L) run_national_early_re()
