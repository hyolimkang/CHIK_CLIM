# ---------------------------------------------------------------------------
# sim_prevacc_age_infections_v12f.R
#
# Draw-level age-structured prevaccination simulator for infections only.
#
# Run after fit_chik_ceara_stan_weekly.R (with USE_DEMOGRAPHIC_FLOW <- TRUE)
# has completed in the same R session.
#
# Objects created:
#   - prevacc_age_inf_array: [age, week, draw] true infections
#   - prevacc_weekly_totals: [week, draw] total infections across ages
#   - prevacc_age_week_summary: median and UI for each age/week
#   - prevacc_weekly_summary: median and UI for weekly totals
#   - prevacc_draw_summary: draw-level rho/S0/R0/total infection diagnostics
# ---------------------------------------------------------------------------

required_pkgs <- c("dplyr", "posterior", "tibble")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(posterior)
  library(tibble)
})

required_objects <- c("fit_v12e", "ce_fit", "years", "UF_CODE", "UF_NAME")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", "),
    "\nRun source('02_Script/fit_chik_ceara_stan_weekly.R') first.",
    call. = FALSE
  )
}

if (!all(c("week_start", "year", "population", "births_weekly") %in% names(ce_fit))) {
  stop("ce_fit must contain week_start, year, population, and births_weekly.", call. = FALSE)
}

# ---- Settings --------------------------------------------------------------
if (!exists("N_PREVACC_DRAWS")) N_PREVACC_DRAWS <- 200L
if (!exists("PREVACC_CI_LEVEL")) PREVACC_CI_LEVEL <- 0.95
if (!exists("MAX_AGE")) MAX_AGE <- 100L
if (!exists("ADD_BIRTHS")) ADD_BIRTHS <- TRUE

mu_death_annual_prevacc <- if (exists("mu_death_annual")) {
  mu_death_annual
} else if (exists("BRAZIL_ANNUAL_DEATH_RATE")) {
  BRAZIL_ANNUAL_DEATH_RATE
} else {
  0.0066
}

gamma_prevacc <- if (exists("GAMMA_WEEK")) GAMMA_WEEK else 1.0
ci_lo <- (1 - PREVACC_CI_LEVEL) / 2
ci_hi <- 1 - ci_lo

# ---- Helpers ---------------------------------------------------------------
extract_fit_draws <- function(fit) {
  if ("CmdStanMCMC" %in% class(fit)) {
    as.data.frame(fit$draws(format = "draws_df"))
  } else {
    as.data.frame(posterior::as_draws_df(fit))
  }
}

extract_index <- function(x) {
  as.integer(sub(".*\\[([0-9]+)\\].*", "\\1", x))
}

get_draw_matrix <- function(draws_df, var, draw_ids) {
  cols <- grep(paste0("^", var, "\\["), names(draws_df), value = TRUE)
  if (length(cols) == 0) stop("No posterior columns found for variable: ", var, call. = FALSE)
  cols <- cols[order(extract_index(cols))]
  as.matrix(draws_df[draw_ids, cols, drop = FALSE])
}

get_draw_scalar <- function(draws_df, var, draw_ids) {
  if (!var %in% names(draws_df)) stop("No posterior column found for variable: ", var, call. = FALSE)
  as.numeric(draws_df[[var]][draw_ids])
}

make_initial_age_population <- function(total_pop, max_age = 100L) {
  ages <- 0:max_age

  if (exists("age_pop_initial", inherits = TRUE)) {
    age_pop <- get("age_pop_initial", inherits = TRUE)
    if (!all(c("age", "population") %in% names(age_pop))) {
      stop("age_pop_initial must have columns: age, population", call. = FALSE)
    }
    out <- rep(0, length(ages))
    keep <- age_pop$age %in% ages
    out[match(age_pop$age[keep], ages)] <- age_pop$population[keep]
    return(out / sum(out) * total_pop)
  }

  # Placeholder age pyramid; replace by age_pop_initial when IBGE age data exist.
  weights <- exp(-0.018 * ages)
  weights <- weights / sum(weights)
  weights * total_pop
}

age_one_week <- function(x) {
  out <- x
  move <- x / 52
  out <- out - move
  out[-1] <- out[-1] + move[-length(move)]
  out[length(out)] <- out[length(out)] + move[length(move)]
  out
}

simulate_prevacc_age_infections_one_draw <- function(pars, ce_fit, years,
                                                     max_age = 100L,
                                                     mu_death_annual = 0.0066,
                                                     add_births = TRUE) {
  ages <- 0:max_age
  n_age <- length(ages)
  n_time <- nrow(ce_fit)
  pmax0 <- function(x) pmax(x, 0)

  p_death <- 1 - exp(-mu_death_annual / 52)
  p_rec <- 1 - exp(-pars$gamma)

  pop_age <- make_initial_age_population(ce_fit$population[1], max_age)
  S <- pop_age * pars$S0_frac
  I <- pop_age / sum(pop_age) * pars$I0
  R <- pmax0(pop_age - S - I)

  age_inf <- matrix(0, nrow = n_age, ncol = n_time)
  pop_balance_resid <- rep(NA_real_, n_time)
  reff <- rep(NA_real_, n_time)

  for (tt in seq_len(n_time)) {
    year_idx <- match(ce_fit$year[tt], years)
    if (is.na(year_idx)) stop("Year not found in years vector: ", ce_fit$year[tt], call. = FALSE)

    beta_t <- pars$beta_t[tt]

    if (tt > 1L) {
      n_prev <- sum(S) + sum(I) + sum(R)
      import_count_week <- pars$import_count_year[year_idx] / 52
      infectious_pressure <- sum(I) / n_prev + import_count_week / n_prev
      lambda_t <- beta_t * infectious_pressure

      new_inf <- S * (1 - exp(-lambda_t))
      rec <- I * p_rec

      S <- S - new_inf
      I <- I + new_inf - rec
      R <- R + rec

      if (add_births) {
        S[1] <- S[1] + ce_fit$births_weekly[tt]
      }

      age_inf[, tt] <- new_inf
    }

    if (tt > 1L) {
      surv <- 1 - p_death
      S <- S * surv
      I <- I * surv
      R <- R * surv
    }

    S <- pmax0(S)
    I <- pmax0(I)
    R <- pmax0(R)

    total_model_pop <- sum(S) + sum(I) + sum(R)
    pop_balance_resid[tt] <- total_model_pop - ce_fit$population[tt]
    reff[tt] <- beta_t / p_rec * sum(S) / total_model_pop

    if (tt < n_time) {
      S <- age_one_week(S)
      I <- age_one_week(I)
      R <- age_one_week(R)
    }
  }

  list(
    age_infections = age_inf,
    weekly_infections = colSums(age_inf),
    pop_balance_resid = pop_balance_resid,
    reff = reff
  )
}

# ---- Draw-level run --------------------------------------------------------
draws_df <- extract_fit_draws(fit_v12e)
n_total_draws <- nrow(draws_df)
n_draws <- min(N_PREVACC_DRAWS, n_total_draws)
set.seed(123)
draw_ids <- sort(sample.int(n_total_draws, n_draws))

message(sprintf(
  "Prevacc age infection run: %d of %d posterior draws, %.0f%% intervals.",
  n_draws, n_total_draws, 100 * PREVACC_CI_LEVEL
))

beta_mat <- get_draw_matrix(draws_df, "beta_t", draw_ids)
import_mat <- get_draw_matrix(draws_df, "import_count_year", draw_ids)
rho_vec <- get_draw_scalar(draws_df, "rho", draw_ids)
S0_vec <- get_draw_scalar(draws_df, "S0_frac", draw_ids)
I0_vec <- get_draw_scalar(draws_df, "I0", draw_ids)

ages <- 0:MAX_AGE
n_age <- length(ages)
n_time <- nrow(ce_fit)

prevacc_age_inf_array <- array(
  NA_real_,
  dim = c(n_age, n_time, n_draws),
  dimnames = list(
    age = as.character(ages),
    week = as.character(ce_fit$week_start),
    draw = as.character(draw_ids)
  )
)
prevacc_reff <- matrix(NA_real_, nrow = n_time, ncol = n_draws)
prevacc_pop_balance_resid <- matrix(NA_real_, nrow = n_time, ncol = n_draws)

for (d in seq_len(n_draws)) {
  pars_d <- list(
    beta_t = beta_mat[d, ],
    import_count_year = import_mat[d, ],
    rho = rho_vec[d],
    S0_frac = S0_vec[d],
    I0 = I0_vec[d],
    gamma = gamma_prevacc
  )

  sim_d <- simulate_prevacc_age_infections_one_draw(
    pars = pars_d,
    ce_fit = ce_fit,
    years = years,
    max_age = MAX_AGE,
    mu_death_annual = mu_death_annual_prevacc,
    add_births = ADD_BIRTHS
  )

  prevacc_age_inf_array[, , d] <- sim_d$age_infections
  prevacc_reff[, d] <- sim_d$reff
  prevacc_pop_balance_resid[, d] <- sim_d$pop_balance_resid

  if (d %% 25 == 0) message(sprintf("  ... %d / %d draws", d, n_draws))
}

# ---- Summaries -------------------------------------------------------------
prevacc_weekly_totals <- apply(prevacc_age_inf_array, c(2, 3), sum)
prevacc_age_totals <- apply(prevacc_age_inf_array, c(1, 3), sum)

q_over_draws <- function(x, prob) {
  apply(x, seq_len(length(dim(x)) - 1), stats::quantile, probs = prob, na.rm = TRUE)
}

prevacc_age_week_median <- apply(prevacc_age_inf_array, c(1, 2), stats::median, na.rm = TRUE)
prevacc_age_week_lwr <- q_over_draws(prevacc_age_inf_array, ci_lo)
prevacc_age_week_upr <- q_over_draws(prevacc_age_inf_array, ci_hi)

prevacc_age_week_summary <- tibble::tibble(
  age = rep(ages, times = n_time),
  week_start = rep(ce_fit$week_start, each = n_age),
  calendar_year = rep(ce_fit$year, each = n_age),
  median = as.vector(prevacc_age_week_median),
  lwr = as.vector(prevacc_age_week_lwr),
  upr = as.vector(prevacc_age_week_upr)
)

prevacc_weekly_summary <- tibble::tibble(
  week_start = ce_fit$week_start,
  calendar_year = ce_fit$year,
  median = apply(prevacc_weekly_totals, 1, stats::median, na.rm = TRUE),
  lwr = apply(prevacc_weekly_totals, 1, stats::quantile, probs = ci_lo, na.rm = TRUE),
  upr = apply(prevacc_weekly_totals, 1, stats::quantile, probs = ci_hi, na.rm = TRUE)
)

prevacc_age_total_summary <- tibble::tibble(
  age = ages,
  median = apply(prevacc_age_totals, 1, stats::median, na.rm = TRUE),
  lwr = apply(prevacc_age_totals, 1, stats::quantile, probs = ci_lo, na.rm = TRUE),
  upr = apply(prevacc_age_totals, 1, stats::quantile, probs = ci_hi, na.rm = TRUE)
)

p_rec_prevacc <- 1 - exp(-gamma_prevacc)
prevacc_draw_summary <- tibble::tibble(
  draw_id = draw_ids,
  rho = rho_vec,
  S0_frac = S0_vec,
  I0 = I0_vec,
  peak_R0 = apply(beta_mat, 1, max) / p_rec_prevacc,
  median_R0 = apply(beta_mat, 1, stats::median) / p_rec_prevacc,
  total_infections = colSums(prevacc_weekly_totals, na.rm = TRUE),
  max_abs_pop_balance_resid = apply(abs(prevacc_pop_balance_resid), 2, max, na.rm = TRUE)
)

message("\nPrevacc infection totals across period:")
print(
  tibble::tibble(
    median = stats::median(prevacc_draw_summary$total_infections),
    lwr = stats::quantile(prevacc_draw_summary$total_infections, ci_lo),
    upr = stats::quantile(prevacc_draw_summary$total_infections, ci_hi)
  )
)

message("\nPosterior draw diagnostics used by prevacc simulator:")
print(
  as.data.frame(dplyr::summarise(
    prevacc_draw_summary,
    rho_median = stats::median(.data$rho),
    rho_lwr = stats::quantile(.data$rho, ci_lo),
    rho_upr = stats::quantile(.data$rho, ci_hi),
    peak_R0_median = stats::median(.data$peak_R0),
    peak_R0_lwr = stats::quantile(.data$peak_R0, ci_lo),
    peak_R0_upr = stats::quantile(.data$peak_R0, ci_hi)
  )),
  digits = 6
)
