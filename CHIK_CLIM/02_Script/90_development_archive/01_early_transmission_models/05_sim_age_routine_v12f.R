# ---------------------------------------------------------------------------
# sim_age_routine_v12f.R
#
# Age-structured historical counterfactual simulator built on the fitted v12f
# transmission trajectory (births in + natural deaths out). This is a scenario
# engine, not a calibration model.
#
# Run after fit_chik_ceara_stan_weekly.R (with USE_DEMOGRAPHIC_FLOW <- TRUE)
# has completed in the same R session.
#
# Demography (matches v12f at the aggregate level)
#   - newborns enter the susceptible age-0 compartment (births_weekly)
#   - a constant per-capita natural mortality (mu_death_annual) removes people
#     from every age and compartment each week
#   - the modelled population therefore evolves freely; it is not pinned to pop[t]
#
# Vaccination
#   - routine_schedule gives the coverage applied to the target-age cohort in
#     each listed calendar year (years not listed get no routine vaccination)
#   - vaccine reduces susceptibility by ve_sus
#
# Objects created:
#   - baseline_sim, routine_sim, cf_weekly, cf_annual
#   - p_cf_weekly_infections, p_cf_annual_averted, p_cf_susceptible_fraction
# ---------------------------------------------------------------------------

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "posterior", "scales", "tibble")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(posterior)
  library(scales)
  library(tibble)
})

required_objects <- c("fit_v12e", "ce_fit", "years", "UF_CODE", "UF_NAME")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", "),
    "\nRun source('02_Script/10_sir_transmission_model/02_fit_chik_ceara_stan_weekly.R') first."
  )
}

# ---- Scenario settings -----------------------------------------------------
target_age <- 12L

# Edit this table to choose which years receive routine vaccination and the
# coverage used in each year. Years not listed receive no routine vaccination.
routine_schedule <- tibble::tibble(
  calendar_year = c(2015L, 2016L, 2017L),
  coverage = c(0.50, 0.50, 0.50)
)

ve_sus <- 0.989
max_age <- 100L

# Natural mortality rate (per capita per year). Prefer the value passed to Stan
# so the simulator demography matches the fitted v12f model.
if (exists("stan_data_v12e") && !is.null(stan_data_v12e$mu_death_annual)) {
  mu_death_annual <- stan_data_v12e$mu_death_annual
} else if (exists("BRAZIL_ANNUAL_DEATH_RATE")) {
  mu_death_annual <- BRAZIL_ANNUAL_DEATH_RATE
} else {
  mu_death_annual <- 0.0066
}

# Use posterior median by default. For uncertainty propagation, set
# POSTERIOR_DRAW_ID before sourcing this file.
if (!exists("POSTERIOR_DRAW_ID")) POSTERIOR_DRAW_ID <- NA_integer_

utils::globalVariables(c(
  "calendar_year", "effective_susceptible_frac", "infections_averted",
  "scenario", "susceptible_frac", "t", "total_infections", "coverage",
  "vaccinated_this_week", "vaccination_coverage", "week_start",
  "ce_fit", "years", "mu_death_annual"
))

extract_fit_draws <- function(fit) {
  if ("CmdStanMCMC" %in% class(fit)) {
    as.data.frame(fit$draws(format = "draws_df"))
  } else {
    as.data.frame(fit)
  }
}

extract_index <- function(x) {
  as.integer(sub(".*\\[([0-9]+)\\].*", "\\1", x))
}

get_vector_from_draw <- function(draws_df, var, draw_id = NA_integer_) {
  cols <- grep(paste0("^", var, "\\["), names(draws_df), value = TRUE)
  if (length(cols) == 0) stop("No posterior columns found for variable: ", var)
  cols <- cols[order(extract_index(cols))]

  if (is.na(draw_id)) {
    apply(draws_df[, cols, drop = FALSE], 2, stats::median, na.rm = TRUE)
  } else {
    as.numeric(draws_df[draw_id, cols, drop = TRUE])
  }
}

get_scalar_from_draw <- function(draws_df, var, draw_id = NA_integer_) {
  if (!var %in% names(draws_df)) stop("No posterior column found for variable: ", var)

  if (is.na(draw_id)) {
    stats::median(draws_df[[var]], na.rm = TRUE)
  } else {
    as.numeric(draws_df[[var]][draw_id])
  }
}

make_initial_age_population <- function(total_pop, max_age = 100L) {
  ages <- 0:max_age

  if (exists("age_pop_initial", inherits = TRUE)) {
    age_pop <- get("age_pop_initial", inherits = TRUE)
    if (!all(c("age", "population") %in% names(age_pop))) {
      stop("age_pop_initial must have columns: age, population")
    }
    out <- rep(0, length(ages))
    out[match(age_pop$age, ages, nomatch = 0)] <- age_pop$population[age_pop$age %in% ages]
    return(out / sum(out) * total_pop)
  }

  # Smooth placeholder age pyramid. Replace with IBGE age-specific population
  # when available.
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

validate_routine_schedule <- function(schedule) {
  if (is.null(schedule)) return(NULL)

  required_cols <- c("calendar_year", "coverage")
  missing_cols <- setdiff(required_cols, names(schedule))
  if (length(missing_cols) > 0) {
    stop(
      "routine_schedule must contain columns: ",
      paste(required_cols, collapse = ", "),
      call. = FALSE
    )
  }

  schedule <- schedule |>
    dplyr::mutate(
      calendar_year = as.integer(.data$calendar_year),
      coverage = as.numeric(.data$coverage)
    ) |>
    dplyr::arrange(.data$calendar_year)

  if (any(is.na(schedule$calendar_year)) || any(is.na(schedule$coverage))) {
    stop("routine_schedule contains missing or non-numeric year/coverage values.", call. = FALSE)
  }

  if (any(schedule$coverage < 0 | schedule$coverage > 1)) {
    stop("routine_schedule$coverage must be between 0 and 1.", call. = FALSE)
  }

  duplicated_years <- schedule$calendar_year[duplicated(schedule$calendar_year)]
  if (length(duplicated_years) > 0) {
    stop(
      "routine_schedule has duplicate years: ",
      paste(unique(duplicated_years), collapse = ", "),
      call. = FALSE
    )
  }

  schedule
}

coverage_for_year <- function(schedule, calendar_year) {
  if (is.null(schedule)) return(0)
  year_match <- schedule$calendar_year == calendar_year
  if (!any(year_match)) return(0)
  schedule$coverage[which(year_match)[1]]
}

simulate_age_sir <- function(pars, ce_fit, years, scenario = "baseline",
                             target_age = 12L,
                             routine_schedule = NULL,
                             ve_sus = 0.80,
                             mu_death_annual = 0.0066,
                             max_age = 100L,
                             vaccinate_full_cohort = TRUE,
                             add_births = TRUE) {

  ages <- 0:max_age
  A <- length(ages)
  N_time <- nrow(ce_fit)

  pmax0 <- function(x) pmax(x, 0)
  routine_schedule <- validate_routine_schedule(routine_schedule)

  # Weekly natural mortality probability (matches v12f).
  p_death <- 1 - exp(-mu_death_annual / 52)
  p_rec <- 1 - exp(-pars$gamma)

  pop_age <- make_initial_age_population(ce_fit$population[1], max_age)

  S_u <- pop_age * pars$S0_frac
  S_v <- rep(0, A)

  I_u <- pop_age / sum(pop_age) * pars$I0
  I_v <- rep(0, A)

  R_u <- pmax0(pop_age - S_u - I_u)
  R_v <- rep(0, A)

  out <- vector("list", N_time)

  for (tt in seq_len(N_time)) {

    year_idx <- match(ce_fit$year[tt], years)
    if (is.na(year_idx)) {
      stop("Year not found in years vector: ", ce_fit$year[tt])
    }

    beta_t <- pars$beta_t[tt]

    vaccinated_this_week <- 0
    vaccinated_s_this_week <- 0
    vaccinated_i_this_week <- 0
    vaccinated_r_this_week <- 0

    # ------------------------------------------------------------
    # 1. Vaccination event (first epi week of a scheduled year)
    # ------------------------------------------------------------
    first_week_start_this_year <- min(ce_fit$week_start[ce_fit$year == ce_fit$year[tt]])

    coverage_this_year <- if (scenario == "routine") {
      coverage_for_year(routine_schedule, ce_fit$year[tt])
    } else {
      0
    }

    vaccination_due <- ce_fit$week_start[tt] == first_week_start_this_year &&
      coverage_this_year > 0

    if (vaccination_due) {
      age_idx <- which(ages %in% target_age)
      if (length(age_idx) == 0) stop("target_age is outside the modelled age range.")

      if (vaccinate_full_cohort) {
        vacc_s <- coverage_this_year * S_u[age_idx]
        vacc_i <- coverage_this_year * I_u[age_idx]
        vacc_r <- coverage_this_year * R_u[age_idx]

        S_u[age_idx] <- S_u[age_idx] - vacc_s
        I_u[age_idx] <- I_u[age_idx] - vacc_i
        R_u[age_idx] <- R_u[age_idx] - vacc_r

        S_v[age_idx] <- S_v[age_idx] + vacc_s
        I_v[age_idx] <- I_v[age_idx] + vacc_i
        R_v[age_idx] <- R_v[age_idx] + vacc_r

        vaccinated_s_this_week <- sum(vacc_s)
        vaccinated_i_this_week <- sum(vacc_i)
        vaccinated_r_this_week <- sum(vacc_r)
        vaccinated_this_week <- vaccinated_s_this_week +
          vaccinated_i_this_week + vaccinated_r_this_week
      } else {
        vacc_s <- coverage_this_year * S_u[age_idx]
        S_u[age_idx] <- S_u[age_idx] - vacc_s
        S_v[age_idx] <- S_v[age_idx] + vacc_s
        vaccinated_s_this_week <- sum(vacc_s)
        vaccinated_this_week <- sum(vacc_s)
      }
    }

    # ------------------------------------------------------------
    # 2. Infection and recovery dynamics
    # ------------------------------------------------------------
    N_prev <- sum(S_u) + sum(S_v) + sum(I_u) + sum(I_v) + sum(R_u) + sum(R_v)

    if (tt == 1L) {
      lambda_t <- 0
      new_inf_u <- rep(0, A)
      new_inf_v <- rep(0, A)
    } else {
      import_count_week <- pars$import_count_year[year_idx] / 52

      infectious_pressure <- (sum(I_u) + sum(I_v)) / N_prev +
        import_count_week / N_prev

      lambda_t <- beta_t * infectious_pressure

      p_inf_u <- 1 - exp(-lambda_t)
      p_inf_v <- 1 - exp(-lambda_t * (1 - ve_sus))

      new_inf_u <- S_u * p_inf_u
      new_inf_v <- S_v * p_inf_v

      rec_u <- I_u * p_rec
      rec_v <- I_v * p_rec

      S_u <- S_u - new_inf_u
      S_v <- S_v - new_inf_v

      I_u <- I_u + new_inf_u - rec_u
      I_v <- I_v + new_inf_v - rec_v

      R_u <- R_u + rec_u
      R_v <- R_v + rec_v

      # Births enter susceptible age 0.
      if (add_births) {
        S_u[1] <- S_u[1] + ce_fit$births_weekly[tt]
      }
    }

    # ------------------------------------------------------------
    # 3. Natural mortality (constant per-capita, all ages/compartments)
    #    Matches v12f: total -> (total + births) * (1 - p_death)
    # ------------------------------------------------------------
    deaths_this_week <- 0
    if (tt > 1L) {
      surv <- 1 - p_death
      pre_death_total <- sum(S_u) + sum(S_v) + sum(I_u) + sum(I_v) + sum(R_u) + sum(R_v)

      S_u <- S_u * surv
      S_v <- S_v * surv
      I_u <- I_u * surv
      I_v <- I_v * surv
      R_u <- R_u * surv
      R_v <- R_v * surv

      deaths_this_week <- pre_death_total * p_death
    }

    S_u <- pmax0(S_u); S_v <- pmax0(S_v)
    I_u <- pmax0(I_u); I_v <- pmax0(I_v)
    R_u <- pmax0(R_u); R_v <- pmax0(R_v)

    # ------------------------------------------------------------
    # 4. Diagnostics and output
    # ------------------------------------------------------------
    total_model_pop <- sum(S_u) + sum(S_v) + sum(I_u) + sum(I_v) + sum(R_u) + sum(R_v)

    effective_S <- sum(S_u) + (1 - ve_sus) * sum(S_v)
    susceptible_total <- sum(S_u) + sum(S_v)
    vaccinated_total <- sum(S_v) + sum(I_v) + sum(R_v)

    out[[tt]] <- tibble::tibble(
      t = tt,
      week_start = ce_fit$week_start[tt],
      calendar_year = ce_fit$year[tt],
      week_of_year = ce_fit$week_of_year[tt],
      scenario = scenario,

      beta_t = beta_t,
      beta_eff = beta_t * effective_S / total_model_pop,

      S_unvaccinated = sum(S_u),
      S_vaccinated = sum(S_v),
      I_unvaccinated = sum(I_u),
      I_vaccinated = sum(I_v),
      R_unvaccinated = sum(R_u),
      R_vaccinated = sum(R_v),

      total_population_modelled = total_model_pop,
      total_population_data = ce_fit$population[tt],
      pop_balance_resid = total_model_pop - ce_fit$population[tt],

      new_inf_unvaccinated = sum(new_inf_u),
      new_inf_vaccinated = sum(new_inf_v),
      total_infections = sum(new_inf_u) + sum(new_inf_v),

      deaths_this_week = deaths_this_week,

      susceptible_frac = susceptible_total / total_model_pop,
      effective_susceptible_frac = effective_S / total_model_pop,

      vaccinated_susceptible_frac = sum(S_v) / total_model_pop,
      vaccinated_total_frac = vaccinated_total / total_model_pop,

      vaccinated_this_week = vaccinated_this_week,
      vaccinated_s_this_week = vaccinated_s_this_week,
      vaccinated_i_this_week = vaccinated_i_this_week,
      vaccinated_r_this_week = vaccinated_r_this_week,
      vaccination_coverage = coverage_this_year,

      lambda = lambda_t,
      Reff = beta_t / p_rec * effective_S / total_model_pop
    )

    # ------------------------------------------------------------
    # 5. Aging
    # ------------------------------------------------------------
    if (tt < N_time) {
      S_u <- age_one_week(S_u)
      S_v <- age_one_week(S_v)
      I_u <- age_one_week(I_u)
      I_v <- age_one_week(I_v)
      R_u <- age_one_week(R_u)
      R_v <- age_one_week(R_v)
    }
  }

  dplyr::bind_rows(out)
}

# ---------------------------------------------------------------------------
# Draw-level prevacc / postvacc runner
#
# Step 1 of the staged plan: instead of a single posterior-median run, loop over
# posterior draws and run BOTH scenarios per draw, then summarise with credible
# intervals.
#   - prevacc  = baseline (no vaccination)
#   - postvacc = routine schedule
# DALY / morbidity layers are intentionally left for a later step.
# ---------------------------------------------------------------------------
draws_df <- extract_fit_draws(fit_v12e)

# Per-draw parameter matrices / vectors.
get_draw_matrix <- function(draws_df, var, draw_ids) {
  cols <- grep(paste0("^", var, "\\["), names(draws_df), value = TRUE)
  if (length(cols) == 0) stop("No posterior columns found for variable: ", var)
  cols <- cols[order(extract_index(cols))]
  as.matrix(draws_df[draw_ids, cols, drop = FALSE])
}

get_draw_scalar <- function(draws_df, var, draw_ids) {
  if (!var %in% names(draws_df)) stop("No posterior column found for variable: ", var)
  as.numeric(draws_df[[var]][draw_ids])
}

# Number of posterior draws to simulate. Subsample for speed; set higher (or to
# nrow(draws_df)) for final runs.
if (!exists("N_SIM_DRAWS")) N_SIM_DRAWS <- 200L
if (!exists("CI_LEVEL")) CI_LEVEL <- 0.95
ci_lo <- (1 - CI_LEVEL) / 2
ci_hi <- 1 - ci_lo

n_total_draws <- nrow(draws_df)
n_draws <- min(N_SIM_DRAWS, n_total_draws)
set.seed(123)
draw_ids <- sort(sample.int(n_total_draws, n_draws))

message(sprintf(
  "Draw-level run: %d of %d posterior draws, %.0f%% credible intervals.",
  n_draws, n_total_draws, 100 * CI_LEVEL
))

beta_mat   <- get_draw_matrix(draws_df, "beta_t", draw_ids)             # [n_draws x N]
import_mat <- get_draw_matrix(draws_df, "import_count_year", draw_ids)  # [n_draws x Y]
rho_vec    <- get_draw_scalar(draws_df, "rho", draw_ids)
S0_vec     <- get_draw_scalar(draws_df, "S0_frac", draw_ids)
I0_vec     <- get_draw_scalar(draws_df, "I0", draw_ids)
gamma_val  <- if (exists("GAMMA_WEEK")) GAMMA_WEEK else 1.0

routine_schedule <- validate_routine_schedule(routine_schedule)
schedule_years_not_in_fit <- setdiff(routine_schedule$calendar_year, unique(ce_fit$year))
if (length(schedule_years_not_in_fit) > 0) {
  warning(
    "routine_schedule includes years outside ce_fit and they will be ignored: ",
    paste(schedule_years_not_in_fit, collapse = ", "),
    call. = FALSE
  )
}

N_time <- nrow(ce_fit)

# Draw-wise weekly storage [N_time x n_draws]
prevacc_inf  <- matrix(NA_real_, N_time, n_draws)
postvacc_inf <- matrix(NA_real_, N_time, n_draws)
prevacc_esf  <- matrix(NA_real_, N_time, n_draws)
postvacc_esf <- matrix(NA_real_, N_time, n_draws)

run_one <- function(d, scenario, ce_fit, years, mu_death_annual) {
  pars_d <- list(
    beta_t = beta_mat[d, ],
    import_count_year = import_mat[d, ],
    rho = rho_vec[d],
    S0_frac = S0_vec[d],
    I0 = I0_vec[d],
    gamma = gamma_val
  )
  simulate_age_sir(
    pars = pars_d,
    ce_fit = ce_fit,
    years = years,
    scenario = scenario,
    target_age = target_age,
    routine_schedule = if (scenario == "routine") routine_schedule else NULL,
    ve_sus = ve_sus,
    mu_death_annual = mu_death_annual,
    max_age = max_age
  )
}

for (d in seq_len(n_draws)) {
  base_d <- run_one(d, "baseline", ce_fit, years, mu_death_annual)
  rout_d <- run_one(d, "routine", ce_fit, years, mu_death_annual)

  prevacc_inf[, d]  <- base_d$total_infections
  postvacc_inf[, d] <- rout_d$total_infections
  prevacc_esf[, d]  <- base_d$effective_susceptible_frac
  postvacc_esf[, d] <- rout_d$effective_susceptible_frac

  if (d %% 25 == 0) message(sprintf("  ... %d / %d draws", d, n_draws))
}

# ---- Weekly summaries (median + credible interval) ------------------------
row_q <- function(mat, p) apply(mat, 1, stats::quantile, probs = p, na.rm = TRUE)
row_med <- function(mat) apply(mat, 1, stats::median, na.rm = TRUE)

averted_inf <- prevacc_inf - postvacc_inf   # per draw, per week

weekly_summary <- dplyr::bind_rows(
  tibble::tibble(
    week_start = ce_fit$week_start, scenario = "prevacc",
    median = row_med(prevacc_inf), lwr = row_q(prevacc_inf, ci_lo), upr = row_q(prevacc_inf, ci_hi)
  ),
  tibble::tibble(
    week_start = ce_fit$week_start, scenario = "postvacc",
    median = row_med(postvacc_inf), lwr = row_q(postvacc_inf, ci_lo), upr = row_q(postvacc_inf, ci_hi)
  )
)

averted_weekly <- tibble::tibble(
  week_start = ce_fit$week_start,
  median = row_med(averted_inf),
  lwr = row_q(averted_inf, ci_lo),
  upr = row_q(averted_inf, ci_hi)
)

esf_summary <- dplyr::bind_rows(
  tibble::tibble(
    week_start = ce_fit$week_start, scenario = "prevacc",
    median = row_med(prevacc_esf), lwr = row_q(prevacc_esf, ci_lo), upr = row_q(prevacc_esf, ci_hi)
  ),
  tibble::tibble(
    week_start = ce_fit$week_start, scenario = "postvacc",
    median = row_med(postvacc_esf), lwr = row_q(postvacc_esf, ci_lo), upr = row_q(postvacc_esf, ci_hi)
  )
)

# ---- Annual summaries (aggregate per draw, then quantiles) ----------------
year_vec <- ce_fit$year
year_levels <- sort(unique(year_vec))

annual_from_weekly <- function(weekly_mat) {
  # returns [n_years x n_draws]
  apply(weekly_mat, 2, function(col) tapply(col, year_vec, sum))
}

prevacc_annual  <- annual_from_weekly(prevacc_inf)
postvacc_annual <- annual_from_weekly(postvacc_inf)
averted_annual_mat <- prevacc_annual - postvacc_annual
pct_reduction_mat  <- averted_annual_mat / pmax(prevacc_annual, 1e-9)

annual_summary <- tibble::tibble(
  calendar_year = as.integer(rownames(prevacc_annual)),
  prevacc_median  = apply(prevacc_annual, 1, stats::median),
  postvacc_median = apply(postvacc_annual, 1, stats::median),
  averted_median  = apply(averted_annual_mat, 1, stats::median),
  averted_lwr     = apply(averted_annual_mat, 1, stats::quantile, probs = ci_lo),
  averted_upr     = apply(averted_annual_mat, 1, stats::quantile, probs = ci_hi),
  pct_reduction_median = apply(pct_reduction_mat, 1, stats::median),
  pct_reduction_lwr    = apply(pct_reduction_mat, 1, stats::quantile, probs = ci_lo),
  pct_reduction_upr    = apply(pct_reduction_mat, 1, stats::quantile, probs = ci_hi)
)

# Total averted across the whole period, per draw -> credible interval
total_averted_by_draw <- colSums(averted_inf, na.rm = TRUE)
total_averted_summary <- tibble::tibble(
  median = stats::median(total_averted_by_draw),
  lwr = stats::quantile(total_averted_by_draw, ci_lo),
  upr = stats::quantile(total_averted_by_draw, ci_hi)
)

# ---- Knife-edge diagnostic: R0 across draws vs averted effect --------------
# R0_t(draw) = beta_t(draw) / (1 - exp(-gamma)). When the per-draw peak R0
# sits just above 1, vaccination tips R_eff below threshold and averts almost
# everything; this block shows how concentrated the posterior is near R0 = 1
# and whether the averted effect is driven by that proximity.
p_rec_draw <- 1 - exp(-gamma_val)
peak_R0_by_draw <- apply(beta_mat, 1, max) / p_rec_draw   # [n_draws]
med_R0_by_draw  <- apply(beta_mat, 1, stats::median) / p_rec_draw

knife_edge_diag <- tibble::tibble(
  peak_R0_median = stats::median(peak_R0_by_draw),
  peak_R0_lwr = stats::quantile(peak_R0_by_draw, ci_lo),
  peak_R0_upr = stats::quantile(peak_R0_by_draw, ci_hi),
  frac_draws_peak_R0_below_1.2 = mean(peak_R0_by_draw < 1.2),
  cor_peakR0_vs_averted = suppressWarnings(
    stats::cor(peak_R0_by_draw, total_averted_by_draw)
  )
)

# ---- Vaccination event check (single representative draw) -----------------
rep_routine <- run_one(
  which.min(abs(rho_vec - stats::median(rho_vec))), "routine",
  ce_fit, years, mu_death_annual
)
vaccination_check <- rep_routine |>
  dplyr::filter(.data$vaccinated_this_week > 0) |>
  dplyr::group_by(calendar_year) |>
  dplyr::summarise(
    scheduled_coverage = max(vaccination_coverage, na.rm = TRUE),
    total_vaccinated = sum(vaccinated_this_week, na.rm = TRUE),
    vaccination_week = week_start[which.max(vaccinated_this_week)],
    .groups = "drop"
  )

message("\nRoutine vaccination schedule:")
print(routine_schedule)
message("\nVaccination event check by year (representative draw):")
print(vaccination_check)

# ---- Baseline vs Stan posterior infection trajectory (medians) ------------
stan_new_inf_med <- get_vector_from_draw(draws_df, "new_inf", NA_integer_)
baseline_check <- tibble::tibble(
  max_abs_diff = max(abs(row_med(prevacc_inf) - stan_new_inf_med), na.rm = TRUE),
  cor_baseline_stan = stats::cor(row_med(prevacc_inf), stan_new_inf_med),
  baseline_peak_week = ce_fit$week_start[which.max(row_med(prevacc_inf))],
  stan_peak_week = ce_fit$week_start[which.max(stan_new_inf_med)]
)

message("\nBaseline (prevacc) median vs Stan posterior median infections:")
print(baseline_check)

message("\nAnnual infections averted (median [CI]):")
print(as.data.frame(annual_summary), digits = 6)

message("\nTotal infections averted over period (median [CI]):")
print(total_averted_summary)

message("\nKnife-edge diagnostic (peak R0 across draws vs averted effect):")
message("  - If peak_R0 sits close to 1, the averted estimate is fragile.")
message("  - A high positive cor_peakR0_vs_averted means draws with larger R0")
message("    avert MORE (bigger baseline outbreak), confirming threshold leverage.")
print(as.data.frame(knife_edge_diag), digits = 6)

# ---- Plots ----------------------------------------------------------------
routine_schedule_label <- paste(
  sprintf("%d: %.0f%%", routine_schedule$calendar_year, 100 * routine_schedule$coverage),
  collapse = ", "
)
scenario_cols <- c(prevacc = "#08519c", postvacc = "#cb181d")

p_cf_weekly_infections <- ggplot(weekly_summary, aes(x = week_start, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) +
  geom_line(aes(y = median), linewidth = 0.8) +
  scale_colour_manual(values = scenario_cols) +
  scale_fill_manual(values = scenario_cols) +
  scale_y_continuous(labels = scales::label_number()) +
  labs(
    x = NULL, y = "Weekly infections", colour = NULL, fill = NULL,
    title = sprintf("Prevacc vs postvacc weekly infections, UF %s", UF_CODE),
    subtitle = sprintf(
      "Draw-level (%d draws); routine age %d (%s), VE_sus = %.1f%%",
      n_draws, target_age, routine_schedule_label, 100 * ve_sus
    )
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_cf_annual_averted <- ggplot(annual_summary, aes(x = calendar_year, y = averted_median)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_col(fill = "#3182bd") +
  geom_errorbar(aes(ymin = averted_lwr, ymax = averted_upr), width = 0.3) +
  scale_x_continuous(breaks = year_levels) +
  scale_y_continuous(labels = scales::label_number(), expand = expansion(mult = c(0.05, 0.10))) +
  labs(
    x = NULL, y = "Infections averted",
    title = sprintf("Annual infections averted by routine age-%d vaccination", target_age),
    subtitle = sprintf("Median and %.0f%% credible interval across %d draws", 100 * CI_LEVEL, n_draws)
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

p_cf_susceptible_fraction <- ggplot(esf_summary, aes(x = week_start, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) +
  geom_line(aes(y = median), linewidth = 0.8) +
  scale_colour_manual(values = scenario_cols) +
  scale_fill_manual(values = scenario_cols) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA)) +
  labs(
    x = NULL, y = "Effective susceptible fraction", colour = NULL, fill = NULL,
    title = sprintf("Effective susceptible fraction, prevacc vs postvacc, UF %s", UF_CODE),
    subtitle = "Vaccinated susceptibles weighted by residual susceptibility, 1 - VE_sus"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

print(p_cf_weekly_infections)
print(p_cf_annual_averted)
print(p_cf_susceptible_fraction)
