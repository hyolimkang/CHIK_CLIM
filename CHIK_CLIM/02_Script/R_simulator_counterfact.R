#-------------------------------------------------------------------------------
#
# purpose:
# reproduce the deterministic Stan dynamics in R and check the R simulator
# against the Stan generated quantities.
#
# Supports BOTH models:
#   - v12e : births only, fractions on observed pop[t]
#   - v12f : births in + natural deaths out, modelled population pop_model
# The script auto-detects which model was fit by looking for v12f-only outputs
# (pop_model_out / deaths_out) in the posterior draws.
#
# run this after fit_chik_ceara_stan_weekly.R has completed in the same R session
#
# required objects
# - fit_v12e
# - ce_fit
# - years
# - GAMMA_WEEK
#
# This file does not fit a model and does not do vaccination.
# It only checks whether the R simulator reproduces Stan generated quantities.
#-------------------------------------------------------------------------------

required_objects <- c("fit_v12e", "ce_fit", "years", "GAMMA_WEEK")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ", paste(missing_objects, collapse = ", "),
    "\nRun source('02_Script/fit_chik_ceara_stan_weekly.R') first."
  )
}

draws_df <- as.data.frame(fit_v12e)

# Detect demographic-flow model (v12f) by its extra generated quantities.
is_v12f <- any(grepl("^pop_model_out\\[", names(draws_df))) ||
  any(grepl("^deaths_out\\[", names(draws_df)))

# Natural mortality rate used by v12f. Prefer the value actually passed to Stan.
if (exists("stan_data_v12e") && !is.null(stan_data_v12e$mu_death_annual)) {
  mu_death_annual <- stan_data_v12e$mu_death_annual
} else if (exists("BRAZIL_ANNUAL_DEATH_RATE")) {
  mu_death_annual <- BRAZIL_ANNUAL_DEATH_RATE
} else {
  mu_death_annual <- 0.0066
}
p_death <- if (is_v12f) 1 - exp(-mu_death_annual / 52) else 0

message(sprintf(
  "Replaying %s dynamics (p_death = %.6g per week).",
  if (is_v12f) "v12f" else "v12e", p_death
))

# posterior median trajectory
get_med_vector <- function(var) {
  cols <- grep(paste0("^", var, "\\["), names(draws_df), value = TRUE)
  if (length(cols) == 0) {
    stop("No columns found for variable: ", var)
  }
  # Make sure columns are ordered as var[1], var[2], ...
  idx <- as.integer(gsub(paste0(var, "\\[|\\]"), "", cols))
  cols <- cols[order(idx)]
  apply(draws_df[, cols, drop = FALSE], 2, median)
}

get_med_scalar <- function(var) {
  if (!var %in% names(draws_df)) {
    stop("No column found for variable: ", var)
  }
  median(draws_df[[var]])
}


# reporting probability
rho <- get_med_scalar("rho")

# initial susceptible fraction
S0_frac <- get_med_scalar("S0_frac")

# initial infectious count
I0 <- get_med_scalar("I0")

# fitted weekly transmission trajectory
beta_t <- get_med_vector("beta_t")

# annual infectious equivalent import count
import_count_year <- get_med_vector("import_count_year")

N <- nrow(ce_fit)

S_frac <- numeric(N)
I_frac <- numeric(N)
R_frac <- numeric(N)
new_inf_count <- numeric(N)
new_inf_frac <- numeric(N)
mu_cases <- numeric(N)
lambda_out <- numeric(N)
pop_model <- numeric(N)   # modelled population (= observed pop[t] under v12e)

# Initial state. Modelled population is anchored to the first observed value.
pop_model[1] <- ce_fit$population[1]
S_frac[1] <- S0_frac
I_frac[1] <- I0 / ce_fit$population[1]
R_frac[1] <- max(1 - S_frac[1] - I_frac[1], 0)

new_inf_count[1] <- 1e-12
new_inf_frac[1] <- new_inf_count[1] / ce_fit$population[1]
mu_cases[1] <- 1e-6
lambda_out[1] <- 0

p_rec <- 1 - exp(-GAMMA_WEEK)

for (t in 2:N) {
  y <- match(ce_fit$year[t], years)

  # Under v12f the previous-week denominator is the modelled population; under
  # v12e the modelled population equals the observed population at every step.
  N_prev <- if (is_v12f) pop_model[t - 1] else ce_fit$population[t - 1]

  import_count_week <- import_count_year[y] / 52
  import_frac       <- import_count_week / N_prev

  S_prev_count <- S_frac[t - 1] * N_prev
  I_prev_count <- I_frac[t - 1] * N_prev
  R_prev_count <- R_frac[t - 1] * N_prev

  infectious_pressure <- I_frac[t - 1] + import_frac

  lambda_t <- beta_t[t] * infectious_pressure
  lambda_out[t] <- lambda_t

  p_inf <- 1 - exp(-lambda_t)

  new_inf_count[t] <- S_prev_count * p_inf
  recov_count <- I_prev_count * p_rec

  # Epidemic transitions + births
  S_post <- S_prev_count - new_inf_count[t] + ce_fit$births_weekly[t]
  I_post <- I_prev_count + new_inf_count[t] - recov_count
  R_post <- R_prev_count + recov_count

  if (is_v12f) {
    # Natural mortality removes a constant per-capita fraction from each
    # compartment; the modelled population evolves freely.
    S_new <- S_post * (1 - p_death)
    I_new <- I_post * (1 - p_death)
    R_new <- R_post * (1 - p_death)
    N_new <- S_new + I_new + R_new

    pop_model[t] <- N_new
    S_frac[t] <- S_new / N_new
    I_frac[t] <- I_new / N_new
    R_frac[t] <- R_new / N_new
    new_inf_frac[t] <- new_inf_count[t] / N_new
  } else {
    # v12e: no deaths, fractions are taken on the observed population.
    pop_model[t] <- ce_fit$population[t]
    S_frac[t] <- S_post / ce_fit$population[t]
    I_frac[t] <- I_post / ce_fit$population[t]
    R_frac[t] <- R_post / ce_fit$population[t]
    new_inf_frac[t] <- new_inf_count[t] / ce_fit$population[t]
  }

  mu_cases[t] <- rho * new_inf_count[t] + 1e-6
}

Reff_t <- beta_t * S_frac / p_rec

stan_S_out <- get_med_vector("S_out")
stan_I_out <- get_med_vector("I_out")
stan_new_inf <- get_med_vector("new_inf")
stan_mu_cases <- get_med_vector("mu_cases")
stan_Reff_t <- get_med_vector("Reff_t")

# Express the R compartments on the same population scale as the Stan output:
#   v12f -> modelled population, v12e -> observed population.
r_pop <- if (is_v12f) pop_model else ce_fit$population
r_S_out <- S_frac * r_pop
r_I_out <- I_frac * r_pop
r_R_out <- R_frac * r_pop
r_new_inf <- new_inf_count
r_mu_cases <- mu_cases
r_Reff_t <- Reff_t


compare_replay <- data.frame(
  quantity = c("S_out", "I_out", "new_inf", "mu_cases", "Reff_t"),
  max_abs_diff = c(
    max(abs(stan_S_out - r_S_out)),
    max(abs(stan_I_out - r_I_out)),
    max(abs(stan_new_inf - r_new_inf)),
    max(abs(stan_mu_cases - r_mu_cases)),
    max(abs(stan_Reff_t - r_Reff_t))
  )
)

print(compare_replay)

par(mfrow = c(2, 2))

plot(stan_S_out, type = "l", col = "black",
     ylab = "S_out", main = "Susceptible")
lines(r_S_out, col = "red")
legend("topright", c("Stan median", "R replay"), col = c("black", "red"), lty = 1)

plot(stan_new_inf, type = "l", col = "black",
     ylab = "new infections", main = "Weekly infections")
lines(r_new_inf, col = "red")

plot(stan_mu_cases, type = "l", col = "black",
     ylab = "mu_cases", main = "Expected reported cases")
lines(r_mu_cases, col = "red")

plot(stan_Reff_t, type = "l", col = "black",
     ylab = "Reff", main = "Effective transmission index")
lines(r_Reff_t, col = "red")
abline(h = 1, lty = 2)

par(mfrow = c(1, 1))


#-------------------------------------------------------------------------------
# Internal mass balance of the R replay
# (S + I + R) divided by the modelled population should be 1 at every week.
#-------------------------------------------------------------------------------
mass_count <- r_S_out + r_I_out + r_R_out
mass_ratio <- mass_count / r_pop
mass_gap <- mass_count - r_pop

mass_check <- data.frame(
  min_mass_ratio = min(mass_ratio, na.rm = TRUE),
  max_mass_ratio = max(mass_ratio, na.rm = TRUE),
  max_abs_mass_error = max(abs(mass_ratio - 1), na.rm = TRUE),
  min_mass_gap = min(mass_gap, na.rm = TRUE),
  max_mass_gap = max(mass_gap, na.rm = TRUE)
)

print(mass_check)

plot(
  ce_fit$week_start,
  mass_ratio,
  type = "l",
  xlab = "",
  ylab = "S + I + R / modelled population",
  main = "Internal mass balance of R replay"
)
abline(h = 1, lty = 2)


#-------------------------------------------------------------------------------
# Demographic drift: modelled population vs observed population.
# Under v12e these are identical by construction. Under v12f the modelled
# population evolves by births - natural deaths, so a non-zero drift is expected
# and shows how far the rough mortality rate pulls the model from the data.
#-------------------------------------------------------------------------------
demo_excess <- pop_model - ce_fit$population
demo_excess_frac <- demo_excess / ce_fit$population

demo_check <- data.frame(
  max_excess_people = max(demo_excess, na.rm = TRUE),
  final_excess_people = tail(demo_excess, 1),
  max_excess_frac = max(demo_excess_frac, na.rm = TRUE),
  final_excess_frac = tail(demo_excess_frac, 1)
)

print(demo_check)

plot(
  ce_fit$week_start,
  demo_excess_frac,
  type = "l",
  xlab = "",
  ylab = "(modelled - observed) / observed population",
  main = "Demographic drift of modelled population"
)
abline(h = 0, lty = 2)
