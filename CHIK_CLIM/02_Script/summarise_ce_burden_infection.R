# summarise_ce_burden_infection.R
#
# Sum grid-cell infection burdens (med_inf, lo_inf, hi_inf) for Brazil and
# Ceará from combined_burden (.RData).
#
# Input (first existing path wins):
#   01_Data/combined_burden.RData
#   01_Data/combined_burden_shrink.RData
#
# Output:
#   03_Output/tables/ce_burden_infection_summary.csv
#   03_Output/tables/bra_burden_infection_summary.csv

# ---- 0. Packages -----------------------------------------------------------
for (p in c("here", "dplyr", "sf", "readr",
            "rnaturalearth", "rnaturalearthdata")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(sf)
  library(readr)
  library(rnaturalearth)
  library(rnaturalearthdata)
})

# ---- 1. Load combined_burden -----------------------------------------------
candidates <- c(
  here("01_Data/combined_burden.RData"),
  here("01_Data/combined_burden_shrink.RData")
)
BURDEN_PATH <- candidates[file.exists(candidates)][1]

if (is.na(BURDEN_PATH)) {
  stop("No combined_burden file found in 01_Data/.")
}

message("Loading: ", BURDEN_PATH)
load(BURDEN_PATH)   # object name: combined_burden

need_cols <- c("x", "y", "country", "tot", "med_inf", "lo_inf", "hi_inf")
miss <- setdiff(need_cols, names(combined_burden))
if (length(miss) > 0) {
  stop("Missing columns in combined_burden: ", paste(miss, collapse = ", "))
}

burden <- combined_burden[, need_cols]
rm(combined_burden)
gc()

# ---- 2. Brazil subset ------------------------------------------------------
bra_burden <- burden |>
  filter(country == "Brazil")

message("Brazil grid cells: ", nrow(bra_burden))

bra_summary <- bra_burden |>
  summarise(
    region   = "Brazil",
    n_cells  = n(),
    pop_tot  = sum(tot, na.rm = TRUE),
    med_inf  = sum(med_inf, na.rm = TRUE),
    lo_inf   = sum(lo_inf, na.rm = TRUE),
    hi_inf   = sum(hi_inf, na.rm = TRUE),
    .groups  = "drop"
  )

# ---- 3. Ceará subset (spatial join) ----------------------------------------
bra_sf <- st_as_sf(bra_burden, coords = c("x", "y"), crs = 4326)

br_states <- rnaturalearth::ne_states(country = "Brazil", returnclass = "sf")
ce_sf <- br_states |>
  filter(grepl("^Cear", name))
ce_sf <- st_transform(ce_sf, crs = st_crs(bra_sf))

ce_burden <- bra_sf |>
  st_join(ce_sf[, "name"], join = st_within) |>
  filter(!is.na(name))

message("Ceará grid cells: ", nrow(ce_burden))

ce_summary <- ce_burden |>
  st_drop_geometry() |>
  summarise(
    region   = "Ceará",
    n_cells  = n(),
    pop_tot  = sum(tot, na.rm = TRUE),
    med_inf  = sum(med_inf, na.rm = TRUE),
    lo_inf   = sum(lo_inf, na.rm = TRUE),
    hi_inf   = sum(hi_inf, na.rm = TRUE),
    .groups  = "drop"
  )

# ---- 4. Save ---------------------------------------------------------------
out_dir <- here("03_Output/tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bra_out <- file.path(out_dir, "bra_burden_infection_summary.csv")
ce_out  <- file.path(out_dir, "ce_burden_infection_summary.csv")


annual_summary <- ce_fit %>%
  group_by(year) %>%
  summarise(
    cases = sum(cases),
    .groups = "drop"
  )

mean_reported_cases <- mean(annual_summary$cases)

ext_med <- 75810.62
ext_lo  <- 53048.57
ext_hi  <- 101390.1

rho_implied_med <- mean_reported_cases / ext_med
rho_implied_lo  <- mean_reported_cases / ext_hi
rho_implied_hi  <- mean_reported_cases / ext_lo

c(
  mean_reported_cases = mean_reported_cases,
  rho_implied_lo = rho_implied_lo,
  rho_implied_med = rho_implied_med,
  rho_implied_hi = rho_implied_hi
)



library(dplyr)

ce_fit <- ce_fit |>
  dplyr::arrange(year, week_of_year)

N <- nrow(ce_fit)

years <- sort(unique(ce_fit$year))
Y <- length(years)
year_id <- match(ce_fit$year, years)

W <- 52
week_id <- as.integer(pmin(ce_fit$week_of_year, 52))

# External model-predicted long-term average annual infections
ext_med <- 75810.62
ext_lo  <- 53048.57
ext_hi  <- 101390.10

# Convert 95% UI to log-scale SD: v9 with informative rho prior
ext_log_sd <- (log(ext_hi) - log(ext_lo)) / (2 * 1.96)

stan_data <- list(
  N = N,
  cases = as.integer(ce_fit$cases),
  
  Y = Y,
  year_id = as.integer(year_id),
  
  W = W,
  week_id = as.integer(week_id),
  
  pop = as.vector(ce_fit$population),
  births_weekly = as.vector(ce_fit$births_weekly),
  
  gamma_fixed = GAMMA_WEEK,
  fit_start = 20,
  
  ext_mean_annual_infections = ext_med,
  ext_log_sd = ext_log_sd,
  
  rho_alpha = 8,
  rho_beta  = 23,
  
  s0_alpha = 20,
  s0_beta  = 3
)