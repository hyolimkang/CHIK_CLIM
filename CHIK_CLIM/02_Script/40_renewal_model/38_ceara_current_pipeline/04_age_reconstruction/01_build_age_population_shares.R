# Age-structured vaccine extension -- Phase 1 output: CE age-population
# SHARES (not absolute counts) by calendar year, single-year ages 0-100.
# Shares are self-normalised (sum to 1 within each year) so they can be
# applied to the ALREADY-FITTED v4.9 weekly N(t) without importing the
# age file's own (slightly different vintage) population totals. 2024
# shares are held constant for 2025 (age file does not cover 2025).

required_packages <- c("here", "dplyr", "readr", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tidyr) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
table_dir <- file.path(root, "03_Output/tables/climate_forced_v4_9")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

CE_UF_CODE <- 23L
MAX_AGE <- 100L

age_pop <- readRDS(file.path(root, "01_Data/ibge_pop_uf_single_age_expanded_2015_2024.rds")) |>
  filter(uf_code == CE_UF_CODE) |>
  arrange(year, age)
stopifnot(all(age_pop$age == rep(0:MAX_AGE, length(unique(age_pop$year)))))

shares_by_year <- age_pop |> group_by(year) |>
  mutate(share = population / sum(population)) |>
  ungroup() |> select(year, age, share)
stopifnot(all(abs(tapply(shares_by_year$share, shares_by_year$year, sum) - 1) < 1e-10))

# Hold 2024 shares constant for 2025 (last available year -- data ends 2024).
shares_2025 <- shares_by_year |> filter(year == max(year)) |> mutate(year = 2025L)
shares_full <- bind_rows(shares_by_year, shares_2025) |> arrange(year, age)

message(sprintf("CE age-population shares: years %d-%d, ages 0-%d, %d rows", min(shares_full$year), max(shares_full$year), MAX_AGE, nrow(shares_full)))
message("2025 shares = 2024 shares (held constant; age file does not cover 2025)")
print(as.data.frame(shares_full |> filter(year %in% c(2015, 2020, 2024, 2025), age %in% c(0, 1, 12, 50, 100))), digits = 4)

write_csv(shares_full, file.path(table_dir, "CE_age_population_shares.csv"))
message("[saved] ", file.path(table_dir, "CE_age_population_shares.csv"))
