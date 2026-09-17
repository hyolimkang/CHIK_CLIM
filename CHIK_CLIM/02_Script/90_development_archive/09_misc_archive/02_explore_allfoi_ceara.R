# explore_allfoi_ceara.R
#
# Open global FOI grid (allfoi_s1.RData) and summarise foi_mid / foi_lo / foi_hi
# for Brazil and Ceará.
#
# NOTE: allfoi_s1.RData is ~3 GB. Run in a fresh R session with >= 8 GB RAM free.

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

# ---- 1. Paths --------------------------------------------------------------
ALLFOI_PATH <- here("01_Data/allfoi_s1.RData")

if (!file.exists(ALLFOI_PATH)) {
  stop("File not found: ", ALLFOI_PATH,
       "\nExpected: CHIK_CLIM/01_Data/allfoi_s1.RData")
}

cat("File:", ALLFOI_PATH, "\n")
cat("Size:", round(file.info(ALLFOI_PATH)$size / 1e9, 2), "GB\n\n")

# ---- 2. Load allfoi (keep key columns only) --------------------------------
message("Loading allfoi ... (may take 1–3 min)")
load(ALLFOI_PATH)   # object name: allfoi

need_cols <- c("x", "y", "country", "tot", "foi_mid", "foi_lo", "foi_hi")
miss <- setdiff(need_cols, names(allfoi))
if (length(miss) > 0) {
  stop("Missing columns in allfoi: ", paste(miss, collapse = ", "))
}

foi <- allfoi[, need_cols]
rm(allfoi)
gc()

# ---- 3. Quick inspect (global) ---------------------------------------------
cat("\n--- global FOI (foi_mid) ---\n")
print(summary(foi$foi_mid))

# ---- 4. Summarise helper ---------------------------------------------------
summarise_foi_region <- function(df, region_label) {
  pos <- df$foi_mid > 0

  tibble(
    region = region_label,
    n_cells = nrow(df),
    n_cells_foi_gt0 = sum(pos, na.rm = TRUE),
    pop_tot = sum(df$tot, na.rm = TRUE),

    # unweighted grid-cell summaries
    foi_mid_mean   = mean(df$foi_mid, na.rm = TRUE),
    foi_mid_median = median(df$foi_mid, na.rm = TRUE),
    foi_lo_mean    = mean(df$foi_lo, na.rm = TRUE),
    foi_hi_mean    = mean(df$foi_hi, na.rm = TRUE),

    # among cells with foi_mid > 0 (CHIK_VIM state-summary style)
    foi_mid_median_pos = median(df$foi_mid[pos], na.rm = TRUE),
    foi_mid_q025_pos   = as.numeric(stats::quantile(df$foi_mid[pos], 0.025, na.rm = TRUE)),
    foi_mid_q975_pos   = as.numeric(stats::quantile(df$foi_mid[pos], 0.975, na.rm = TRUE)),

    # population-weighted mean FOI (cells with valid pop > 0)
    foi_mid_pop_wmean = {
      ok <- !is.na(df$foi_mid) & !is.na(df$tot) & df$tot > 0
      stats::weighted.mean(df$foi_mid[ok], df$tot[ok])
    },
    foi_lo_pop_wmean = {
      ok <- !is.na(df$foi_lo) & !is.na(df$tot) & df$tot > 0
      stats::weighted.mean(df$foi_lo[ok], df$tot[ok])
    },
    foi_hi_pop_wmean = {
      ok <- !is.na(df$foi_hi) & !is.na(df$tot) & df$tot > 0
      stats::weighted.mean(df$foi_hi[ok], df$tot[ok])
    },

    # expected ever-infected counts (FOI x population per cell, summed)
    expected_inf_mid = sum(df$foi_mid * df$tot, na.rm = TRUE),
    expected_inf_lo  = sum(df$foi_lo * df$tot, na.rm = TRUE),
    expected_inf_hi  = sum(df$foi_hi * df$tot, na.rm = TRUE)
  )
}

# ---- 5. Brazil subset ------------------------------------------------------
bra_foi <- foi |>
  filter(country == "Brazil")

message("Brazil grid cells: ", nrow(bra_foi))

bra_summary <- summarise_foi_region(bra_foi, "Brazil")

# ---- 6. Ceará subset (spatial join) ----------------------------------------
bra_sf <- st_as_sf(bra_foi, coords = c("x", "y"), crs = 4326)

br_states <- ne_states(country = "Brazil", returnclass = "sf")
ce_sf <- br_states |>
  filter(grepl("^Cear", name)) |>
  st_transform(crs = st_crs(bra_sf))

ce_foi_sf <- bra_sf |>
  st_join(ce_sf[, "name"], join = st_within) |>
  filter(!is.na(name))

message("Ceará grid cells: ", nrow(ce_foi_sf))

ce_foi <- st_drop_geometry(ce_foi_sf)
ce_summary <- summarise_foi_region(ce_foi, "Ceará")

foi_summary <- bind_rows(bra_summary, ce_summary)

# ---- 7. Print + save -------------------------------------------------------
cat("\n--- Brazil / Ceará FOI summary ---\n")
print(as.data.frame(foi_summary), row.names = FALSE)

out_dir <- here("03_Output/tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_csv <- file.path(out_dir, "ce_foi_summary_brazil_ceara.csv")
out_rds <- here("01_Data/ce_foi_grid_brazil_ceara.rds")

write_csv(foi_summary, out_csv)
saveRDS(list(bra_foi = bra_foi, ce_foi = ce_foi, summary = foi_summary),
        out_rds)

cat("\nSaved summary: ", out_csv, "\n", sep = "")
cat("Saved grids:   ", out_rds, "\n", sep = "")
