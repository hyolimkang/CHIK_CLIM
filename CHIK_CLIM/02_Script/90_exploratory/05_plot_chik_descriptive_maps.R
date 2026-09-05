# ---------------------------------------------------------------------------
# plot_chik_descriptive_maps.R
#
# Goal : First descriptive maps of reported chikungunya cases in Brazilian
#        municipalities, 2015-2024.
#
# Inputs:
#   01_Data/chik_brazil_muni_month_2015_2024.rds   (from clean_chik_sinan_brazil.R)
#
# Outputs (written to 03_Output/figures/):
#   chik_map_yearly_2015_2024.png      yearly cumulative cases, 10-panel facet
#   chik_map_seasonality.png           mean monthly cases by calendar month
#   chik_map_cumulative_total.png      total cases over the full period
#   chik_national_monthly_curve.png    companion national time series
#
# Notes:
#   * We plot RAW notified cases on a log scale. Incidence rates per
#     100k require IBGE population projections per municipality-year,
#     which can be added later (e.g. via the {sidrar} package).
#   * Municipality boundary file comes from {geobr}. It returns 7-digit
#     IBGE codes; SINAN uses 6-digit codes -> we map by
#     code_muni6 = code_muni %/% 10.
# ---------------------------------------------------------------------------

# ---- 0. Packages ----------------------------------------------------------

for (p in c("here", "dplyr", "tidyr", "ggplot2", "sf", "geobr",
            "lubridate", "scales", "patchwork")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(ggplot2)
  library(sf);   library(geobr); library(lubridate)
  library(scales); library(patchwork)
})

# ---- 1. Config ------------------------------------------------------------

YEAR_START <- 2015
YEAR_END   <- 2024

DATA_DIR <- here::here("01_Data")
FIG_DIR  <- here::here("03_Output", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

PANEL_PATH <- file.path(DATA_DIR,
                        sprintf("chik_brazil_muni_month_%d_%d.rds",
                                YEAR_START, YEAR_END))

# Which case definition to plot.
#  - "cases_notified"      = every suspected notification (most permissive)
#  - "cases_confirmed"     = CLASSI_FIN == 13 (chikungunya confirmed)
#  - "cases_lab_confirmed" = above + CRITERIO == 1
CASE_VAR <- "cases_notified"

# ---- 2. Load case panel ---------------------------------------------------

if (!file.exists(PANEL_PATH))
  stop("Panel not found: ", PANEL_PATH,
       "\nRun 02_Script/00_data_prep/02_clean_chik_sinan_brazil.R first.")

panel <- readRDS(PANEL_PATH) %>%
  rename(cases = !!sym(CASE_VAR)) %>%
  mutate(
    year  = lubridate::year(year_month),
    month = lubridate::month(year_month)
  )

message(sprintf("[load]  %s rows, %d muni, %d months, total cases = %s",
                format(nrow(panel), big.mark = ","),
                dplyr::n_distinct(panel$muni6),
                dplyr::n_distinct(panel$year_month),
                format(sum(panel$cases, na.rm = TRUE), big.mark = ",")))

# ---- 3. Load municipality polygons (cached by geobr) ---------------------

message("[geobr] downloading / loading municipality polygons (cached)")
muni_sf <- geobr::read_municipality(year = 2020, simplified = TRUE,
                                    showProgress = FALSE) %>%
  mutate(muni6 = sprintf("%06d", code_muni %/% 10L)) %>%
  select(muni6, code_muni, abbrev_state, geom)

states_sf <- geobr::read_state(year = 2020, simplified = TRUE,
                               showProgress = FALSE)

# Sanity check on the join
join_diag <- panel %>% distinct(muni6) %>%
  mutate(in_sf = muni6 %in% muni_sf$muni6) %>% pull(in_sf) %>% mean()
message(sprintf("[join]  %.1f%% of muni codes in panel match the shapefile",
                100 * join_diag))

# ---- 4. Helper for a clean map theme -------------------------------------

theme_map <- function() {
  theme_void(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(colour = "grey30"),
      legend.position  = "right",
      strip.text       = element_text(face = "bold")
    )
}

# Log-pseudo scale so 0-case municipalities appear in light grey.
case_fill <- function(name = "Cases") {
  scale_fill_viridis_c(
    name        = name,
    option      = "magma",
    direction   = -1,
    trans       = "log1p",
    breaks      = c(0, 10, 100, 1000, 10000),
    labels      = c("0", "10", "100", "1k", "10k"),
    na.value    = "grey90"
  )
}

# ---- 5. Map A: yearly cumulative, faceted --------------------------------

yearly <- panel %>%
  group_by(muni6, year) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop")

map_yearly_sf <- muni_sf %>%
  left_join(yearly, by = "muni6")

p_yearly <- ggplot(map_yearly_sf) +
  geom_sf(aes(fill = cases), colour = NA) +
  geom_sf(data = states_sf, fill = NA, colour = "grey20", linewidth = 0.15) +
  case_fill("Reported\ncases") +
  facet_wrap(~ year, ncol = 4) +
  coord_sf(crs = 4674, datum = NA) +
  labs(
    title    = "Chikungunya reported cases per municipality, Brazil",
    subtitle = sprintf("Yearly cumulative, %d-%d  (variable: %s)",
                       YEAR_START, YEAR_END, CASE_VAR),
    caption  = "Source: SINAN (Ministério da Saúde). Polygons: geobr 2020."
  ) +
  theme_map()

ggsave(file.path(FIG_DIR, "chik_map_yearly_2015_2024.png"),
       p_yearly, width = 12, height = 8.5, dpi = 200, bg = "white")
message("[plot]  ", file.path(FIG_DIR, "chik_map_yearly_2015_2024.png"))

# ---- 6. Map B: seasonality (mean monthly cases by calendar month) --------

month_lbl <- c("Jan","Feb","Mar","Apr","May","Jun",
               "Jul","Aug","Sep","Oct","Nov","Dec")

seasonal <- panel %>%
  group_by(muni6, month) %>%
  summarise(mean_cases = mean(cases, na.rm = TRUE), .groups = "drop") %>%
  mutate(month_lbl = factor(month_lbl[month], levels = month_lbl))

map_season_sf <- muni_sf %>%
  left_join(seasonal, by = "muni6")

p_season <- ggplot(map_season_sf) +
  geom_sf(aes(fill = mean_cases), colour = NA) +
  geom_sf(data = states_sf, fill = NA, colour = "grey20", linewidth = 0.15) +
  scale_fill_viridis_c(
    name      = "Mean monthly\ncases",
    option    = "magma",
    direction = -1,
    trans     = "log1p",
    breaks    = c(0, 1, 10, 100, 1000),
    labels    = c("0", "1", "10", "100", "1k"),
    na.value  = "grey90"
  ) +
  facet_wrap(~ month_lbl, ncol = 4) +
  coord_sf(crs = 4674, datum = NA) +
  labs(
    title    = "Chikungunya seasonality by municipality, Brazil",
    subtitle = sprintf("Mean reported cases per calendar month, %d-%d",
                       YEAR_START, YEAR_END),
    caption  = "Source: SINAN (Ministério da Saúde). Polygons: geobr 2020."
  ) +
  theme_map()

ggsave(file.path(FIG_DIR, "chik_map_seasonality.png"),
       p_season, width = 12, height = 9, dpi = 200, bg = "white")
message("[plot]  ", file.path(FIG_DIR, "chik_map_seasonality.png"))

# ---- 7. Map C: cumulative total ------------------------------------------

cumulative <- panel %>%
  group_by(muni6) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop")

map_cum_sf <- muni_sf %>% left_join(cumulative, by = "muni6")

p_cum <- ggplot(map_cum_sf) +
  geom_sf(aes(fill = cases), colour = NA) +
  geom_sf(data = states_sf, fill = NA, colour = "grey20", linewidth = 0.2) +
  case_fill("Reported\ncases") +
  coord_sf(crs = 4674, datum = NA) +
  labs(
    title    = "Chikungunya cumulative reported cases, Brazil",
    subtitle = sprintf("All months, %d-%d", YEAR_START, YEAR_END),
    caption  = "Source: SINAN (Ministério da Saúde). Polygons: geobr 2020."
  ) +
  theme_map()

ggsave(file.path(FIG_DIR, "chik_map_cumulative_total.png"),
       p_cum, width = 7.5, height = 8, dpi = 200, bg = "white")
message("[plot]  ", file.path(FIG_DIR, "chik_map_cumulative_total.png"))

# ---- 8. Companion: national monthly time series --------------------------

national <- panel %>%
  group_by(year_month) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop")

p_ts <- ggplot(national, aes(x = year_month, y = cases)) +
  geom_line(colour = "firebrick", linewidth = 0.6) +
  geom_area(fill = "firebrick", alpha = 0.15) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "Brazil-wide monthly reported chikungunya cases",
    subtitle = sprintf("%d-%d  (variable: %s)",
                       YEAR_START, YEAR_END, CASE_VAR),
    x = NULL, y = "Cases",
    caption = "Source: SINAN."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "chik_national_monthly_curve.png"),
       p_ts, width = 9, height = 4, dpi = 200, bg = "white")
message("[plot]  ", file.path(FIG_DIR, "chik_national_monthly_curve.png"))

message("\n[done] all figures written to: ", FIG_DIR)
