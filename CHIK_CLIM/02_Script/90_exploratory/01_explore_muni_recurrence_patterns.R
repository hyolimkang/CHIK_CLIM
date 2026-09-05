# ---------------------------------------------------------------------------
# explore_muni_recurrence_patterns.R
#
# Goal : Descriptive exploration of municipality-level chikungunya dynamics:
#        recurrent vs sporadic reporting patterns, maps, and climate contrast.
#
# Inputs:
#   01_Data/chik_brazil_muni_month_2015_2024.rds
#
# Outputs (03_Output/figures/ and 03_Output/tables/):
#   chik_recurrence_map.png              municipality pattern class map
#   chik_recurrence_metrics_scatter.png    2D classification space
#   chik_example_muni_curves.png           exemplar time series (recurrent vs sporadic)
#   chik_climate_recurrent_vs_sporadic.png climate boxplots by pattern class
#   chik_muni_recurrence_metrics.csv     per-municipality metrics table
# ---------------------------------------------------------------------------

for (p in c("here", "dplyr", "tidyr", "ggplot2", "sf", "geobr",
            "lubridate", "scales", "patchwork", "purrr", "ggh4x", "geofacet")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(ggplot2)
  library(sf);   library(geobr); library(lubridate)
  library(scales); library(patchwork); library(purrr);library(ggh4x); library(stringr)
  library(geofacet)
})

# ---- 1. Config ------------------------------------------------------------

YEAR_START <- 2015
YEAR_END   <- 2024
N_YEARS    <- YEAR_END - YEAR_START + 1L
N_MONTHS   <- N_YEARS * 12L

DATA_DIR  <- here::here("01_Data")
FIG_DIR   <- here::here("03_Output", "figures")
TAB_DIR   <- here::here("03_Output", "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

PANEL_PATH <- file.path(DATA_DIR,
                        sprintf("chik_brazil_muni_month_%d_%d.rds",
                                YEAR_START, YEAR_END))
CASE_VAR <- "cases_notified"

# Minimum monthly cases to count as "active" (reduces single-notification noise)
ACTIVE_MONTH_MIN <- 1L
# Year is an "outbreak year" if annual cases exceed this fraction of the muni's
# own peak-year total (or at least ACTIVE_MONTH_MIN * 3 months worth)
OUTBREAK_YEAR_FRAC <- 0.15

# ---- 2. Load panel --------------------------------------------------------

if (!file.exists(PANEL_PATH))
  stop("Panel not found: ", PANEL_PATH,
       "\nRun 02_Script/00_data_prep/02_clean_chik_sinan_brazil.R first.")

panel <- readRDS(PANEL_PATH) %>%
  rename(cases = !!sym(CASE_VAR)) %>%
  mutate(
    year  = year(year_month),
    month = month(year_month)
  )

message(sprintf("[load] %d rows, %d municipalities, %d months",
                nrow(panel), n_distinct(panel$muni6), n_distinct(panel$year_month)))

# # ---- 2. State / region lookup -------------------------------------------
# IBGE: first 2 digits of muni6 = state (UF) code.

uf_lookup <- tibble::tribble(
  ~uf_code, ~uf, ~uf_name,           ~region,
  "11",   "RO", "Rondônia",          "Norte",
  "12",   "AC", "Acre",              "Norte",
  "13",   "AM", "Amazonas",          "Norte",
  "14",   "RR", "Roraima",           "Norte",
  "15",   "PA", "Pará",              "Norte",
  "16",   "AP", "Amapá",             "Norte",
  "17",   "TO", "Tocantins",         "Norte",
  "21",   "MA", "Maranhão",          "Nordeste",
  "22",   "PI", "Piauí",             "Nordeste",
  "23",   "CE", "Ceará",             "Nordeste",
  "24",   "RN", "Rio Grande do Norte","Nordeste",
  "25",   "PB", "Paraíba",           "Nordeste",
  "26",   "PE", "Pernambuco",        "Nordeste",
  "27",   "AL", "Alagoas",           "Nordeste",
  "28",   "SE", "Sergipe",           "Nordeste",
  "29",   "BA", "Bahia",             "Nordeste",
  "31",   "MG", "Minas Gerais",      "Sudeste",
  "32",   "ES", "Espírito Santo",    "Sudeste",
  "33",   "RJ", "Rio de Janeiro",    "Sudeste",
  "35",   "SP", "São Paulo",         "Sudeste",
  "41",   "PR", "Paraná",            "Sul",
  "42",   "SC", "Santa Catarina",    "Sul",
  "43",   "RS", "Rio Grande do Sul", "Sul",
  "50",   "MS", "Mato Grosso do Sul","Centro-Oeste",
  "51",   "MT", "Mato Grosso",       "Centro-Oeste",
  "52",   "GO", "Goiás",             "Centro-Oeste",
  "53",   "DF", "Distrito Federal",  "Centro-Oeste"
)

# Order: by region, then alphabetic UF inside each region.
uf_lookup <- uf_lookup %>%
  mutate(region = factor(region,
                         levels = c("Norte","Nordeste","Sudeste",
                                    "Sul","Centro-Oeste"))) %>%
  arrange(region, uf) %>%
  mutate(uf = factor(uf, levels = uf))

region_palette <- c(
  "Norte"        = "#1b9e77",
  "Nordeste"     = "#d95f02",
  "Sudeste"      = "#7570b3",
  "Sul"          = "#e7298a",
  "Centro-Oeste" = "#66a61e"
)


# ---- 3. Explore panel data -------------------------------
glimpse(panel)
n_distinct(panel$muni6) # 5152
n_distinct(panel$year_month)
mean(panel$cases_confirmed > 0) 

# 3-1. national
national <- panel |>
  group_by(year_month) |>
  summarise(cases = sum(cases_confirmed), .groups = "drop") |>
  ggplot(aes(year_month, cases)) +
  geom_col(fill = "steelblue")


# 3-2. state
state <- panel |>
  mutate(uf_code = substr(muni6, 1, 2)) |>
  group_by(year_month, uf_code) |>
  summarise(cases = sum(cases_confirmed), .groups = "drop") |>
  left_join(uf_lookup, by = "uf_code")  |>
  mutate(uf_name = factor(uf_name, levels = uf_lookup$uf_name))

p_state <- state |>
  ggplot(aes(year_month, cases)) +
  geom_col(fill = "steelblue") + 
  facet_nested_wrap(~ region + uf_name, scales = "free_y", ncol = 5)+
  theme_bw()

## observations: 
# Two peak outbreaks (early and late): Pernambuco: Mato Grosso, Ceara, Alagoas, Paraiba, Pernambuco, Piaui
# Very sporadic (one huge or moderate): Roraima (early), Parana (late), Rio Grande do Sul (late)
# Small number of cases (recurrent): Santa Catarina / Rondonia / Amazonas / Acre / Distrito Fedral
# Early year one peak: Para, Maranhao
# Later year one peak: Tocantins
# Endemic/ regular:  Rio Grande de norte, Bahia
# Early outbreak + interepidemic + later endemic: Sergipe
# Later year outbreaks: Sao Paulo/Goias/Minas Gerais/Espirito Santo

# 4. phenotype classification
source("02_Script/00_data_prep/03_fetch_ibge_population.R")
pop <- readRDS("01_Data/ibge_pop_muni_year_2015_2024.rds")

state_pop <- pop |>
  group_by(uf, year) |>
  summarise(population = sum(population), .groups = "drop")

state <- state |>
  mutate(year = lubridate::year(year_month)) |>
  left_join(state_pop, by = c("uf", "year"))

state <- state |>
  mutate(
    incidence_per_100k = cases / population * 100000
  )

ggplot(state, aes(year_month, incidence_per_100k)) +
  geom_col(fill = "steelblue") +
  facet_nested_wrap(~ region + uf_name, scales = "free_y", ncol = 5)+
  theme_bw()

MONTHLY_OB_CUTOFF <- 10 # monthly incidence per 100k
ANNUAL_OB_CUTOFF <- 100 # annual incidence per 100k

state2 <- state |>
  mutate(
    year = year(year_month),
    month = month(year_month),
    time_num = as.numeric(year_month)
  )

# annual incidence table
state_annual <- state |>
  mutate(year = lubridate::year(year_month)) |>
  group_by(uf, uf_name, region, year) |>
  summarise(
    annual_cases = sum(cases, na.rm = TRUE),
    population = mean(population, na.rm = TRUE),
    annual_inc = annual_cases / population * 1e5,
    .groups = "drop"
  ) |>
  mutate(
    outbreak_year = annual_inc >= ANNUAL_OB_CUTOFF
  )

annual_metrics <- state_annual |>
  group_by(uf) |>
  summarise(
    n_outbreak_years = sum(outbreak_year, na.rm = TRUE),
    first_outbreak_year = ifelse(
      any(outbreak_year, na.rm = TRUE),
      min(year[outbreak_year], na.rm = TRUE),
      NA_integer_
    ),
    last_outbreak_year = ifelse(
      any(outbreak_year, na.rm = TRUE),
      max(year[outbreak_year], na.rm = TRUE),
      NA_integer_
    ),
    year_of_max_inc = year[which.max(annual_inc)],
    max_annual_inc = max(annual_inc, na.rm = TRUE),
    total_cases_10y = sum(annual_cases, na.rm = TRUE),
    mean_population = mean(population, na.rm = TRUE),
    cumulative_inc_10y = total_cases_10y / mean_population * 1e5,
    mean_annual_inc = mean(annual_inc, na.rm = TRUE),
    .groups = "drop"
  )


monthly_metrics <- state |>
  mutate(
    time_num = as.numeric(year_month)
  ) |>
  group_by(uf, uf_name, region) |>
  summarise(
    max_monthly_inc = max(incidence_per_100k, na.rm = TRUE),
    mean_monthly_inc = mean(incidence_per_100k, na.rm = TRUE),
    median_monthly_inc = median(incidence_per_100k, na.rm = TRUE),
    
    n_outbreak_months = sum(incidence_per_100k >= MONTHLY_OB_CUTOFF, na.rm = TRUE),
    prop_outbreak_months = mean(incidence_per_100k >= MONTHLY_OB_CUTOFF, na.rm = TRUE),
    prop_zero_months = mean(cases == 0, na.rm = TRUE),
    
    peak_month = year_month[which.max(incidence_per_100k)],
    peak_to_median_ratio = max_monthly_inc / (median_monthly_inc + 0.01),
    
    cog_num = ifelse(
      sum(incidence_per_100k, na.rm = TRUE) > 0,
      sum(time_num * incidence_per_100k, na.rm = TRUE) /
        sum(incidence_per_100k, na.rm = TRUE),
      NA_real_
    ),
    cog_date = as.Date(cog_num, origin = "1970-01-01"),
    .groups = "drop"
  )

state_metrics <- monthly_metrics |>
  left_join(annual_metrics, by = "uf")

state_metrics <- state_metrics |>
  mutate(
    phenotype = case_when(
      n_outbreak_years == 0 ~ "Low detected activity",
      
      n_outbreak_years == 1 & year_of_max_inc <= 2018 ~ 
        "Early single outbreak",
      
      n_outbreak_years == 1 & year_of_max_inc >= 2021 ~ 
        "Late single outbreak",
      
      n_outbreak_years == 1 ~ 
        "Single outbreak, mid-period",
      
      n_outbreak_years == 2 & first_outbreak_year <= 2018 & last_outbreak_year >= 2021 ~ 
        "Early and late two-wave",
      
      n_outbreak_years == 2 & year_of_max_inc <= 2018 ~ 
        "Two-wave, early-dominant",
      
      n_outbreak_years == 2 & year_of_max_inc >= 2021 ~ 
        "Two-wave, late-dominant",
      
      n_outbreak_years == 2 ~ 
        "Two-wave, mixed timing",
      
      n_outbreak_years >= 3 & prop_zero_months < 0.30 ~ 
        "Persistent reported activity",
      
      n_outbreak_years >= 3 ~ 
        "Intermittent recurrent outbreaks",
      
      TRUE ~ "Mixed / uncertain"
    )
  )

state_pheno <- state2 |>
  left_join(
    state_metrics |> select(uf, phenotype),
    by = "uf"
  )

ggplot(state_pheno, aes(year_month, incidence_per_100k)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~ phenotype + uf_name, scales = "free_y", ncol = 5) +
  theme_bw() +
  labs(
    title = "State-level chikungunya incidence by preliminary outbreak phenotype",
    x = NULL,
    y = "Incidence per 100,000"
  )

# Sensitivity grid
thresholds <- expand.grid(
  annual_cutoff = c(50, 100, 200, 500),
  monthly_cutoff = c(5, 10, 20)
)


# metric 
cluster_features <- state_metrics |>
  select(uf,
         max_monthly_inc, mean_monthly_inc,
         n_outbreak_months, prop_zero_months,
         n_outbreak_years, cumulative_inc_10y,
         peak_to_median_ratio) |>
  tibble::column_to_rownames("uf") |>
  scale()

# 2. Hierarchical clustering
hc <- hclust(dist(cluster_features), method = "ward.D2")
plot(hc, hang = -1, cex = 0.7)

# 3. Cut to k clusters
state_metrics$cluster <- as.factor(cutree(hc, k = 4))
# 4. Compare rule-based phenotype vs data-driven cluster
comparison_table <- table(
  rule_based = state_metrics$phenotype,
  cluster    = state_metrics$cluster
)
print(comparison_table)

# 5. Characterize each cluster
state_metrics |>
  group_by(cluster) |>
  summarise(
    states = paste(sort(uf), collapse = ", "),
    mean_burden = round(mean(cumulative_inc_10y), 0),
    n = n()
  ) |>
  arrange(cluster)

state_metrics <- state_metrics |>
  mutate(
    cluster_name = case_when(
      cluster == "1" ~ "D: Low / sparse",
      cluster == "2" ~ "B: Recurrent high burden",
      cluster == "3" ~ "A: Hyperendemic outlier (CE)",
      cluster == "4" ~ "C: Sporadic single outbreak",
      TRUE ~ NA_character_
    ),
    cluster_name = factor(cluster_name, levels = c(
      "A: Hyperendemic outlier (CE)",
      "B: Recurrent high burden",
      "C: Sporadic single outbreak",
      "D: Low / sparse"
    )),
    
    # first/last outbreak year NA 처리
    first_outbreak_year_clean = case_when(
      is.na(first_outbreak_year) ~ "No outbreak detected",
      TRUE ~ as.character(first_outbreak_year)
    ),
    last_outbreak_year_clean = case_when(
      is.na(last_outbreak_year) ~ "No outbreak detected",
      TRUE ~ as.character(last_outbreak_year)
    )
  )

state_metrics |>
  select(uf, cluster, cluster_name, cumulative_inc_10y,
         n_outbreak_years, first_outbreak_year_clean) |>
  arrange(cluster_name) |>
  print(n = 27)


### state profile
state_profile <- state_metrics |>
  mutate(intro_year = first_outbreak_year) |>
  select(uf, uf_name, region, phenotype, cluster_name,
         intro_year, cumulative_inc_10y, n_outbreak_years,
         mean_population)


states_sf <- geobr::read_state(year = 2020, showProgress = FALSE)
extract_state_mean <- function(raster_path, var_name) {
  r <- terra::rast(raster_path)
  v <- terra::extract(r, states_sf, fun = mean, na.rm = TRUE)[, 2]
  tibble::tibble(
    abbrev_state = states_sf$abbrev_state,
    !!var_name := as.numeric(v)
  )
}

climate_uf <- extract_state_mean(
  "01_Data/final/PRCP_TerraClim_2010_2020_005dg_masked_.tif",
  "PRCP"
) |>
  left_join(
    extract_state_mean(
      "01_Data/final/Dengue_temperature_suitaiblity_masked_.tif",
      "Tsuit"
    ),
    by = "abbrev_state"
  ) |>
  left_join(
    extract_state_mean(
      "01_Data/final/Tmax_TerraClim_2010_2020_005dg_masked_.tif",
      "Tmax"
    ),
    by = "abbrev_state"
  ) |>
  left_join(
    extract_state_mean(
      "01_Data/final/Tmean_TerraClim_2010_2020_005dg_masked_.tif",
      "Tmean"
    ),
    by = "abbrev_state"
  ) |>
  left_join(
    extract_state_mean(
      "01_Data/final/Tmin_TerraClim_2010_2020_005dg_masked_.tif",
      "Tmin"
    ),
    by = "abbrev_state"
  ) |>
  rename(uf = abbrev_state)

state_climate <- state_profile |>
  left_join(climate_uf, by = "uf")
state_climate |> filter(is.na(PRCP))

ggplot(state_climate, aes(cluster_name, Tsuit)) +
  geom_boxplot(aes(fill = cluster_name), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(colour = region), width = 0.15, size = 3) +
  ggrepel::geom_text_repel(aes(label = uf, colour = region), 
                           size = 2.8, max.overlaps = 20) +
  scale_fill_brewer(palette = "Set2") +
  scale_colour_manual(values = region_palette) +
  labs(title = "Thermal suitability by outbreak phenotype cluster",
       x = NULL, y = "Mean Tsuit (state)",
       fill = "Cluster", colour = "Region") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggplot(state_climate, aes(cluster_name, Tmean)) +
  geom_boxplot(aes(fill = cluster_name), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(colour = region), width = 0.15, size = 3) +
  ggrepel::geom_text_repel(aes(label = uf, colour = region), 
                           size = 2.8, max.overlaps = 20) +
  scale_fill_brewer(palette = "Set2") +
  scale_colour_manual(values = region_palette) +
  labs(title = "Mean temp by outbreak phenotype cluster",
       x = NULL, y = "Mean temp (state)",
       fill = "Cluster", colour = "Region") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggplot(state_climate, aes(cluster_name, Tmin)) +
  geom_boxplot(aes(fill = cluster_name), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(colour = region), width = 0.15, size = 3) +
  ggrepel::geom_text_repel(aes(label = uf, colour = region), 
                           size = 2.8, max.overlaps = 20) +
  scale_fill_brewer(palette = "Set2") +
  scale_colour_manual(values = region_palette) +
  labs(title = "Min temp by outbreak phenotype cluster",
       x = NULL, y = "Min temp (state)",
       fill = "Cluster", colour = "Region") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggplot(state_climate, aes(cluster_name, PRCP)) +
  geom_boxplot(aes(fill = cluster_name), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(colour = region), width = 0.15, size = 3) +
  ggrepel::geom_text_repel(aes(label = uf, colour = region), 
                           size = 2.8, max.overlaps = 20) +
  scale_fill_brewer(palette = "Set2") +
  scale_colour_manual(values = region_palette) +
  labs(title = "PRCP by outbreak phenotype cluster",
       x = NULL, y = "PRCP (state)",
       fill = "Cluster", colour = "Region") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))


ggplot(state_climate, aes(x = PRCP, y = Tsuit, colour = cluster_name)) +
  geom_point(aes(shape = cluster_name), size = 4, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = uf), size = 3, max.overlaps = 20) +
  
  # "Sweet spot" 영역 강조
  annotate("rect", 
           xmin = 600, xmax = 1300, 
           ymin = 2.8e6, ymax = 3.6e6,
           alpha = 0.1, fill = "orange") +
  annotate("text", x = 950, y = 3.65e6, 
           label = "Recurrent zone\n(dry + warm)", 
           size = 3, color = "darkorange") +
  
  scale_colour_brewer(palette = "Set2") +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(
    x = "Mean annual precipitation (mm)",
    y = "Thermal suitability (Tsuit)",
    colour = "Phenotype cluster",
    shape = "Phenotype cluster",
    title = "Chikungunya outbreak phenotype in climate space",
    subtitle = "Orange box = apparent 'sweet spot' for recurrent transmission"
  ) +
  theme_bw()


# Tmin으로 바꾼 2D scatter
ggplot(state_climate, aes(x = PRCP, y = Tmin, 
                          colour = cluster_name, shape = cluster_name)) +
  geom_point(size = 4, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = uf), size = 3, max.overlaps = 20) +
  
  # Thermal threshold line (Aedes aegypti minimum ~15-16°C)
  geom_hline(yintercept = 16, linetype = "dashed", 
             colour = "grey40", linewidth = 0.8) +
  annotate("text", x = 2400, y = 16.3, 
           label = "Ae. aegypti\nthermal threshold (~16°C)", 
           size = 2.8, colour = "grey40", hjust = 1) +
  
  # Recurrent zone
  annotate("rect",
           xmin = 600, xmax = 1300,
           ymin = 17, ymax = 22,
           alpha = 0.1, fill = "orange") +
  annotate("text", x = 950, y = 22.5,
           label = "Recurrent zone\n(dry + warm)",
           size = 3, colour = "darkorange") +
  
  scale_colour_brewer(palette = "Set2") +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(
    x = "Mean annual precipitation (mm)",
    y = "Mean minimum temperature (°C)",
    colour = "Phenotype cluster",
    shape = "Phenotype cluster",
    title = "Chikungunya outbreak phenotype in climate space",
    subtitle = "Dashed line = approximate Aedes aegypti thermal threshold"
  ) +
  theme_bw()


### quadrat
vars_to_test <- c("Tmin", "PRCP", "Tsuit", 
                  "Tmean", "Tmax",
                  "mean_population",  
                  "cumulative_inc_10y") 

auc_results <- sapply(vars_to_test, function(v) {
  if (v %in% names(state_climate)) {
    tryCatch({
      r <- roc(state_climate$is_recurrent, 
               state_climate[[v]], quiet = TRUE)
      auc(r)
    }, error = function(e) NA)
  } else NA
})

tibble(variable = vars_to_test, AUC = round(auc_results, 3)) |>
  arrange(desc(AUC)) |>
  print()

# 그림 1: PRCP × n_outbreak_years
p1 <- ggplot(state_climate,
             aes(x = PRCP, y = n_outbreak_years,
                 colour = cluster_name, shape = cluster_name)) +
  
  # 배경 색칠
  annotate("rect", xmin = -Inf, xmax = 1308,
           ymin = -Inf, ymax = Inf,
           fill = "#FFF3CD", alpha = 0.4) +
  annotate("rect", xmin = 1308, xmax = Inf,
           ymin = -Inf, ymax = Inf,
           fill = "#D1ECF1", alpha = 0.4) +
  
  geom_vline(xintercept = 1308, linetype = "dashed",
             colour = "grey40", linewidth = 0.8) +
  
  # Zone labels
  annotate("text", x = 700, y = 4.3,
           label = "Dry zone (< 1308mm)\nRecurrent endemic",
           size = 3.2, colour = "darkorange", fontface = "bold", hjust = 0) +
  annotate("text", x = 1400, y = 4.3,
           label = "Wet zone (≥ 1308mm)\nLow / sporadic",
           size = 3.2, colour = "steelblue", fontface = "bold", hjust = 0) +
  
  # Sul annotation
  annotate("text", x = 1550, y = 0.3,
           label = "Sul states: thermally\nexcluded (Tmin < 16°C)",
           size = 2.8, colour = "grey40", fontface = "italic") +
  annotate("segment",
           x = 1530, xend = 1480, y = 0.2, yend = 0.05,
           arrow = arrow(length = unit(0.2, "cm")),
           colour = "grey50") +
  
  geom_jitter(aes(size = log10(mean_population)),
              height = 0.1, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = uf), size = 3,
                           max.overlaps = 20, box.padding = 0.4) +
  
  scale_colour_brewer(palette = "Set2") +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_size_continuous(name = "Population",
                        range = c(3, 9),
                        breaks = c(6, 6.5, 7, 7.5),
                        labels = c("1M", "3M", "10M", "30M")) +
  labs(x = "Mean annual precipitation (mm)",
       y = "Number of outbreak years (2015-2024)",
       colour = "Phenotype", shape = "Phenotype",
       title = "Precipitation predicts outbreak recurrence",
       subtitle = "AUC = 0.88, optimal threshold = 1308mm") +
  theme_bw()


x_range <- range(state_climate$PRCP, na.rm = TRUE)
y_range <- range(state_climate$cumulative_inc_10y, na.rm = TRUE)
cat("PRCP range:", x_range, "\n")
cat("Incidence range:", y_range, "\n")

# 수정된 그림 2
p2 <- ggplot(state_climate,
             aes(x = PRCP, y = cumulative_inc_10y,
                 colour = cluster_name, shape = cluster_name)) +
  
  # 배경 색칠 — y 범위를 log scale이라 명시적으로 지정
  annotate("rect",
           xmin = min(state_climate$PRCP, na.rm = TRUE) * 0.9,
           xmax = 1308,
           ymin = min(state_climate$cumulative_inc_10y, na.rm = TRUE) * 0.5,
           ymax = max(state_climate$cumulative_inc_10y, na.rm = TRUE) * 2,
           fill = "#FFF3CD", alpha = 0.4) +
  annotate("rect",
           xmin = 1308,
           xmax = max(state_climate$PRCP, na.rm = TRUE) * 1.1,
           ymin = min(state_climate$cumulative_inc_10y, na.rm = TRUE) * 0.5,
           ymax = max(state_climate$cumulative_inc_10y, na.rm = TRUE) * 2,
           fill = "#D1ECF1", alpha = 0.4) +
  
  geom_vline(xintercept = 1308, linetype = "dashed",
             colour = "grey40", linewidth = 0.8) +
  
  # Zone labels — 데이터 분포 보고 빈 공간에
  annotate("text",
           x = 800,
           y = max(state_climate$cumulative_inc_10y, na.rm = TRUE) * 1.5,
           label = "Dry zone (< 1308mm)\nRecurrent endemic",
           size = 3.2, colour = "darkorange", fontface = "bold") +
  annotate("text",
           x = 2000,
           y = max(state_climate$cumulative_inc_10y, na.rm = TRUE) * 1.5,
           label = "Wet zone (≥ 1308mm)\nLow / sporadic",
           size = 3.2, colour = "steelblue", fontface = "bold") +
  
  # Sul annotation — Sul states 실제 위치 근처
  annotate("text",
           x = 1650,
           y = min(state_climate$cumulative_inc_10y, na.rm = TRUE) * 0.7,
           label = "Sul: Tmin < 16°C\n(thermal exclusion)",
           size = 2.8, colour = "grey40", fontface = "italic") +
  
  geom_point(aes(size = log10(mean_population)), alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = uf), size = 3,
                           max.overlaps = 20, box.padding = 0.4) +
  
  scale_colour_brewer(palette = "Set2") +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_y_log10(labels = scales::comma) +
  scale_size_continuous(name = "Population",
                        range = c(3, 9),
                        breaks = c(6, 6.5, 7, 7.5),
                        labels = c("1M", "3M", "10M", "30M")) +
  labs(x = "Mean annual precipitation (mm)",
       y = "Cumulative incidence per 100,000 (log scale)",
       colour = "Phenotype", shape = "Phenotype",
       title = "Precipitation threshold separates transmission regimes",
       subtitle = "Dashed line = ROC-optimal threshold (1308mm, AUC = 0.88)") +
  theme_bw()

# 두 그림 합치기
p1 / p2 +
  plot_annotation(
    title = "PRCP as primary climate determinant of chikungunya dynamics"
  )

ggsave("03_Output/figures/prcp_transmission_regimes.png",
       width = 10, height = 14, dpi = 300)