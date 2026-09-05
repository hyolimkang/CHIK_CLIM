# ---------------------------------------------------------------------------
# plot_chik_state_comparison.R
#
# Goal : Compare Brazil-wide chikungunya reporting pattern with state-
#        level patterns (27 UFs), grouped by macro-region (N/NE/SE/S/CO).
#
# Inputs:
#   01_Data/chik_brazil_muni_month_2015_2024.rds
#
# Outputs (03_Output/figures/):
#   chik_state_curves_facet.png        27 small-multiples + national overlay
#   chik_state_heatmap.png             states x year-month heatmap
#   chik_region_curves_overlay.png     5 region lines + national, log y
#   chik_state_ranking_bar.png         total cases per state, by region
#   chik_state_share_stacked.png       monthly share of national cases per region
# ---------------------------------------------------------------------------

# ---- 0. Packages ----------------------------------------------------------

for (p in c("here", "dplyr", "tidyr", "ggplot2", "lubridate",
            "scales", "forcats", "RColorBrewer", "patchwork")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(ggplot2)
  library(lubridate); library(scales); library(forcats)
  library(RColorBrewer); library(patchwork)
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

# Case definition: cases_notified | cases_confirmed | cases_lab_confirmed
CASE_VAR <- "cases_notified"

# ---- 2. State / region lookup -------------------------------------------
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

# ---- 3. Load panel and roll up to state x month --------------------------

panel <- readRDS(PANEL_PATH) %>%
  rename(cases = !!sym(CASE_VAR)) %>%
  mutate(uf_code = substr(muni6, 1, 2)) %>%
  left_join(uf_lookup, by = "uf_code")

# Sanity
miss <- panel %>% filter(is.na(uf)) %>% nrow()
if (miss > 0)
  warning(sprintf("[warn] %d rows have muni6 with no UF match", miss))

state_month <- panel %>%
  group_by(region, uf, year_month) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop")

region_month <- state_month %>%
  group_by(region, year_month) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop")

national_month <- panel %>%
  group_by(year_month) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop")

message(sprintf("[load]  national total cases = %s",
                format(sum(national_month$cases), big.mark = ",")))

# ---- 4. Plot A: 27 state facets, national curve as background reference --

# Re-scale national curve so it fits visually inside each facet.
nat_max  <- max(national_month$cases, na.rm = TRUE)
nat_norm <- national_month %>%
  mutate(cases_norm = cases / nat_max)

facet_data <- state_month %>%
  group_by(uf) %>%
  mutate(state_max = max(cases, na.rm = TRUE)) %>%
  ungroup()

p_state_facet <- ggplot(facet_data,
                        aes(x = year_month, y = cases)) +
  # National pattern as light-grey background reference, scaled per panel.
  geom_area(
    data = function(d) d %>% distinct(uf, state_max) %>%
      tidyr::crossing(nat_norm) %>%
      mutate(cases = cases_norm * state_max),
    fill = "grey80", alpha = 0.7
  ) +
  geom_line(aes(colour = region), linewidth = 0.55) +
  facet_wrap(~ uf, scales = "free_y", ncol = 5) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = scales::comma) +
  scale_colour_manual(values = region_palette, name = "Region") +
  labs(
    title    = "Monthly reported chikungunya cases by Brazilian state",
    subtitle = sprintf("Coloured line = state; grey area = national pattern (rescaled). %d-%d",
                       YEAR_START, YEAR_END),
    x = NULL, y = "Cases (per-state y-axis)",
    caption  = "Source: SINAN."
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title       = element_text(face = "bold", size = 13),
        strip.text       = element_text(face = "bold"),
        legend.position  = "top",
        panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "chik_state_curves_facet.png"),
       p_state_facet, width = 13, height = 10, dpi = 200, bg = "white")
message("[plot]  ", file.path(FIG_DIR, "chik_state_curves_facet.png"))

# ---- 5. Plot B: state x year-month heatmap -------------------------------

heat <- state_month %>%
  mutate(uf = fct_rev(uf))   # rows: regions go top-down N -> CO

p_heat <- ggplot(heat, aes(x = year_month, y = uf, fill = cases)) +
  geom_tile(colour = NA) +
  scale_fill_viridis_c(
    name      = "Cases",
    option    = "magma",
    direction = -1,
    trans     = "log1p",
    breaks    = c(0, 10, 100, 1000, 10000, 100000),
    labels    = c("0","10","100","1k","10k","100k"),
    na.value  = "grey95"
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
  labs(
    title    = "Chikungunya monthly cases per Brazilian state",
    subtitle = sprintf("State rows ordered N -> NE -> SE -> S -> CO. %d-%d",
                       YEAR_START, YEAR_END),
    x = NULL, y = NULL,
    caption  = "Source: SINAN."
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title       = element_text(face = "bold"),
        panel.grid       = element_blank(),
        axis.ticks       = element_blank(),
        axis.text.y      = element_text(family = "mono", size = 9))

ggsave(file.path(FIG_DIR, "chik_state_heatmap.png"),
       p_heat, width = 11, height = 8, dpi = 200, bg = "white")
message("[plot]  ", file.path(FIG_DIR, "chik_state_heatmap.png"))

# ---- 6. Plot C: region overlay (5 regions + national, log y) -------------

p_region <- ggplot() +
  geom_line(data = region_month,
            aes(x = year_month, y = cases + 1, colour = region),
            linewidth = 0.8) +
  geom_line(data = national_month,
            aes(x = year_month, y = cases + 1),
            colour = "black", linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = max(national_month$year_month), y = max(national_month$cases),
           label = "Brazil (national)", hjust = 1, vjust = -0.5,
           size = 3, fontface = "italic") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_log10(labels = scales::comma) +
  scale_colour_manual(values = region_palette, name = "Region") +
  labs(
    title    = "Chikungunya monthly cases: national vs. macro-regions",
    subtitle = sprintf("Log y-axis. %d-%d", YEAR_START, YEAR_END),
    x = NULL, y = "Cases (+1, log scale)",
    caption  = "Source: SINAN."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")

ggsave(file.path(FIG_DIR, "chik_region_curves_overlay.png"),
       p_region, width = 10, height = 5, dpi = 200, bg = "white")
message("[plot]  ", file.path(FIG_DIR, "chik_region_curves_overlay.png"))

# ---- 7. Plot D: total cases per state, ranked, coloured by region --------

state_total <- state_month %>%
  group_by(region, uf) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(cases)) %>%
  mutate(uf = factor(uf, levels = uf))

p_rank <- ggplot(state_total,
                 aes(x = cases, y = fct_rev(uf), fill = region)) +
  geom_col() +
  scale_x_continuous(labels = scales::comma) +
  scale_fill_manual(values = region_palette, name = "Region") +
  labs(
    title    = "Total reported chikungunya cases per state",
    subtitle = sprintf("Cumulative, %d-%d", YEAR_START, YEAR_END),
    x = "Cases", y = NULL,
    caption  = "Source: SINAN."
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank())

ggsave(file.path(FIG_DIR, "chik_state_ranking_bar.png"),
       p_rank, width = 8, height = 7, dpi = 200, bg = "white")
message("[plot]  ", file.path(FIG_DIR, "chik_state_ranking_bar.png"))

# ---- 8. Plot E: monthly share of national cases per region ---------------

share <- region_month %>%
  group_by(year_month) %>%
  mutate(share = ifelse(sum(cases) == 0, NA_real_, cases / sum(cases))) %>%
  ungroup()

p_share <- ggplot(share,
                  aes(x = year_month, y = share, fill = region)) +
  geom_area(position = "stack", alpha = 0.9) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = region_palette, name = "Region") +
  labs(
    title    = "Region's share of Brazil's monthly chikungunya cases",
    subtitle = sprintf("Each month sums to 100%%. %d-%d", YEAR_START, YEAR_END),
    x = NULL, y = "Share of national cases",
    caption  = "Source: SINAN."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")

ggsave(file.path(FIG_DIR, "chik_state_share_stacked.png"),
       p_share, width = 10, height = 4.5, dpi = 200, bg = "white")
message("[plot]  ", file.path(FIG_DIR, "chik_state_share_stacked.png"))

message("\n[done] all figures written to: ", FIG_DIR)
