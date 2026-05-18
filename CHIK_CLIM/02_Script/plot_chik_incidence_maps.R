# ===========================================================================
# plot_chik_incidence_maps.R
#
# 목적
# ----
# raw cases가 아니라 **인구 보정된 incidence rate (per 100,000)** 기준으로
# chikungunya burden을 다시 본다. raw cases map은 인구 많은 도시(SP,RJ 등)이
# 항상 진해 보이기 때문에, 진짜 risk pattern을 보려면 incidence를 봐야 한다.
#
# 입력
# ----
#   01_Data/chik_brazil_muni_month_2015_2024.rds   (clean_chik_sinan_brazil.R 결과)
#   01_Data/ibge_pop_muni_year_2015_2024.rds       (fetch_ibge_population.R 결과)
#
# 출력 (03_Output/figures/)
#   chik_map_incidence_yearly.png       연도별 incidence rate map (10-panel)
#   chik_map_incidence_cumulative.png   전체 기간 누적 incidence (1 panel)
#   chik_incidence_region_overlay.png   region별 incidence 시계열
#   chik_incidence_state_heatmap.png    state × month heatmap (incidence)
# ===========================================================================

# ---- 0. Packages ---------------------------------------------------------

for (p in c("here", "dplyr", "tidyr", "ggplot2", "sf", "geobr",
            "lubridate", "scales", "forcats")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(ggplot2)
  library(sf); library(geobr); library(lubridate)
  library(scales); library(forcats)
})

# ---- 1. Config -----------------------------------------------------------

YEAR_START <- 2015
YEAR_END   <- 2024

DATA_DIR <- here::here("01_Data")
FIG_DIR  <- here::here("03_Output", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

CASES_PATH <- file.path(DATA_DIR,
                        sprintf("chik_brazil_muni_month_%d_%d.rds",
                                YEAR_START, YEAR_END))
POP_PATH   <- file.path(DATA_DIR,
                        sprintf("ibge_pop_muni_year_%d_%d.rds",
                                YEAR_START, YEAR_END))

# 케이스 정의 (cases_notified | cases_confirmed | cases_lab_confirmed)
CASE_VAR <- "cases_notified"

# ---- 2. UF / region lookup -----------------------------------------------

uf_lookup <- tibble::tribble(
  ~uf_code, ~uf, ~region,
  "11","RO","Norte","12","AC","Norte","13","AM","Norte",
  "14","RR","Norte","15","PA","Norte","16","AP","Norte","17","TO","Norte",
  "21","MA","Nordeste","22","PI","Nordeste","23","CE","Nordeste",
  "24","RN","Nordeste","25","PB","Nordeste","26","PE","Nordeste",
  "27","AL","Nordeste","28","SE","Nordeste","29","BA","Nordeste",
  "31","MG","Sudeste","32","ES","Sudeste","33","RJ","Sudeste","35","SP","Sudeste",
  "41","PR","Sul","42","SC","Sul","43","RS","Sul",
  "50","MS","Centro-Oeste","51","MT","Centro-Oeste",
  "52","GO","Centro-Oeste","53","DF","Centro-Oeste"
) %>%
  mutate(region = factor(region,
                         levels = c("Norte","Nordeste","Sudeste",
                                    "Sul","Centro-Oeste")))

region_palette <- c(
  "Norte"="#1b9e77","Nordeste"="#d95f02","Sudeste"="#7570b3",
  "Sul"="#e7298a","Centro-Oeste"="#66a61e"
)

# ---- 3. Load case panel + population panel -------------------------------

stopifnot(file.exists(CASES_PATH), file.exists(POP_PATH))
panel <- readRDS(CASES_PATH) %>%
  rename(cases = !!sym(CASE_VAR))
# --- NSE 해설 ---------------------------------------------------------
# `!!sym(CASE_VAR)`는 "지금 CASE_VAR라는 R 변수에 들어있는 문자열
# ('cases_notified')을 컬럼 이름으로 풀어 써라"라는 뜻.
# dplyr는 보통 컬럼 이름을 그냥 적으면 받는데(`rename(cases = cases_notified)`),
# 문자열 변수에서 컬럼 이름을 가져오고 싶을 때 sym() + !! 가 필요함.
# 즉 위 코드는 본질적으로:
#   rename(cases = cases_notified)
# 과 같다.

pop <- readRDS(POP_PATH)
message(sprintf("[load]  cases panel %s rows | pop panel %s rows",
                format(nrow(panel), big.mark = ","),
                format(nrow(pop),   big.mark = ",")))

# ---- 4. Cases + population join → incidence ------------------------------
# muni × year 단위로 인구를 합쳐서, 월 패널에 join한다.

panel2 <- panel %>%
  mutate(year = lubridate::year(year_month)) %>%
  left_join(pop %>% select(muni6, year, population),
            by = c("muni6", "year")) %>%
  mutate(
    # incidence = cases per 100,000 (월별)
    inc_per100k = ifelse(is.na(population) | population == 0,
                         NA_real_, cases / population * 1e5)
  )

# 인구를 못 찾은 muni 비율 점검 (대개 1% 미만이어야 정상)
miss_pop <- panel2 %>%
  distinct(muni6, .keep_all = TRUE) %>%
  summarise(n_total  = n(),
            n_no_pop = sum(is.na(population)),
            pct      = 100 * n_no_pop / n_total)
message(sprintf("[diag]  munis with no IBGE pop match: %d/%d (%.2f%%)",
                miss_pop$n_no_pop, miss_pop$n_total, miss_pop$pct))

# ---- 5. Polygons ---------------------------------------------------------

message("[geobr] loading municipality + state polygons (cached)")
muni_sf <- geobr::read_municipality(year = 2020, simplified = TRUE,
                                    showProgress = FALSE) %>%
  mutate(muni6 = sprintf("%06d", code_muni %/% 10L)) %>%
  select(muni6, geom)

states_sf <- geobr::read_state(year = 2020, simplified = TRUE,
                               showProgress = FALSE)

# ---- 6. 공통 theme / color scale -----------------------------------------

theme_map <- function() {
  theme_void(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(colour = "grey30"),
          legend.position = "right",
          strip.text    = element_text(face = "bold"))
}

# log1p scale은 0 case도 보여주면서 100 vs. 10000을 시각적으로 구분.
incidence_fill <- function(name = "Incidence\n(per 100k)") {
  scale_fill_viridis_c(
    name      = name,
    option    = "magma",
    direction = -1,
    trans     = "log1p",
    breaks    = c(0, 1, 10, 100, 1000, 10000),
    labels    = c("0","1","10","100","1k","10k"),
    na.value  = "grey90"
  )
}

# ---- 7. Map A : 연도별 cumulative incidence rate -------------------------

# 각 (muni × year) 의 cases와 그 해 인구로 incidence 계산.
yearly_inc <- panel2 %>%
  group_by(muni6, year) %>%
  summarise(cases     = sum(cases, na.rm = TRUE),
            population = mean(population, na.rm = TRUE),  # 그 해 인구
            .groups   = "drop") %>%
  mutate(inc_per100k = ifelse(is.na(population) | population == 0,
                              NA_real_, cases / population * 1e5))

# muni_sf에는 5,570개 muni가 다 있지만 yearly_inc에는 그 해 case도 0이고
# population도 NA인 muni가 있을 수 있다. left_join 후 year가 NA인 행은
# facet에서 "NA"라는 panel을 만들어버리니까 미리 제외한다.
map_year_sf <- muni_sf %>%
  left_join(yearly_inc, by = "muni6") %>%
  filter(!is.na(year))

p_year <- ggplot(map_year_sf) +
  geom_sf(aes(fill = inc_per100k), colour = NA) +
  geom_sf(data = states_sf, fill = NA, colour = "grey20", linewidth = 0.15) +
  incidence_fill() +
  facet_wrap(~ year, ncol = 4) +
  coord_sf(crs = 4674, datum = NA) +
  labs(
    title    = "Chikungunya annual incidence per municipality (per 100,000)",
    subtitle = sprintf("Yearly, %d-%d  (case variable: %s)",
                       YEAR_START, YEAR_END, CASE_VAR),
    caption  = "Source: SINAN cases / IBGE estimated population. Polygons: geobr 2020."
  ) +
  theme_map()

ggsave(file.path(FIG_DIR, "chik_map_incidence_yearly.png"),
       p_year, width = 12, height = 8.5, dpi = 200, bg = "white")
message("[plot]  ", file.path(FIG_DIR, "chik_map_incidence_yearly.png"))

# ---- 8. Map B : 전체 기간 cumulative incidence ---------------------------

cum_inc <- panel2 %>%
  group_by(muni6) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE),
            # 인구는 매년 변하므로 평균값을 person-years 근사로 사용:
            # cumulative incidence ≈ total_cases / mean(pop) / n_years * 1e5
            mean_pop   = mean(population, na.rm = TRUE),
            n_years    = dplyr::n_distinct(year),
            .groups = "drop") %>%
  mutate(inc_per100k_per_yr = ifelse(is.na(mean_pop) | mean_pop == 0,
                                     NA_real_,
                                     total_cases / mean_pop * 1e5 / n_years))

map_cum_sf <- muni_sf %>% left_join(cum_inc, by = "muni6")

p_cum <- ggplot(map_cum_sf) +
  geom_sf(aes(fill = inc_per100k_per_yr), colour = NA) +
  geom_sf(data = states_sf, fill = NA, colour = "grey20", linewidth = 0.2) +
  incidence_fill("Avg annual\nincidence\n(per 100k)") +
  coord_sf(crs = 4674, datum = NA) +
  labs(
    title    = "Chikungunya average annual incidence, Brazil",
    subtitle = sprintf("Mean per-year incidence, %d-%d", YEAR_START, YEAR_END),
    caption  = "Source: SINAN / IBGE."
  ) +
  theme_map()

ggsave(file.path(FIG_DIR, "chik_map_incidence_cumulative.png"),
       p_cum, width = 7.5, height = 8, dpi = 200, bg = "white")
message("[plot]  ", file.path(FIG_DIR, "chik_map_incidence_cumulative.png"))

# ---- 9. Region 시계열 (incidence) ----------------------------------------

panel_uf <- panel2 %>%
  mutate(uf_code = substr(muni6, 1, 2)) %>%
  left_join(uf_lookup, by = "uf_code")

region_month <- panel_uf %>%
  filter(!is.na(region), !is.na(population)) %>%
  group_by(region, year_month) %>%
  summarise(cases      = sum(cases, na.rm = TRUE),
            population = sum(population, na.rm = TRUE),
            .groups    = "drop") %>%
  mutate(inc_per100k = cases / population * 1e5)

p_region <- ggplot(region_month,
                   aes(x = year_month, y = inc_per100k, colour = region)) +
  geom_line(linewidth = 0.8) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_log10(labels = scales::comma) +
  scale_colour_manual(values = region_palette, name = "Region") +
  labs(
    title    = "Monthly chikungunya incidence by macro-region (per 100k)",
    subtitle = sprintf("Log y-axis. %d-%d", YEAR_START, YEAR_END),
    x = NULL, y = "Cases per 100,000 (log scale)",
    caption  = "Source: SINAN / IBGE."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")

ggsave(file.path(FIG_DIR, "chik_incidence_region_overlay.png"),
       p_region, width = 10, height = 5, dpi = 200, bg = "white")
message("[plot]  ", file.path(FIG_DIR, "chik_incidence_region_overlay.png"))

# ---- 10. State x month heatmap (incidence) -------------------------------

state_month_inc <- panel_uf %>%
  filter(!is.na(region), !is.na(population)) %>%
  group_by(region, uf, year_month) %>%
  summarise(cases      = sum(cases, na.rm = TRUE),
            population = sum(population, na.rm = TRUE),
            .groups    = "drop") %>%
  mutate(inc_per100k = cases / population * 1e5)

# state ordering: region 순서 → alphabetic UF
uf_order <- uf_lookup %>% arrange(region, uf) %>% pull(uf)
state_month_inc <- state_month_inc %>%
  mutate(uf = factor(uf, levels = uf_order))

p_heat <- ggplot(state_month_inc,
                 aes(x = year_month, y = forcats::fct_rev(uf),
                     fill = inc_per100k)) +
  geom_tile() +
  scale_fill_viridis_c(
    name      = "Incidence\n(per 100k)",
    option    = "magma",
    direction = -1,
    trans     = "log1p",
    breaks    = c(0, 1, 10, 100, 1000),
    labels    = c("0","1","10","100","1k"),
    na.value  = "grey95"
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
  labs(
    title    = "Chikungunya monthly incidence per state (per 100k)",
    subtitle = sprintf("State rows N -> NE -> SE -> S -> CO. %d-%d",
                       YEAR_START, YEAR_END),
    x = NULL, y = NULL,
    caption  = "Source: SINAN / IBGE."
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(family = "mono", size = 9))

ggsave(file.path(FIG_DIR, "chik_incidence_state_heatmap.png"),
       p_heat, width = 11, height = 8, dpi = 200, bg = "white")
message("[plot]  ", file.path(FIG_DIR, "chik_incidence_state_heatmap.png"))

message("\n[done] figures written to: ", FIG_DIR)
