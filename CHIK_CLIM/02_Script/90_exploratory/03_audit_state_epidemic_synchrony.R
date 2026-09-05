# ---------------------------------------------------------------------------
# audit_state_epidemic_synchrony.R
#
# Goal : Audit, for all 27 Brazilian states (UFs), (1) when the first major
#        chikungunya epidemic occurred, and (2) whether municipalities
#        inside that state moved together (one synchronous state-wide wave)
#        or asynchronously (many small, staggered municipality-level
#        outbreaks that happen to sum to a state "wave" but are not really
#        one epidemic). Then track how that within-state spatial pattern
#        changed 2015-2024 (e.g. was Bahia's 2022 outbreak one synchronous
#        event, or a patchwork of local outbreaks -- and does that differ
#        from Bahia's earlier waves?).
#
# Input :
#   01_Data/chik_brazil_muni_week_2015_2024.rds
#     -> produced by 02_Script/00_data_prep/02_clean_chik_sinan_brazil.R (section 6,
#        "Weekly" panel). NOTE: the .rds on disk currently carries a few
#        extra comorbidity columns (cases_w_comorb / cases_diabetic / ...)
#        that a newer, trimmed version of clean_chik_sinan_brazil.R no
#        longer writes -- harmless, we simply don't use them here. The
#        columns we rely on (muni6, week_start, epi_year, epi_week,
#        cases_notified, cases_confirmed, ...) are unchanged. It is a
#        complete muni6 x week_start panel (zero-filled), Sunday-start
#        epi-weeks, 2015-01-04 .. 2024-12-29 (plus a partial 2014-12-28 wk).
#   01_Data/ibge_pop_muni_year_2015_2024.rds
#     -> produced by 02_Script/00_data_prep/03_fetch_ibge_population.R. Muni6 x year
#        population, used to turn state-level case counts into incidence.
#
# Method (all thresholds are named constants in Section 1, easy to vary):
#   A. Aggregate muni-week cases to state-week, express as incidence per
#      100k, lightly smooth (3-wk centered mean) to damp reporting noise.
#   B. Detect "epidemic waves" per state = runs of weeks with smoothed
#      incidence >= WAVE_INC_THRESHOLD, short dips bridged, short/small
#      runs dropped. The first surviving wave = that state's first major
#      epidemic.
#   C. For each wave, look back at the underlying municipalities within
#      the wave's time window (+/- a buffer): which munis were "active"
#      (>= ACTIVE_MUNI_MIN_CASES), when did each active muni peak, and how
#      correlated are their weekly curves? Two numbers summarise this:
#        - peak_time_sd_weeks : weighted SD of each active muni's peak
#          week (weighted by muni's case load). Small = everyone peaked
#          together (synchronous). Large = staggered.
#        - cor_sync_index     : mean pairwise correlation of active munis'
#          weekly case curves within the window. High = curves co-move.
#      These feed a simple rule-based label per wave (Section 6).
#
# Outputs:
#   03_Output/tables/chik_state_epidemic_wave_audit.csv   (state x wave)
#   03_Output/tables/chik_state_first_epidemic_summary.csv (state x 1 row)
#   03_Output/figures/chik_state_wave_diagnostic_facets.png  (QA: 27 states)
#   03_Output/figures/chik_state_epidemic_timeline.png       (Gantt overview)
#   03_Output/figures/chik_bahia_muni_week_heatmap.png       (BA deep dive)
#   03_Output/figures/chik_synchrony_trend_over_time.png     (system-wide)
#   03_Output/figures/chik_state_synchrony_by_wave_facet.png (per-state)
#   03_Output/figures/chik_state_first_epidemic_ranked.png   (ranked bar)
# ---------------------------------------------------------------------------

# ---- 0. Packages -----------------------------------------------------------

for (p in c("here", "dplyr", "tidyr", "ggplot2", "lubridate", "scales",
            "forcats", "purrr", "tibble", "stringr", "patchwork", "ggh4x")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(ggplot2)
  library(lubridate); library(scales); library(forcats); library(purrr)
  library(tibble); library(stringr); library(patchwork); library(ggh4x)
})

# geobr/sf are used ONLY to order Bahia's municipalities north->south in the
# deep-dive heatmap. Optional: if unavailable or offline, we fall back to
# ordering by each muni's first active week instead (script still runs).
HAVE_GEOBR <- requireNamespace("geobr", quietly = TRUE) &&
              requireNamespace("sf",    quietly = TRUE)

# ---- 1. Config --------------------------------------------------------------

YEAR_START <- 2015
YEAR_END   <- 2024

DATA_DIR <- here::here("01_Data")
FIG_DIR  <- here::here("03_Output", "figures")
TAB_DIR  <- here::here("03_Output", "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

PANEL_PATH <- file.path(DATA_DIR,
                        sprintf("chik_brazil_muni_week_%d_%d.rds",
                                YEAR_START, YEAR_END))
POP_PATH   <- file.path(DATA_DIR,
                        sprintf("ibge_pop_muni_year_%d_%d.rds",
                                YEAR_START, YEAR_END))

# Case definition: cases_notified | cases_confirmed | cases_lab_confirmed
# Notified is used by default (matches the rest of the repo's convention) --
# it avoids biasing "synchrony" by places/periods with different lab
# confirmation capacity.
CASE_VAR <- "cases_notified"

# --- Wave-detection thresholds (state level) ---
SMOOTH_WINDOW      <- 3    # weeks, centered moving average
WAVE_INC_THRESHOLD <- 1    # smoothed cases per 100k per week to count as "active"
MIN_GAP_WEEKS      <- 3    # bridge dips of up to this many weeks below threshold
MIN_WAVE_DURATION  <- 4    # weeks; shorter runs are dropped as noise
MIN_WAVE_CASES     <- 50   # state-wide cases within the run; drops tiny blips

# --- Municipality-level synchrony thresholds ---
ACTIVE_MUNI_MIN_CASES <- 20  # muni must report >= this many cases within the
                             # wave window to be treated as "active" for the
                             # synchrony calculation (avoids 1-2 case noise)
WAVE_BUFFER_WEEKS     <- 4   # pad this many weeks before/after the state wave
                             # window when looking at muni-level curves, so a
                             # muni that peaks slightly before/after the
                             # state-level window is still captured

message("[config] case_var = ", CASE_VAR,
        " | wave_inc_threshold = ", WAVE_INC_THRESHOLD, "/100k/wk",
        " | active_muni_min_cases = ", ACTIVE_MUNI_MIN_CASES)

# ---- 2. UF / region lookup (same convention as other 02_Script files) ------

uf_lookup <- tibble::tribble(
  ~uf_code, ~uf, ~uf_name,            ~region,
  "11",   "RO", "Rondônia",           "Norte",
  "12",   "AC", "Acre",               "Norte",
  "13",   "AM", "Amazonas",           "Norte",
  "14",   "RR", "Roraima",            "Norte",
  "15",   "PA", "Pará",               "Norte",
  "16",   "AP", "Amapá",              "Norte",
  "17",   "TO", "Tocantins",          "Norte",
  "21",   "MA", "Maranhão",           "Nordeste",
  "22",   "PI", "Piauí",              "Nordeste",
  "23",   "CE", "Ceará",              "Nordeste",
  "24",   "RN", "Rio Grande do Norte","Nordeste",
  "25",   "PB", "Paraíba",            "Nordeste",
  "26",   "PE", "Pernambuco",         "Nordeste",
  "27",   "AL", "Alagoas",            "Nordeste",
  "28",   "SE", "Sergipe",            "Nordeste",
  "29",   "BA", "Bahia",              "Nordeste",
  "31",   "MG", "Minas Gerais",       "Sudeste",
  "32",   "ES", "Espírito Santo",     "Sudeste",
  "33",   "RJ", "Rio de Janeiro",     "Sudeste",
  "35",   "SP", "São Paulo",          "Sudeste",
  "41",   "PR", "Paraná",             "Sul",
  "42",   "SC", "Santa Catarina",     "Sul",
  "43",   "RS", "Rio Grande do Sul",  "Sul",
  "50",   "MS", "Mato Grosso do Sul", "Centro-Oeste",
  "51",   "MT", "Mato Grosso",        "Centro-Oeste",
  "52",   "GO", "Goiás",              "Centro-Oeste",
  "53",   "DF", "Distrito Federal",   "Centro-Oeste"
) %>%
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

class_palette <- c(
  "Synchronous (state-wide)"          = "#d7191c",
  "Mixed / intermediate"              = "#fdae61",
  "Asynchronous (localized/staggered)"= "#2c7bb6",
  "Single-municipality dominated"     = "#7570b3",
  "No major epidemic detected"        = "grey70"
)

# ---- 3. Load data ------------------------------------------------------------

stopifnot(file.exists(PANEL_PATH), file.exists(POP_PATH))

panel_raw <- readRDS(PANEL_PATH) %>%
  as_tibble() %>%
  transmute(muni6, week_start,
            cases = .data[[CASE_VAR]]) %>%
  mutate(uf_code = substr(muni6, 1, 2)) %>%
  left_join(uf_lookup, by = "uf_code") %>%
  filter(!is.na(uf))   # drops handful of foreign/invalid codes, if any

pop <- readRDS(POP_PATH)

state_pop_year <- pop %>%
  group_by(uf, year) %>%
  summarise(population = sum(population, na.rm = TRUE), .groups = "drop")

message(sprintf(
  "[load] %s muni-week rows | %d munis | %d weeks (%s .. %s)",
  format(nrow(panel_raw), big.mark = ","),
  n_distinct(panel_raw$muni6), n_distinct(panel_raw$week_start),
  min(panel_raw$week_start), max(panel_raw$week_start)))

# ---- 4. State-week incidence series (smoothed) -----------------------------

roll_mean_centered <- function(x, k) {
  # simple centered moving average, edges fall back to the raw value
  n <- length(x)
  out <- x
  half <- k %/% 2
  for (i in seq_len(n)) {
    lo <- max(1, i - half); hi <- min(n, i + half)
    out[i] <- mean(x[lo:hi], na.rm = TRUE)
  }
  out
}

state_week <- panel_raw %>%
  group_by(region, uf, uf_name, week_start) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  mutate(year = lubridate::year(week_start)) %>%
  left_join(state_pop_year, by = c("uf", "year")) %>%
  arrange(uf, week_start) %>%
  group_by(uf) %>%
  mutate(
    inc_per100k = ifelse(is.na(population) | population == 0,
                         NA_real_, cases / population * 1e5),
    inc_smooth  = roll_mean_centered(inc_per100k, SMOOTH_WINDOW)
  ) %>%
  ungroup()

# ---- 5. Wave detection per state --------------------------------------------

detect_state_waves <- function(df) {
  # df: one state's weekly series, already ordered by week_start
  above <- df$inc_smooth >= WAVE_INC_THRESHOLD
  above[is.na(above)] <- FALSE
  if (!any(above)) return(tibble())

  r      <- rle(above)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  runs <- tibble(value = r$values, start_idx = starts, end_idx = ends) %>%
    filter(value) %>%
    arrange(start_idx)

  # merge runs separated by a gap <= MIN_GAP_WEEKS
  merged <- vector("list", nrow(runs))
  cur <- runs[1, ]
  k <- 1
  if (nrow(runs) > 1) {
    for (i in 2:nrow(runs)) {
      gap <- runs$start_idx[i] - cur$end_idx - 1
      if (gap <= MIN_GAP_WEEKS) {
        cur$end_idx <- runs$end_idx[i]
      } else {
        merged[[k]] <- cur; k <- k + 1
        cur <- runs[i, ]
      }
    }
  }
  merged[[k]] <- cur
  merged <- bind_rows(merged[seq_len(k)])

  merged %>%
    rowwise() %>%
    mutate(
      start_week     = df$week_start[start_idx],
      end_week       = df$week_start[end_idx],
      duration_weeks = end_idx - start_idx + 1L,
      total_cases    = sum(df$cases[start_idx:end_idx]),
      peak_idx       = start_idx - 1L + which.max(df$cases[start_idx:end_idx]),
      peak_week      = df$week_start[peak_idx],
      peak_cases     = df$cases[peak_idx]
    ) %>%
    ungroup() %>%
    filter(duration_weeks >= MIN_WAVE_DURATION,
           total_cases    >= MIN_WAVE_CASES) %>%
    arrange(start_idx) %>%
    mutate(wave_id = row_number()) %>%
    select(wave_id, start_week, end_week, peak_week,
           duration_weeks, total_cases, peak_cases)
}

waves_state <- state_week %>%
  group_by(region, uf, uf_name) %>%
  group_modify(~ detect_state_waves(.x)) %>%
  ungroup()

message(sprintf("[waves] %d candidate epidemic waves detected across %d states",
                nrow(waves_state), n_distinct(waves_state$uf)))

# ---- 6. Municipality-level synchrony metrics per wave ----------------------

n_munis_by_state <- panel_raw %>%
  distinct(uf, muni6) %>%
  count(uf, name = "n_total_munis")

compute_wave_synchrony <- function(uf_i, start_week, end_week) {
  win_start <- start_week - lubridate::weeks(WAVE_BUFFER_WEEKS)
  win_end   <- end_week   + lubridate::weeks(WAVE_BUFFER_WEEKS)

  sub <- panel_raw %>%
    filter(uf == uf_i, week_start >= win_start, week_start <= win_end)

  muni_tot <- sub %>%
    group_by(muni6) %>%
    summarise(total = sum(cases, na.rm = TRUE), .groups = "drop") %>%
    filter(total >= ACTIVE_MUNI_MIN_CASES)

  n_active <- nrow(muni_tot)
  if (n_active == 0) {
    return(tibble(n_active_munis = 0L, peak_time_sd_weeks = NA_real_,
                  cor_sync_index = NA_real_, top3_share = NA_real_))
  }

  active_sub <- sub %>% filter(muni6 %in% muni_tot$muni6)

  peak_tbl <- active_sub %>%
    group_by(muni6) %>%
    slice_max(cases, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(muni_tot, by = "muni6") %>%
    mutate(t_weeks = as.numeric(week_start - start_week) / 7)

  w <- peak_tbl$total
  peak_time_sd_weeks <- if (n_active >= 2) {
    wm <- weighted.mean(peak_tbl$t_weeks, w)
    sqrt(sum(w * (peak_tbl$t_weeks - wm)^2) / sum(w))
  } else NA_real_

  cor_sync_index <- NA_real_
  if (n_active >= 3) {
    wide <- active_sub %>%
      select(muni6, week_start, cases) %>%
      pivot_wider(names_from = muni6, values_from = cases) %>%
      arrange(week_start) %>%
      select(-week_start) %>%
      as.matrix()
    cm <- suppressWarnings(cor(wide, use = "pairwise.complete.obs"))
    upper <- cm[upper.tri(cm)]
    cor_sync_index <- mean(upper, na.rm = TRUE)
  }

  top3_share <- sum(sort(muni_tot$total, decreasing = TRUE)[1:min(3, n_active)]) /
                sum(muni_tot$total)

  tibble(n_active_munis = n_active,
         peak_time_sd_weeks = peak_time_sd_weeks,
         cor_sync_index = cor_sync_index,
         top3_share = top3_share)
}

wave_synchrony <- waves_state %>%
  mutate(uf_chr = as.character(uf)) %>%
  rowwise() %>%
  mutate(sync = list(compute_wave_synchrony(uf_chr, start_week, end_week))) %>%
  ungroup() %>%
  unnest(sync) %>%
  left_join(n_munis_by_state %>% mutate(uf = as.character(uf)),
            by = c("uf_chr" = "uf")) %>%
  mutate(pct_active_munis = n_active_munis / n_total_munis) %>%
  select(-uf_chr)

message("[sync] municipality-level synchrony metrics computed for each wave")

# ---- 7. Classification ------------------------------------------------------
#
# Heuristic, transparent thresholds (not a statistical test) -- meant to
# separate three qualitatively different spatial signatures of a "state
# wave": (i) it really is one shared epidemic across many munis moving
# together; (ii) it is many local outbreaks staggered in time that only
# look like one wave once summed to the state level; (iii) a single
# municipality (usually the state capital) is doing almost all of it.

state_epidemic_audit <- wave_synchrony %>%
  mutate(
    synchrony_class = case_when(
      n_active_munis <= 1 ~ "Single-municipality dominated",
      peak_time_sd_weeks <= 2 & cor_sync_index >= 0.5 ~
        "Synchronous (state-wide)",
      peak_time_sd_weeks >= 4 | cor_sync_index < 0.3 ~
        "Asynchronous (localized/staggered)",
      TRUE ~ "Mixed / intermediate"
    ),
    is_first_major_epidemic = wave_id == 1,
    # The panel's first observed week is 2014-12-28 (a partial epi-week
    # right before YEAR_START). If a wave's start is pinned to that exact
    # week, cases were already substantial on day 1 of the data -- the
    # true start almost certainly predates the data window (left-censored),
    # not an artefact of the smoothing. Flag it rather than reporting a
    # false-precision start date.
    left_censored = start_week == min(state_week$week_start)
  ) %>%
  arrange(region, uf, wave_id)

# States with zero detected waves still need a row so the audit covers all 27.
states_no_wave <- uf_lookup %>%
  filter(!uf %in% state_epidemic_audit$uf) %>%
  mutate(wave_id = NA_integer_, synchrony_class = "No major epidemic detected",
         is_first_major_epidemic = FALSE, left_censored = FALSE)

state_epidemic_audit <- bind_rows(state_epidemic_audit, states_no_wave) %>%
  mutate(uf = factor(uf, levels = levels(uf_lookup$uf))) %>%
  arrange(region, uf, wave_id)

readr::write_csv(state_epidemic_audit,
                 file.path(TAB_DIR, "chik_state_epidemic_wave_audit.csv"))
message("[save] ", file.path(TAB_DIR, "chik_state_epidemic_wave_audit.csv"))

first_epidemic_summary <- state_epidemic_audit %>%
  filter(is_first_major_epidemic | is.na(wave_id)) %>%
  select(region, uf, uf_name, start_week, peak_week, end_week,
         duration_weeks, total_cases, n_active_munis, n_total_munis,
         pct_active_munis, peak_time_sd_weeks, cor_sync_index,
         synchrony_class, left_censored) %>%
  arrange(start_week)

readr::write_csv(first_epidemic_summary,
                 file.path(TAB_DIR, "chik_state_first_epidemic_summary.csv"))
message("[save] ", file.path(TAB_DIR, "chik_state_first_epidemic_summary.csv"))

n_censored <- sum(first_epidemic_summary$left_censored, na.rm = TRUE)
if (n_censored > 0) {
  message(sprintf(
    "\n[caveat] %d state(s) already had substantial cases in the panel's very
first observed week (%s) -- their true first-epidemic start almost
certainly predates the data window (consistent with chikungunya's known
2014 introduction into Brazil via Oiapoque/AP and Feira de Santana/BA).
These are flagged left_censored = TRUE; treat their 'start_week' as an
upper bound, not the true onset: %s",
    n_censored, min(state_week$week_start),
    paste(first_epidemic_summary$uf[first_epidemic_summary$left_censored],
          collapse = ", ")))
}

message("\n[summary] first major epidemic per state:")
print(first_epidemic_summary %>%
        select(uf, start_week, left_censored, synchrony_class), n = 30)

message("\n[summary] Bahia (BA) -- all detected waves, oldest first:")
print(state_epidemic_audit %>% filter(uf == "BA") %>%
        select(wave_id, start_week, peak_week, end_week, total_cases,
               n_active_munis, pct_active_munis, peak_time_sd_weeks,
               cor_sync_index, synchrony_class))

# ---- 8. Plot A: per-state diagnostic facets (QA on wave detection) --------

state_week_labelled <- state_week %>%
  mutate(uf_name = factor(uf_name, levels = uf_lookup$uf_name))

wave_rects <- state_epidemic_audit %>%
  filter(!is.na(wave_id)) %>%
  mutate(uf_name = factor(uf_name, levels = uf_lookup$uf_name))

p_diag <- ggplot() +
  geom_rect(data = wave_rects,
           aes(xmin = start_week, xmax = end_week,
               ymin = -Inf, ymax = Inf, fill = synchrony_class),
           alpha = 0.35) +
  geom_line(data = state_week_labelled,
           aes(x = week_start, y = inc_per100k),
           colour = "grey20", linewidth = 0.3) +
  geom_hline(yintercept = WAVE_INC_THRESHOLD, linetype = "dashed",
            colour = "grey40", linewidth = 0.3) +
  facet_wrap(~ uf_name, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = class_palette, name = "Wave classification") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title    = "State-level weekly incidence with detected epidemic waves",
    subtitle = sprintf(
      "Dashed line = wave threshold (%.0f/100k/wk, %d-wk smoothed). Shaded = detected wave, coloured by synchrony class. %d-%d",
      WAVE_INC_THRESHOLD, SMOOTH_WINDOW, YEAR_START, YEAR_END),
    x = NULL, y = "Incidence per 100,000",
    caption = "Source: SINAN / IBGE."
  ) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 7),
        plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "chik_state_wave_diagnostic_facets.png"),
       p_diag, width = 15, height = 11, dpi = 200, bg = "white")
message("[plot] ", file.path(FIG_DIR, "chik_state_wave_diagnostic_facets.png"))

# ---- 9. Plot B: Gantt-style overview of all waves, all states --------------

gantt_data <- state_epidemic_audit %>%
  filter(!is.na(wave_id)) %>%
  mutate(uf = factor(uf, levels = rev(levels(uf_lookup$uf))))

p_timeline <- ggplot(gantt_data) +
  geom_segment(aes(x = start_week, xend = end_week, y = uf, yend = uf,
                   colour = synchrony_class,
                   linewidth = log10(total_cases + 1)),
              lineend = "round") +
  geom_point(data = gantt_data %>% filter(is_first_major_epidemic),
            aes(x = peak_week, y = uf), shape = 8, size = 2, colour = "black") +
  scale_colour_manual(values = class_palette, name = "Synchrony class") +
  scale_linewidth_continuous(name = "Total cases\n(log scale)", range = c(0.6, 4),
                            breaks = log10(c(100, 1000, 10000, 100000) + 1),
                            labels = c("100","1k","10k","100k")) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title    = "Detected chikungunya epidemic waves per Brazilian state",
    subtitle = "Star = peak week of each state's first major epidemic. States ordered N -> NE -> SE -> S -> CO.",
    x = NULL, y = NULL,
    caption  = "Source: SINAN / IBGE."
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(family = "mono", size = 9),
        legend.position = "right",
        panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "chik_state_epidemic_timeline.png"),
       p_timeline, width = 12, height = 9, dpi = 200, bg = "white")
message("[plot] ", file.path(FIG_DIR, "chik_state_epidemic_timeline.png"))

# ---- 10. Plot C: Bahia deep dive -- municipality x week heatmap -----------

ba_panel <- panel_raw %>% filter(uf == "BA")
ba_first_active <- ba_panel %>%
  filter(cases > 0) %>%
  group_by(muni6) %>%
  summarise(first_active_week = min(week_start), .groups = "drop")

muni_order <- ba_first_active$muni6[order(ba_first_active$first_active_week)]
order_used <- "first active week (chronological onset)"

if (HAVE_GEOBR) {
  ba_centroids <- tryCatch({
    ba_sf <- geobr::read_municipality(code_muni = "BA", year = 2020,
                                      simplified = TRUE, showProgress = FALSE)
    ba_sf %>%
      mutate(muni6 = sprintf("%06d", code_muni %/% 10L),
            lat = sf::st_coordinates(sf::st_centroid(sf::st_geometry(ba_sf)))[, 2]) %>%
      sf::st_drop_geometry() %>%
      select(muni6, lat)
  }, error = function(e) NULL)
  if (!is.null(ba_centroids) && nrow(ba_centroids) > 0) {
    ord <- ba_centroids %>%
      filter(muni6 %in% ba_panel$muni6) %>%
      arrange(desc(lat)) %>%   # north (higher lat) at top
      pull(muni6)
    if (length(ord) == n_distinct(ba_panel$muni6)) {
      muni_order <- ord
      order_used <- "latitude, north -> south"
    }
  }
}
message(sprintf("[bahia] municipality rows ordered by: %s%s", order_used,
                if (order_used != "latitude, north -> south")
                  " (geobr polygons unavailable/offline -- fell back)" else ""))

ba_panel <- ba_panel %>%
  mutate(muni6 = factor(muni6, levels = muni_order))

ba_wave_lines <- state_epidemic_audit %>%
  filter(uf == "BA", !is.na(wave_id))

p_ba_heat <- ggplot(ba_panel, aes(x = week_start, y = muni6, fill = cases)) +
  geom_tile() +
  geom_vline(data = ba_wave_lines, aes(xintercept = start_week),
            linetype = "dashed", colour = "white", linewidth = 0.25) +
  geom_vline(data = ba_wave_lines, aes(xintercept = end_week),
            linetype = "dashed", colour = "white", linewidth = 0.25) +
  scale_fill_viridis_c(name = CASE_VAR, option = "magma", direction = -1,
                      trans = "log1p",
                      breaks = c(0, 1, 10, 100, 1000),
                      labels = c("0","1","10","100","1k"),
                      na.value = "grey95") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
  labs(
    title    = "Bahia: municipality x week chikungunya cases, 2015-2024",
    subtitle = sprintf(
      "Rows ordered by %s; dashed lines = detected state-wide wave boundaries (%d waves). Contrast row-spread within each wave to judge synchrony.",
      order_used, nrow(ba_wave_lines)),
    x = NULL, y = NULL,
    caption  = "Source: SINAN."
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())

ggsave(file.path(FIG_DIR, "chik_bahia_muni_week_heatmap.png"),
       p_ba_heat, width = 11, height = 9, dpi = 200, bg = "white")
message("[plot] ", file.path(FIG_DIR, "chik_bahia_muni_week_heatmap.png"))

# ---- 11. Plot D: system-wide trend in synchrony, 2015-2024 -----------------

trend_data <- state_epidemic_audit %>%
  filter(!is.na(wave_id), n_active_munis >= 2)

p_trend <- ggplot(trend_data, aes(x = start_week, y = cor_sync_index)) +
  geom_point(aes(colour = region, size = total_cases), alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE, colour = "black",
             linewidth = 0.6, span = 1) +
  scale_colour_manual(values = region_palette, name = "Region") +
  scale_size_continuous(name = "Total cases", range = c(1.5, 8),
                        labels = scales::comma) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title    = "Has within-state municipality synchrony changed over time?",
    subtitle = "Each point = one detected state wave. Higher = municipalities moved together; lower = staggered local outbreaks.",
    x = NULL, y = "Mean pairwise correlation of active munis' weekly cases",
    caption  = "Source: SINAN. Waves with < 2 active municipalities excluded."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

ggsave(file.path(FIG_DIR, "chik_synchrony_trend_over_time.png"),
       p_trend, width = 10, height = 6, dpi = 200, bg = "white")
message("[plot] ", file.path(FIG_DIR, "chik_synchrony_trend_over_time.png"))

# ---- 12. Plot E: per-state synchrony trajectory across successive waves --

p_state_traj <- state_epidemic_audit %>%
  filter(!is.na(wave_id)) %>%
  mutate(uf_name = factor(uf_name, levels = uf_lookup$uf_name)) %>%
  ggplot(aes(x = wave_id, y = cor_sync_index)) +
  geom_line(colour = "grey50", linewidth = 0.4) +
  geom_point(aes(colour = synchrony_class, size = total_cases)) +
  facet_wrap(~ uf_name, ncol = 5) +
  scale_colour_manual(values = class_palette, name = "Synchrony class") +
  scale_size_continuous(name = "Total cases", range = c(1, 5),
                        labels = scales::comma) +
  scale_x_continuous(breaks = scales::breaks_pretty(3)) +
  labs(
    title    = "Synchrony of successive waves within each state",
    subtitle = "Wave 1 = state's first major epidemic. Falling line = later waves more staggered/local than the first.",
    x = "Wave order (chronological)", y = "Mean pairwise correlation (active munis)",
    caption  = "Source: SINAN."
  ) +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 7),
        legend.position = "bottom")

ggsave(file.path(FIG_DIR, "chik_state_synchrony_by_wave_facet.png"),
       p_state_traj, width = 14, height = 11, dpi = 200, bg = "white")
message("[plot] ", file.path(FIG_DIR, "chik_state_synchrony_by_wave_facet.png"))

# ---- 13. Plot F: states ranked by timing of first major epidemic ---------

ranked <- first_epidemic_summary %>%
  filter(!is.na(start_week)) %>%
  arrange(start_week) %>%
  mutate(uf_name = factor(uf_name, levels = rev(uf_name)))

p_ranked <- ggplot(ranked, aes(x = start_week, y = uf_name, colour = region)) +
  geom_segment(aes(x = as.Date(sprintf("%d-01-01", YEAR_START)),
                  xend = start_week, yend = uf_name), linewidth = 0.4,
              colour = "grey80") +
  geom_point(aes(size = total_cases, shape = left_censored)) +
  scale_colour_manual(values = region_palette, name = "Region") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     name = "Start date",
                     labels = c(`FALSE` = "Observed onset",
                               `TRUE` = "Left-censored\n(true start predates data)")) +
  scale_size_continuous(name = "Total cases\n(first wave)", range = c(2, 8),
                        labels = scales::comma) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title    = "When did each state's first major chikungunya epidemic start?",
    subtitle = "Triangle = wave already underway at the very first observed week (2014-12-28); true onset predates the SINAN extract.",
    x = NULL, y = NULL,
    caption  = "Source: SINAN. States with no detected wave are omitted."
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank())

ggsave(file.path(FIG_DIR, "chik_state_first_epidemic_ranked.png"),
       p_ranked, width = 9, height = 8, dpi = 200, bg = "white")
message("[plot] ", file.path(FIG_DIR, "chik_state_first_epidemic_ranked.png"))

message("\n[done] all tables written to: ", TAB_DIR)
message("[done] all figures written to: ", FIG_DIR)
