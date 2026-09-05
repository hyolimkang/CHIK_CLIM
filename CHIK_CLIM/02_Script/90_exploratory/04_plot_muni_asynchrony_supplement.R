# ---------------------------------------------------------------------------
# plot_muni_asynchrony_supplement.R
#
# Goal : ONE supplementary figure showing, for each state's largest detected
#        chikungunya wave, how spread out (asynchronous) each active
#        municipality's peak week is relative to the state-level peak.
#        This is the evidence figure behind "the fitted state-level curve is
#        a sum of staggered local epidemics, not one shared epidemic" --
#        used to justify not over-interpreting fitted state-level S_t /
#        beta_t as a single well-mixed biological process (disease-only
#        direct-effect framing as primary; spatial audit as diagnostic
#        support, not a new transmission model).
#
# Inputs:
#   01_Data/chik_brazil_muni_week_2015_2024.rds
#   03_Output/tables/chik_state_epidemic_wave_audit.csv
#     -> produced by 02_Script/90_exploratory/03_audit_state_epidemic_synchrony.R.
#        Run that script first if this table doesn't exist yet.
#        Re-uses its ACTIVE_MUNI_MIN_CASES / WAVE_BUFFER_WEEKS thresholds so
#        the SD/correlation annotated here match that table exactly.
#
# Output:
#   03_Output/figures/chik_supplement_muni_peak_asynchrony.png
# ---------------------------------------------------------------------------

# ---- 0. Packages ------------------------------------------------------------

for (p in c("here", "dplyr", "tidyr", "ggplot2", "lubridate", "readr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(ggplot2)
  library(lubridate); library(readr)
})

# ---- 1. Config ---------------------------------------------------------------

DATA_DIR <- here::here("01_Data")
FIG_DIR  <- here::here("03_Output", "figures")
TAB_DIR  <- here::here("03_Output", "tables")

PANEL_PATH <- file.path(DATA_DIR, "chik_brazil_muni_week_2015_2024.rds")
AUDIT_PATH <- file.path(TAB_DIR, "chik_state_epidemic_wave_audit.csv")

CASE_VAR <- "cases_notified"

# Which wave illustrates each state: "largest" (by total_cases -- usually
# the wave a state-level fit would target, most data-rich) or "first" (that
# state's first major epidemic, wave_id == 1). Default: largest.
WAVE_SELECTION <- "largest"

# Restrict to the states actually used in the paper's fitting, e.g.:
#   c("BA","CE","PE","RN","PB","MA","PI","SE","AL","MG","RJ")
# Leave NULL to include every state with a detected wave.
STATES_TO_PLOT <- NULL

# Same thresholds as audit_state_epidemic_synchrony.R -- kept identical on
# purpose so this figure's numbers agree with the main audit table.
ACTIVE_MUNI_MIN_CASES <- 20
WAVE_BUFFER_WEEKS     <- 4

stopifnot(file.exists(PANEL_PATH))
if (!file.exists(AUDIT_PATH))
  stop("Run 02_Script/90_exploratory/03_audit_state_epidemic_synchrony.R first -- ",
       AUDIT_PATH, " not found.")

# ---- 2. UF lookup (same convention as the rest of 02_Script) ---------------

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
                                    "Sul","Centro-Oeste")))

# ---- 3. Pick one illustrative wave per state -------------------------------

audit <- readr::read_csv(AUDIT_PATH, show_col_types = FALSE) %>%
  filter(!is.na(wave_id))

target_wave <- if (WAVE_SELECTION == "largest") {
  audit %>% group_by(uf) %>%
    slice_max(total_cases, n = 1, with_ties = FALSE) %>% ungroup()
} else {
  audit %>% filter(is_first_major_epidemic)
}

if (!is.null(STATES_TO_PLOT)) {
  target_wave <- target_wave %>% filter(uf %in% STATES_TO_PLOT)
}

message(sprintf("[select] illustrating the '%s' wave for %d states",
                WAVE_SELECTION, nrow(target_wave)))

# ---- 4. Load muni-week panel -------------------------------------------------

panel_raw <- readRDS(PANEL_PATH) %>%
  as_tibble() %>%
  transmute(muni6, week_start, cases = .data[[CASE_VAR]]) %>%
  mutate(uf_code = substr(muni6, 1, 2)) %>%
  left_join(uf_lookup, by = "uf_code") %>%
  filter(!is.na(uf))

# ---- 5. Per-municipality peak week within each selected wave's window -----

get_muni_peaks <- function(uf_i, start_week, end_week, peak_week) {
  win_start <- start_week - lubridate::weeks(WAVE_BUFFER_WEEKS)
  win_end   <- end_week   + lubridate::weeks(WAVE_BUFFER_WEEKS)

  sub <- panel_raw %>%
    filter(uf == uf_i, week_start >= win_start, week_start <= win_end)

  muni_tot <- sub %>%
    group_by(muni6) %>%
    summarise(total = sum(cases, na.rm = TRUE), .groups = "drop") %>%
    filter(total >= ACTIVE_MUNI_MIN_CASES)

  if (nrow(muni_tot) == 0) return(tibble())

  sub %>%
    filter(muni6 %in% muni_tot$muni6) %>%
    group_by(muni6) %>%
    slice_max(cases, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(muni_tot, by = "muni6") %>%
    mutate(offset_weeks = as.numeric(week_start - peak_week) / 7) %>%
    select(muni6, week_start, cases, total, offset_weeks)
}

muni_peaks <- target_wave %>%
  rowwise() %>%
  mutate(peaks = list(get_muni_peaks(uf, start_week, end_week, peak_week))) %>%
  ungroup() %>%
  select(uf, uf_name, region, wave_id, total_cases, n_active_munis,
         peak_time_sd_weeks, cor_sync_index, peaks) %>%
  unnest(peaks)

message(sprintf(
  "[compute] municipality-level peak weeks extracted for %d states, %d municipality-waves",
  n_distinct(muni_peaks$uf), nrow(muni_peaks)))

# ---- 6. Plot ------------------------------------------------------------------

panel_labels <- muni_peaks %>%
  distinct(uf, uf_name, n_active_munis, peak_time_sd_weeks, cor_sync_index) %>%
  mutate(
    label = sprintf("N=%d munis\nSD=%.1f wk\nr=%s",
                    n_active_munis, peak_time_sd_weeks,
                    ifelse(is.na(cor_sync_index), "NA",
                          sprintf("%.2f", cor_sync_index)))
  )

uf_order <- muni_peaks %>%
  distinct(uf, uf_name, region) %>%
  arrange(region, uf) %>%
  pull(uf_name)

muni_peaks   <- muni_peaks   %>% mutate(uf_name = factor(uf_name, levels = uf_order))
panel_labels <- panel_labels %>% mutate(uf_name = factor(uf_name, levels = uf_order))

weighted_mean_offset <- muni_peaks %>%
  group_by(uf_name) %>%
  summarise(wmean = weighted.mean(offset_weeks, total), .groups = "drop")

p <- ggplot(muni_peaks, aes(x = offset_weeks)) +
  geom_histogram(binwidth = 1, fill = "steelblue", colour = "white",
                linewidth = 0.1) +
  geom_vline(xintercept = 0, colour = "firebrick", linetype = "dashed",
            linewidth = 0.4) +
  geom_vline(data = weighted_mean_offset, aes(xintercept = wmean),
            colour = "grey30", linetype = "dotted", linewidth = 0.4) +
  geom_text(data = panel_labels, aes(label = label),
           x = Inf, y = Inf, hjust = 1.05, vjust = 1.1, size = 2.4,
           colour = "grey20", lineheight = 0.9) +
  facet_wrap(~ uf_name, scales = "free", ncol = 5) +
  labs(
    title    = "Within-state asynchrony: municipality peak timing relative to the state-level peak",
    subtitle = sprintf(
      paste("Each state's %s detected epidemic wave. Red dashed = state peak week (0).",
            "Grey dotted = case-weighted mean municipality peak."),
      WAVE_SELECTION),
    x = "Municipality peak week, relative to state peak (weeks)",
    y = "Number of municipalities",
    caption  = "Source: SINAN. Municipalities with < 20 cases within the wave window (±4 wk buffer) excluded."
  ) +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11),
        strip.text = element_text(face = "bold", size = 7.5),
        panel.grid.minor = element_blank())

n_states <- n_distinct(muni_peaks$uf)
ncol_fig <- min(5, n_states)
nrow_fig <- ceiling(n_states / ncol_fig)

out_path <- file.path(FIG_DIR, "chik_supplement_muni_peak_asynchrony.png")
ggsave(out_path, p, width = 3 * ncol_fig, height = 2.7 * nrow_fig,
       dpi = 250, bg = "white", limitsize = FALSE)
message("[plot] ", out_path)
