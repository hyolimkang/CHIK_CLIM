# Step 6 pilot: centre of gravity by year (CE, BA, RR)
library(dplyr)
library(ggplot2)
library(lubridate)

panel <- readRDS("01_Data/chik_brazil_muni_month_2015_2024.rds")
pop   <- readRDS("01_Data/ibge_pop_muni_year_2015_2024.rds")

uf_lookup <- tibble::tribble(
  ~uf_code, ~uf,
  "11", "RO", "12", "AC", "13", "AM", "14", "RR", "15", "PA", "16", "AP", "17", "TO",
  "21", "MA", "22", "PI", "23", "CE", "24", "RN", "25", "PB", "26", "PE", "27", "AL",
  "28", "SE", "29", "BA", "31", "MG", "32", "ES", "33", "RJ", "35", "SP",
  "41", "PR", "42", "SC", "43", "RS", "50", "MS", "51", "MT", "52", "GO", "53", "DF"
)

state_pop <- pop |>
  group_by(uf, year) |>
  summarise(population = sum(population), .groups = "drop")

state <- panel |>
  mutate(uf_code = substr(muni6, 1, 2)) |>
  group_by(year_month, uf_code) |>
  summarise(cases = sum(cases_confirmed), .groups = "drop") |>
  left_join(uf_lookup, by = "uf_code") |>
  mutate(
    year  = year(year_month),
    month = month(year_month)
  ) |>
  left_join(state_pop, by = c("uf", "year")) |>
  mutate(incidence_per_100k = cases / population * 1e5)

state_year_cog <- state |>
  group_by(uf, year) |>
  summarise(
    cog_month = {
      w <- incidence_per_100k
      if (sum(w, na.rm = TRUE) > 0) sum(month * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
      else NA_real_
    },
    .groups = "drop"
  )

focus <- state_year_cog |>
  filter(uf %in% c("CE", "BA", "RR")) |>
  arrange(uf, year)

print(focus)

dir.create("03_Output/figures", recursive = TRUE, showWarnings = FALSE)

p <- ggplot(focus, aes(year, cog_month, colour = uf)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 2015:2024) +
  scale_y_continuous(limits = c(1, 12), breaks = 1:12) +
  labs(
    title = "Centre of gravity (month) by year — pilot",
    subtitle = "CE = two-wave; BA = persistent; RR = early sporadic",
    x = "Year", y = "Mean epidemic month (1=Jan, 12=Dec)",
    colour = "UF"
  ) +
  theme_minimal()

ggsave("03_Output/figures/step6_cog_by_year_CE_BA_RR.png", p, width = 8, height = 5, dpi = 150)
message("[saved] 03_Output/figures/step6_cog_by_year_CE_BA_RR.png")
