# Bahia global-q local-consistency audit -- Section 5.
#
# CRUDE case-implied cumulative infection fraction per municipality,
# assuming constant q and ignoring transmission dynamics entirely. This is
# NOT a Stan estimate -- explicitly labelled throughout.

required_packages <- c("dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tibble); library(tidyr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir <- file.path(root, "03_Output/tables/bahia_global_q_local_consistency_audit")

inputs <- readRDS(file.path(table_dir, "audit_inputs_bundle.rds"))
bahia_panel <- inputs$bahia_panel

checkpoints <- as.Date(c("2016-12-31", "2017-12-31", "2018-12-31", "2020-12-31", "2022-12-31", "2025-12-21"))
q_grid <- c(0.05, 0.10, 0.136, 0.15, 0.20, 0.25, 0.30, 0.359)

# Cumulative reported incidence R_i(t) per municipality at each checkpoint,
# using each municipality's population AT that checkpoint week.
cum_incidence <- bind_rows(lapply(checkpoints, function(cp) {
  up_to <- bahia_panel |> filter(week_start <= cp)
  pop_at_cp <- bahia_panel |> filter(week_start == max(week_start[week_start <= cp])) |>
    select(muni6, population_at_checkpoint = population)
  up_to |>
    group_by(muni6, name_muni) |>
    summarise(cumulative_cases = sum(cases_confirmed, na.rm = TRUE), .groups = "drop") |>
    left_join(pop_at_cp, by = "muni6") |>
    mutate(checkpoint = cp, cumulative_reported_incidence = cumulative_cases / population_at_checkpoint)
}))

# Crude case-implied cumulative infection fraction under each candidate q.
implied <- cum_incidence |>
  crossing(q = q_grid) |>
  mutate(A_i = cumulative_reported_incidence / q)

write_csv(implied, file.path(table_dir, "bahia_municipality_implied_attack_by_q.csv"))
message("[saved] ", file.path(table_dir, "bahia_municipality_implied_attack_by_q.csv"), " (", nrow(implied), " rows)")

# ---- Summary distribution + flagging, per checkpoint x q -------------------
summary_tbl <- implied |>
  group_by(checkpoint, q) |>
  summarise(
    n_munis = n(),
    median_A = median(A_i), pop_weighted_mean_A = weighted.mean(A_i, w = population_at_checkpoint),
    p90_A = quantile(A_i, .90), p95_A = quantile(A_i, .95), max_A = max(A_i),
    n_gt_025 = sum(A_i > 0.25), pct_gt_025 = 100 * mean(A_i > 0.25),
    n_gt_050 = sum(A_i > 0.50), pct_gt_050 = 100 * mean(A_i > 0.50),
    n_gt_075 = sum(A_i > 0.75), pct_gt_075 = 100 * mean(A_i > 0.75),
    n_gt_100 = sum(A_i > 1.00), pct_gt_100 = 100 * mean(A_i > 1.00),
    .groups = "drop"
  )

message("\n=== Section 5: municipality-level crude implied attack fraction summary (2025 checkpoint) ===")
print(as.data.frame(summary_tbl |> filter(checkpoint == max(checkpoints))), digits = 3)

write_csv(summary_tbl, file.path(table_dir, "bahia_municipality_implied_attack_summary.csv"))
message("\n[saved] ", file.path(table_dir, "bahia_municipality_implied_attack_summary.csv"))

message("\n=== Full checkpoint x q summary ===")
print(as.data.frame(summary_tbl |> select(checkpoint, q, median_A, pop_weighted_mean_A, max_A, pct_gt_050, pct_gt_100)), digits = 3)

# ---- Flag audit for A_i > 1: are these plausibly duplicate/reporting issues?
flagged <- implied |> filter(A_i > 1) |>
  left_join(bahia_panel |> group_by(muni6) |> summarise(muni_max_pop = max(population), .groups = "drop"), by = "muni6")
message(sprintf("\nTotal municipality x checkpoint x q rows with A_i > 1: %d (out of %d)", nrow(flagged), nrow(implied)))
if (nrow(flagged) > 0) {
  flagged_summary <- flagged |> filter(q == 0.136) |> arrange(desc(A_i)) |> select(muni6, name_muni, checkpoint, cumulative_cases, population_at_checkpoint, cumulative_reported_incidence, A_i)
  message("\nAt q=0.136, municipality x checkpoint rows with A_i > 1 (largest first):")
  print(as.data.frame(head(flagged_summary, 20)), digits = 3)
  write_csv(flagged, file.path(table_dir, "bahia_municipality_flagged_A_gt1.csv"))
  message("[saved] ", file.path(table_dir, "bahia_municipality_flagged_A_gt1.csv"))
} else {
  message("No municipality x checkpoint x q combination exceeds A_i = 1.")
}
