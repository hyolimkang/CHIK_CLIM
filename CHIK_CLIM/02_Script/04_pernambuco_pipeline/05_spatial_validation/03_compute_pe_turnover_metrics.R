# Pernambuco spatial-turnover diagnostic -- Sections 4 & 7.
#
# Section 4: for every pair of successive major waves -- Spearman rank
# correlation of municipality incidence, top-10%/20% Jaccard overlap,
# 50%/80% case-contributor overlap, and the share of wave w+1 cases coming
# from municipalities that were low/zero-incidence in wave w.
#
# Section 7: the same core metrics (Spearman, top-10% Jaccard, top-20%
# Jaccard, 80%-contributor overlap) for EVERY pairwise combination of major
# waves, not just adjacent ones, for the turnover-matrix heatmaps.
#
# State-agnostic logic, ported near-verbatim from
# 31_bahia_spatial_turnover_diagnostic/scripts/03_compute_turnover_metrics.R.

required_packages <- c("dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tibble); library(tidyr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/tables/pernambuco_v4_9_replication/spatial_turnover_diagnostic")

burden <- read_csv(file.path(table_dir, "pe_municipality_wave_burden.csv"), show_col_types = FALSE)
waves <- read_csv(file.path(table_dir, "pe_wave_definitions.csv"), show_col_types = FALSE) |> arrange(start_week)
wave_ids <- waves$wave_id

top_set <- function(df, prop) {
  n <- ceiling(nrow(df) * prop)
  df |> arrange(desc(incidence_iw)) |> slice_head(n = n) |> pull(muni6)
}
contributor_set <- function(df, prop) {
  df <- df |> arrange(desc(reported_cases_iw)) |> mutate(cum_share = cumsum(reported_cases_iw) / sum(reported_cases_iw))
  df |> dplyr::filter(lag(cum_share, default = 0) < prop) |> pull(muni6)
}
jaccard <- function(a, b) {
  if (length(a) == 0L && length(b) == 0L) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

pair_metrics <- function(w1, w2) {
  d1 <- burden |> dplyr::filter(wave_id == w1)
  d2 <- burden |> dplyr::filter(wave_id == w2)
  joined <- inner_join(d1 |> select(muni6, incidence_iw, reported_cases_iw, cumulative_previous_incidence_iw),
                        d2 |> select(muni6, incidence_iw, reported_cases_iw),
                        by = "muni6", suffix = c("_w1", "_w2"))

  spearman_rho <- suppressWarnings(cor(joined$incidence_iw_w1, joined$incidence_iw_w2, method = "spearman", use = "complete.obs"))

  top10_w1 <- top_set(d1, 0.10); top10_w2 <- top_set(d2, 0.10)
  top20_w1 <- top_set(d1, 0.20); top20_w2 <- top_set(d2, 0.20)
  contrib50_w1 <- contributor_set(d1, 0.50); contrib50_w2 <- contributor_set(d2, 0.50)
  contrib80_w1 <- contributor_set(d1, 0.80); contrib80_w2 <- contributor_set(d2, 0.80)

  # Section 4D: contribution to wave w2 from municipalities that were low/zero in wave w1
  d1_ranked <- d1 |> mutate(incidence_rank = rank(incidence_iw, ties.method = "average") / n())
  bottom50_w1 <- d1_ranked |> dplyr::filter(incidence_rank <= 0.50) |> pull(muni6)
  bottom25_w1 <- d1_ranked |> dplyr::filter(incidence_rank <= 0.25) |> pull(muni6)
  zero_w1 <- d1 |> dplyr::filter(reported_cases_iw == 0) |> pull(muni6)

  total_w2_cases <- sum(d2$reported_cases_iw)
  share_from <- function(muni_set) {
    if (total_w2_cases == 0) return(NA_real_)
    sum(d2$reported_cases_iw[d2$muni6 %in% muni_set]) / total_w2_cases
  }

  tibble(
    wave_from = w1, wave_to = w2,
    spearman_incidence_corr = spearman_rho,
    top10pct_jaccard = jaccard(top10_w1, top10_w2),
    top20pct_jaccard = jaccard(top20_w1, top20_w2),
    contrib50pct_overlap_jaccard = jaccard(contrib50_w1, contrib50_w2),
    contrib80pct_overlap_jaccard = jaccard(contrib80_w1, contrib80_w2),
    share_w2_cases_from_bottom50_w1 = share_from(bottom50_w1),
    share_w2_cases_from_bottom25_w1 = share_from(bottom25_w1),
    share_w2_cases_from_zero_w1 = share_from(zero_w1)
  )
}

# ---- Adjacent-wave pairs (Section 4) --------------------------------------
adjacent_pairs <- tibble(w1 = wave_ids[-length(wave_ids)], w2 = wave_ids[-1])
adjacent_metrics <- bind_rows(Map(pair_metrics, adjacent_pairs$w1, adjacent_pairs$w2))
adjacent_metrics$comparison_type <- "adjacent"

# ---- All pairwise combinations (Section 7) ---------------------------------
all_pairs <- as_tibble(t(combn(wave_ids, 2))) |> setNames(c("w1", "w2"))
all_pairwise_metrics <- bind_rows(Map(pair_metrics, all_pairs$w1, all_pairs$w2))
all_pairwise_metrics$comparison_type <- "all_pairwise"

turnover_metrics <- bind_rows(adjacent_metrics, all_pairwise_metrics)

message("=== Adjacent-wave turnover metrics (Section 4) ===")
print(as.data.frame(adjacent_metrics |> select(wave_from, wave_to, spearman_incidence_corr, top10pct_jaccard, top20pct_jaccard,
                                                 contrib80pct_overlap_jaccard, share_w2_cases_from_bottom50_w1,
                                                 share_w2_cases_from_zero_w1)))

write_csv(turnover_metrics, file.path(table_dir, "pe_wave_turnover_metrics.csv"))
message("\n[saved] ", file.path(table_dir, "pe_wave_turnover_metrics.csv"))

message("\n=== Summary across adjacent-wave transitions ===")
message(sprintf("Median Spearman incidence correlation: %.3f (range %.3f-%.3f)",
                 median(adjacent_metrics$spearman_incidence_corr), min(adjacent_metrics$spearman_incidence_corr), max(adjacent_metrics$spearman_incidence_corr)))
message(sprintf("Median top-10%% Jaccard: %.3f | top-20%% Jaccard: %.3f | 80%%-contributor Jaccard: %.3f",
                 median(adjacent_metrics$top10pct_jaccard), median(adjacent_metrics$top20pct_jaccard), median(adjacent_metrics$contrib80pct_overlap_jaccard)))
message(sprintf("Median share of next-wave cases from previously bottom-50%% municipalities: %.3f | from zero-case municipalities: %.3f",
                 median(adjacent_metrics$share_w2_cases_from_bottom50_w1), median(adjacent_metrics$share_w2_cases_from_zero_w1)))
