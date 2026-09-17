# Bahia global-q local-consistency audit -- Sections 7, 8, 9, 10.
#
# Aligns municipality cumulative reported incidence to each survey's actual
# collection window (using the SAME survey-week weights already used in
# the Stan serology observation model -- uniform weights over the window),
# then computes crude case-implied local infection under candidate q's and
# the crude "implied q" diagnostic with Wilson binomial CIs.

required_packages <- c("dplyr", "readr", "tibble", "tidyr", "binom")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tibble); library(tidyr); library(binom) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
table_dir <- file.path(root, "03_Output/tables/bahia_global_q_local_consistency_audit")

inputs <- readRDS(file.path(table_dir, "audit_inputs_bundle.rds"))
bahia_panel <- inputs$bahia_panel
sero_audit <- inputs$sero_windows$audit
dates_all <- sort(unique(bahia_panel$week_start))

# Map each survey to its municipality by IBGE muni6 code (robust to
# accented-character encoding issues with regex name matching):
# 292630 = Riachao do Jacuipe, BA; 292740 = Salvador, BA.
survey_muni_map <- tribble(
  ~sero_id, ~muni6,
  "U03", "292630",
  "U04", "292740",
  "U19", "292740",
  "U20", "292740",
  "U21", "292740",
  "U22", "292740"
)

muni_lookup <- bahia_panel |> distinct(muni6, name_muni)
resolved_muni <- survey_muni_map |> left_join(muni_lookup, by = "muni6")

message("=== Section 7: survey -> municipality resolution ===")
print(as.data.frame(resolved_muni |> select(sero_id, name_muni, muni6)))
if (any(is.na(resolved_muni$muni6))) stop("STOP: could not resolve a survey's municipality by name match -- check name_muni spellings.")

# For each survey, build that municipality's weekly cumulative incidence
# trajectory, then average over the SAME window (uniform weights, matching
# the Stan serology model's window-average convention) used already.
compute_R_survey <- function(sero_id_val) {
  audit_row <- sero_audit |> dplyr::filter(sero_id == sero_id_val)
  muni6_val <- resolved_muni$muni6[resolved_muni$sero_id == sero_id_val]
  window_dates <- dates_all[audit_row$start_idx:(audit_row$start_idx + audit_row$n_weeks - 1)]

  muni_series <- bahia_panel |> dplyr::filter(muni6 == muni6_val) |> arrange(week_start) |>
    mutate(cumulative_cases = cumsum(cases_confirmed), cumulative_incidence = cumulative_cases / population)

  window_incidence <- muni_series |> dplyr::filter(week_start %in% window_dates) |> pull(cumulative_incidence)
  tibble(sero_id = sero_id_val, muni6 = muni6_val, name_muni = resolved_muni$name_muni[resolved_muni$sero_id == sero_id_val],
         R_survey = mean(window_incidence), window_start = min(window_dates), window_end = max(window_dates))
}

R_survey_tbl <- bind_rows(lapply(sero_audit$sero_id, compute_R_survey)) |>
  left_join(sero_audit |> select(sero_id, n_positive, n_tested, observed_prevalence), by = "sero_id")

message("\n=== Section 7: municipality cumulative reported incidence aligned to survey windows ===")
print(as.data.frame(R_survey_tbl |> select(sero_id, name_muni, window_start, window_end, R_survey, observed_prevalence)), digits = 4)
message("\nNOTE: Chapada, Pau da Lima, Alto do Cabrito, Marechal Rondon, Nova Constituinte, and Rio Sena are SUB-municipality survey sites. Municipality-level reported incidence is CONTEXT only, not a site-specific numerator (Section 7 caveat).")

# ---- Section 8: case-implied local infection under candidate q -----------
q_grid <- c(0.05, 0.10, 0.136, 0.15, 0.20, 0.25, 0.30, 0.359)
case_implied <- R_survey_tbl |>
  crossing(q = q_grid) |>
  mutate(p_case_implied = R_survey / q)

case_implied_wide <- case_implied |>
  select(sero_id, name_muni, observed_prevalence, R_survey, q, p_case_implied) |>
  pivot_wider(names_from = q, values_from = p_case_implied, names_prefix = "p_case_implied_q")

message("\n=== Section 8: case-implied local infection fraction under candidate q ===")
print(as.data.frame(case_implied_wide), digits = 3)
write_csv(case_implied_wide, file.path(table_dir, "bahia_serology_case_implied_attack_by_q.csv"))
message("[saved] ", file.path(table_dir, "bahia_serology_case_implied_attack_by_q.csv"))

# ---- Section 9: serology-implied q (crude), with Wilson CI on prevalence --
wilson_ci <- binom::binom.confint(R_survey_tbl$n_positive, R_survey_tbl$n_tested, method = "wilson")
implied_q_tbl <- R_survey_tbl |>
  mutate(
    sero_wilson_lo = wilson_ci$lower, sero_wilson_hi = wilson_ci$upper,
    q_implied = R_survey / observed_prevalence,
    q_implied_low_untruncated = R_survey / sero_wilson_hi,   # untruncated
    q_implied_high_untruncated = R_survey / sero_wilson_lo,  # untruncated
    q_implied_low_clamped = pmin(pmax(q_implied_low_untruncated, 0), 1),
    q_implied_high_clamped = pmin(pmax(q_implied_high_untruncated, 0), 1)
  )

message("\n=== Section 9: serology-implied q (crude), untruncated then clamped to [0,1] for plotting ===")
print(as.data.frame(implied_q_tbl |> select(sero_id, name_muni, R_survey, observed_prevalence, sero_wilson_lo, sero_wilson_hi,
                                              q_implied, q_implied_low_untruncated, q_implied_high_untruncated)), digits = 3)
write_csv(implied_q_tbl, file.path(table_dir, "bahia_serology_implied_q_audit.csv"))
message("[saved] ", file.path(table_dir, "bahia_serology_implied_q_audit.csv"))

# ---- Section 10: U19-U22 within-Salvador variation ------------------------
u_block <- implied_q_tbl |> dplyr::filter(sero_id %in% c("U19", "U20", "U21", "U22"))
message("\n=== Section 10: U19-U22 (same Salvador municipality-level R_survey, different community seroprevalence) ===")
print(as.data.frame(u_block |> select(sero_id, R_survey, observed_prevalence, q_implied)), digits = 3)
message(sprintf("\nAll four share R_survey = %.5f (identical Salvador municipality context). Observed community seroprevalence ranges %.1f%%-%.1f%% (%.1fx variation), so q_implied ranges %.3f-%.3f purely from WITHIN-municipality community heterogeneity in seroprevalence, not from any difference in q, reporting, or municipality-level case burden.",
                 unique(u_block$R_survey), 100*min(u_block$observed_prevalence), 100*max(u_block$observed_prevalence),
                 max(u_block$observed_prevalence)/min(u_block$observed_prevalence), min(u_block$q_implied), max(u_block$q_implied)))
