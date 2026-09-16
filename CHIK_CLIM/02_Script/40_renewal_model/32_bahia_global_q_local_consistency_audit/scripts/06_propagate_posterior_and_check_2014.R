# Bahia global-q local-consistency audit -- Sections 11 & 12.

required_packages <- c("rstan", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir <- file.path(root, "03_Output/tables/bahia_global_q_local_consistency_audit")

# ---- Section 11: pre-2015 exposure check -----------------------------------
message("=== Section 11: pre-2015 Bahia chikungunya data check ===")
raw_sinan_files <- list.files(file.path(root, "01_Data/sinan_chik_csv"))
message("Raw SINAN files present: ", paste(raw_sinan_files, collapse = ", "))
message("Earliest raw SINAN file: CHIKBR15 (2015). No 2014 Bahia chikungunya surveillance file exists in this repository.")
message("CONCLUSION (Section 11): 2014 data are UNAVAILABLE. This diagnostic does NOT alter the Stan analysis (2015 start unchanged). ",
        "We state clearly, WITHOUT estimating a magnitude: observed 2016-2018 seroprevalence may partly reflect pre-2015 chikungunya exposure ",
        "(chikungunya was first confirmed circulating in Brazil in 2014), so q_implied values computed here from 2015+ reported cases only ",
        "could be biased DOWNWARD (i.e. true q could be somewhat higher than the crude q_implied suggests) to an unknown degree.")

# ---- Section 12: propagate the full q posterior -----------------------------
inputs <- readRDS(file.path(table_dir, "audit_inputs_bundle.rds"))
q_draws <- inputs$q_draws
R_survey_tbl <- read_csv(file.path(table_dir, "bahia_serology_implied_q_audit.csv"), show_col_types = FALSE)

posterior_case_implied <- lapply(seq_len(nrow(R_survey_tbl)), function(i) {
  draws <- R_survey_tbl$R_survey[i] / q_draws
  tibble(sero_id = R_survey_tbl$sero_id[i], name_muni = R_survey_tbl$name_muni[i],
         observed_prevalence = R_survey_tbl$observed_prevalence[i], R_survey = R_survey_tbl$R_survey[i],
         median = median(draws), lo50 = quantile(draws, .25), hi50 = quantile(draws, .75),
         lo95 = quantile(draws, .025), hi95 = quantile(draws, .975))
}) |> bind_rows()

message("\n=== Section 12: posterior-propagated case-implied local infection fraction (using full q posterior) ===")
print(as.data.frame(posterior_case_implied), digits = 4)
write_csv(posterior_case_implied, file.path(table_dir, "bahia_global_q_posterior_local_attack_summary.csv"))
message("[saved] ", file.path(table_dir, "bahia_global_q_posterior_local_attack_summary.csv"), " (survey-level)")

# Municipality-level A_i(2025) posterior distribution (all 414 munis).
implied <- read_csv(file.path(table_dir, "bahia_municipality_implied_attack_by_q.csv"), show_col_types = FALSE,
                     col_types = cols(muni6 = col_character()))
r2025 <- implied |> filter(checkpoint == as.Date("2025-12-21"), q == 0.136) |>
  distinct(muni6, name_muni, cumulative_reported_incidence)

muni_posterior_A <- lapply(seq_len(nrow(r2025)), function(i) {
  draws <- r2025$cumulative_reported_incidence[i] / q_draws
  tibble(muni6 = r2025$muni6[i], name_muni = r2025$name_muni[i],
         R_2025 = r2025$cumulative_reported_incidence[i],
         median_A = median(draws), lo95_A = quantile(draws, .025), hi95_A = quantile(draws, .975))
}) |> bind_rows()

write_csv(muni_posterior_A, file.path(table_dir, "bahia_global_q_posterior_municipality_attack_2025.csv"))
message("[saved] ", file.path(table_dir, "bahia_global_q_posterior_municipality_attack_2025.csv"), " (", nrow(muni_posterior_A), " municipalities, full q-posterior propagated)")
message(sprintf("\nAcross all 414 municipalities (q posterior propagated), median A_i(2025) ranges %.5f to %.3f; %d municipalities have posterior median A_i > 0.5; %d have 95%% CrI upper bound > 1.0.",
                 min(muni_posterior_A$median_A), max(muni_posterior_A$median_A),
                 sum(muni_posterior_A$median_A > 0.5), sum(muni_posterior_A$hi95_A > 1.0)))
