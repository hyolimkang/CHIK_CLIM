# Serology audit -- PREPARATION ONLY (task Section 15). No serology
# likelihood is used anywhere in v4.1 q-only. This documents every
# Ceará-relevant candidate serology observation found in the project, for a
# POSSIBLE future model version, and explicitly flags provenance gaps and
# double-counting risk against the existing long-term FOI product.
#
# Source trace (full detail in this script's accompanying commit message /
# conversation record; summarised here):
#
# 1. `01_Data/serology/CountryModel.xlsx` (sheet `inclusion`, 602 rows) is
#    this project's own global chikungunya serosurvey systematic-review
#    line list (age-stratified N/N.pos per study). Two Brazil/Ceará
#    entries: study_no=165 "Brazil (Juazeiro do Norte)" (Barreto et al.,
#    2018, ELISA IgG+IgM, population-based, 4 age strata summing to
#    N=404, N.pos=103 -- EXACTLY matches the number hardcoded in every
#    renewal-model diagnose script since 05_v3_0_spatial); and
#    study_no=124 "Brazil (Ceará)" (Braga et al., 2018-2019, ELISA
#    IgG+IgM, general Ceará state -- NOT municipality-specific, N=41,
#    N.pos=4).
# 2. `01_Data/chik_foi.csv` (76 rows, per-study catalytic-model FOI fit
#    derived FROM CountryModel.xlsx) confirms study_no=165 (NumTested=404)
#    and mislabels study_no=124 as "Quixadá, Ceará" with NumTested=41 --
#    this is NOT the same record as the "Quixada" N=409/N.pos=289 used in
#    the renewal models.
# 3. The renewal-model "Quixada" N=409/N.pos=289 value (first hardcoded in
#    05_v3_0_spatial/prepare_renewal_v3_0_spatial_data.R) does NOT appear
#    anywhere in CountryModel.xlsx, chik_foi.csv, or markdownLAC.Rmd. The
#    project's own `serology_metadata.csv`
#    (03_Output/tables/episode_susceptibility/ceara_pilot_v1/) already
#    flags both records' source as "publication provenance not stored in
#    this repository" -- but only Juazeiro do Norte's numbers independently
#    reconcile against the systematic-review workbook; Quixada's do not.
# 4. The long-term gridded FOI surface (`01_Data/allfoi_s1.RData`, a global
#    100-member ensemble) is built upstream (outside this repo) FROM
#    CountryModel.xlsx's point/study data. Any future serology likelihood
#    that also uses Juazeiro do Norte's study 165 (or other
#    CountryModel.xlsx studies) would double-count information already
#    baked into allfoi_s1.RData's absolute FOI scale -- this must be
#    avoided (Section 15's explicit instruction).

required_packages <- c("here", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(tibble); library(readr) })

serology_audit_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

run_prepare_serology_audit <- function() {
  root <- serology_audit_root()
  out_dir <- file.path(root, "02_Script", "40_renewal_model", "20_v4_1_q_only")

  audit <- tribble(
    ~location, ~state, ~municipality, ~survey_start, ~survey_end, ~target_population, ~age_range,
    ~sampling_design, ~n_tested, ~n_positive, ~assay, ~repeated_cohort_or_independent, ~population_representative,
    ~geographic_scale, ~already_used_in_longterm_FOI_product, ~citation_source,

    "Juazeiro do Norte", "Ceara", "Juazeiro do Norte", "2018-06-03", "2018-12-30", "general population",
    "4 age strata (pooled all ages)", "population-based (per CountryModel.xlsx pop_group)", 404L, 103L,
    "ELISA IgG+IgM", "independent (cross-sectional)", "yes (population-based design, per source workbook)",
    "municipality-specific", "YES -- study_no=165 (Barreto et al. 2018) is part of CountryModel.xlsx, the point-data input to the global allfoi_s1.RData FOI surface. Using this survey again in a v4.x serology likelihood would double-count information already in the model's fixed long-term FOI scale.",
    "Traced exactly to 01_Data/serology/CountryModel.xlsx, sheet 'inclusion', study_no=165, author 'FKA, Barreto', year 2018 (sum of 4 age-stratified rows: N=69+114+44+177=404, N.pos=15+38+9+41=103). Also appears in 01_Data/chik_foi.csv (NumTested=404).",

    "Ceara (general, Braga et al.)", "Ceara", NA_character_, "2018", "2019", "general population",
    "4 age strata (pooled all ages)", "unclear (per CountryModel.xlsx pop_group)", 41L, 4L,
    "ELISA IgG+IgM", "independent (cross-sectional)", "unclear",
    "state-level (not municipality-specific)", "YES -- study_no=124 (Braga et al.) is part of CountryModel.xlsx and thus already in allfoi_s1.RData's input data.",
    "01_Data/serology/CountryModel.xlsx, sheet 'inclusion', study_no=124, author 'DAO, Braga'. Mislabelled 'Quixadá, Ceará' in 01_Data/chik_foi.csv (NumTested=41) -- NOT the same record as the renewal model's 'Quixada' N=409/N.pos=289.",

    "Quixada (as used in renewal models)", "Ceara", "Quixada", "2018-06-03", "2019-12-29", "unknown",
    "unknown", "unknown", 409L, 289L, "unknown",
    "unknown", "unknown",
    "municipality-specific (claimed)", "UNKNOWN -- cannot be checked against allfoi_s1.RData's input list because its source study cannot be identified.",
    "*** PROVENANCE NOT FOUND. *** Hardcoded since 05_v3_0_spatial/prepare_renewal_v3_0_spatial_data.R and propagated unchanged through v2.1/v2.2/v2.3/v3.1/v3.2/v4.0/v4.1 diagnose scripts. Does NOT appear in CountryModel.xlsx, chik_foi.csv, or markdownLAC.Rmd. The repo's own serology_metadata.csv already states 'publication provenance not stored in this repository.' DO NOT use in any future likelihood until a citable source is located -- treat as an external-consistency display item only, exactly as it is used today (never fitted)."
  )

  write_csv(audit, file.path(out_dir, "serology_audit_candidates.csv"))
  message("[serology audit] ", nrow(audit), " candidate record(s) documented (PREPARATION ONLY -- no likelihood use in v4.1 q-only).")
  message("[serology audit] *** Quixada (409/289) has NO traceable source anywhere in this repository -- flagged, not resolved. ***")
  message("[serology audit] Juazeiro do Norte (404/103) IS traceable and IS already part of the long-term FOI surface's input data -- any future serology likelihood must not double-count it.")
  invisible(audit)
}

if (sys.nframe() == 0L) run_prepare_serology_audit()
