# Mato Grosso v4.9 Input Audit

Frozen v4.9 architecture, reused EXACTLY (renewal transmission, lifelong
immunity, demographic S/U bookkeeping, common generation interval, NB2
observation model, harmonic seasonal R0(t), year effects, A_year,
existing importation convention). No MT-specific model changes.

## Key input dimensions

- N (weeks): 573
- Model years: 11 (2015-2025)
- First / last date: 2015-01-04 / 2025-12-21
- Total reported cases (2015-2025): 84,803
- Population range: 3,301,457 to 3,921,606
- Generation interval: G=8 (generation_weights(), Gamma(shape=4,rate=2) discretised, sums to 1) -- UNCHANGED
- Importation: imports_per_week = 1 -- UNCHANGED frozen default
- Demographic accounting max abs error: 0.00e+00 (must be ~0)
- Week continuity: TRUE
- NA values: 0
- Max inter-major-wave gap: 184 weeks (vs 399-week Ceara-2022 seeding standard) -- NO SEED USED

## Sources reused (unmodified)

- `02_Script/40_renewal_model/00_shared/01_build_ceara_state_weekly.R` (`build_state_weekly()`, UF-agnostic)
- `02_Script/stan/renewal_bahia_v4_9_global_q_multisite_serology.stan` (global-q, J_sero>=0, state-agnostic)
- `02_Script/40_renewal_model/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R` (`generation_weights()`, `compute_hmc_gate()`)
- `02_Script/40_renewal_model/22_v4_3_short_2015_2019_q_calibration/scripts/02_fit_v4_3.R` (`make_v4_3_data()`, `load_td14_scalar_inits()`)
- `03_Output/tables/national_wave_census_v2/brazil_chik_wave_census_v2.csv` (episode definitions, MT rows, unmodified)

## MT-specific preprocessing (unavoidable, mechanical only)

- `02_Script/00_data_prep/09f_fetch_mato_grosso_sinasc_weekly_births.R` (UF_PREFIX="51", verbatim copy of the RJ script)
- `02_Script/00_data_prep/10f_build_mato_grosso_weekly_demography.R` (UF_CODE=51L, verbatim copy)

CE/BA/PE/RJ models/results/scripts are UNMODIFIED by this work.
