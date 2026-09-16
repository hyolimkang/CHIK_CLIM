# Rio de Janeiro v4.9 Replication -- Model/Provenance Audit

RJ is fit as an INDEPENDENT external state, exactly as Ceara, Bahia, and
Pernambuco were. No Pernambuco-specific spatial structures, regional q,
regional transmission effects, ridge reparameterisations, or mobility are
imported. CE/BA/PE models, scripts, and results are unmodified.

## Version / provenance

- Git commit at start of RJ work: `133226a2e8642e41ac549bf9297ad053fd500311` (2026-09-12)
- Analysis directory: `02_Script/40_renewal_model/34_rio_de_janeiro_v4_9_replication/` (new, does not touch `30_bahia_v4_9_replication/` or `33_pernambuco_v4_9_replication/`)
- Random seed: `20260915L` (Stan `seed=`), `td14_scalar_inits` default seed `20260912L` (unchanged shared default)

## Frozen specification reused EXACTLY (no edits)

| Component | Source | Notes |
|---|---|---|
| Stan model | `02_Script/stan/renewal_bahia_v4_9_global_q_multisite_serology.stan` | Global-q, `J_sero>=0` (verified state-agnostic; only used with `J_sero=0` for Phase 1) |
| State-weekly case aggregation | `02_Script/40_renewal_model/00_shared/01_build_ceara_state_weekly.R` -- `build_state_weekly(panel, "33")` | UF-agnostic, unedited |
| Generation interval | `generation_weights()` in `17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R` -- G=8, Gamma(shape=4,rate=2) discretised | Unedited |
| HMC gate | `compute_hmc_gate()`, same file | divergences==0, max_treedepth_hits==0, Rhat<=1.01, min ESS>=100, BFMI>=0.3 |
| Data assembly helper | `make_v4_3_data()` in `22_v4_3_short_2015_2019_q_calibration/scripts/02_fit_v4_3.R` | Unedited |
| Dispersed inits | `load_td14_scalar_inits()`, same file | Unedited |
| Episode/wave definitions | `03_Output/tables/national_wave_census_v2/brazil_chik_wave_census_v2.csv`, `state=="RJ"` rows | National, state-agnostic, unmodified table -- SAME algorithm already used for CE/BA/PE |
| Importation convention | `imports_per_week = 1` | Frozen default, unchanged |
| Seasonality / year effects | Harmonic (`seasonal_sin/cos`, second harmonic `sin2/cos2`), `year_effect_prior_sd=0.40`, hierarchical `A_year` | Unchanged priors/structure |
| q prior | `logit_q ~ Normal(logit(0.10), 1.0)` | SAME broad prior as CE/BA/PE -- NOT tuned for RJ |

## RJ-specific preprocessing (new, mechanical only -- UF code/prefix substitution, no logic changes)

- `02_Script/00_data_prep/09d_fetch_rio_de_janeiro_sinasc_weekly_births.R` (`UF_PREFIX="33"`, verbatim copy of `09c` (Pernambuco), reuses the existing national SINASC cache -- no new downloads needed)
- `02_Script/00_data_prep/10d_build_rio_de_janeiro_weekly_demography.R` (`UF_CODE=33L`, verbatim copy of `10c`)
- `02_Script/40_renewal_model/34_rio_de_janeiro_v4_9_replication/scripts/01_rj_data_audit.R` -- `03_rj_v49_input_audit.R` (RJ-specific wrapper scripts, following the PE `01_pe_data_audit.R` template line-for-line with UF substitution)
- `scripts/04_fit_rj_global_q_case_only.R` (RJ-specific wrapper around the frozen shared helpers and Stan file, `J_sero=0` for Phase 1)

## Explicitly NOT done (per instruction)

- No fixed-q sweep before the global-q case-only fit
- No serology in the Phase-1 likelihood
- No Pernambuco affine/curved q-alpha_R ridge reparameterisation
- No spatial strata, regional q, regional R0, mobility, weekly AR1/dynamic R
- No modification to CE/BA/PE code, data, or results
