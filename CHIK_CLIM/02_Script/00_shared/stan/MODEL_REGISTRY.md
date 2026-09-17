# Stan Model Registry

Every `.stan` file that existed under the old `02_Script/stan/` directory,
classified by whether any current (non-archived) pipeline script still
compiles/fits it. "Used by" lists the current consumer script(s) only --
see `git log` or the corresponding `90_development_archive` folder for
historical usage. "Fit output" gives the current location of that model's
posterior draws (`03_Output/...`), its compiled-model cache location if
different, or "n/a" if no fit is retained.

Status values: `current` (actively fit/read by a production pipeline),
`initialization_only` (fit once and reused as a fixed starting point, never
refit as part of normal operation), `sensitivity` (an alternative/what-if
variant explored alongside the accepted model), `legacy` (superseded,
frozen history), `failed` (an attempted extension that hit a computational
or identifiability blocker and was abandoned).

## current/initialization/

| Stan file | Status | Used by | Region(s) | Scientific role | Fit output |
|---|---|---|---|---|---|
| `renewal_ceara_v4_0_minimal_no_vaccine.stan` | initialization_only | `00_shared/legacy_model_functions/17_v4_0_minimal_no_vaccine/` (sourced by every state's Stage 02 for `generation_weights()`/`compute_hmc_gate()`) | All states | All-age, no-vaccination renewal benchmark; its `td14`-tagged fit is the fixed initial-value source for every state's v4.9 fit | `03_Output/02_ceara_pipeline/model_fits/initialization/` |

## current/ceara/

| Stan file | Status | Used by | Region(s) | Scientific role | Fit output |
|---|---|---|---|---|---|
| `renewal_ceara_v4_3_legacy_prior.stan` | sensitivity | `00_shared/legacy_model_functions/22_v4_3.../scripts/02_fit_v4_3.R` | Ceara | Legacy (pre-weak-prior) q calibration variant, kept for comparison | `03_Output/02_ceara_pipeline/model_fits/baseline/v4_3_short_2015_2019_q_calibration/` |
| `renewal_ceara_v4_3_weak_q.stan` | current | `.../scripts/02_fit_v4_3.R`, `03_fit_mode_targeted.R`, `04_mode_audit_analysis.R` | Ceara | Short-period (2015-2019) q calibration with a weak q prior -- part of the accepted baseline lineage | `03_Output/02_ceara_pipeline/model_fits/baseline/v4_3_short_2015_2019_q_calibration/` |
| `renewal_ceara_v4_3_weak_q_truncated.stan` | sensitivity | `.../scripts/05_fit_truncated_partition.R`, `08_partition_bridge_sampling.R` | Ceara | Truncated-partition sensitivity check on the v4.3 q calibration (low-incidence half) | `03_Output/02_ceara_pipeline/model_fits/baseline/v4_3_short_2015_2019_q_calibration/` |
| `renewal_ceara_v4_3_weak_q_truncated_high.stan` | sensitivity | `.../scripts/05_fit_truncated_partition.R`, `08_partition_bridge_sampling.R` | Ceara | Truncated-partition sensitivity check (high-incidence half) | `03_Output/02_ceara_pipeline/model_fits/baseline/v4_3_short_2015_2019_q_calibration/` |
| `renewal_ceara_v4_8_seeded_recurrence.stan` | current | `00_shared/legacy_model_functions/27_v4_8.../scripts/01_fit_v4_8.R` | Ceara | Adds a conditioned recurrence seed after the long post-2015-epidemic quiet gap -- next step in the accepted baseline lineage after v4.3 | `03_Output/02_ceara_pipeline/model_fits/baseline/v4_8_seeded_recurrence/` |
| `renewal_ceara_v4_9_climate_forced_canary.stan` | current | `02_ceara_pipeline/02_transmission_fitting/01_fit_ceara_climate_v49.R` | Ceara | Extends the accepted v4.9 baseline with temperature/precipitation anomaly covariates in log R0(t) -- the accepted climate-forced model | `03_Output/02_ceara_pipeline/model_fits/climate_forced/` |
| `renewal_ceara_v4_9_hierarchical_seasonality.stan` | current | `00_shared/legacy_model_functions/28_v4_9.../scripts/01_fit_v4_9.R`; hard prerequisite of `02_ceara_pipeline/01_input_preparation/00_run_stage_01.R` | Ceara | Adds hierarchical (year-varying) harmonic seasonality on top of v4.8 -- **the accepted Ceara baseline** that the climate-forced and vaccine-counterfactual work extend | `03_Output/02_ceara_pipeline/model_fits/baseline/v4_9_hierarchical_seasonality/` |
| `renewal_ceara_v4_9_serology_ablation.stan` | sensitivity | `00_shared/legacy_model_functions/28_v4_9.../scripts/04_fit_v4_9_serology_ablation.R` | Ceara | Ablates the serology likelihood term from v4.9 to test its influence on posterior q | `03_Output/02_ceara_pipeline/model_fits/baseline/v4_9_hierarchical_seasonality/` |

## current/multistate/

Generic models fit **verbatim** against a different state's data -- no
per-state duplicate copies of the same Stan code.

| Stan file | Status | Used by | Region(s) | Scientific role | Fit output |
|---|---|---|---|---|---|
| `renewal_bahia_v4_9_global_q_multisite_serology.stan` | current | Bahia Stage 02 (`04_fit_bahia_global_q_multisite.R`), Pernambuco Stage 02 (`02_fit_pe_global_q.R`), RJ Stage 01/02 (state+city), MT Stage 01/02 | Bahia, Pernambuco, Rio de Janeiro, Mato Grosso | State-agnostic global-q renewal model with an optional multisite-serology likelihood term (`J_sero>=0`) -- the shared v4.9-replication template for every non-Ceara state | `03_Output/{03_bahia,04_pernambuco,05_rio_de_janeiro,06_mato_grosso}_pipeline/model_fits/v4_9_replication/` |
| `renewal_ceara_v4_9_bahia_multisite_serology.stan` | current | Bahia Stage 02 (`03_fit_bahia_multisite_serology.R`, `04_plot_bahia_multisite_six_panel.R`), Pernambuco Stage 02 (`03_fit_pe_fixed_q_sweep.R`) | Bahia, Pernambuco | Ceara-origin v4.9 model extended with a multisite serology term, predating the fully state-agnostic `renewal_bahia_v4_9_global_q_multisite_serology.stan` template | `03_Output/{03_bahia,04_pernambuco}_pipeline/model_fits/v4_9_replication/` |
| `renewal_ceara_v4_9_hierarchical_seasonality_optional_serology.stan` | current | Bahia Stage 02 (`01_fit_bahia_v4_9.R`) | Bahia | The accepted Ceara v4.9 structure with an optional (togglable) serology term, used for Bahia's own v4.9 replication fit | `03_Output/03_bahia_pipeline/model_fits/v4_9_replication/` |

## current/national/

| Stan file | Status | Used by | Region(s) | Scientific role | Fit output |
|---|---|---|---|---|---|
| `renewal_ceara_episode_re.stan` | current | `07_national_pipeline/01_episode_renewal/`, `02_major_episode_renewal/`, `06_national_wave_analysis/01_fit_brazil_major_wave_early_re.R`, `08_first_epidemic_climate_transmission/` | National (Ceara episodes + Brazil-wide major-wave early Re) | Episode-level renewal model for early-phase effective-R estimation, fit both per-episode (4/6/8-week windows) and per-major-wave | `03_Output/07_national_pipeline/model_fits/{episode_renewal,major_episode_renewal,national_wave_analysis}/` |
| `episode_sequential_susceptibility.stan` | current | `07_national_pipeline/03_episode_susceptibility_reconstruction/02_fit_ceara_episode_susceptibility.R` | National (Ceara episodes) | Period-level (not weekly) susceptibility reconstruction across sequential Ceara epidemic episodes | `03_Output/07_national_pipeline/model_fits/episode_susceptibility/` |
| `chik_dynamic_annual_foi_shape_v2.stan` | current | `07_national_pipeline/06_national_wave_analysis/02_fit_brazil_annual_shape_susceptibility.R`, `11_annual_foi_shape_v2/` | National | Annual-resolution dynamic FOI model with a flexible within-year shape term -- current annual susceptibility-reconstruction model, fit under 6 sensitivity configurations (`M1v2_*`) | `03_Output/07_national_pipeline/model_fits/annual_foi_shape_v2/` |
| `chik_dynamic_annual_foi_nb_v2.stan` | current | `07_national_pipeline/11_annual_foi_shape_v2/02_fit_annual_foi_shape_v2.R` (NB-observation comparison arm) | National | NB2-observation-model variant of the annual dynamic-FOI model, fit alongside the SHAPE variants for model comparison | `03_Output/07_national_pipeline/model_fits/annual_foi_shape_v2/` |

## current/pernambuco/

| Stan file | Status | Used by | Region(s) | Scientific role | Fit output |
|---|---|---|---|---|---|
| `renewal_pernambuco_v4_9_spatial0.stan` | current | `04_pernambuco_pipeline/05_spatial_validation/18_fit_pe_spatial0.R`, `21_pe_q2_feasibility_diagnostic.R` | Pernambuco | Accepted, simpler single-stratum spatial extension of the PE v4.9 replication | `00_shared/stan/current/pernambuco/` (compiled cache); fit under `04_pernambuco_pipeline`'s own sensitivity-identification outputs |
| `renewal_pernambuco_v4_9_5strata_spatial.stan` | sensitivity (halted) | `04_pernambuco_pipeline/05_spatial_validation/12_fit_pe_5strata_spatial.R` | Pernambuco | 5-regional-stratum spatial disaggregation, a more detailed alternative to `spatial0` -- **halted before fitting**; see `PE_SPATIAL1_GEOMETRY_DIAGNOSTIC.md` and Stage 05's graceful skip of steps 13-15 | n/a (no accepted fit; pilot input audit only) |
| `renewal_pernambuco_v4_9_global_q_curved_reparam.stan` | sensitivity | `04_pernambuco_pipeline/04_sensitivity_identification/07_fit_pe_curved_ridge_reparam.R`, `09_fit_pe_curved_longwarmup_diagnostic.R` | Pernambuco | Curved-ridge reparameterisation explored for q identifiability | `04_pernambuco_pipeline`'s sensitivity-identification tables |
| `renewal_pernambuco_v4_9_global_q_ridge_reparam.stan` | sensitivity | `04_pernambuco_pipeline/04_sensitivity_identification/01_fit_pe_global_q_ridge_reparam.R`, `02_diagnose_pe_global_q_ridge_reparam.R` | Pernambuco | Ridge reparameterisation explored for q identifiability | `04_pernambuco_pipeline`'s sensitivity-identification tables |

## current/rio_de_janeiro/

| Stan file | Status | Used by | Region(s) | Scientific role | Fit output |
|---|---|---|---|---|---|
| `renewal_riodejaneiro_city_v4_9_global_q_logitnormal_serology.stan` | current | `05_rio_de_janeiro_pipeline/02_transmission_fitting/05_fit_rj_city_u10.R`, `06_fit_rj_city_u10_ad099.R` | Rio de Janeiro (city) | City-level global-q model with a logit-normal serology likelihood, for the RJ-city arm of the v4.9 replication | `03_Output/05_rio_de_janeiro_pipeline/model_fits/v4_9_replication/` |

## Legacy (moved to `02_Script/90_development_archive/08_old_stan_models/`)

All frozen/superseded Stan sources -- and their compiled-model caches --
below. Fitted results (where they exceed the small compiled-cache size)
moved separately to `03_Output/90_development_archive/model_fits/<version>/`.

| Stan file | Status | Historical role |
|---|---|---|
| `ceara_weekly.stan`, `ceara_weekly_v6.stan`, `ceara_weekly_v8.stan`, `ceara_weekly_v9.stan`, `ceara_weekly_v12.stan`, `ceara_weekly_v12d.stan`, `ceara_weekly_v12e.stan`, `ceara_weekly_v12f.stan`, `ceara_weekly_v12f_demography.stan`, `ceara_weekly_v13.stan` | legacy | Earliest weekly-transmission prototypes, predating the `renewal_*` version-numbered lineage (`90_development_archive/01_early_transmission_models/`) |
| `chik_dynamic_annual_foi.stan` | legacy | Original (pre-shape-v2) dynamic annual FOI model, including its stationary-AR1 (M1) susceptibility track -- superseded by `chik_dynamic_annual_foi_shape_v2.stan` (`90_development_archive/06_legacy_susceptibility/`) |
| `renewal_ceara_v1.stan` through `renewal_ceara_v3_2_coarse_regional.stan` (v1, v2, v2_1, v2_2, v2_3_state, v3_0_spatial, v3_1_minimal_spatial, v3_2_coarse_regional) | legacy | Sequential early renewal-model development, including the spatial-disaggregation exploratory tracks |
| `renewal_ceara_v4_1_minimal_no_vaccine.stan`, `renewal_ceara_v4_1_modular_q.stan`, `renewal_ceara_v4_1_q_only.stan` | legacy | v4.1 modular-q identifiability exploration |
| `renewal_ceara_v4_2_q_juazeiro_serology.stan`, `renewal_ceara_v4_2_q_juazeiro_serology_binomial.stan` | legacy | v4.2 Juazeiro-serology-constrained q calibration |
| `renewal_ceara_v4_3_4week_diagnostic.stan` | legacy | v4.3-era 4-week-observation diagnostic, superseded within the v4.3 lineage; fit only by the now-archived v4.4 fitting script |
| `renewal_ceara_v4_5_burden_shape.stan` | legacy | v4.5 burden-shape observation model |
| `renewal_ceara_v4_6_fixed_q_sweep.stan` | legacy | v4.6/v4.7 fixed-q sweep (full-period variant fit by the archived v4.7 script reuses this same source) |
| `renewal_ceara_v4_9_climate_canary_primary.stan` | legacy | Early draft of the climate-forced canary model, superseded by `renewal_ceara_v4_9_climate_forced_canary.stan` before ever being wired to a fitting script |
| `renewal_ceara_v5_0a_dynamic_R_direct_serology.stan`, `renewal_ceara_v5_0b_lowrank_dynamic_R.stan` | failed | Dynamic-R (weekly AR1 / low-rank) extension -- computationally infeasible; see `V5_0_DYNAMIC_R_ASSESSMENT.md` and `FAILED_v5_0a_weekly_AR1_computationally_infeasible/` |

## Notes on classification method

Every file above was classified by grepping every current pipeline folder
(`00_shared` through `07_national_pipeline`) for a literal reference to its
filename, then verifying the surviving ambiguous cases by reading the
producing/consuming script directly. `.rds` files were **not** classified
by extension or by proximity to a `.stan` file of the same name alone --
file size was used to separate genuine posterior fits (10MB-500MB,
depending on chain/iteration settings) from `rstan_options(auto_write =
TRUE)` compiled-model caches (a consistent ~1.7-1.9MB regardless of model
complexity, since that reflects compiled C++ template code, not posterior
draws). Two files (`ce_foi_summary_brazil_ceara.csv`-adjacent
`renewal_ceara_v2_2_fit.rds`, read optionally by a current diagnostic
script purely to quote a documentation snippet) kept their legacy
classification since their only *producer* is an archived script; the
current script's optional read path was updated to the new archive
location.
