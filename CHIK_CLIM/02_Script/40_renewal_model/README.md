# Renewal-model development order

## CURRENT CÉARA PRODUCTION WORKFLOW

**For the current Ceará analysis (climate-forced v4.9 transmission model,
age reconstruction, and the historical age-12 vaccine counterfactual +
NNV), use `38_ceara_current_pipeline/`.** Start with its
`00_README_EXECUTION_ORDER.md`.

Use folders `01` -> `09` there for current analyses. Folders `01_v1`
through `37_ceara_historical_age12_vaccine_counterfactual` below primarily
document model development, replication across other states, or prior
analysis structure, and should **not** be manually executed sequentially
for the current Ceará analysis — `38_ceara_current_pipeline/` sources the
specific frozen functions it still depends on directly, so you should not
need to open those folders yourself.

The rest of this document (below) is retained as the historical
development-order record and is unchanged by that refactor.

---

The R scripts are organised in the order in which the renewal models were developed. Stan source files and saved posterior objects remain in `02_Script/stan/` because fitters refer to that central, versioned location.

1. `00_shared/` — state-weekly panel builder shared by the state-level models.
2. `01_v1/` — baseline state renewal model.
3. `02_v2/` — susceptible-depletion state model.
4. `03_v2_1/` — weekly-population, serology, reporting, comparison, and full-period stress-test workflow.
5. `04_v2_2/` — subsequent state-level demographic version.
6. `05_v3_0_spatial/` — municipality-level spatial prototype.
7. `06_v3_1_minimal_spatial/` — minimal spatial pilot.
8. `07_v3_2_coarse_regional/` — eight-unit coarse-regional pilot.
9. `09_v2_3_state/` — state-level stabilization model: unchanged explicit S/U demography, direct informative overall detection, stationary seasonal AR(1) R0 process, and tempered serology calibration. Recurrent introduction is intentionally absent from this initial version.

10. `10_episode_renewal/` — outbreak-level early-phase renewal analysis: one effective reproduction number per audit-defined Ceará outbreak, without a long-term weekly transmission process.

11. `11_major_episode_renewal/` — fixed six-week early-phase Re analysis for audit-defined Ceará candidate waves with at least 4,000 reported cases.

12. `12_episode_susceptibility_reconstruction/` - period-level sequential S/U reconstruction pilot. It estimates one cumulative force of infection per epidemic or inter-episode period, reports filtered onset susceptibility separately from smoothed sensitivity estimates, and includes no weekly transmission, R0, Reff, introductions, or climate component. The first Ceara pilot output is isolated under `episode_susceptibility/ceara_pilot_v1/`; its current HMC diagnostic gate status must be checked before use.

Version-specific outputs are under `03_Output/figures/renewal_<version>/` and `03_Output/tables/renewal_<version>/`. The v2-versus-v2.1 comparison and v2.1 full-period stress test are subfolders of `renewal_v2_1/`.

`03_Output/tables/renewal_model_development/renewal_convergence_comparison.csv` is an execution-based convergence ledger. It compares sampler performance, not scientific adequacy or predictive validity.
