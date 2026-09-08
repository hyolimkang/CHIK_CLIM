# Renewal-model development order

The R scripts are organised in the order in which the renewal models were developed. Stan source files and saved posterior objects remain in `02_Script/stan/` because fitters refer to that central, versioned location.

1. `00_shared/` — state-weekly panel builder shared by the state-level models.
2. `01_v1/` — baseline state renewal model.
3. `02_v2/` — susceptible-depletion state model.
4. `03_v2_1/` — weekly-population, serology, reporting, comparison, and full-period stress-test workflow.
5. `04_v2_2/` — subsequent state-level demographic version.
6. `05_v3_0_spatial/` — municipality-level spatial prototype.
7. `06_v3_1_minimal_spatial/` — minimal spatial pilot.
8. `07_v3_2_coarse_regional/` — eight-unit coarse-regional pilot.

Version-specific outputs are under `03_Output/figures/renewal_<version>/` and `03_Output/tables/renewal_<version>/`. The v2-versus-v2.1 comparison and v2.1 full-period stress test are subfolders of `renewal_v2_1/`.

`03_Output/tables/renewal_model_development/renewal_convergence_comparison.csv` is an execution-based convergence ledger. It compares sampler performance, not scientific adequacy or predictive validity.
