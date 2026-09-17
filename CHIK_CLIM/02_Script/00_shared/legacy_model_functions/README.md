# Legacy model functions -- still-live shared dependencies

These 4 folders look like frozen version history (and are numbered the
same as when they lived under `40_renewal_model/`), but they are **not**
dead code: every state's Stage 02 transmission-fitting script (Ceará,
Bahia, Pernambuco, Rio de Janeiro, Mato Grosso) sources one or more of
them to get shared functions (`generation_weights()`, `compute_hmc_gate()`,
`make_*_data()`/`make_*_init_fn()` helpers, etc.) that the frozen v4.9
model structure is built on.

| Folder | Provides |
|---|---|
| `17_v4_0_minimal_no_vaccine/` | Base generation-interval/HMC-gate machinery every later version sources |
| `22_v4_3_short_2015_2019_q_calibration/` | q-calibration fitting machinery |
| `27_v4_8_seeded_recurrence/` | Seeded-recurrence fitting machinery |
| `28_v4_9_hierarchical_seasonality/` | The current frozen v4.9 structure every state's "v4.9 replication" is built on |

Their own accepted Stan fit outputs (`.rds`, chain CSVs) live under
`03_Output/02_ceara_pipeline/model_fits/{initialization,baseline}/<version>/`
per the project's extension-based placement rule (`.R` -> `02_Script`,
`.rds` -> `03_Output/`) -- only the scripts and their audit/design docs
stayed here.

Do not move or rename these without first re-checking every live pipeline
that sources them (`grep -rl "legacy_model_functions" 02_Script/`).
