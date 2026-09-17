# DATA PIPELINE

Shared, state-spanning data acquisition and preparation. Runs BEFORE any
per-state or national analysis pipeline (`02_ceara_pipeline/`, ...,
`07_national_pipeline/`). Consolidated from the former `00_data_prep/` and
the fetch/build portion of `30_climate_covariates_dlnm/` -- no script logic
changed, only location (and a few outdated cross-reference comments/messages
updated to the new paths).

## Stages

| Step | Folder | Contents | Notes |
|---|---|---|---|
| 01 | `01_surveillance/` | Fetch + clean national SINAN chikungunya case data | `01_fetch_sinan.R` -> `02_clean_sinan.R` |
| 02 | `02_population/` | IBGE population (muni-year), population projections (UF), single-year age population (UF) | Three independent fetches, no inter-dependency |
| 03 | `03_demography/` | Per-state SINASC weekly births (Ceara, Bahia, Pernambuco, Rio de Janeiro state + city, Mato Grosso) -> per-state weekly demography build | Kept as separate per-state scripts (`01_09*_fetch_*_sinasc_weekly_births.R` -> `03_10*_build_*_weekly_demography.R`), NOT force-genericised into one script, to avoid any risk of altering per-state logic. There is no separate "fetch deaths" step in the current pipeline -- death counts are obtained as part of each state's weekly-demography build script, not fetched independently. |
| 04 | `04_geography/` | Municipality centroids (used for climate spatial join) | |
| 05 | `05_climate/` | ERA5-Land temperature/precipitation fetch -> national weekly climate table -> national municipality-week DLNM panel | `04_build_climate_panel.R` produces `chik_dlnm_panel_muni_week_2015_2025.rds`, the core panel read by every downstream state/national pipeline |

## Deferred / not yet resolved (see parent PR conversation)

Two files from `30_climate_covariates_dlnm/` were NOT moved here because their
canonical status is still undecided (left in place at
`02_Script/30_climate_covariates_dlnm/` pending a decision, not archived):
- `00_current_covariates_baseline.R` -- a near-duplicate of `01_Data/current_covariates_baseline.R`
  exists with materially different (more developed) imputation logic; which
  one is canonical has not been decided.
- `05_dlnm_model.R` -- an actual DLNM statistical model, possibly superseded
  by the renewal-model framework (v1-v4.9); current-vs-archive status not
  yet decided.

Superseded: `00_data_prep/05_fetch_ce_births_ibge.R` (an earlier, Ceara-only,
IBGE-based births source) was archived to
`90_development_archive/09_misc_archive/05_fetch_ce_births_ibge_LEGACY.R`,
superseded by the SINASC-based per-state approach in `03_demography/`.

## Running

```r
source("02_Script/01_data_pipeline/99_run_data_pipeline.R")
```

Each stage's scripts are independent fetch/build steps (no source()
cross-references between them); the orchestrator simply runs them in the
order above. Most involve external API calls (SINAN, IBGE, ERA5-Land) and
are NOT re-run automatically if their output already exists -- check each
script's own output-existence guard before assuming a call will hit the
network.
