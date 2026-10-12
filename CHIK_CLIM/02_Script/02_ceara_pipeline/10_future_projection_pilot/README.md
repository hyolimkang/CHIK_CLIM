# Future projection pilot (Steps 1-2 only)

Isolated, additive analysis module. Does **not** modify any accepted
Stan model, posterior fit, or existing output. Re-analyses the already-
fitted climate-forced v4.9 canary posterior (Step 1) and builds a
standalone, pure-R stochastic future epidemic generator for 2026-2035
under stationary/resampled historical climate and **no vaccination**
(Step 2). No SSP/CMIP6 data are used anywhere in this module -- see
`STEP1_CLIMATE_INFERENCE_ASSESSMENT.md` and
`STEP2_FUTURE_GENERATOR_ASSESSMENT.md` in
`03_Output/02_ceara_pipeline/reports/future_projection_pilot/` for the
full findings and what is still required before Step 3 (actual
SSP/CMIP6 integration) or a vaccination extension.

## Execution order

Open `CHIK_CLIM.Rproj` first (sources `00_project_setup.R` via
`.Rprofile`, defining `ROOT`/`project_path()`). Then:

```r
# Step 1: climate-transmission inference audit (no Stan refit -- re-analyses
# the two already-saved posteriors on disk)
source("02_Script/02_ceara_pipeline/10_future_projection_pilot/00_run_step1.R")

# Step 2: mandatory validation, then the simulation ensemble, then figures.
# Set FUTURE_PILOT_ENSEMBLE_MODE to "smoke" (10 draws x 5 replicates,
# ~1 minute total) or "pilot" (100 draws x 20 replicates for the reference
# scenario, ~2 minutes total; default "smoke" if unset). The structural
# sensitivity grid (26 scenarios) always runs at smoke scale regardless
# of this setting -- see STEP2_FUTURE_GENERATOR_ASSESSMENT.md Section 6
# for why only the reference scenario is run at full pilot scale here.
Sys.setenv(FUTURE_PILOT_ENSEMBLE_MODE = "pilot")
source("02_Script/02_ceara_pipeline/10_future_projection_pilot/00_run_step2.R")
```

Each of those two runners is also independently re-runnable; the numbered
scripts underneath can be sourced individually if you only need one piece:

| Script | Role |
|---|---|
| `00_config_future_pilot.R` | Central configuration: paths, future horizon (2026-01-04 to 2035-12-30, 522 weeks), structural arms, importation grid, ensemble sizes, wave-detection settings. Sourced by everything else. |
| `01_step1_climate_posterior_reanalysis.R` | Step 1 numeric analysis: equation reconstruction/verification, beta_T/beta_P posterior vs. prior, extended correlation table, climate multiplier, full-series S(t) sensitivity, variance decomposition. Writes tables to `03_Output/.../tables/future_projection_pilot/`. |
| `02_step1_figures.R` | FIG01-04. |
| `03_future_renewal_functions.R` | Reusable, dependency-free simulation functions (no file I/O): future seasonal regressors, year-effect/climate-block/demography builders, the core renewal recursion (deterministic AND stochastic modes), the observation layer, and the adapted wave detector. Sourced by both the unit tests and the ensemble runner. |
| `04_run_unit_tests.R` | **Run this before trusting any ensemble output.** Section 2.10's 15 mandatory checks (A-J plus 2.1/2.2); stops with an error if anything fails. |
| `05_run_simulation_ensemble.R` | The simulation ensemble: one reference scenario (Y1 x C1 x D1, lambda_import=1) at the requested `FUTURE_PILOT_ENSEMBLE_MODE` scale, plus a 26-scenario structural sensitivity grid (year-effect x climate x importation, plus a demography-arm pair) always at smoke scale. Writes full trajectories (not just summaries) to `03_Output/.../results/future_projection_pilot/`. |
| `06_step2_figures.R` | FIG05-09. |

## Outputs

```
03_Output/02_ceara_pipeline/
  reports/future_projection_pilot/   STEP1_CLIMATE_INFERENCE_ASSESSMENT.md, STEP2_FUTURE_GENERATOR_ASSESSMENT.md
  tables/future_projection_pilot/    STEP1_*.csv, STEP2_unit_test_QA_summary.csv, STEP2_future_demography_*.csv
  figures/future_projection_pilot/   FIG01-09.{pdf,png@600dpi,svg}
  results/future_projection_pilot/   step1_beta_and_multiplier_draws.rds, STEP1_variance_decomposition_per_draw.csv,
                                      STEP2_reference_scenario_{trajectories.rds,annual_summary.csv,cumulative_summary.csv,waves.csv},
                                      STEP2_sensitivity_grid_cumulative_summary.csv
```

## Scope reminder

- Step 1 draws conclusions ONLY from the two posteriors already on disk
  (climate-forced q=0.05 canary; non-climate baseline q=0.05). No new
  Stan fit, no lag-spec sensitivity fit, no multi-state pooling.
- Step 2's climate arms (C0/C1/C2) use stationary or resampled
  **historical** weather only. No SSP/CMIP6 forcing. No vaccination.
- Every Step 2 figure carries the caption "Illustrative scenario
  simulation; no SSP forcing or vaccination applied." -- the trajectories
  are a scientifically auditable demonstration of a validated engine
  under labelled structural assumptions, not a forecast.
