# 90_development_archive

This directory records model development, superseded analyses, and legacy
workflows.

**Current production pipelines must NOT depend directly on files in this
directory.**

```
90_development_archive/
├── 01_early_transmission_models/    <- earliest weekly-transmission prototypes (pre-"renewal_*" numbering)
├── 02_historical_renewal_versions/  <- v1 through v5.0, the frozen renewal-model version lineage
├── 05_legacy_vaccine/               <- superseded vaccine-simulator prototypes
├── 06_legacy_susceptibility/        <- original (pre-shape-v2) dynamic annual-FOI/susceptibility track
├── 07_old_climate_experiments/      <- superseded future-climate-projection pilots (incl. WorldClim.R)
├── 08_old_stan_models/              <- legacy .stan sources + their compiled-model caches (see MODEL_REGISTRY.md)
└── 09_misc_archive/                 <- one-off exploratory/descriptive scripts, superseded by current pipelines
```

Corresponding outputs live under `03_Output/90_development_archive/`
(figures/tables by version name, `model_fits/` by version name).

## Known exceptions -- current pipelines that DO source archived code

The rule above is violated in three places, all pre-dating this
reorganisation and left as-is (not a new issue introduced by this
cleanup, and fixing it -- extracting shared logic into `00_shared/` --
would be a functional restructuring beyond a structural relocation, so it
is reported here rather than silently changed):

- `07_national_pipeline/10_grid_susceptibility_reconstruction/01_reconstruct_brazil_grid_susceptibility.R`
  sources `09_misc_archive/15_build_brazil_uf_foi_ensemble.R` directly.
- `07_national_pipeline/11_annual_foi_shape_v2/01_prepare_annual_foi_shape_v2.R`
  and `02_fit_annual_foi_shape_v2.R` both source
  `06_legacy_susceptibility/prepare_dynamic_annual_foi_data.R` directly
  (for its `build_dynamic_annual_foi_data()` helper).
- `00_shared/legacy_model_functions/22_v4_3.../scripts/04_mode_audit_analysis.R`
  and `00_shared/legacy_model_functions/27_v4_8.../scripts/02_plot_v4_8_sweep.R`
  each read one archived fit/output (`21_v4_2_q_juazeiro_serology/` and
  `26_v4_7_fixed_q_full_period/` respectively) purely for external
  cross-version comparison plots -- not required for the current model to
  fit, only for these specific comparison figures.

If you are extending any of these three scripts, be aware they will break
if `06_legacy_susceptibility/`, `09_misc_archive/`, or the two named
version folders are ever deleted (as opposed to merely reorganised within
the archive).

## Stan-path staleness inside archived scripts

Archived model-development scripts that reference a Stan source file by
its old `02_Script/stan/<name>.stan` path were **not** rewritten as part
of separating Stan sources from outputs (see
`00_shared/stan/MODEL_REGISTRY.md`) -- that directory no longer exists.
Their own `.stan` source now lives at
`08_old_stan_models/<name>.stan`. This does not affect any current
pipeline (none of these scripts are on a current execution path); it only
means an archived script cannot be re-run verbatim without first updating
that one path by hand.

## Why some version-numbered folders live in `00_shared/` instead

Four folders that look like they belong here by numbering
(`17_v4_0_minimal_no_vaccine`, `22_v4_3_short_2015_2019_q_calibration`,
`27_v4_8_seeded_recurrence`, `28_v4_9_hierarchical_seasonality`) are
**not** here -- they live under
`02_Script/00_shared/legacy_model_functions/` because every state's
Stage 02 transmission-fitting script still sources them. See that
folder's own `README.md`.
