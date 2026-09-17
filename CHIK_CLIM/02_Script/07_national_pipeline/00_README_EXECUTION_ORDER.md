# NATIONAL PIPELINE

Consolidated from 7 formerly-separate folders under `40_renewal_model/`
(`10_episode_renewal`, `11_major_episode_renewal`,
`12_episode_susceptibility_reconstruction`, `13_national_episode_census`,
`14_national_wave_census_v2`, `15_national_wave_analysis`,
`16_grid_susceptibility_reconstruction`) plus `40_susceptibility/02_annual_foi_shape_v2`
-- no scientific logic changed, only location (git mv) and internal path
references. `15_national_wave_analysis` (26 files) was itself split into 4
thematic folders (06-09 below) by explicit decision, since it bundled 4
distinct sub-studies.

This pipeline does **not** use the state pipelines' 01-09 stage template --
it bundles several thematically distinct national-scope analyses rather
than one sequential fit/diagnose chain. Folder numbers reflect topic
grouping inherited from the source material, not execution order; see
`99_run_national_pipeline.R`'s `EXECUTION_ORDER` for the actual
dependency-respecting run order.

| Folder | Was | Content |
|---|---|---|
| `01_episode_renewal` | `10_episode_renewal` | Ceara early-phase episode-level Re (frozen K-episode renewal Stan program) |
| `02_major_episode_renewal` | `11_major_episode_renewal` | Ceara major-outbreak (6-week window) episode-level Re |
| `03_episode_susceptibility_reconstruction` | `12_episode_susceptibility_reconstruction` | Ceara episode-sequential susceptibility pilot (prepare/fit/diagnose) |
| `04_national_episode_census` | `13_national_episode_census` | Frozen Brazil-wide chikungunya episode census |
| `05_national_wave_census_v2` | `14_national_wave_census_v2` | Frozen Brazil-wide chikungunya wave census (v2) -- everything downstream reads this |
| `06_national_wave_analysis` | `15_national_wave_analysis` (00, 01-08) | Nationwide early-Re + annual-SHAPE retrospective susceptibility: fit, weekly reconstruction, master table, diagnostics |
| `07_provisional_climate_susceptibility_checkpoint` | `15_national_wave_analysis` (09-14) | Exploratory: does FOI-informed susceptibility add explanatory power for early Re beyond climate alone? (GAM + Monte Carlo) |
| `08_first_epidemic_climate_transmission` | `15_national_wave_analysis` (15-21) | First-epidemic climate-transmission pilot across all 27 UFs (see its own `FIRST_EPIDEMIC_CLIMATE_TRANSMISSION_REPORT.md`) |
| `09_recurrent_climate_susceptibility_phase` | `15_national_wave_analysis` (22-25) | Recurrent-epidemic climate x susceptibility phase analysis for BA/RJ/MT |
| `10_grid_susceptibility_reconstruction` | `16_grid_susceptibility_reconstruction` | Bottom-up (grid-cell) reconstruction of Brazil-wide susceptibility |
| `11_annual_foi_shape_v2` | `40_susceptibility/02_annual_foi_shape_v2` | State-specific batch extension of `chik_dynamic_annual_foi_shape_v2.stan` (SHAPE/NB2 sensitivity variants) |

## Pre-existing, not-to-be-fixed finding

`03_episode_susceptibility_reconstruction`'s HMC pass/fail gate reports
only 3 of 7 episodes passing the standard gate (divergences/treedepth/
Rhat/ESS/BFMI). This is a property of the frozen pilot-v1 fit itself, not
introduced by or fixed in this refactor -- it is reported as-is by
`00_run_folder.R`.

```r
# every output already exists on disk -- this just reuses everything and reports so
source("02_Script/07_national_pipeline/99_run_national_pipeline.R")

# run only specific folders (comma-separated, matching EXECUTION_ORDER names)
Sys.setenv(NATIONAL_RUN_ONLY = "06_national_wave_analysis,07_provisional_climate_susceptibility_checkpoint")
source("02_Script/07_national_pipeline/99_run_national_pipeline.R")

# allow rebuilding/refitting whatever is actually missing (single project-wide gate)
Sys.setenv(NATIONAL_ALLOW_REFIT = "TRUE")
source("02_Script/07_national_pipeline/99_run_national_pipeline.R")
```

Outputs remain at their original locations under `03_Output/tables|figures/
{national_episode_census,national_wave_census_v2,national_wave_analysis,
national_grid_susceptibility,annual_foi_shape_v2}/`, `02_Script/stan/`
(episode/major-episode/episode-susceptibility/annual-foi-shape-v2 fits), and
`03_Output/{tables,figures}/{episode_susceptibility,renewal_episode_re,
renewal_major_episode_re}/` -- this refactor did not move or rename any
existing result.

See `REFACTOR_MIGRATION_MAP.csv` (in `02_ceara_pipeline/`, extended with a
National section) for the full old -> new script mapping.
