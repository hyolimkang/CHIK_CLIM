# Dynamic annual FOI prototype: provenance audit

## Decision

**STOP: do not fit `chik_dynamic_annual_foi.stan` using
`brazil_uf_foi_weighted.csv` as a fixed annual infection hazard.**

The available file is useful as a spatially smoothed chikungunya exposure or
transmission-burden index.  The material preserved in this repository does
not support interpreting every value in that surface as a realized,
state-specific *long-term mean annual infection hazard*.  The requested
dynamic model would require exactly that interpretation in order for the
constraint

`mean_y(lambda[i, y]) = baseline_foi[i]`

to have a defensible epidemiological meaning.

This is a scientific/provenance stop, not a Stan or data-engineering failure.

## Local provenance chain

1. `02_Script/90_exploratory/14_build_uf_foi_weighted.R` creates
   `01_Data/brazil_uf_foi_weighted.rds` and
   `03_Output/tables/brazil_uf_foi_weighted.csv`.
2. That script does **not** estimate FOI. It loads the external gridded object
   `01_Data/allfoi_s1.RData`, spatially joins grid-cell points to Natural
   Earth Brazilian UF polygons, and calculates
   `weighted.mean(foi_mid, weight = tot)` within each UF.
3. The requested Ceará baseline is therefore exactly the population-weighted
   average of that grid surface: `foi_mid_wmean = 0.01216802` for `uf = CE`.
   `foi_lo_wmean` and `foi_hi_wmean` are weighted averages of grid-cell point
   bounds, not a propagated state-level interval; the aggregation script
   explicitly documents this limitation.

## What the upstream FOI inputs represent

The locally retained upstream construction material is:

- `01_Data/chik_foi.csv`, a table of location-level `mfoi`, `mfoi_lo`,
  `mfoi_hi`, model type, seroprevalence category, and number tested; and
- `01_Data/serology/markdownLAC.Rmd` (dated 2023-07-06), which estimates
  those quantities from age-seroprevalence data using JAGS catalytic models;
- `01_Data/current_covariates_baseline.R` and archived spatial scripts, which
  use `mfoi` as a response for spatial/covariate prediction.

For a constant catalytic model, the retained code uses
`1 - exp(-lambda * age)`. With age measured in years, `lambda` has an annual
hazard interpretation for that particular model.

However, the same input table also contains epidemic models. For example, the
Ceará/Quixadá model in `markdownLAC.Rmd` represents seropositivity through
terms such as `1 - exp(-lambda1 - lambda2 - lambda3)` for successive epidemic
exposures. Those `lambda` values are event-specific cumulative forces, not
annual rates. The `chik_foi.csv` input mixes `constant` and `epidemic` model
labels and multiple `epimodel` values. The retained spatial workflow then
uses the resulting `mfoi` quantities together as a single response surface.

Consequently, the resulting grid value is not documented as a common-period,
annual hazard. It is a spatially predicted quantity informed by heterogeneous
serology-derived exposure estimates and environmental/demographic covariates.

## Years, cases, reporting, and population weighting

- The direct serology models use study-specific historical intervals. The
  Ceará example includes breakpoints through 2019. There is no common calendar
  period for all grid predictions in the retained materials.
- The spatial covariate workflow uses static or multi-year covariates (for
  example, several 2010--2020 climate summaries and a 2020 vector layer), not
  an annual 2014--2025 surveillance time series.
- No retained upstream FOI-construction script reads SINAN, `cases_notified`,
  `cases_confirmed`, or a Brazilian state-year surveillance panel. Thus use of
  the project's 2014--2025 surveillance data is **not evidenced locally**.
  The complete original global mapping work package and its data dictionary
  are not present, so absence from the retained scripts cannot prove that no
  external surveillance input was used.
- The retained catalytic code conditions on seropositivity and sample size;
  it does not contain an infection-to-case reporting/detection model.
- `14_build_uf_foi_weighted.R` applies population weighting only at the final
  Brazilian UF aggregation step, with the grid population column `tot`.
  This does not turn the spatial prediction into an observed UF-level annual
  hazard.

## Information reuse and Ceará serology

The upstream input table contains Brazilian location-level serology inputs,
including Juazeiro do Norte (`study_no = 165`, 404 tested) and Quixadá/Ceará
(`study_no = 124`). The gridded baseline can therefore indirectly reuse local
serological information also held in this repository. This would make a
dynamic fit exploratory even if the annual-hazard interpretation were
resolved. Municipality observations must not be treated as independent,
state-representative validation after conditioning on a surface informed by
them.

## Required resolution before an annual dynamic FOI fit

Provide one of the following before implementation resumes:

1. the original global FOI model documentation/posterior defining
   `allfoi_s1$foi_mid` as a common-period annual hazard, including its data
   sources and reporting treatment; or
2. a separately estimated state-level baseline with an explicit annual-hazard
   definition and non-overlapping calibration evidence; or
3. approval to reframe the proposed model so this surface is used only as a
   covariate/prior scale index, **not** as a constraint on the arithmetic mean
   annual infection hazard.

Until then, no Stan model, annual fit, weekly susceptibility post-processing,
or national expansion has been created for this proposal.
