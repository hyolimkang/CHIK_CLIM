# Brazilian UF long-term annual chikungunya FOI aggregation

## Purpose and scope

This workflow reconstructs a Brazilian state (UF)-level **long-term average
annual chikungunya force of infection (FOI)** from the original 5-km `allfoi`
ensemble.  It is an aggregation workflow only.  It does not fit a dynamic FOI,
susceptibility, surveillance, climate, renewal, reporting, or reproduction
number model.

Run:

```r
source("02_Script/90_exploratory/15_build_brazil_uf_foi_ensemble.R")
run_brazil_uf_foi_aggregation()
```

The script is project-root aware, so it can be run from either the outer Git
repository or the inner `CHIK_CLIM` R project.

## Inputs and geographic assignment

- Original grid predictions: `01_Data/allfoi_s1.RData`, object `allfoi`.
- Brazil selection: exactly `iso3 == "BRA"`; no country-name fallback is used.
- Population weight: grid-cell `tot`.
- FOI ensemble: `foi1` through `foi100`; cell-level `foi_mid`, `foi_lo`, and
  `foi_hi` are retained only for a comparison audit.
- Official geography: `01_Data/ibge_muni_polygons.rds`, created from official
  IBGE municipality boundaries. Municipal polygons are dissolved to the 27
  IBGE UF codes before the point-in-polygon assignment.

Grid centroids are assigned to a UF by an `st_within()` spatial join. A centroid
that is outside every official UF polygon is retained and assigned to its
nearest UF, with both the number of such assignments and the maximum nearest
distance recorded. Ambiguous assignments cause the workflow to stop. The
`near_uf_boundary` audit flag denotes a cell with an eight-directionally
adjacent 5-km grid cell assigned to a different UF; it does not affect any
estimate or coverage calculation.

## Ensemble-first aggregation

For ensemble member `m` and UF `u`, the primary arithmetic annual-FOI summary
is calculated before summarising across ensemble members:

\[
\bar\lambda_{u,m} =
  \frac{\sum_{g \in u} N_g\lambda_{g,m}}
       {\sum_{g \in u} N_g},
\]

where `N_g = tot`. A cell with a missing/non-finite member prediction is not
treated as zero: it is excluded from that member's denominator, and the
member-specific population coverage is reported.

An alternative population-equivalent annual FOI is also computed for each
member:

\[
p_{g,m} = 1 - \exp(-\lambda_{g,m}),\qquad
\bar p_{u,m} = \frac{\sum_g N_gp_{g,m}}{\sum_gN_g},\qquad
\lambda^{eq}_{u,m}=-\log(1-\bar p_{u,m}).
\]

The final summaries provide median, mean, SD, 2.5%, 25%, 75%, and 97.5%
quantiles across the 100 state-level member estimates for both definitions.
They are not derived by aggregating `foi_mid`, `foi_lo`, or `foi_hi` cell
summaries.

## Outputs

- `01_Data/brazil_uf_foi_ensemble.csv`: 2,700 rows (27 UFs × 100 members),
  containing `foi_pwmean` and `foi_equiv`.
- `01_Data/brazil_uf_foi_longterm.csv`: one row per UF, including both
  ensemble distributions and population-coverage fields.
- `03_Output/tables/brazil_uf_foi_aggregation/aggregation_audit.csv`: Brazil
  cell, prediction availability, fallback, boundary-adjacency, and coverage
  totals.
- `03_Output/tables/brazil_uf_foi_aggregation/population_coverage_by_member.csv`:
  member- and UF-specific valid population denominators.
- `03_Output/tables/brazil_uf_foi_aggregation/uf_cell_assignment_audit.csv`:
  direct versus nearest-fallback assignment counts by UF.
- `03_Output/tables/brazil_uf_foi_aggregation/old_vs_ensemble_first_aggregation.csv`:
  a clearly labelled audit comparison with the old, non-primary aggregation of
  cell summary fields.
- `03_Output/tables/brazil_uf_foi_aggregation/ceara_foi_aggregation_audit.csv`:
  the Ceará row from that comparison.
- `03_Output/figures/brazil_uf_foi_aggregation/brazil_uf_foi_aggregation_diagnostics.pdf`:
  coverage, arithmetic-versus-equivalent FOI, uncertainty, geographic, and
  old-versus-new audit figures.

## Interpretation boundary

These outputs describe long-term average annual infection hazards inherited
from the original `allfoi` prediction ensemble. They must not be interpreted
as a dynamic annual FOI series, current incidence, attack rate, susceptibility,
or a renewal-model transmission parameter.
