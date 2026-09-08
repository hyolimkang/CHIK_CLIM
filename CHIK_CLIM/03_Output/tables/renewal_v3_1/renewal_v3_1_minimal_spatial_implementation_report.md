# v3.1 minimal spatial renewal implementation report

## Model structure

- Six units: Fortaleza, Juazeiro do Norte, Quixada, and fixed high/moderate/low burden pooled strata.
- Fixed thresholds: high >= 250 total cases; moderate 25 to < 250; low < 25.
- The named municipalities are excluded from pooled strata before threshold assignment.
- Unit-specific S/U/X states use the v2.2 demographic recursion; statewide susceptibility is population-weighted and derived.
- Serology likelihood includes only Juazeiro and Quixada. Fortaleza is external only.
- Reporting is common across units: p_symp * rho_sym.

## Checks

- 184 municipalities collapsed to 6 units across 261 weeks.
- Maximum deterministic S + U accounting error: 8.38e-09.
- Maximum municipality-sum versus v2.2 state population gap: 1.415%.
- Unit assignment: C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM/03_Output/tables/renewal_v3_1_minimal_spatial_unit_assignment_high250_moderate25.csv.
- Stan parameter dimension: 312.
- Stan compilation: successful.

## Pilot diagnostics

- Four-chain pilot elapsed time: 123.6 seconds.
- Divergences: 27; treedepth hits: 0; maximum Rhat: 2.677.
- Minimum bulk ESS: 7.7; minimum tail ESS: 8.9; minimum BFMI: 0.331.
- This short pilot is a computational diagnostic, not a converged production analysis.

## Remaining limitations

- Municipality deaths remain an explicit population-share allocation of the Ceará annual total before aggregation, not observed municipality mortality.
- No serology anchor is created for Fortaleza or pooled units; their uncertainty should remain data-driven.
- Reporting, seeding, and susceptibility can remain confounded. The intended future test is robustness of the Ceará trajectory to refined pooling, not a preferred point estimate.
