# v3.0 spatial renewal implementation report

## Files created

- 02_Script/stan/renewal_ceara_v3_0_spatial.stan
- 02_Script/40_renewal_model/prepare_renewal_v3_0_spatial_data.R
- 02_Script/40_renewal_model/fit_renewal_v3_0_spatial.R
- 02_Script/40_renewal_model/diagnose_renewal_v3_0_spatial.R

## Differences from v2.2

- Replaces the single state S/U/X trajectory with one trajectory per municipality.
- Applies the v2.2 demographic recursion independently to every municipality.
- Replaces state-to-site serology offsets with direct municipality-window immune proportions.
- Uses one shared weekly log-R0 random walk plus non-centred municipality R0 offsets.
- Uses one common reporting level (p_symp * rho_sym_global); v2.2's Ceará reporting offset and trend are omitted.

## Data and deterministic checks

- Fitting period: 2015-01-04 to 2019-12-29; 261 complete weeks.
- Ceará panel: 48024 municipality-week rows, 184 municipalities. Stan data: 13 municipalities.
- Maximum demographic identity error |S + U - N_start| in deterministic recursion: 8.38e-09.
- Maximum absolute difference between summed municipality population and v2.2 state stock: 1.415%.
- Births: SINASC individual records by maternal residence, weekly aggregation.
- Deaths: IBGE Ceara annual all-cause deaths, day-weighted then allocated by municipality start-week population share.

## Serology mapping

- Juazeiro do Norte (IBGE municipality code 230730): 103/404, 2018-06-03 to 2018-12-30.
- Quixada (IBGE municipality code 231130): 289/409, 2018-06-03 to 2019-12-29.
- Fortaleza is retained only for an external posterior predictive comparison and is not in the likelihood.

## Model size and compilation

- Estimated parameter dimension: 375 (104 municipality seed hazards, 253 shared-week R0 states, 13 municipality R0 standard-normal effects, and 5 scalar parameters).
- Stan compilation: successful.

## Pilot HMC

- Pilot HMC (debug subset only): divergences 34; maximum-treedepth hits 5; minimum E-BFMI 0.691.
- One short chain is a computational smoke test and has no convergence claim.

## Remaining warnings

- Municipal all-cause deaths are a transparent population-share allocation of the documented Ceará annual total, not municipality-specific mortality observations.
- Municipality annual population stocks are linearly interpolated between 1-January annual anchors; their summed stock is checked against the v2.2 state demographic series.
- Reporting, initial infections, and susceptibility remain potentially confounded. Compilation or a short pilot is not scientific validation.
- A future all-Ceará 2015-2019 production fit must be assessed with multi-chain convergence, posterior predictive, serology, and sensitivity diagnostics before interpretation.
