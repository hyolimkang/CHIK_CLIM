# Recurrent epidemic climate-susceptibility analysis

## Design

For each recurrent episode, the mechanistic prediction is
`Re_pred = R_climate_hat * (S_pre / N)`.

- `R_climate_hat` is the population-level prediction from the already fitted
  first-epidemic climate GAM, evaluated at the recurrent episode's climate.
  The UF random effect is excluded. It is a climate-defined transmission
  potential, not an observed state-specific R0.
- `S_pre / N` comes from the existing long-term state reconstruction in the
  week immediately before the early-growth window.
- `Re_obs` is the independently estimated 6-week early reproduction number
  from the national case-growth pipeline. It is not derived from the long-term
  susceptibility reconstruction.

The figures do not refit the climate GAM or any long-term model.

## Episode groups

Primary figures and validation use 16 recurrent episodes: Bahia (8), Rio de
Janeiro (5), and Mato Grosso (3). Bahia is serology anchored; Rio de Janeiro
and Mato Grosso remain partially identified case-based reconstructions.

Ceara is shown only as a conditional fixed-q sensitivity analysis (six
episodes across q = 0.05, 0.10, 0.15, and 0.20). Pernambuco remains excluded
because its absolute susceptibility is not robustly identified.

## Figure interpretation

`RECURRENT_CLIMATE_SUSCEPTIBILITY_PHASE_PLANE.png` is the primary phase plane.
Point area is **reported cumulative incidence per 100,000**, not raw case
counts, so the displayed epidemic-size measure is comparable across states.
Point fill is `Re_obs` and point shape denotes state. The threshold is
`R_climate_hat * S_pre/N = 1`.

`RECURRENT_CLIMATE_SUSCEPTIBILITY_PHASE_OUTCOMES.png` shows the two outcome
relationships directly: observed early Re versus pre-epidemic susceptibility,
and reported cumulative incidence versus pre-epidemic susceptibility. The
linear summaries are descriptive and unadjusted; they are not causal effect
estimates. The outcome panel carries the relevant 95% intervals for Re and S.

`RECURRENT_CE_FIXED_Q_SENSITIVITY_PANEL.png` keeps all Ceara q scenarios
separate. It must not be interpreted as a uniquely reconstructed susceptibility
trajectory.

Across the selected episodes, `R_climate_hat` is tightly clustered (about
2.0-2.6), while `S_pre / N` covers a much wider range (about 55%-100%). Every
selected major recurrent episode lies above the theoretical growth threshold.
That pattern is compatible with the mechanism, but is not a prospective test:
the episodes were selected because they grew into major waves.

## Mechanistic validation

`RECURRENT_RE_PREDICTED_VS_OBSERVED.png` uses posterior medians for the 16
primary episodes and an equal-scale y = x reference line.

| model | episodes | correlation with Re_obs | MAE | RMSE | Re > 1 agreement |
|---|---:|---:|---:|---:|---:|
| climate only (`R_climate_hat`) | 16 | 0.411 | 0.783 | 0.869 | 100% |
| climate x susceptibility (`Re_pred`) | 16 | 0.370 | 0.733 | 0.827 | 100% |

The calibration regression for `Re_obs ~ Re_pred` has slope 0.958 and
intercept 0.074. Adding susceptibility improves point accuracy modestly (MAE
by 0.050 and RMSE by 0.042), but does not improve rank correlation in this
small sample. The perfect Re > 1 agreement is uninformative because all
episodes were selected major waves.

This remains mechanistically consistent predictive evidence, not causal proof
of climate effects or a validated forecasting model. Reported incidence is a
case-based outcome and can still reflect reporting and surveillance differences
between states and years.

## Outputs

Tables: `RECURRENT_EPISODES.csv`,
`RECURRENT_CLIMATE_SUSCEPTIBILITY_EPISODES.csv`, and
`RECURRENT_RE_VALIDATION_METRICS.csv`.

Figures: `RECURRENT_CLIMATE_SUSCEPTIBILITY_PHASE_PLANE.png`,
`RECURRENT_CLIMATE_SUSCEPTIBILITY_PHASE_OUTCOMES.png`,
`RECURRENT_CE_FIXED_Q_SENSITIVITY_PANEL.png`, and
`RECURRENT_RE_PREDICTED_VS_OBSERVED.png`.
