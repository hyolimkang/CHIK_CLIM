# All-opportunities climate-susceptibility analysis

## Purpose

This analysis expands the recurrent-episode phase plane to every eligible
post-first-wave, climate-supported sliding six-week opportunity in Bahia (BA),
Rio de Janeiro (RJ), and Mato Grosso (MT).  It is descriptive retrospective
evidence about whether a favourable climate-defined transmission potential is
expressed differently at different reconstructed susceptible fractions.

The quantities are kept distinct:

- `R_climate_hat` is the population-level estimate from the pre-existing
  first-epidemic GAM, with its UF random effect excluded.  It is climate
  transmission potential, not an observed state-specific reproduction number.
- `S_pre/N` is the posterior median susceptible proportion from the existing
  state reconstruction in the week immediately before the six-week growth
  window.
- `R_e_pred = R_climate_hat * S_pre/N` is a mechanistic derived quantity, not
  a separately fitted or observed `R_e`.
- When a future major wave occurs, `R_e_obs` is its independently estimated
  six-week early reproduction number from the frozen case-growth pipeline.

## Opportunity and outcome definitions

An opportunity is retained when all of the following hold:

1. It is a weekly sliding six-week window after the first frozen major wave
   has ended in the state.
2. Its GAM climate summaries and reconstructed `S_pre/N` are available.
3. Its temperature and precipitation combination is within the first-episode
   GAM training support.

For the primary outcome, follow-up starts in the week after the opportunity
anchor and lasts eight weeks.  A positive outcome is an onset in that interval
of a `major_epidemic_primary` wave from the already frozen national census.
That census requires both at least 100 reported cases and cumulative incidence
of at least 5 per 100,000; its onset rule also requires a two-week sustained
crossing of the rise threshold in the smoothed series.  The six-week version
of the same onset outcome is retained in the output table as a shorter-horizon
sensitivity.

`Non-outbreak` means a climate-supported opportunity with complete eight-week
follow-up and no such frozen major-wave onset.  It does **not** mean no
transmission, no imported infections, or no minor/local rise in reported
cases.  Windows with incomplete eight-week follow-up are right-censored and
are excluded from outcome-coloured phase planes rather than counted as
non-outbreaks.

## Figures and tables

- `RECURRENT_ALL_OPPORTUNITIES_PHASE_PLANE_FINAL.png`: all complete opportunity
  windows.  Blue bins retain the count of every non-outbreak window.  Because
  weekly sliding windows can all point to the same future wave, each distinct
  future major wave is shown once: the latest eligible opportunity before its
  onset.  Its outline colour is observed early `R_e`, and its size is the
  linked future wave's reported cumulative incidence per 100,000.
- `RECURRENT_CLIMATE_SUSCEPTIBILITY_TRAJECTORIES.png`: for each state, climate
  potential, pre-window susceptibility, and their product over time.  Grey
  bands are frozen major-wave intervals and the lower-panel dashed line is
  `R_e_pred = 1`.
- `RECURRENT_CLIMATE_SUSCEPTIBILITY_PHASE_PLANE.png`: recurrent episodes only,
  upgraded with predicted-`R_e` contours at 1, 1.5, 2, and 2.5.
- `RECURRENT_CLIMATE_SUSCEPTIBILITY_OPPORTUNITIES_WEEKLY.csv`: analysis rows,
  including six- and eight-week outcome indicators, censoring status, and the
  linked future-wave severity fields.
- `RECURRENT_ALL_OPPORTUNITIES_OUTCOME_SUMMARY.csv`: descriptive count and
  predictor summaries by state and eight-week outcome category.

## Interpretation limits

The windows overlap weekly and therefore are not independent observations;
the figure is not a regression test or a calibrated forecast evaluation.
Climate predictors use the lags established in the first-episode GAM and
should not be interpreted as a new causal climate model.  A major wave needs
seeding/importation and appropriate surveillance in addition to a favourable
climate and susceptible pool; neither importation nor a detection mechanism is
modelled here.  Finally, the existing long-period susceptibility posterior is
a full-history, retrospective (smoothed) reconstruction, so `S_pre/N` has
future-information leakage for a prospective-prediction interpretation.  BA
is serology anchored, whereas RJ and MT remain partially identified
case-based reconstructions; results must retain those differences in
evidential weight.
