# Brazil First-Epidemic Climate–Transmission Pilot (Temperature + Precipitation Only)

## Scope and design

Population: the chronologically **first** qualifying chikungunya epidemic
in each of the 27 Brazilian UFs (`wave_order == 1` among
`major_epidemic_primary` waves in the frozen national wave census).
Transmission-potential proxy: **early-wave Re**, reused from the
already-frozen, unmodified early-Re Stan pipeline
(`01_fit_brazil_major_wave_early_re.R`) — no refitting, no new
estimation method. Susceptibility (S(t) from the v4.9 long-term
reconstruction) is **not used anywhere** in this analysis.

**Implementation note on "weekly" Re(t)** (documented per instruction to
flag any necessary interpretation): the frozen early-Re method estimates
one constant Re per (wave, window-length), for window lengths 4, 6 and 8
weeks — it does not produce a continuously time-varying weekly
trajectory, and building one would require a new estimation method,
which was out of scope. To honor "reuse the existing method exactly"
while still producing more than one time-anchored observation per
episode, each of the three already-fitted window-length Re estimates is
treated as one early-growth observation, anchored to the **last week of
its window** (`onset_week + (window_weeks - 1)` weeks). This gives up to
3 weekly-anchored Re observations per eligible UF, spanning roughly a
month within the same early-growth phase — all drawn from the exact
frozen posteriors, with zero re-estimation.

## 1. Eligibility (Section 1)

27/27 UFs assessed; **18 eligible**. `RS` has no qualifying major
epidemic at all in the national census. The other 8 exclusions (AC, AM,
AP, DF, MS, RO, SE, TO) all fail the pre-existing, unmodified
`re6_analysis_usable` gate (`re6_prefit_eligible & re6_fit_pass`),
overwhelmingly for "fewer than 30 raw cases" in the early window. See
`FIRST_EPISODE_CLIMATE_ELIGIBILITY.csv`.

## 2–4. Weekly Re + climate covariates

54 weekly-Re observations (18 UFs x 3 window-anchors). Climate:
population-weighted mean Tmean (deg C) and population-weighted total
weekly PRCP (mm), from the already-built, unmodified
`brazil_chik_uf_weekly_climate.csv` (built by looping the existing
`build_state_weekly()` over all 27 UFs against the national ERA5-Land-
derived DLNM panel). Three pre-specified lag windows (fixed before
inspecting results): primary (temp t-2:t, precip t-6:t-2), sensitivity A
(temp t-1:t, precip t-4:t-1), sensitivity B (temp t-3:t-1, precip
t-8:t-3). See `first_epidemic_climate_training_data.csv` (primary spec)
and `FIRST_EPISODE_CLIMATE_ALL_SPECS.csv` (all three).

## 5–7. Primary GAM and prediction grid

`log(Re) ~ s(temperature, k=4) + s(precipitation, k=4) + s(UF, bs="re")`,
`mgcv::gam(method="REML")`, N=54, 18 UF clusters.

| term | edf | p-value |
|---|---|---|
| s(temperature) | 1.00 | 0.618 (not significant) |
| s(precipitation) | 1.00 | **0.012** |
| s(UF) random effect | 15.37 (of 17 max) | <2e-16 |

Deviance explained = 88.1%, adj. R² = 0.823 — but this is dominated by
the UF random effect (see concurvity below), not by climate. Both
climate smooths reduced to k=1 (linear) — no evidence of nonlinearity
supported at this sample size. Residuals: SD=0.091, Shapiro-Wilk
p=0.577 (no departure from normality detected).

**Concurvity is severe**: UF vs. temperature ≈ 0.93-0.95, UF vs.
precipitation ≈ 0.81-0.84 (worst-case measure). Because each UF
contributes only one first epidemic, "UF identity" and "this UF's
typical early-epidemic climate" are almost the same information — the
UF random effect and the climate smooths compete for the same
variance. The reported deviance explained/adj. R² should **not** be read
as climate explaining 88% of Re variation; most of that is the UF
intercept itself, which is collinear with climate.

Prediction grid restricted to the observed range (temperature
21.6-28.4°C, precipitation 12.7-341.2 mm): `climate_R0_prediction_grid.csv`.
Median-Re point-estimate model and Monte Carlo (Section 6) versions both
saved (`first_epidemic_climate_gam.rds`, `climate_R0_prediction_grid_mc.csv`).

## 6. Re uncertainty propagation

200/200 Monte Carlo iterations succeeded, resampling one posterior Re
draw per (UF, window-anchor) from the already-saved early-Re posterior
draws and refitting the GAM each time (no new Re estimation). Across
iterations: temperature and precipitation smooths stayed linear (edf
medians =1, up to ~2.1/1.9 at the extremes), median deviance explained
across the resampled fits = 66.6% (lower than the single median-Re fit's
88.1%, as expected once full Re uncertainty is included).

## 8. Figures

`FIRST_EPISODE_RE_BY_UF.png`, `CLIMATE_R0_TEMPERATURE_RESPONSE.png`,
`CLIMATE_R0_PRECIPITATION_RESPONSE.png`, `CLIMATE_TRANSMISSION_SURFACE.png`,
`FIRST_EPISODE_RE_PREDICTED_VS_OBSERVED.png`, plus a diagnostic
`FIRST_EPISODE_LOO_UF_SENSITIVITY.png`.

## 9. Diagnostics

- N = 54 observations, 18 eligible UFs, 3.0 obs/UF.
- Temperature range 21.58-28.42°C; precipitation range 12.68-341.24 mm (5-week cumulative).
- **Lag sensitivity** (`first_epidemic_climate_lag_sensitivity.csv`):
  precipitation is significant under primary (p=0.012) and sensitivity B
  (p=0.0012) but **not** under sensitivity A (p=0.103); temperature is
  never significant under any of the three specs (p=0.21-0.62).
- **Leave-one-UF-out** (`first_epidemic_climate_loo_uf_sensitivity.csv`,
  figure above): temperature-effect significance is stable (never
  significant in any of the 18 refits). Precipitation stays significant
  (p<0.05) in **17 of 18** refits — it loses significance only when PI is
  excluded (p=0.083). The predicted temperature-response curve itself
  shifts by a median of 6.3% and a maximum of 10.3% of its value when any
  single UF is dropped — a real but modest sensitivity, not evidence that
  one or two UFs are driving the whole result.

## Answers to the nine required questions

**1. How many UFs provided usable first-epidemic early-Re information?**
18 of 27.

**2. What temperature range is represented?**
21.58°C to 28.42°C (population-weighted mean, weeks t-2:t relative to
each Re anchor week).

**3. What precipitation range is represented?**
12.68 mm to 341.24 mm (population-weighted total, weeks t-6:t-2, 5-week
cumulative).

**4. Is there evidence of a nonlinear temperature-transmission relationship?**
No. The temperature smooth reduces to edf=1 (linear) under all three lag
specifications, and the linear temperature effect is not statistically
significant in any of them (p=0.21-0.62). Direction is consistently
positive (higher temperature -> modestly higher predicted Re) across the
full model and every leave-one-UF-out refit, but the data cannot
distinguish this from no effect.

**5. Is there evidence of a precipitation-transmission relationship?**
Yes, but with a caveat. Precipitation is significantly, positively
associated with early Re in the primary specification (p=0.012) and
robust to dropping all but one UF (PI). However it loses significance
under the sensitivity-A lag window (p=0.103), so the finding is
lag-window-sensitive, not airtight.

**6. How stable are these relationships to plausible lag windows?**
Moderately stable for precipitation (significant in 2/3 pre-specified
windows, always positive in sign, effect size varies), and consistently
null for temperature (never significant, always small positive edf=1
slope) across all 3 windows.

**7. Are results dominated by one or two UFs?**
No single UF dominates by itself (precipitation significance survives
17/18 leave-one-out refits; curve shifts are modest, median 6.3%/max
10.3%). But the model as a whole is dominated by the **UF random
effect**, not by any one UF — 15.4 of a possible 17 degrees of freedom go
to UF identity, and that identity is highly concurvous with climate
itself. This is a structural feature of having only 18 independent
clusters with essentially one climate "value" each, not a one-or-two-UF
outlier problem.

**8. Over which climate range is R0_climate_hat supported by observed data?**
Only the rectangle actually spanned by the 54 training observations
(see the white points overlaid on `CLIMATE_TRANSMISSION_SURFACE.png`):
approximately 21.6-28.4°C x 13-341 mm, and even within that rectangle the
data are unevenly distributed (dense in the 25-27°C / 100-300 mm region,
sparse at the temperature extremes). The prediction grid is deliberately
restricted to this observed range; extrapolation beyond it is not
supported.

**9. Is the model sufficiently stable to proceed to recurrent-epidemic prediction?**
Conditionally, with an important caveat that must travel with any later
use. The model is *computationally* stable (clean REML fit, no
convergence issues, modest leave-one-out sensitivity) and gives a
*plausible, sign-consistent* precipitation association, but:
(a) severe UF-climate concurvity means the climate effect and the
"which-UF" effect are not cleanly separable at this sample size,
(b) temperature shows no detectable effect, and
(c) precipitation's significance is not robust to one of the two
plausible sensitivity lag windows.
This supports treating the fitted mapping as a **provisional, precipitation-
led empirical association** suitable for later comparison against
recurrent-epidemic weeks (where within-UF, within-season climate variation
will help separate climate from UF identity) — but not yet as a
validated predictive R0(climate) function.

## Outputs

Tables: `FIRST_EPISODE_CLIMATE_ELIGIBILITY.csv`,
`FIRST_EPISODE_WEEKLY_RE.csv`, `first_epidemic_climate_training_data.csv`,
`FIRST_EPISODE_CLIMATE_ALL_SPECS.csv`,
`first_epidemic_climate_model_diagnostics.csv`,
`climate_R0_prediction_grid.csv`, `climate_R0_prediction_grid_mc.csv`,
`first_epidemic_climate_gam_mc_summary.csv`,
`first_epidemic_climate_lag_sensitivity.csv`,
`first_epidemic_climate_loo_uf_sensitivity.csv`,
`first_epidemic_climate_diagnostics_overview.csv`.

Model objects: `first_epidemic_climate_gam.rds`.

Figures: `FIRST_EPISODE_RE_BY_UF.png`,
`CLIMATE_R0_TEMPERATURE_RESPONSE.png`,
`CLIMATE_R0_PRECIPITATION_RESPONSE.png`, `CLIMATE_TRANSMISSION_SURFACE.png`,
`FIRST_EPISODE_RE_PREDICTED_VS_OBSERVED.png`,
`FIRST_EPISODE_LOO_UF_SENSITIVITY.png`.

## STOP HERE

Per instruction: recurrent epidemics are not analysed, reconstructed
susceptibility is not combined with this model, and no vaccination
simulation is run in this pilot.
