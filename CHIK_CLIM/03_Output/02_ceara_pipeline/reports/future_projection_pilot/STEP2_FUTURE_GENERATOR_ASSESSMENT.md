# STEP 2 -- Stochastic future epidemic generator assessment (Ceara, 2026-2035, no vaccination)

Produced by `04_run_unit_tests.R`, `05_run_simulation_ensemble.R` and
`06_step2_figures.R` under `02_Script/02_ceara_pipeline/10_future_projection_pilot/`.
**Illustrative scenario simulation; no SSP forcing or vaccination applied.**
Pure-R simulation throughout -- no Stan model was fit or refit.

## 1. Exact equations implemented

**Future transmission potential** (Section 2.2):
```
log_R0_future(t) = alpha_R + year_effect_future(t) + A_year_future(t) *
                    (beta_sin1*sin1(t) + beta_cos1*cos1(t)) + beta_sin2*sin2(t) + beta_cos2*cos2(t)
                  + beta_T_eff * z_T_future(t) + beta_P_eff * z_P_future(t)
```
`t` continues the exact historical week index (574, 575, ... after the
historical 1..573), so the harmonic phase is unbroken across the
2025/2026 boundary -- verified numerically (unit test "2.2", 0.0 error)
to reproduce the historical R0(t) exactly when fed historical
coefficients/climate/year-effects.

**Future year structure** -- two arms:
- **Y0 (reference)**: `year_effect_future = 0`, `A_year_future = 1` for every future year (neutral/flat).
- **Y1 (resampled)**: `z_year ~ Normal(0,1)` and `z_season ~ Normal(0,1)` freshly drawn per future year (NOT copied from 2015-2025), using the model's own documented construction: `year_effect = 0.40*(z_year - mean(z_year))`, `A_year = exp(sigma_season_year_draw * z_season)`, where `sigma_season_year_draw` is carried from the SAME posterior draw as every other parameter used in that simulation.

**Future climate arms** (Section 2.3): C0 (`z_T=z_P=0`), C1 (whole historical-calendar-year blocks of the fitted model's own `z_T_anom`/`z_P_anom` resampled with replacement, seasonally aligned, never shuffled week-by-week -- preserves joint T-P structure since both come from the same original week), C2 (identical resampled block to C1, `beta_T_eff=beta_P_eff=0`).

**Stochastic renewal** (Section 2.5):
```
infectiousness(t) = sum_g w_g * X(t-g)            [SAME G=8 Gamma(4,2) kernel as the accepted fit]
lambda(t) = R0_future(t) * infectiousness(t) / N(t) + M(t) / N(t)
M(t) ~ Poisson(lambda_import)                      [stochastic mode; = lambda_import in deterministic mode]
p(t) = 1 - exp(-lambda(t))
X(t) ~ Binomial(round(S(t)), p(t))                 [stochastic mode; = S(t)*p(t) in deterministic mode]
```
Demographic update identical in form to the accepted model:
`S(t+1) = S(t)-X(t) + births(t) - deaths(t)*sfrac + reconciliation(t)*sfrac`, `sfrac = (S(t)-X(t)) / ((S(t)-X(t))+(U(t)+X(t)))`.

**Observation layer** (Section 2.7): `Y(t) ~ NegativeBinomial(mean = q*X(t), size = phi_obs_draw)`, q=0.05 fixed (same as the accepted fit). Latent `X` and reported `Y` are stored separately throughout; only `X` is treated as the primary biological output.

## 2. Assumptions and their justification

| Assumption | Justification / disclosure |
|---|---|
| Initial state = SAME posterior draw's terminal `S`, `U`, last 8 weeks of `X`, and every transmission/climate parameter | Preserves joint posterior correlation (2.1); verified by unit test "E" (`cor(alpha_R, S_terminal)` identical on a sample vs. the full posterior: -0.951 both). |
| Terminal-state bridge | The accepted Stan fit's `S[N]`/`U[N]` are PRE-update (the model's `if (t<N)` guard never applies the final week's demography). A dedicated `advance_terminal_state()` step completes that one missing week so the future horizon starts exactly at `N_end(N_hist)` -- verified to land within 3.9e-08 of the data value across 300 draws (unit test "2.1"). |
| Lifelong natural immunity (unchanged from the accepted model) | Inherited, not re-examined here. Over a 10-year horizon this is an increasingly consequential structural assumption -- chikungunya long-term immune duration is not perfectly established in the literature, and this pilot cannot speak to what happens if immunity actually wanes. |
| `S(t)` rounded to the nearest integer for the Binomial trial count | Population counts are in the millions; rounding error is <1 individual against counts of 10^6-10^7, i.e. negligible. Not separately bias-tested beyond this scale argument. |
| Importation `M(t)` enters only as a hazard addition (`M(t)/N(t)`), never as separately counted cases | Matches the accepted model's treatment of `imports_per_week` exactly; avoids double-counting (explicit design requirement). |
| `lambda_import` grid {0, 0.25, 1, 4}/week | A prespecified STRUCTURAL range, not fitted to any data. 1/week matches the accepted fit's fixed `imports_per_week`. |
| D1 demography: IBGE UF (2024 revision) annual growth/birth/death RATES applied to this model's OWN continuing population level | The only future demographic source in the repo. Only RATES are borrowed (never IBGE's absolute population level, which uses a different base/definition and would otherwise create an artificial break at the 2025/2026 boundary). |
| D0 demography: final historical year's own observed rates, frozen and repeated | Required sensitivity alternative; also the fallback if D1 were judged unsuitable. |
| Outbreak/wave detection | A simplified, single-series adaptation of the accepted national wave-census "primary" peak-trough definition (3-week smoothing, 30% trough, 5% onset threshold, major = >=100 reported cases AND >=5/100k), applied to simulated REPORTED cases only, never to latent infections directly. |

## 3. Mandatory validation results (Section 2.10)

All 15 checks **PASS** (`03_Output/02_ceara_pipeline/tables/future_projection_pilot/STEP2_unit_test_QA_summary.csv`):

| Check | Result |
|---|---|
| A. Historical replay (deterministic, historical forcing/seed) reproduces accepted S/U/X | **Exact**: 0.00e+00 relative error |
| B. Mass conservation S(t)+U(t)=N_start(t) | max abs error 2.4e-08 |
| C. Zero introduction + zero history -> zero infections (deterministic AND stochastic) | Exact: max X = 0 |
| F. Forcing identity: beta_T=beta_P=0 == z_T=z_P=0 | Exact: 0.00e+00 difference |
| 2.2. `compute_R0_future()` reproduces historical R0(t) exactly | Exact: 0.00e+00 |
| D. Stochastic replicates: distinct trajectories, fadeout permitted | 16/20 unique cumulative outcomes; 20/20 replicates show >=1 zero-infection week |
| E. Posterior parameters/terminal states propagate jointly | cor(alpha_R, S_terminal) identical on sample vs. full posterior (-0.951 both) |
| G. No future data leakage | Future horizon starts 2026-01-04, strictly after the 2025-12-21 historical end; `simulate_renewal_core()` takes no observed-future-data argument |
| H. Demographic plausibility | S,U,N all non-negative/positive across the historical replay |
| I. Reproducibility | Same seed -> identical trajectory; different seed -> different trajectory |
| 2.1. Terminal-state bridge lands exactly on N_end(N_hist) | max abs error 3.9e-08 across 300 draws |

**Code correctness is strongly supported by this battery**, in particular the EXACT (floating-point-only) historical replay -- the future engine is not an approximation of the accepted mechanism, it is the same mechanism continued forward in time plus the documented stochastic/structural extensions.

## 4. Simulation results (reference scenario: Y1 x C1 x D1, lambda_import=1, pilot scale 100 draws x 20 replicates = 2,000 simulations)

- 10-year cumulative infection incidence (per 100,000): min 126, p10 4,544, median **31,923**, p90 52,993, max 66,247 (FIG09).
- Wide, markedly right-skewed and multimodal outcome distribution (FIG06, FIG09) -- NOT a tight band around a single "expected" trajectory.
- 82,302 candidate waves detected by the adapted peak-trough rule across the 2,000 simulated reported-case series, of which 13,912 (16.9%) meet the "major epidemic" threshold (>=100 cases AND >=5/100k); the remainder are small stochastic fluctuations in low-incidence periods that the simplified single-series detector does not suppress as aggressively as the full national wave-census algorithm (outbreak-definition sensitivity, Section 2.9's explicit caveat -- this inflates the raw wave COUNT but the `major_epidemic` flag is the metric that should be used for recurrence-interval interpretation, not raw wave count).
- FIG07's three representative realisations illustrate the mechanism directly: a low-incidence realisation that never depletes susceptibles below ~60%; a median realisation with one sharp depletion event (~65%->~37% S/N) around 2032-2033 coinciding with a ~100,000/week peak; a high-incidence realisation with a single very large epidemic (~480,000/week peak) crashing S/N from ~65% to ~15% around 2034.

## 5. Sensitivity analyses

### 5.1 Future year effects (Y0 vs. Y1) -- the DOMINANT driver

| | median 10-year cumulative infections (lambda_import=1, C1, D1) |
|---|---|
| Y0 (reference/neutral) | **~8,100** |
| Y1 (resampled) | **~3,256,000** |

A ~400-fold difference in the MEDIAN outcome between the two year-effect
treatments, for the same climate arm and importation level. Under Y0,
R0(t) repeats an identical, bounded seasonal cycle every year (peak
~1.9-2.0; FIG05 top panel) and R_eff rarely sustains above 1 for long,
so the system stays close to a low endemic/fade-out equilibrium. Under
Y1, independently resampled year-to-year variation occasionally pushes
R0 to 4-8 in a given year (FIG05 bottom panel), which is sufficient to
ignite a large epidemic that depletes a large share of the susceptible
pool in a few months (FIG07). **This is the single most consequential
structural choice in the whole generator** -- far more consequential than
the climate arm or the importation-pressure choice below -- and should be
reported as such whenever these results are used or extended.

### 5.2 Climate coefficients (C0 vs. C1 vs. C2)

| | median 10-year cumulative infections (Y1, lambda_import=1, D1) |
|---|---|
| C0 (climate-neutral) | 3,478,942 |
| C1 (resampled weather variability) | 3,255,705 |
| C2 (resampled weather, beta forced to 0) | 1,938,566 |

Differences exist but are second-order relative to Section 5.1, consistent
with Step 1's finding that climate explains only ~1% of historical log
R0(t) variance against seasonality's ~80% and year effects' ~19%. The
C0-vs-C2 gap (both have no climate EFFECT, by construction, since C0 has
z=0 and C2 has beta=0) reflects only year-effect resampling noise between
independently drawn ensembles, not a climate signal -- a reminder that
differences of this size are within the structural noise floor of a
smoke-scale (10 draws x 5 replicates = 50 simulations per cell)
comparison and should not be over-interpreted.

### 5.3 Introduction pressure (lambda_import)

| lambda_import | Y1 median cum_10y | Y1 mean cum_10y |
|---|---|---|
| 0 | 10,198 | 806,373 |
| 0.25 | 2,579,253 | 2,451,602 |
| 1 (reference) | 3,255,705 | 2,729,799 |
| 4 | 2,869,432 | 2,523,434 |

**Important interpretive note**: `lambda_import=0` does NOT mean "no
possible future transmission" -- per the task's own explicit instruction
(2.1), each simulation's recent infection history (the last 8 historical
weeks of that posterior draw's `X`) is carried forward, not reset to
zero. So even with zero NEW weekly introductions after 2025, a handful
of draws still ignite a large epidemic from that historical carry-over
once a high-R0 year is resampled, while most draws instead see the
carry-over decay to extinction first -- producing the large mean/median
gap at lambda_import=0 (a genuinely bimodal outcome, not an error; this
is distinct from, and not contradicted by, unit test C's zero-
introduction check, which deliberately starts from a FRESH, zero
infection history to isolate the mechanism). Once any ongoing
introduction pressure exists (lambda_import>=0.25), large epidemics
become far more reliable. The non-monotonicity between lambda_import=1
and 4 (median lower at 4 than at 1) is very likely Monte Carlo noise from
the smoke-scale ensemble (n=50 per cell against a right-skewed,
multimodal outcome distribution) rather than a genuine effect -- **this
comparison has NOT been empirically resolved at adequate precision** and
should not be read as evidence that higher importation pressure reduces
cumulative infections.

### 5.4 Demography (D1 vs. D0)

Median 10-year cumulative infections: D1 (IBGE-rate-based) 3,222,289 vs.
D0 (stationary 2025 rates) 3,093,276 -- within ~4% of each other,
consistent with the two arms' population trajectories differing by only
~1% by 2035. The generator's conclusions are **not materially sensitive**
to this particular choice at the scale of the other structural
sensitivities above.

## 6. Failure modes and unresolved issues

1. **Y1's iid-per-year resampling is itself an untested structural
   choice.** Treating each future year's deviation as an independent draw
   from the same distribution the historical fit used is a defensible
   "uncertainty-driven sensitivity" construction, but there is no
   empirical check here that real year-to-year chikungunya transmission
   variability in Ceara is well-described by 11 independent historical
   draws extrapolated this way (vs., e.g., some persistence/autocorrelation,
   or a trend).
2. **Outbreak/wave detector over-counts minor fluctuations** relative to
   the full national wave-census algorithm it was adapted from (no
   `trough_persistence` requirement was carried over in the simplified
   version). Use `major_epidemic==TRUE` rows for any recurrence-interval
   claim, not the raw wave count.
3. **The lambda_import=1-vs-4 comparison, and most of the climate-arm
   comparison, were run only at smoke scale** (10 draws x 5 replicates
   per cell) and are not precisely resolved -- only the single reference
   scenario was run at the recommended pilot scale (100x20).
4. **Lifelong immunity** is carried over unexamined; this is increasingly
   load-bearing the further the horizon extends.
5. **No age structure or vaccination** in this pilot (by design/instruction).
6. No out-of-sample or external validation of the SIMULATED future
   trajectories is possible or attempted -- by construction, nothing
   about 2026-2035 can be checked against reality yet.

## 7. Requirements before SSP projections or vaccination analyses

- Resolve the lag-spec and multi-state climate-coefficient sensitivity
  gaps identified in Step 1 before trusting beta_T/beta_P as a portable
  climate-transmission relationship to force with SSP data.
- Re-run the full sensitivity grid (not just the reference scenario) at
  pilot scale (100x20) or larger to resolve the importation-pressure and
  climate-arm comparisons with adequate Monte Carlo precision.
- Decide, and empirically motivate if possible, a better-justified
  distributional model for future year-to-year transmission variability
  than i.i.d. resampling of the historical hierarchical prior, given how
  dominant this one choice is (Section 5.1).
- Strengthen the wave-detection algorithm (carry over `trough_persistence`
  and the other national-census refinements) before using recurrence-
  interval outputs for any downstream burden or NNV-style calculation.
- Only then introduce actual SSP/CMIP6 climate forcing (replacing arms
  C0-C2) and, separately, layer the existing validated age-cohort/
  vaccination architecture on top of this aggregate engine.

## 8. Conclusion

**Classification: B -- partially validated, requiring the specified
improvements above.**

The simulation ENGINE itself is rigorously validated: exact historical
replay, exact forcing identity, exact mass conservation, correct
zero-introduction behaviour, correct joint posterior propagation, and
full reproducibility all pass without qualification (Section 3). This is
strong evidence of **code correctness and mechanistic faithfulness** to
the accepted historical model. What is NOT yet validated is the
SCENARIO layer built on top of that engine: the magnitude of future
year-to-year variability (by far the dominant driver of outcomes),
precise quantification of the importation and some climate-arm
sensitivities (currently smoke-scale only), and the simplified
outbreak-detection rule. None of the Step 2 outputs -- including FIG05-09
-- should be read as validated forecasts or policy-ready projections;
they are a scientifically auditable demonstration that the generator
produces epidemiologically coherent, internally consistent, mechanism-
faithful trajectories under clearly labelled structural assumptions, and
a map of exactly which of those assumptions matter most for the eventual
SSP/vaccination work.
