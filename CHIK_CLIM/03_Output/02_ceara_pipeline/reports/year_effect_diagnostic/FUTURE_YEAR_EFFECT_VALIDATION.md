# Future year-effect validation diagnostic (Ceara, follow-up to the future_projection_pilot)

Produced by `02_Script/02_ceara_pipeline/11_year_effect_diagnostic/`.
**ONE targeted diagnostic step.** Re-analyses the already-fitted accepted
climate-forced v4.9 canary posterior (no refit) and reuses, unmodified,
the Step 2 simulation functions (`03_future_renewal_functions.R`). Does
**not** proceed to SSP projection or vaccine simulation.

## 0. The bug being fixed

`05_run_simulation_ensemble.R`'s `run_scenario()` seeds both posterior-draw
selection (`sub_seed("drawselect", scenario_id)`) and future year-effect
innovations (`sub_seed("yeareffect", scenario_id, d)`) by `scenario_id`.
Since the Y0 and Y1 sensitivity-grid scenarios have different
`scenario_id` strings, **comparing Y0 vs. Y1 in the original sensitivity
grid used different posterior draws, different initial conditions
(because the draw differs), and different stochastic innovations** --
every reported contrast was confounded with ordinary Monte Carlo sampling
noise on top of the structural difference. This diagnostic uses a PAIRED
design instead: one common draw-selection seed, one common future
z_year/z_season innovation per draw (shared across all 4 arms), one
common climate input (C0, zero), one common demography, and one common
starting RNG seed per (draw, replicate) reset before each arm's
simulation (common random numbers, exact for the first draws consumed
each week; it decoheres after trajectories diverge enough to change
`rbinom()`'s internal algorithm path -- disclosed, not hidden, and does
not affect the validity of the paired summary statistics below, which
use the realised outcomes directly).

## 1. Historical posterior vs. future Y1 sampling distribution (Section 2)

**year_effect** (`year_effect = 0.40*(z_year - mean(z_year))`, `year_effect_prior_sd = 0.40` is a FIXED data value, NOT an estimated parameter):

| | median | 95% range (across draws/realisations) |
|---|---|---|
| Historical within-draw SD (genuine interannual variability implied by the fit) | **0.138** | [0.120, 0.157] |
| Historical within-draw range (max-min across the 11 years) | 0.439 | [0.374, 0.511] |
| Prior-implied SD if z_year were freely iid N(0,1) (finite-sample correction for Y=11) | 0.381 | -- |
| **Future Y1 within-realisation SD (10 fresh years)** | **0.385** | [0.218, 0.583] |
| Future Y1 within-realisation range | 1.205 | [0.665, 1.924] |
| **Ratio, future/historical (SD)** | **2.79x** | |
| Ratio, future/historical (range) | 2.74x | |

The historical FITTED between-year spread (0.138) is only **36% of what
the fixed prior scale (0.40) implies for freely-sampled years** (0.381) --
the likelihood pulls the 11 actual years much closer together than the
hyperparameter alone would allow. Because Y1 resamples fresh
`z_year ~ Normal(0,1)` for each future year using that SAME fixed 0.40
scale (not the fitted, data-informed spread), **future year_effect
sampling injects about 2.8x more year-to-year spread than the historical
posterior actually exhibits.**

**A_year** (`A_year = exp(sigma_season_year * z_season_year)`, `sigma_season_year` **IS** a genuine posterior-estimated parameter, posterior median 0.476, 95% CrI [0.342, 0.651] against a Normal(0,0.20)-truncated prior -- the data pushed this parameter well above its prior scale, i.e. real season-to-season amplitude variation is empirically supported):

| | median | 95% range |
|---|---|---|
| Historical within-draw SD | 0.614 | [0.449, 0.841] |
| Historical within-draw range | 2.082 | [1.509, 2.926] |
| **Future Y1 within-realisation SD** (using the SAME draw's posterior `sigma_season_year`) | **0.491** | [0.230, 1.170] |
| Future Y1 within-realisation range | 1.535 | [0.700, 3.857] |
| **Ratio, future/historical (SD)** | **0.80x** | |
| Ratio, future/historical (range) | 0.74x | |

Future A_year sampling is, if anything, **slightly more conservative**
than the historical fitted spread -- because it correctly propagates the
posterior-estimated `sigma_season_year` per draw, this mechanism IS
empirically calibrated.

**Frequency of extreme future peak R0** (2,000 realisations, C0
climate-neutral, using real posterior draws' alpha_R/seasonal
coefficients): the historical fitted annual peak R0 across all 1,600
draws x 11 years ranges up to **3.80** (99th pct 3.33, 95th pct 3.12,
median 1.92). Future Y1-sampled annual peaks exceed this historical
**maximum** in **49.9% of realisations**, exceed the historical 99th
percentile in 72.1%, and exceed the historical median in 100%. **Future
R0 extremes routinely exceed anything the historical fit ever produced.**
This is an extrapolation beyond the fitted historical range, not merely
"more uncertainty" -- flagged explicitly, not treated as a hard
biological impossibility (the historical record is only 11 years and
could itself under-represent rare extremes), but as **unsupported by
the fitted model as currently specified.**

**Is Y1 relying on the prior scale rather than an empirically estimated
future variance?** For year_effect: **yes** -- `year_effect_prior_sd=0.40`
is a fixed hyperparameter never updated by the likelihood; nothing in the
Stan model estimates "how much years should differ," only which
direction each specific year deviates, regularised by that fixed scale.
For A_year: **no** -- `sigma_season_year` is a true free parameter with
its own prior, informed by the data, and correctly propagated per-draw
into the future sampling.

## 2. Paired 4-arm decomposition (Sections 3-4)

Full-scale paired run: 100 posterior draws x 20 stochastic replicates x 4
arms = 8,000 simulations (smoke test at 10x5 confirmed the same pattern
first). C0 climate-neutral, lambda_import=1, D1 (IBGE-rate) demography --
all four arms held identical except the year_effect/A_year construction.

| Arm | median 10y cum. infections | S/N min (median) | major-epidemic proportion |
|---|---|---|---|
| A: neutral (year_effect=0, A_year=1) | 7,676 | 0.636 | 0.0% |
| **B: year-level only** (year_effect=Y1, A_year=1) | **2,359,031** | 0.407 | 99.4% |
| C: amplitude only (year_effect=0, A_year=Y1) | 9,603 | 0.636 | 19.0% |
| D: both (= the original Y1 arm) | 2,616,051 | 0.373 | 99.7% |

(Arm A's ~7,676 and arm D's ~2,616,051 closely reproduce the original
non-paired sensitivity grid's ~8,100 and ~3.26M -- confirming the
original finding's existence was real, while this paired design removes
the sampling-noise confound from the explanation.)

**Paired contrasts** (same draw x replicate across arms, n=2,000 pairs):

| Contrast | median difference | 95% interval | P(difference > 0) |
|---|---|---|---|
| B - A (year-level effect) | **+2,351,463** | [24,856, 5,295,314] | **100%** |
| C - A (amplitude effect) | +1,926 | [-2,315, 116,828] | 75.5% |
| D - A (both, total) | +2,607,990 | [34,770, 5,757,273] | 100% |
| Interaction (D-B-C+A) | +41,857 | [-724,665, 2,969,886] | 58.6% (~null) |

On `log1p(cum_10y)` (more appropriate given the extreme right-skew and
population bound): B-A = 5.70, C-A = 0.23, D-A = 5.78, interaction =
-0.20 (small, slightly negative -- sub-additive, consistent with B alone
already driving most realisations to near-saturating depletion, leaving
little additional room for C to add).

**The year-level effect (B) alone accounts for ~90% of the total (D-A)
effect on the raw cumulative-infection scale, and ~99% on the log
scale.** The amplitude-only effect (C) is small and not even reliably
positive (CI crosses zero; only 75.5% of paired draws show C>A). The
interaction is small and statistically indistinguishable from zero.
**This fully persists under the paired design** -- it is not an artifact
of the original ensemble's confounded sampling.

## 3. Figures

- **FIG10** (`FIG10_historical_vs_future_distributions.{pdf,png,svg}`): historical posterior vs. future Y1 pooled densities for year_effect and A_year. Visually: year_effect's future density is far flatter/wider than the historical (multi-modal, tight) density; A_year's two densities nearly overlap.
- **FIG11** (`FIG11_annual_peak_R0_by_arm.*`): annual peak R0 violin+box by arm plus historical reference. Historical, A and C are visually similar and tight; B and D show long upper tails reaching R0~10.
- **FIG12** (`FIG12_paired_cumulative_infections.*`): (A) log-scale cumulative-infection distributions by arm; (B) forest plot of the paired contrasts above with 95% intervals.

## 4. Answers to the five questions

**1. Is future Y1 variability larger than the historical posterior supports?**
For the year-level mechanism: **yes, substantially** (~2.8x the fitted
historical within-draw spread, and future peak R0 exceeds the historical
maximum in ~50% of realisations). For the seasonal-amplitude mechanism:
**no** -- it is well-calibrated to the posterior (ratio 0.74-0.80x,
slightly conservative).

**2. Which of the two year-varying components drives the extreme R0 / outcomes?**
**The year-level effect (`year_effect`, mechanism B), overwhelmingly.**
It alone reproduces ~90-99% of the full Y1 (arm D) effect; the
seasonal-amplitude mechanism (C) contributes a small, only marginally
positive effect; the interaction between them is small and not
distinguishable from zero.

**3. Does the large infection difference persist under paired simulations?**
**Yes.** Fixing the posterior-draw/innovation confound with a paired,
common-random-number design, the B-A and D-A contrasts remain large,
strictly positive in 100% of 2,000 paired comparisons, and of the same
order of magnitude as the original (confounded) sensitivity-grid finding.

**4. Is Y1 justified as a central future scenario, or only as a structural sensitivity scenario?**
**Only as a structural sensitivity scenario**, and specifically because
of its year-level component. `year_effect_prior_sd=0.40` is a
researcher-chosen regularisation constant, not an estimated
between-year-variance parameter; resampling fresh years at that fixed
scale demonstrably overstates the variability the fitted historical
process itself exhibits. Y1 remains legitimate and useful as an
upper-bound / stress-test structural scenario (exactly as the Step 2
report already labelled it), but should NOT be read as a best-estimate
or central future expectation without first addressing point 5. The
amplitude-only mechanism (C), by contrast, is reasonably well-calibrated
and could be used closer to a central assumption on its own.

**5. What remains necessary to calibrate future interannual variability?**
- Replace the fixed `year_effect_prior_sd=0.40` with an actually
  ESTIMATED between-year variance hyperparameter (a hierarchical
  structure analogous to how `sigma_season_year` already has its own
  prior and posterior), so future sampling can legitimately draw from a
  data-informed scale instead of an arbitrary fixed one. This requires a
  Stan model change and refit -- explicitly NOT done here, per
  instruction.
- Alternatively, as a non-refitting interim fix, calibrate the future
  year-effect sampling scale directly to the empirical within-draw
  spread measured here (median SD 0.138, i.e. an effective scale of
  about 0.138/sqrt((11-1)/11) ≈ 0.145 instead of 0.40) and document this
  explicitly as a post-hoc calibration, not a refit.
- With only 11 historical years, the fitted between-year spread is
  itself estimated from a small sample and could be an underestimate of
  the true long-run interannual variability -- this diagnostic reduces,
  but does not eliminate, uncertainty about the right future scale.
- The extreme-R0 extrapolation finding (Section 1) means any future
  scenario relying on the current Y1 construction should be explicitly
  labelled as testing values beyond the historically fitted range, not
  as a validated plausible scenario.

## 5. Files changed / created

No existing file modified. New, isolated module only:

```
02_Script/02_ceara_pipeline/11_year_effect_diagnostic/
  00_config_year_effect_diagnostic.R
  00_run_year_effect_diagnostic.R
  01_historical_vs_future_variability.R
  02_paired_arm_simulation.R
  03_make_diagnostic_figures.R

03_Output/02_ceara_pipeline/
  reports/year_effect_diagnostic/FUTURE_YEAR_EFFECT_VALIDATION.md   (this file)
  tables/year_effect_diagnostic/DIAG_*.csv                          (11 tables)
  figures/year_effect_diagnostic/FIG10-12.{pdf,png@600dpi,svg}
  results/year_effect_diagnostic/DIAG_*.{csv,rds}                   (realisation-level draws, full trajectories)
```

## 6. Unresolved uncertainties

- The common-random-number pairing is exact only for the first random
  draws consumed each week; it decoheres once arms' (size, p) for the
  Binomial step diverge enough to change R's internal sampling
  algorithm. The paired summary statistics (computed on realised
  outcomes, not on raw RNG alignment) are unaffected by this, but it
  means the reported 95% intervals are not maximally variance-reduced.
- The historical-vs-future "pooled" comparison in FIG10 mixes posterior
  uncertainty and genuine interannual variability on the historical side
  (disclosed in the figure caption); the within-draw-only metrics in
  Section 1's tables are the methodologically cleaner comparison and
  should be preferred for any quantitative claim.
- 11 historical years is a short record for estimating any between-year
  variance parameter, fixed or hierarchical.
- This diagnostic used C0 (climate-neutral) throughout to isolate the
  year-effect mechanisms, per instruction; it does not re-examine how
  these findings interact with the climate arms (C1/C2) from Step 2.
