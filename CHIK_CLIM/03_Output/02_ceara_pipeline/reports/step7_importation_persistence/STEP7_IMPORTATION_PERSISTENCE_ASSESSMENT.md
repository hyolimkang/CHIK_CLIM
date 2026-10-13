# STEP 7 -- Importation / Persistence Sensitivity Assessment

Controlled mechanistic sensitivity analysis, NOT model calibration. No
Stan refitting; reuses the H1/H2 training-only posteriors unmodified.
No lambda_import value is selected as a new reference. **PROVISIONAL**:
H1/H2 training-only posteriors carry low but nonzero HMC divergences
and incomplete within-chain mixing (tau_year upper tail,
`HMC_RELIABILITY_AUDIT.md`) -- every probability/interval below is
conditional on this. Not committed. Scripts:
`02_Script/02_ceara_pipeline/16_step7_importation_persistence/`.

## Phase 0 -- Baseline reproduction and audit

Reproduced the lambda_import=1 baseline for H1/G2 and H2/G2 using the
SAME training-only posteriors (`PHASE0_posterior_identity_check.csv`:
N_train=261/365, 1600 post-warmup draws each, matching Step 6 exactly).

| Metric | H1 new | H1 Step6 | rel.diff | H2 new | H2 Step6 | rel.diff |
|---|---|---|---|---|---|---|
| Entry-2022 weighted infectiousness (median) | 4.879 | 4.77 | 2.3% (agree) | 871.3 | 886.4 | 1.7% (agree) |
| S/N, Jan 2022 (median) | 0.7417 | 0.741 | 0.09% (agree) | 0.7439 | 0.744 | 0.01% (agree) |
| Peak Reff, 2022 (median) | 1.327 | 1.41 | 5.9% | 1.346 | 1.40 | 3.9% |
| P(exceed 2022 observed) | 0.0% | 0.3% | n/a (floor) | 14.35% | 11.9% | within MC range |

**Cause of the (modest) disagreement, investigated explicitly rather
than assumed**: this module's simulation engine draws Poisson
introductions via a shared-uniform inverse-CDF method,
`qpois(runif(1), lambda)`, instead of Step 6's direct `rpois(1, lambda)`
-- REQUIRED to enable the common-random-number lambda comparisons in
Phases 1-3 (documented in `01_simulation_engine.R`). This consumes a
different number of random draws per week than `rpois()`, so even at
an IDENTICAL lambda and an identical starting seed, the two engines'
per-week stochastic infection draws diverge from week 1 onward --
functionally equivalent to drawing an independent Monte Carlo
realisation, on top of an independently sampled 100-of-1600 posterior
draw subset (different seed stream than Step 6's `BASE_SEED_YV+777`).
A cluster-bootstrap check on this run's own data
(`PHASE0_addendum_prob_major_discrepancy_check.csv`) found P(major
epidemic, H1/G2, lambda=1) = 19.9% [12.7%, 27.75%] -- Step 6's 10.55%
falls OUTSIDE this interval, and the per-draw SD of P(major) is large
(0.378), confirming substantial between-posterior-draw heterogeneity
in this statistic specifically. **This is disclosed, not hidden**: the
baseline was NOT silently replaced, and every Phase 1-6 comparison
below uses THIS module's own, internally-consistent, shared-uniform-
coupled baseline throughout (never mixing in Step 6's numbers for a
quantitative contrast).

**Established/major-epidemic definition audit**
(`PHASE0_established_major_crosstab_H2.csv`): cross-tabulating the two
flags for H2/G2 2022 found 43/2000 trajectories (2.15%) are
`major_epidemic=TRUE` but `established=FALSE`. Cause: `established`
requires `M_import>0` coincident with `R_eff>1` WITHIN the evaluation
year itself; a trajectory already carrying substantial infections from
EARLIER periods can grow into a major epidemic via pure local
transmission with zero NEW imports in that specific year, and such
cases are never flagged `established` by this rule. This is a genuine
definitional precision issue (not a bug) that recurs throughout this
analysis -- see Phase 2's Experiment B below.

**Known limitations carried forward** (quantified in Step 6, not
re-derived or corrected here): future year-effect recentering (shrinks
marginal SD by sqrt((K-1)/K), K=holdout years) and climate-anomaly
standardisation using full-2015-2025 scaling constants.

## Phase 1 -- Importation intensity grid (FIG24)

Four prespecified levels, common random numbers within each holdout
(same seed, same sampled posterior draws, same z_year/z_season/climate
per draw -- only the shared-uniform-coupled M_t differs):

| Arm | H1 entry infect. | H1 P(established) | H1 P(major) | H1 P(exceed obs.) | H1 median 2022 cases | H2 entry infect. | H2 P(major) | H2 P(exceed obs.) |
|---|---|---|---|---|---|---|---|---|
| I0 (λ=0) | 0.000 | 0.0% | 12.4% | 0.0% | 0 | 871.3* | 97.15% | 14.55% |
| I1 (λ=0.25) | 1.688 | 79.25% | 14.95% | 0.0% | 12 | 871.3* | 97.20% | 14.35% |
| I2 (λ=1, baseline) | 4.879 | 89.9% | 19.9% | 0.0% | 40 | 871.3* | 97.15% | 14.35% |
| I3 (λ=4) | 17.438 | 89.8% | 31.1% | 1.55% | 138 | 871.3* | 97.35% | 14.65% |

\* H2's entry infectiousness is determined ENTIRELY by the training
posterior (fitted to real 2015-2021 surveillance) and is, by
construction, invariant to the HOLDOUT-period lambda -- confirmed
numerically identical to 1 decimal place across all four H2 rows.

**Findings**: H1's entry infectiousness and P(major) rise monotonically
with lambda, but even at 4x baseline H1 remains ~50x below H2's entry
level and its median 2022 case count (138) remains 369x below the
observed 50,986. P(establishment) saturates almost immediately
(79%→90%→90% from λ=0.25 to λ=4) -- establishing SOME transmission is
not the limiting step at any realistic lambda; reaching MAJOR-epidemic
scale is. H2's outcomes are essentially flat across the entire grid,
confirming its advantage is structural (starting state), not an
artefact of holdout-period importation intensity.

## Phase 2 -- Pre-2022 vs. 2022-onward importation (FIG25)

H1/G2 only. Experiment A varies 2020-2021 lambda (2022+ fixed at 1);
Experiment B varies 2022+ lambda (2020-2021 fixed at 1).
(`PHASE2_pre_post_2022_split_summary.csv`)

| Varied lambda | Exp A: P(major), pre-2022 varied | Exp B: P(major), 2022+ varied |
|---|---|---|
| 0 | 14.2% | 18.65% |
| 0.25 | 16.45% | 18.80% |
| 1 | 19.9% (=baseline) | 19.9% (=baseline) |
| 4 | **29.0%** | **24.10%** |

**Range across the grid: Experiment A spans 14.8 percentage points,
Experiment B spans 5.45 -- pre-2022 pressure moves P(major epidemic)
~2.7x more than an equally-sized change to 2022-onward pressure.**
Experiment B's P(established) jumps from 0% (post=0) to 85.1%
(post=0.25) purely because of the definitional artefact above (zero
NEW imports in 2022 makes `established` tautologically FALSE) -- but
P(major epidemic) at post=0 (18.65%) is barely different from the
post=1 baseline (19.9%), confirming that whether or not FRESH
importation occurs during 2022 itself matters far less than how much
standing infection already existed beforehand.

## Phase 3 -- Introduction timing at fixed total pressure (FIG26)

H1/G2, lambda=1. T1 redistributes each (draw, replicate)'s REALISED
weekly M_t within each calendar year (yearly totals preserved exactly:
max abs diff = 0 across 12,000 draw x replicate x year cells,
`PHASE3_timing_permutation_total_preservation_check.csv`), replayed
with the same RNG-stream alignment as T0.

| | T0 (original timing) | T1 (permuted timing) |
|---|---|---|
| Entry-2022 weighted infectiousness (median) | 4.879 | 5.087 |
| P(established) | 89.9% | 89.65% |
| P(major epidemic) | 19.9% | 19.85% |
| P(exceed observed 2022) | 0.0% | 0.05% |

Paired per-(draw,replicate) delta in 2022 reported cases: median=0
[-1,083.6, 803.0] (n=2,000) -- centred exactly on zero.
**Timing, independent of total magnitude, has essentially NO
detectable effect on any 2022 outcome.** This explicitly does not
imply real introductions are uniformly distributed across the year --
it means that, WITHIN this model's own mechanics, redistributing a
fixed total amount of pressure across weeks of a year does not change
the downstream epidemic outcome in any way this analysis can detect.

## Phase 4 -- Epidemic amplification (descriptive)

H1/G2, lambda grid I0/I2/I3 (`PHASE4_growth_dynamics_summary.csv`).
Among trajectories that DO cross the major-epidemic threshold, growth
is fast once underway regardless of lambda (median week-to-100-latent-
infections = 1-2.5 weeks from the start of 2022); among those that
don't, most either never reach 100 cumulative latent infections
(79% at λ=0, dropping to 1.8% at λ=4) or go fully extinct after their
one favourable window (66.2% at λ=0, **0% at λ=1 and λ=4** -- with any
ongoing background importation, trajectories essentially never go
fully extinct once a favourable period begins, yet most STILL fail to
reach major-epidemic scale). The ratio of local to introduction hazard
during favourable weeks is large (10^2-10^6, depending on lambda) for
BOTH major and non-major trajectories once transmission is under way --
**local transmission, once triggered, overwhelmingly dominates the
force of infection regardless of lambda; lambda's role is in whether
and how reliably a chain gets triggered and sustained early on, not in
how strongly it amplifies once going.** This is consistent with Phase
1-3's finding that neither susceptibility nor peak Reff (which are
essentially identical to H2's) is the limiting factor -- the bottleneck
is specifically in reaching a large enough STANDING infection base
before the 2022 season's favourable window closes.

## Phase 5 -- Generator-definition sensitivity (recentering, OPTIONAL arm)

H1/G2, lambda=1, single isolated change (never combined with an
importation change in the same arm). Pooled SD of year_effect across
ALL draws x holdout years: centred (accepted) = 0.180, uncentred
(diagnostic) = 0.194 -- ratio 0.929, closely matching the theoretical
sqrt(5/6)=0.913 for H1's K=6 holdout years, directly confirming Step
6's recentering finding via an actual counterfactual simulation (not
just a formula). Downstream effect: entry infectiousness 4.879→4.382,
P(established) 89.9%→89.25%, P(major) 19.9%→18.4%, **P(exceed
observed) 0.0%→1.2%** -- the uncentred arm's heavier, unconstrained
tail modestly RAISES the probability of an extreme outcome even though
its central tendency is similar-to-lower. A real, disclosed, but
second-order effect relative to the lambda/pre-post findings above
(≤1.2 percentage points vs. Phase 1's up to 1.55-point and Phase 2's
up to 14.8-point swings).

## Phase 6 -- Predictive performance and uncertainty (FIG27)

Full ensembles (100 draws x 20 replicates) throughout; the engine's
correctness was validated at Phase 0 before any full-scale arm was
launched. `PHASE6_annual_predictive_performance_grid.csv`,
`PHASE6_cumulative_predictive_performance_grid.csv`,
`PHASE6_exceedance_probability_mc_grid.csv` (posterior-draw-cluster
bootstrap, separating parameter from stochastic uncertainty).

- Multi-year cumulative H1 WIS improves monotonically with lambda
  (34,548→31,315→26,338→22,201) -- but the predictive MEDIAN remains
  >=10x below observed (628→1,379→2,353→5,404 vs. 56,903) at every
  level, and 95% coverage of the observed 2022 value is **never**
  achieved at ANY tested lambda, including λ=4.
- By year: raising lambda is NOT uniformly beneficial across the
  holdout. 2022 WIS improves steadily with lambda; 2023-2025 WIS is
  actually BEST at the CURRENT baseline (λ=1) and modestly worsens at
  λ=4 -- a real, if small, trade-off, not a free improvement.
  2020-2021 (pre-resurgence years, low observed counts) are
  essentially insensitive to lambda.
- H2's exceedance probability (14.35-14.65%) and cumulative WIS
  (15,237-15,484) are flat across the ENTIRE grid, with substantial
  between-draw parameter uncertainty (SD~0.34-0.35) dwarfing the
  lambda-driven variation -- confirming H2's performance is governed by
  its starting state, not holdout-period importation assumptions.
- H1's exceedance-probability parameter uncertainty is near-zero at
  low lambda (a floor effect: almost no sampled draw ever exceeds the
  threshold) and only becomes meaningfully estimable at λ=4 (SD=0.10).

## Answers to the ten questions

**1. Is the 2020-2021 near-fadeout behaviour sensitive to lambda_import?**
Yes, clearly: entry infectiousness rises monotonically and substantially
(0→1.7→4.9→17.4) across the grid (Phase 1). But "sensitive" is not
"sufficient" -- see Q2.

**2. Does increased pre-2022 introduction pressure materially increase
standing infectiousness at the start of 2022?**
Yes, materially (Phase 1-2), but not SUFFICIENTLY: even at 4x baseline,
H1's entry infectiousness (17.4) remains ~50x below H2's real-data-
informed level (871.3).

**3. Is the 2022 epidemic magnitude more sensitive to introductions
before 2022 or during 2022?**
Before 2022, by a factor of ~2.7x in P(major epidemic) range across
the same prespecified grid (Phase 2) -- a clear, controlled,
quantified result.

**4. Does introduction timing matter when total introduction pressure
is unchanged?**
No. Phase 3 found no detectable effect of within-year timing on any
outcome (paired delta centred exactly at zero).

**5. Does the existing constant-rate Poisson structure generate
epidemiologically plausible infection persistence?**
Partially: it does not go fully extinct once established under
realistic lambda (0% extinction after a favourable window begins, at
λ>=1), but the MAGNITUDE of persistence it generates even at 4x
baseline (median entry infectiousness 17.4) remains an order of
magnitude below what the real 2020-2021 epidemic (as captured by H2's
training posterior) evidently sustained (871.3).

**6. Can introduction pressure explain H1's poor 2022 retrospective
predictive performance?**
Partially. It explains a meaningful FRACTION of the gap (WIS improves
~36% from λ=0 to λ=4) and the controlled pre/post-2022 contrast
cleanly identifies WHERE the sensitivity is concentrated (before 2022).
It does NOT explain the full gap: even the most favourable tested
level leaves the median 369x below observed and never achieves 95%
coverage of the observed 2022 value.

**7. Does an increase in introductions improve the full holdout
distribution, or does it create unrealistically high incidence in
other years?**
Mostly improves, with one disclosed exception: WIS improves or is
flat for 2020-2021 and 2022 at every tested lambda, but 2023-2025 WIS
is best at the CURRENT baseline (λ=1) and modestly worse at λ=4 --
raising lambda for the whole holdout is not a free win across all years.

**8. Are other mechanisms, such as seasonal transmission amplification,
more important?**
No evidence of that here. Susceptibility (Phase 2 of Step 6) and peak
Reff (this module's Phase 0-1 headline metrics) are both essentially
identical between H1 and H2, and Phase 4's hazard decomposition shows
local transmission dominates by orders of magnitude once a chain is
established, regardless of lambda -- amplification capacity is not the
bottleneck; getting enough STANDING infection to amplify from is.

**9. Is a structural revision scientifically justified?**
A revision to the STRUCTURE of the introduction/persistence mechanism
(not merely its rate) is justified by the evidence above: pre-2022
pressure matters most, timing within a constant-rate process does not
matter at all, and even a 4x rate increase leaves a ~50x (entry
infectiousness) to 369x (median 2022 cases) gap. A constant, smoothly-
distributed Poisson rate -- at ANY magnitude plausible as a "background"
rate -- appears structurally unable to reproduce the scale of standing
infection H2's real-data-informed posterior implies for late 2021.
**Proposed minimum structural change** (NOT implemented here): replace
the single constant `lambda_import` with an episodic/hurdle
reintroduction process -- most weeks `M_t ~ Poisson(lambda_background)`
(small, as now) but with a per-period probability `p_event` of a
larger reintroduction episode, `M_t ~ Poisson(lambda_event)` with
`lambda_event >> lambda_background` -- epidemiologically motivated by
the documented 2022 ECSA-lineage introduction into Ceará (consistent
with Task D/E's independent finding that a single new lineage, not a
gradually-accumulating background rate, likely drove the 2022
resurgence). This needs >=2 new parameters (`p_event`, `lambda_event`),
almost certainly requires JOINT historical refitting (not a post-hoc
simulation tweak) since tau_year's current posterior may already be
implicitly absorbing some of what is actually introduction-driven
variability, and introduces a DISCRETE latent state per period that
Stan can only handle via marginalisation -- a materially more complex
HMC geometry than the current model, likely to interact unfavourably
with the ALREADY-disclosed tau_year mixing difficulties
(`HMC_RELIABILITY_AUDIT.md`) and warranting its own dedicated
development and HMC-reliability validation cycle before any
scientific use.

**10. What external information would be necessary to estimate
introduction intensity rather than treating it as an assumption?**
Genomic/phylogenetic surveillance data capable of dating and counting
distinct lineage-introduction events into Ceará (to directly inform
`p_event`/`lambda_event` rather than inferring them purely from case
counts, which this analysis's own identifiability concerns -- shared
with the tau_year HMC caveats -- suggest may not cleanly separate
"many small introductions" from "rare large ones" using surveillance
data alone); alternatively, regional case data from neighbouring
states or travel/mobility data as a time-varying covariate proxy for
introduction risk.

## Classification

**B -- Introduction timing/persistence is important and a time-varying
or episodic process should be investigated.**

Not selected merely because one lambda reproduces 2022 -- NO tested
lambda (up to 4x baseline) comes remotely close to reproducing the
observed magnitude (369x gap persists even at λ=4), so this is
explicitly NOT "A" with an asterisk. The qualifier on B, directly
supported by Phase 3's null timing result, is important: the needed
revision is not "introductions need better within-year timing" (ruled
out) but "the STRUCTURAL FORM of the introduction process -- constant
small-rate vs. episodic/lumpy -- appears inadequate," since the
PERSISTENCE/magnitude dimension of B (Phases 1-2, 4) is strongly and
specifically supported while the TIMING dimension (Phase 3) is not.
Not "C": susceptibility and transmission/amplification capacity were
directly tested and found essentially equivalent between H1 and H2,
and local transmission dominates by orders of magnitude once
established regardless of lambda -- there is no positive evidence
pointing at transmission, immunity, or observation structure as the
limiting mechanism. Not "D": the evidence discriminates cleanly among
the candidate mechanisms tested here (pre-2022 > 2022-onward; timing
null; transmission/susceptibility not implicated) -- the remaining
uncertainty is about WHICH introduction-process structure and
parameters to use, not about which broad mechanism is responsible.

## Files, scripts, and reused fits

```
02_Script/02_ceara_pipeline/16_step7_importation_persistence/
  00_config_step7.R, 01_simulation_engine.R,
  02_phase0_baseline_audit.R, 02b_phase0_addendum_prob_major_check.R,
  03_phase1_importation_grid.R, 04_phase2_pre_post_2022_split.R,
  05_phase3_timing_permutation.R, 06_phase5_recentering_sensitivity.R,
  07_phase4_hazard_decomposition.R, 08_phase6_predictive_performance.R,
  09_make_figures.R

Reused (NOT refit, NOT modified):
  model_fits/future_year_effect_validation/ce_tau_year_training_{H1,H2}.rds

New outputs (isolated):
  03_Output/02_ceara_pipeline/{tables,results}/step7_importation_persistence/*.csv,*.rds
  03_Output/02_ceara_pipeline/figures/step7_importation_persistence/FIG24-27.{pdf,png@600dpi}
  03_Output/02_ceara_pipeline/reports/step7_importation_persistence/STEP7_IMPORTATION_PERSISTENCE_ASSESSMENT.md (this file)
```

**No new Stan fitting was performed or is required for the findings
above.** If the proposed episodic-reintroduction structural revision
(Q9) is pursued, it WOULD require new Stan fitting (joint, not a
simulation-only change) -- proposed, not launched, per the task's
explicit instruction.

## Summary for the user

**1. Analyses completed**: Phase 0 baseline reproduction + definitional
audit; Phase 1 four-level importation grid (H1 primary, H2 secondary,
full 100x20 ensembles); Phase 2 pre-/post-2022 controlled split
(2 experiments x 4 levels); Phase 3 timing-permutation diagnostic with
an explicit total-preservation check; Phase 4 hazard decomposition and
growth-dynamics analysis (descriptive); Phase 5 isolated recentering
counterfactual; Phase 6 full predictive-performance grid with
cluster-bootstrap uncertainty; four figures (FIG24-27).

**2. Key numerical findings**: pre-2022 importation pressure moves
P(major epidemic) ~2.7x more than equally-sized 2022-onward changes;
within-year timing has no detectable effect; even 4x baseline leaves
H1's entry infectiousness ~50x below H2's real-data-informed level and
median 2022 cases 369x below observed; local transmission dominates
the hazard by orders of magnitude once established, at every lambda
tested; raising lambda improves 2022 and cumulative WIS but mildly
worsens 2023-2025 fit.

**3. Is constant importation adequate?** No -- not even at 4x the
current baseline rate. It is directionally correct (more pressure
before 2022 helps) but insufficient in magnitude by roughly two orders
of magnitude on the metric that matters most (entry infectiousness).

**4. Is a structural change justified?** Yes, on the introduction-
process STRUCTURE (episodic/time-varying vs. constant-rate), not on
transmission, susceptibility, or observation structure, which are not
implicated by this analysis. Proposed minimum change specified above
(Q9); not implemented.

**5. Remaining identifiability/HMC limitations**: every quantitative
result here is provisional on the tau_year HMC reliability caveats
already documented. Separately, this module's own baseline differs
non-trivially from Step 6's point estimates for probability-type
metrics (P(major), P(exceed)) due to RNG-method and posterior-draw-
subsampling differences -- disclosed, investigated, and shown not to
reflect an implementation bug, but a reminder that these particular
statistics carry substantial Monte Carlo uncertainty that should
always be reported with a cluster-bootstrap interval (as done
throughout Phase 6), never as a bare point estimate.

**6. One recommended next action**: before any further importation-
structure development, resolve the tau_year HMC reliability caveats
(the proposed episodic-reintroduction model would add a harder,
discrete-mixture geometry on top of an already-imperfect hierarchical
fit) -- then, if pursued, design the episodic-reintroduction model as
a dedicated, separately-validated development effort (its own HMC
audit cycle), explicitly informed by any available genomic/lineage
surveillance data on the 2022 introduction (Q10), rather than as an
extension bolted onto the current tau_year/constant-Poisson structure.

STOP after Step 6/7's scope as instructed -- Step 8 / SSP work not started.
