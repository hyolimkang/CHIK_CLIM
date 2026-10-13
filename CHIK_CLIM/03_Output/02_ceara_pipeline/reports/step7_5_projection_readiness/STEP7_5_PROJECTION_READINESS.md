# STEP 7.5 -- Structural Discrimination and Projection Readiness Assessment

Diagnostic, not calibration. No Stan refitting; reuses H1/H2 training-
only posteriors unmodified; no holdout cases used for calibration; no
SSP/ISIMIP or vaccination implemented. **PROVISIONAL**: every result
below is conditional on the tau_year HMC reliability caveats
(`HMC_RELIABILITY_AUDIT.md`: mild divergences concentrated at larger
tau_year, incomplete within-chain mixing, upper-tail width uncertain).
Scripts: `02_Script/02_ceara_pipeline/17_step7_5_projection_readiness/`.
Not committed.

## Phase 0 -- Reproducibility and implementation audit

Full 8-item audit: `PHASE0_audit_synthesis.csv`. Summary (citing the
exact prior source for each, not re-derived):

| # | Item | Status |
|---|---|---|
| 1 | Step6 vs Step7 baseline discrepancy | Resolved (disclosed -- caused by Step 7's shared-uniform Poisson coupling consuming RNG differently than Step 6's `rpois()`, confirmed via cluster-bootstrap, not a bug) |
| 2 | Posterior draw selection / replicate matching | PASS |
| 3 | Future year-effect recentering | Disclosed, quantified (SD ratio 0.929 vs. theoretical 0.913 for H1), not corrected in main arms |
| 4 | Historical climate-standardisation leakage | Disclosed, quantified (7.6-13.1% SD inflation), not corrected |
| 5 | Establishment/persistence/major-epidemic definitions | Audited; same-year-import tautology confirmed (43/2000 H2 cases major-but-not-established); THIS module uses a revised, year-agnostic chain-activity proxy in Phase 3 |
| 6 | Importation double-counting | PASS (R0 reconstruction to 4.4e-16; q reproduced to 4 s.f.) |
| 7 | Training-only latent history transfer | PASS |
| 8 | HMC reliability | Unresolved, explicitly carried forward |

**Smoke test** (episodic engine, 10x5, before any full-scale run):
realised mean M_t=1.04 (target 1.0), burst frequency=10.29% (target
p=0.10), S+U=N exactly, aggregate reported/latent ratio=0.04996
(target q=0.05) -- passed; full-scale runs proceeded.
(`PHASE0_episodic_engine_smoketest.csv`)

## Phase 1 -- Pre-resurgence infection persistence, 2020-2021 (FIG_S75_1)

H1/G2, lambda=1, full ensemble (100 draws x 20 replicates), reusing
Step 7's already-validated baseline. Annual coverage of OBSERVED
2020-2021 cases: `PHASE1_annual_2020_2021_coverage.csv` -- 2020
covered at 50%/95% (pred. median 484 vs. observed 1,013); 2021 covered
at 95% but NOT 50% (pred. median 99 vs. observed 1,318) -- the
predictive median is already notably low relative to observed by 2021,
one full year before the resurgence.

**Year-by-year persistence** (`PHASE1_annual_persistence_summary.csv`):
at the MEDIAN trajectory, zero weeks have zero latent infections in
EITHER 2020 or 2021 (out of ~52), and the median longest run of
locally-dominant transmission spans the full year -- the typical
trajectory sustains persistently LOW, not fully extinct, transmission.
P(near-total fadeout, >=45 zero weeks) = 0.0% in both years. Reff(t)
crosses 1 twice per year with a clear seasonal pattern (median 17
weeks/year with Reff>1) -- FIG_S75_1C shows transmission potential is
present and adequate.

**Entry-2022 comparison** (`PHASE1_entry2022_H1_vs_H2reference.csv`):
H1 (simulated) median = 4.88 [0.42, 2,070]; H2 (training-posterior,
diagnostic reference ONLY, never used to initialise H1) median = 871.3
-- a ~179x gap.

**A/B/C/D decomposition** (`PHASE1_ABCD_classification.csv`):
- **A (insufficient local transmission growth): NOT primarily** --
  Reff>1 occurs in most trajectories most years; local hazard
  dominates import hazard by orders of magnitude once underway.
- **B (repeated stochastic fadeout): NOT primarily, but not
  irrelevant to the tail** -- the median trajectory never fully fades
  out; fadeout is a real but secondary/tail phenomenon.
- **C (insufficient introduction pressure): PRIMARY** -- a ~179x gap
  between H1's simulated entry infectiousness and H2's real-data-
  informed reference, consistent with Step 7's independently-derived
  ~170-190x estimate.
- **D (state-propagation/observation errors): NOT SUPPORTED** -- all
  numerical-consistency checks pass to floating-point precision.

**The failure is not sudden at the start of 2022** -- it develops
gradually: 2020 is reasonably well covered, 2021's predictive median
is already substantially below observed, and by January 2022 the
simulated entry state is ~179x below what the real data (via H2)
imply it should be.

## Phase 2 -- Constant vs. episodic introduction, equal mean pressure (FIG_S75_2)

Model A (existing, constant Poisson(lambda=1), REUSED from Step 7's
baseline) vs. Model B (episodic: `B_t~Bernoulli(p)`,
`M_t|B_t=1~Poisson(lambda/p)`, `E[M_t]=1` exactly for p=0.10 and
p=0.25 -- prespecified illustrative values, NOT fit). Realised-process
check (`PHASE2_process_realisation_check.csv`): mean M_t = 1.000 /
1.008 / 1.002 for A / B(p=0.10) / B(p=0.25) respectively; burst
frequency 10.08%/24.55% for the two B arms (vs. A's 63.2% of weeks
with any M_t>0, as expected for Poisson(1)) -- the processes are
correctly realised at their intended means.

## Phase 3 -- Separating introduction from local amplification (FIG_S75_3)

Full structural comparison: `PHASE3_structural_comparison_summary.csv`
(10 requested outcomes); `PHASE3_WIS_by_year_structural_arms.csv`;
`PHASE3_other_year_peaks_artificial_outbreak_check.csv`.

| | Model A (constant) | Model B (p=0.10) | Model B (p=0.25) |
|---|---|---|---|
| Median entry-2022 infectiousness | 4.879 | 5.550 | 5.220 |
| P(chain active 2020-2021, revised year-agnostic definition) | 100% | 100% | 100% |
| P(major epidemic, 2022) | 19.90% | 19.55% | 19.50% |
| P(reaches observed 2022 burden) | 0.00% | 0.15% | 0.05% |
| Median cumulative reported, full period | 2,353 | 2,242 | 2,204 |
| 95% interval width, cumulative | 184,591 | 175,547 | 176,674 |
| Full-period cumulative WIS | 26,338 | 26,164 | 26,133 |

**At EQUAL mean introduction pressure, episodic clustering (within the
tested p=0.10-0.25 range) changes outcomes only marginally** -- entry
infectiousness, P(major epidemic), cumulative cases, interval width
and WIS are all close across the three arms; none of the differences
look structurally decisive. The "chain active" revised definition
(a run of >=4 consecutive locally-dominant weeks ANYWHERE in
2020-2021, NOT requiring a same-year fresh import) is 100% in all
three arms -- confirming Phase 0 item 5's definitional concern was
real: the ORIGINAL establishment definition's near-0%/90%-type swings
seen in Step 7 were partly an artefact of the same-year-import
requirement, not evidence of a different underlying transmission
reality. **The non-2022-year peak check found no evidence of
"artificial outbreaks" elsewhere in the holdout from episodic
clustering** -- median and p95/p99 annual peaks in 2020, 2021, 2023-
2025 are similar or slightly LOWER under Model B than Model A.

**Interpretation, stated carefully per the task's explicit caution**:
this does NOT mean introduction-process structure is irrelevant in
general -- it means that, specifically, REDISTRIBUTING the SAME small
mean pressure into occasional MODEST bursts (burst-week mean = lambda/p
= 10 or 4 infections) does not materially change outcomes. The burst
magnitudes tested here are themselves still small relative to the
~179x gap. **Whether a rarer, much LARGER-magnitude event (more
consistent with a genuine lineage-introduction hypothesis, e.g.
tens-to-hundreds of simultaneous introductions) would behave
differently is a DISTINCT, untested hypothesis** -- this module
deliberately did not test it, since no such magnitude was prespecified
and doing so post hoc would risk exactly the "select p to reproduce
2022" behaviour the task explicitly prohibits.

## Phase 4 -- Structural and identifiability assessment (synthesis, no new computation)

**Can the data distinguish interannual transmission variability
(tau_year), external introductions, persistence of existing chains,
seasonal amplification, and reporting uncertainty?** Partially, and
unevenly. Susceptibility and seasonal/transmission-potential components
(peak R0/Reff) are well-constrained and consistent across H1/H2 (Step
6, Step 7, this module's Phase 1) -- these are NOT confounded with the
introduction question. tau_year itself is estimated with a stable
LOCATION across priors and training windows (`HMC_RELIABILITY_AUDIT.md`)
but a provisional upper tail. The genuinely UNDER-IDENTIFIED piece is
introduction-process MAGNITUDE and STRUCTURE: Step 7 showed scaling the
mean (0 to 4x) helps but falls ~2 orders of magnitude short; this
module's Phase 2-3 showed that, at a FIXED small mean, redistributing
timing/clustering does not help either. **Would a new stochastic
importation model duplicate variance currently absorbed by tau_year?**
Plausibly, in part -- both mechanisms currently compete to explain
"how much a given year deviates from the seasonal/climate-driven
baseline," and historical case data alone, during QUIET pre-resurgence
periods, may not cleanly separate "an unusually favourable tau_year
draw" from "an unusually large introduction event" -- this is an
identifiability concern structurally analogous to (and likely
entangled with) the tau_year HMC upper-tail caveat already documented,
not a new, independent problem. **Has the 2022 observation-conditioned
seed (used in the ACCEPTED, non-holdout fit) masked an inability to
explain the outbreak's initiation?** Likely yes, in the specific sense
that the accepted fit's seed mechanism supplies the initiating
infections by construction rather than requiring the model to generate
them -- this module's holdout experiments (which never use that seed)
are precisely what exposes the gap the seeded fit does not have to
confront. No single unobserved mechanism is claimed as identified here
without supporting data (per the task's explicit instruction) --
Phase 1's A/B/C/D decomposition and Phase 2-3's structural null result
jointly narrow the explanation to "introduction magnitude," without
claiming a specific biological mechanism (new lineage vs. travel
surge vs. other) is confirmed.

## Phase 5 -- Projection readiness decision

Full table: `PHASE5_projection_readiness_decision.csv`. Summary:

**A. Absolute incidence forecasting (2026-2045 baseline).** Supported:
internally consistent mechanics; well-estimated susceptibility/
transmission potential. NOT supported: reliable absolute median/point
forecasts -- the one real historical resurgence is underpredicted by
369x even under the most favourable tested importation assumption.
**Not currently defensible** without a structurally revised,
independently-justified introduction process and resolved HMC caveats.

**B. Conditional climate-scenario analysis (reference vs. SSP, labelled).**
Supported, WITH EXPLICIT LABELLING: relative contrasts computed from
the same internally-consistent, mechanistically-verified machinery.
NOT supported: that the relative contrast is free of the introduction-
process limitation, or that absolute magnitudes under any scenario are
reliable. **Restricted, clearly-labelled conditional experiments are
defensible**; unlabelled or absolute-magnitude claims are not.

**C. Vaccination strategy sensitivity (range of epidemic assumptions).**
Supported: comparing strategies WITHIN a transparent range of
epidemic-assumption scenarios (the Step 7/7.5 grid) is methodologically
sound. NOT supported: that relative vaccination effects are
automatically reliable merely because they are relative -- vaccination-
relevant quantities (timing, age/spatial structure) may be
disproportionately sensitive to the SAME introduction-process gap in
ways not tested here. **Defensible only with an explicit, disclosed
range of scenarios and no claim that any one scenario is the true
future.**

No single PASS/FAIL verdict is given, per instruction -- the generator
is simultaneously inadequate for (A) and conditionally usable for (B)
and (C), under the stated labelling and scope restrictions.

## Final scientific decision

**1. Is the H1 failure mainly missing introductions, inadequate
persistence, or insufficient amplification?**
Mainly introduction MAGNITUDE (explanation C), not amplification
capacity (A, not supported) and not primarily persistence/fadeout in
the classic sense (B, secondary at the median). The gap is ~179x in
entry infectiousness and is not closed by scaling the existing
constant rate up to 4x (Step 7) nor by redistributing the same mean
into modest bursts (this module).

**2. Does episodic introduction produce meaningful differences at
equal expected pressure?**
No, not at the tested burst magnitudes (p=0.10/0.25, burst-week means
of 10/4). P(major epidemic), entry infectiousness, cumulative cases
and WIS are all close across Model A and both Model B arms.

**3. Is there independent justification for using an episodic process?**
Not yet demonstrated by this analysis for the TESTED magnitudes.
Independent justification would require either genomic/phylogenetic
evidence of discrete reintroduction events, or testing substantially
LARGER (not merely more clustered) burst magnitudes than prespecified
here -- neither was in scope for this diagnostic pass.

**4. Can the existing posterior be retained, or would a revised joint
historical model be necessary?**
The existing tau_year/transmission posterior can be RETAINED for the
components it is well-identified for (susceptibility, seasonal/
climate transmission potential) -- these are not implicated. A revised
model would be necessary specifically for the introduction-process
component, and per Phase 4, this would need JOINT historical refitting
(not a simulation-only patch) since tau_year and introduction magnitude
likely compete to explain the same historical variance.

**5. Is absolute 2026-2045 baseline incidence projection currently
defensible?**
No.

**6. Are restricted conditional SSP experiments defensible?**
Yes, IF explicitly labelled as conditional on the current, not-yet-
validated epidemic-generator assumptions, and restricted to RELATIVE
(reference-vs-scenario) contrasts rather than absolute magnitudes.

**7. What single model improvement or additional dataset would
provide the greatest scientific value?**
Genomic/phylogenetic surveillance data capable of dating and sizing
actual introduction/reintroduction events into Ceará around 2020-2022
-- this would directly resolve the introduction-magnitude
identifiability gap (Phase 4) that is currently the dominant,
best-supported explanation for H1's retrospective failure, and would
allow a properly justified (not retrofitted) structural revision.

## Parameter assumptions, provenance, and reproducibility

- Posteriors: `ce_tau_year_training_{H1,H2}.rds`
  (`model_fits/future_year_effect_validation/`), REUSED unmodified,
  N_train=261/365, 1,600 post-warmup draws each (4 chains x 400).
- Draws/replicates: 100 posterior draws x 20 stochastic replicates per
  main arm (Phases 1-3); 10x5 smoke test (Phase 0) run first.
- Seeds: `SEED_H1_STEP75 = BASE_SEED_STEP7 + 1` (IDENTICAL to Step 7's
  own H1 baseline seed -- same sampled draws/z_year/z_season/climate/
  initial state across Model A and all Model B arms); episodic Phase 2
  uses a documented 2-uniforms-per-week shared coupling (burst
  indicator + conditional magnitude), common WITHIN Model B's own
  p-variants, NOT claimed paired against Model A (disclosed in
  `01_episodic_engine.R`).
- HMC caveats: carried forward in full from `HMC_RELIABILITY_AUDIT.md`
  and Step 6/7's own disclosures; not re-audited here.
- Monte Carlo uncertainty: 95% intervals reported throughout
  (`lo95`/`hi95` columns); Step 7's cluster-bootstrap methodology is
  the reference approach for probability-type metrics and was not
  re-run per-arm here to keep this module's compute modest, per the
  "do not automatically run extremely large simulations" instruction --
  flagged explicitly as a limitation below.

## Files produced

```
02_Script/02_ceara_pipeline/17_step7_5_projection_readiness/
  00_config_step7_5.R, 01_episodic_engine.R,
  02_phase0_audit_and_smoketest.R, 03_phase2_episodic_vs_constant.R,
  04_phase1_pre_resurgence_persistence.R, 05_phase3_structural_comparison.R,
  06_make_figures.R, 07_make_readiness_table.R

03_Output/02_ceara_pipeline/
  tables/step7_5_projection_readiness/PHASE*.csv (9 tables, including
    PHASE5_projection_readiness_decision.csv)
  results/step7_5_projection_readiness/PHASE2_structural_arms.rds,
    PHASE1_weekly_2020_2021_predictive_vs_observed.csv
  figures/step7_5_projection_readiness/FIG_S75_1-3.{pdf,png@600dpi}
  reports/step7_5_projection_readiness/STEP7_5_PROJECTION_READINESS.md (this file)
```

**No new Stan fitting was performed.** Reused: H1/H2 training-only
posteriors; Step 7's lambda=1 baseline simulation (Model A, identical
seed). New: 2 full-scale episodic-process simulations (Model B,
p=0.10/0.25), both validated by a prior smoke test.

## Limitations not to lose track of

- Probability-type outcomes (P(major), P(reaches observed)) in this
  module's own tables do not carry a per-arm cluster-bootstrap interval
  (unlike Step 7's Phase 6) -- Step 7 already demonstrated these can
  carry substantial between-posterior-draw-subsample uncertainty; the
  close A-vs-B similarity reported in Phase 3 should be read as "no
  large, structurally obvious difference," not as a precisely-bounded
  null result.
- The episodic-process test is necessarily limited to the two
  prespecified, modest burst magnitudes; it does not and cannot rule
  out a qualitatively different, much-larger-magnitude reintroduction
  hypothesis.
- All findings remain PROVISIONAL pending tau_year HMC resolution.

STOP after reporting, per instruction. No revised reference generator
implemented; no SSP projection started; nothing committed.
