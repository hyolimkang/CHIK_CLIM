# Rio de Janeiro -- U10 Case-Only External Validation

**No refit performed.** Uses the existing HMC-PASS global-q case-only
posterior (`outputs/caseonly/rj_global_q_case_only.rds`, adapt_delta=0.98)
and the 5 already-fitted targeted fixed-q posteriors. STOP for review
before any serology-augmented HMC, per instruction.

## Section 1: case-only statewide immune fraction over the U10 window

U10 (Perisse 2020) collection window mapped to the model's weekly index:
**2018-07-01 to 2018-10-28** (18 model weeks, matching the reported
July-October 2018 window).

| Quantity | Value |
|---|---|
| p_state_U10 posterior median | **0.1216** |
| 50% CrI | [0.1141, 0.1299] |
| 95% CrI | [0.0990, 0.1458] |
| U10 observed weighted prevalence | 0.180 |
| U10 observed 95% CI | [0.148, 0.212] |
| (U10 raw prevalence, for reference) | 371/2120 = 0.175 |

## Section 2: compatibility quantification

| Quantity | Value |
|---|---|
| Posterior probability p_state_U10 in [0.148, 0.212] | **0.0177** |
| Difference (observed - posterior), median [95% CrI] | 0.0584 [0.0342, 0.0810] |
| 95% CrI overlap (posterior vs observed) | **FALSE** (posterior upper 0.1458 vs observed lower 0.148 -- misses by 0.002, i.e. essentially touching) |
| Ratio (observed / posterior median) | 1.48 |

### Fixed-q comparison at the same U10 window

| q (fixed) | p_state_U10 median | 95% CrI | Posterior probability in observed CI |
|---|---|---|---|
| 0.0100 | 0.1667 | [0.1482, 0.1893] | **0.976** |
| 0.0125 | 0.1419 | [0.1256, 0.1605] | 0.256 |
| 0.0150 (~global-q median) | 0.1215 | [0.1078, 0.1392] | 0.002 |
| 0.0175 | 0.1054 | [0.0935, 0.1201] | 0.000 |
| 0.0200 | 0.0925 | [0.0823, 0.1057] | 0.000 |

**U10 independently favours a LOWER q than the global-q posterior's own
median (0.0150).** q=0.0100 -- notably, the one fixed-q value in the
targeted sweep that achieved a completely clean HMC PASS on its own --
matches the observed U10 CI almost exactly (97.6% posterior probability
of falling inside it). This is an internally coherent signal, not a
coincidence confined to one diagnostic: the same q value is favoured by
both the HMC-cleanliness check (Section on fixed-q sweep) and this
independent serology comparison.

Full table: `RJ_U10_fixedq_comparison.csv`.

## Section 3: city vs. state (not conflated)

U10 is a Rio de Janeiro **CITY** survey; the case-only posterior evaluated
here is the **STATEWIDE** immune fraction. No state=city equivalence is
assumed. The observed point estimate sits ABOVE the statewide posterior
(18.0% vs 12.2%) -- exactly the direction expected if the capital city
(dense, historically the epicentre of RJ's chikungunya introduction) has
higher cumulative transmission than the state average, which includes
many lower-density interior/coastal municipalities. This is a
biologically coherent direction of discrepancy, not an unexplained
anomaly. **q has NOT been updated based on this comparison.**

## Section 4: HEMORIO (U35) -- NOT used

Recorded only as auxiliary/unusable evidence: blood-donor selection bias,
2019-2022 pooled sampling (no within-period resolution), and the primary
CHIKV endpoint definition and exact numerator are not recoverable from the
source (`HOLD_ENDPOINT_DEFINITION` flag, per `RJ_SEROLOGY_AUDIT.md`). Not
used in this or any quantitative comparison.

## Section 5: Decision

**B. MODERATE GEOGRAPHIC DISCREPANCY.**

- The 95% CrIs are extremely close (posterior upper bound 0.1458 vs
  observed lower bound 0.148 -- a gap of only 0.002, i.e. they very
  nearly touch) rather than being separated by an order of magnitude.
- The ratio (1.48x) is well within what city-vs-state urban/rural
  heterogeneity could plausibly explain -- not classification C (major
  conflict), which would require a difference too large to attribute to
  geography.
- But the posterior probability of falling inside the observed CI is low
  (1.8%) and the point estimates differ meaningfully (18.0% vs 12.2%) --
  not classification A (compatible), which would require substantial
  probability mass overlap.

This also is NOT a symmetric "case-only q is fine" result: the SAME
external check reveals that the case-only model's own preferred q
(~0.015) is somewhat too high relative to what even the CITY-level
(upper-bound-biased) serology would suggest is plausible statewide --
reinforcing, from a completely independent source, the q/S identification
concern already flagged from `cor(logit_q, alpha_R) ~= -0.84`.

## Section 6/7 status

Per instruction, Section 6 (compatible -> retain case-only as leading
candidate) and Section 7 (not compatible -> build a Rio CITY auxiliary
model first) are both written for the two extreme outcomes; this result
(B, moderate discrepancy) sits between them. **No new serology-augmented
HMC has been run.** Reporting for review, as instructed ("STOP for review
before any new serology-augmented HMC").

## Outputs

- `RJ_U10_CASEONLY_VALIDATION.md` (this file)
- `03_Output/figures/rio_de_janeiro_v4_9_replication/RJ_U10_caseonly_overlay.png`
- `03_Output/tables/rio_de_janeiro_v4_9_replication/RJ_U10_fixedq_comparison.csv`
