# Pernambuco Spatial-0-Q2: Two-q Feasibility Diagnostic

**Result: STOP before full HMC. The Recife case/serology conflict is NOT
solvable by differential ascertainment (q_Recife, q_rest) alone under the
current shared-R0(t) structure.** No Q2 Stan model was built or fit. No
production run, no vaccine projection.

## Baseline (archived, not modified)

The prior model is archived as:

**SPATIAL-0-COMMON-Q -- MODEL-CRITICISM BASELINE -- NOT FOR VACCINE
PROJECTIONS.** Its files (`outputs/spatial0/`,
`03_Output/figures/pernambuco_v4_9_replication/spatial0/`,
`03_Output/tables/.../spatial_model/PE_spatial0_*`) are unmodified by this
diagnostic.

## Method (no Stan re-fit)

This model's structure has one exploitable property: `q` is a pure
observation-scaling parameter (`expected_reported_cases = q * X`) that
does **not** feed back into the deterministic S/U/X renewal recursion at
all. So for any FIXED shared transmission intercept `alpha_R` (which fixes
R0(t), and hence X, S, U, `immune_prop` for every stratum, since R0(t) is
common to all 5 strata by design), each region's optimal ascertainment
`q` can be solved in closed form (`q_hat = sum(observed cases) /
sum(latent infections)`, the NB2/Poisson-optimal pure scale factor) --
and the U14 fit can be checked directly against `immune_prop`, which does
not depend on q at all. This allows a full feasibility scan over `alpha_R`
without any MCMC: for every candidate `alpha_R`, a single deterministic
forward pass of the exact `renewal_pernambuco_v4_9_spatial0.stan`
transformed-parameters block (reimplemented in R) gives `q_Recife`-optimal,
`q_rest`-optimal, and the implied Recife U14-window prevalence
simultaneously.

Fixed inputs: posterior-median seasonal/year parameters from the already-
fitted SPATIAL-0-COMMON-Q model (chains 2, 4, 5): `beta_sin1=0.170,
beta_cos1=0.101, beta_sin2=0.122, beta_cos2=0.060`, `year_effect` and
`A_year` at their per-year posterior medians. Fitted alpha_R was 0.063.
Grid: `alpha_R` in [-0.24, 1.56] (19 points, step 0.1), spanning well below
and above the fitted value.

## Result

| alpha_R | q_Recife (optimal) | q_rest (optimal) | Implied Recife U14 prevalence | max R0 | min S/N |
|---|---|---|---|---|---|
| 0.063 (fitted) | 0.153 | 0.044 | 2.5% | 1.89 | 0.996 |
| 0.163 | 0.076 | 0.022 | 19.5% | 2.09 | 0.906 |
| **0.20-0.22 (interpolated: crosses observed U14)** | **approx 0.04-0.05** | **approx 0.010-0.015** | **approx 37%** | approx 2.2 | approx 0.85 |
| 0.263 | 0.028 | 0.008 | 53.8% | 2.31 | 0.446 |
| 0.363 | 0.021 | 0.006 | 71.9% | 2.55 | 0.256 |

Full grid: `PE_Q2_feasibility_grid.csv`. Figure:
`PE_Q2_feasibility_diagnostic.png` (3 panels: q_Recife-optimal, q_rest-optimal,
implied U14 prevalence, all vs alpha_R, with the diagnostic bounds and the
reported U14 95% CI shaded).

**The three curves do not have a common feasible zone.** To bring the
implied Recife U14 prevalence up to the observed 37.2% requires
`alpha_R` approx 0.20-0.22 (interpolated from 19.5% at alpha_R=0.163 to
53.8% at alpha_R=0.263). But AT that same alpha_R, `q_rest`-optimal has
already fallen to approximately 0.010-0015 -- **well below the diagnostic
floor of 0.03**, and less than half of even the lower bound. `q_rest`-optimal
crosses below 0.03 already at alpha_R approx 0.09-0.10, long before U14 is
anywhere close to satisfied (implied prevalence only approx 5-10% there).

Mechanistically: R0(t) is shared across all 5 strata (Section 6, by
design). Raising the common transmission scale enough to let Recife's own
susceptible pool deplete to 37% inflates latent infections **in every
stratum**, not just Recife. Since the other 4 regions' observed case
counts do not change, the only way to keep matching them is for
`q_rest` to shrink in lock-step -- and it shrinks to implausible,
single-digit-percent-or-lower ascertainment long before Recife's serology
target is reached. `q_Recife` alone (unconstrained, allowed up to 36 at
the low end of the grid, which is itself already nonsensical since q must
be <=1) cannot compensate for this, because `q_Recife` and `q_rest` are
mechanically decoupled from each other but NOT decoupled from the shared
`alpha_R` that governs both.

No biologically extreme R0 or impossible depletion is even reached before
the ascertainment bound is violated -- the failure occurs first and
entirely through the ascertainment-plausibility channel, not through
susceptible-depletion exhaustion. (At alpha_R approx 0.26, where U14 first
exceeds its target, max R0 approx 2.3 and min S/N approx 0.45 -- both still
unremarkable; the diagnostic fails on `q_rest` well before either of
those become extreme.)

## Decision

**STOP before full HMC**, per the pre-registered Section 1 rule: "If no
such [feasible] region exists: STOP before full HMC. The problem is not
solved by differential ascertainment alone."

No Q2 Stan model (`renewal_pernambuco_v4_9_q2.stan`) was written. No 4-chain
HMC run was launched. Sections 2-13 of the Q2 specification (model
construction, survey-design-aware serology likelihood, HMC test,
computational/scientific/identification gates) are **not executed**,
since Section 1 explicitly gates all of them on this feasibility result.

**Final classification: Q2 FAIL -- case/serology conflict remains** (assessed
at the feasibility-diagnostic stage, prior to any HMC fit).

## Interpretation

The conflict identified in `PE_U14_LIKELIHOOD_AUDIT.md` is not an artefact
of the single shared `q_PE` specifically -- it is a structural consequence
of pairing (a) a SHARED transmission trajectory R0(t) across all 5 strata
with (b) Recife's serology requiring a transmission intensity that the
other 4 regions' surveillance data flatly contradict under any
ascertainment level in a defensible range. Splitting q into two
parameters relieves the SYMPTOM (a single q having to serve two
irreconcilable masters) but not the CAUSE (the shared R0(t) coupling
Recife's required transmission intensity to everyone else's case counts).
Resolving this would require either allowing regional transmission
intensity to differ (which Spatial-1's sigma_region attempt already
failed computationally, per `PE_SPATIAL1_GEOMETRY_DIAGNOSTIC.md`) or
reconsidering whether Recife's serology and the statewide surveillance
system are measuring compatible quantities at all -- both are scientific
questions beyond the scope of "try another q parameterisation," and per
the standing stop-rule discipline of this project, are not pursued
further here without explicit direction.

## Outputs

- `PE_SPATIAL0_Q2_FEASIBILITY.md` (this file)
- `03_Output/tables/pernambuco_v4_9_replication/spatial_model/PE_Q2_feasibility_grid.csv`
- `03_Output/tables/pernambuco_v4_9_replication/spatial_model/PE_Q2_feasibility_assessment.csv`
- `03_Output/figures/pernambuco_v4_9_replication/spatial0_q2/PE_Q2_feasibility_diagnostic.png`

No `PE_SPATIAL0_Q2_REPORT.md` and none of the 5 Q2-model-result figures
were created, since no Q2 model was built or fit (Section 1's STOP gate
fired before Section 2).
