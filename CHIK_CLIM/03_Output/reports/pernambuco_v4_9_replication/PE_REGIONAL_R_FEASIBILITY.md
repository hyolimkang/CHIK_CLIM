# Pernambuco Regional Transmission-Intercept Feasibility Test

## REGIONAL-R FEASIBILITY FAIL

**"STOP PE absolute susceptibility inference."**

No Stan/HMC model was built or fit. This is a deterministic,
profile/optimisation-based feasibility test only (no MCMC), per the
pre-registered Section 1 gate.

## Method

Model tested: `log R0[r,t] = alpha_global + delta_region[r] +
common_temporal[t]`, `sum(delta_region)=0`, ONE common `q`, common
temporal shape (year effects, harmonic seasonality, A_year) fixed at the
Spatial-0-COMMON-Q representative posterior draw. Equivalent
reparameterisation used for optimisation: 5 free per-region "net levels"
`L[r] = alpha_global + delta_region[r]` (any 5 real numbers are reachable
this way, so the zero-sum constraint restricts interpretation, not the
achievable parameter space -- confirmed algebraically and empirically
below).

**First attempt (superseded):** an unconstrained multi-start
Nelder-Mead over `(alpha_global, 4 free contrasts)` jointly minimising
case deviance + weighted serology log-likelihood found only a tiny
best-fit perturbation (`exp(delta_region)` in [0.97, 1.05]) and still
badly failed U14 (implied prevalence 5.4%). Diagnosing why revealed the
core issue: **when only asked to minimise total case deviance, the
optimiser can "cancel out" any imposed increase in Recife's transmission
level by lowering the shared `alpha_global` (and by sacrificing the
smallest-caseload region, Sertao, as a compensating sink) -- Recife's own
NET level barely moved (0.075 to 0.075 across a d_Recife sweep of 0 to
0.35), because the case likelihood is scale-invariant along the
`q`-vs-transmission-level ridge that has affected every model in this
project.**

**Corrected profile (reported below):** `L_Recife` is fixed directly (not
`delta_region[Recife]` alone), removing the cancellation route, while the
OTHER 4 regions' levels are **fully and independently** re-optimised at
every grid point (their own best case fit, unconstrained by Recife's
level). This is the most favourable possible test of the fixed-regional-
intercept hypothesis: it asks whether Recife can be pushed into a
high-transmission/high-depletion regime while every other region is
completely free to stay wherever fits its own cases best.

## Result

| L_Recife | Total case deviance | dev_Recife (Recife only) | q_common | Implied U14 prevalence | max R0 | min S/N |
|---|---|---|---|---|---|---|
| 0.063 (= original Spatial-0-COMMON-Q fit) | 68,387 | 22,872 | 0.146 | 2.5% | 1.89 | 0.906 |
| 0.163 | 233,616 | 189,232 | 0.073 | 19.5% | 2.09 | 0.802 |
| **0.213 (U14 satisfied: 36.9%)** | **355,058** | **309,852** | 0.046 | **36.9%** | 2.19 | 0.620 |
| 0.263 | 432,627 | 385,232 | 0.032 | 53.8% | 2.31 | 0.446 |
| 0.413 | 526,347 | 476,387 | 0.026 | 74.8% | 2.68 | 0.226 |

Full grid: `PE_regional_R_profile.csv` (15 points). Figures:
`PE_regional_R_U14_vs_cases.png` (U14 fit and per-region deviance vs
L_Recife, log scale), `PE_regional_R_susceptibility.png` (regional S/N at
the U14-matching point).

**At the L_Recife value needed to satisfy U14 (approx 0.213, implied
prevalence 36.9%, within the reported 95% CI [0.340, 0.404]), Recife's OWN
case-fit deviance is 309,852 -- a 13.5x deterioration from the original
fit's 22,872.** This happens even though the other 4 regions were left
completely free to re-optimise independently, and even though `q_common`
(0.046) is close to what would separately be Recife-optimal at that
transmission level (consistent with the earlier Q2 diagnostic). The
failure is not an ascertainment problem and not a spillover problem from
the other regions -- it is that **the epidemic trajectory implied by
Recife's own case series is fundamentally incompatible in scale/shape with
a transmission level high enough to produce 37% cumulative infection**,
regardless of how every other free parameter in this model class is set.

Max R0 (approx 2.2) and min S/N (approx 0.62) at the U14-matching point are
not yet extreme by this project's prior standards (R0 up to approx 5 was
tolerated as "not yet pathological" in earlier Bahia/PE diagnostics) -- the
case-fit collapse is decisive on its own, well before biological
implausibility becomes the binding constraint.

## Section 6 check: are fixed intercepts even the right kind of mechanism?

Independent of the above failure, the pre-model spatial audit
(`PE_5strata_input_audit.md`) already showed that the DOMINANT stratum
switches between eras: Recife/Metropolitana dominate the 2016 and 2021
epidemics, but **Vale do Sao Francisco e Araripe dominates the 2022
epidemic** (63% of that wave's cases) -- 2 of 5 strata each dominate at
least one wave, with the identity of the dominant region changing across
time. This is independent evidence that **even if the Recife-specific
scale problem above did not exist, a single FIXED per-region multiplier
would not represent which region drives which epidemic year** -- the data
pattern already visible in the turnover diagnostics is more consistent
with region-BY-TIME heterogeneity than with fixed regional offsets. Per
instruction, this is reported, not pursued (no region-specific year
effects were built).

## Decision

**REGIONAL-R FEASIBILITY FAIL.** No coherent region of (alpha_global,
delta_region) parameter space exists where Recife's serology, Recife's own
case series, and the other 4 regions' case series are all reasonably
reproduced under one common q and shared temporal shape. The failure is
severe (>10x case-deviance deterioration) and occurs through a case-fit
channel, not merely a mild biological-plausibility tradeoff.

Per the pre-registered stop rule:

- Do NOT proceed to region-specific year effects.
- Do NOT build the `PE Spatial-Rfixed` candidate model (Section 9's
  construction is gated on this feasibility test passing; it did not).
- STOP all further Pernambuco absolute-susceptibility model expansion:
  no five q parameters, no hierarchical q, no mobility matrices, no
  municipality-level latent models, no dynamic reporting.

**Recorded conclusion: the available Pernambuco surveillance and serology
data cannot be reconciled by a parsimonious spatial renewal model with
either regional ascertainment heterogeneity (Q2, already failed) or fixed
regional transmission heterogeneity (this test, failed). The downstream
vaccine analysis for Pernambuco must use susceptibility scenarios /
threshold analysis rather than a PE-specific absolute susceptibility
posterior.**

## Outputs

- `PE_REGIONAL_R_FEASIBILITY.md` (this file)
- `03_Output/tables/pernambuco_v4_9_replication/spatial_model/PE_regional_R_profile.csv`
  (the corrected, decisive profile; supersedes the exploratory multi-start
  result, retained as `PE_regional_R_recife_profile.csv` /
  `PE_regional_R_feasibility_result.rds` for transparency)
- `03_Output/figures/pernambuco_v4_9_replication/spatial_Rfixed_feasibility/PE_regional_R_U14_vs_cases.png`
- `03_Output/figures/pernambuco_v4_9_replication/spatial_Rfixed_feasibility/PE_regional_R_susceptibility.png`
