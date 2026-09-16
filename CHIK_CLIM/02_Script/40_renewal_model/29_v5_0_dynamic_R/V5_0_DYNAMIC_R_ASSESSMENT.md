# V5.0 Dynamic-R Assessment — Final Decision: Freeze v4.9

Scope: (0) confirms the marginal Rhat~1.02 seen in v4.9's q=0.05 canary
was Monte Carlo noise, not a structural problem; (1) tests whether the
Juazeiro serology likelihood is responsible for v4.9's 2017 epidemic
underprediction (it is not); (2) tests two dynamic-R0 parameterisations
(v5.0a weekly AR(1), v5.0b low-rank cyclic basis) intended to let a
regularized short-timescale departure from v4.9's smooth harmonic R0(t)
sharpen the 2017 peak. Both were computationally/geometrically
non-viable. **Decision: freeze v4.9 as the final parsimonious Ceará
long-horizon renewal model.**

## 1. Phase 0 — v4.9 longer-run HMC confirmation (q=0.05)

The original q=0.05 canary (400 warmup / 200 sampling) showed: divergences
0, max-treedepth hits 0, max Rhat ≈1.019, min bulk ESS ≈291, BFMI
0.83–0.96 across chains. Rerun with doubled iterations (800 warmup / 400
sampling), same model, same priors, same data:

| | 400/200 (original) | 800/400 (longer) |
|---|---|---|
| divergences | 0 | 0 |
| max treedepth hits | 0 | 0 |
| max Rhat | 1.019 | **1.0008 (rounds to 1.00)** |
| min bulk ESS | 291 | 570 |
| BFMI (per chain) | 0.955, 0.961, 0.825, 0.840 | 0.887, 0.947, 0.908, 0.722 |

The parameters carrying the highest Rhat in the longer run were
`z_season_year[7]` (1.0029), `beta_sin2` (1.0029), `z_season_year[10]`
(1.0029), `z_year[10]` (1.0022), `z_season_year[4]` (1.0020) — the
non-centred raw parameters of the hierarchical year-effect / seasonal-
amplitude structure, which is exactly where slow-mixing funnel-adjacent
geometry is expected in non-centred hierarchical models and is not a
concern at Rhat≈1.001–1.003.

Headline posterior summaries were stable across the iteration increase:

| | 400/200 | 800/400 |
|---|---|---|
| 2017 total (median) | 89,159 | 88,710 |
| 2017 peak (median) | 5,420 | 5,336 |
| 2022 total (median, post-seed) | 51,178 | 50,628 |
| 2022 peak (median) | 3,422 | 3,364 |
| 2018 immune fraction (median) | 0.264 | 0.2628 |

**Conclusion: the original 1.019 was Monte Carlo noise, not a structural
HMC failure.** v4.9 at q=0.05 is a clean, well-behaved fit.

## 2–3. Phase 1 — serology-likelihood ablation (diagnostic only)

`renewal_ceara_v4_9_serology_ablation.stan`: byte-identical to v4.9 with
only the `sero_pos ~ beta_binomial(...)` line removed from the model
block (all other data, priors, the 2022 seed mechanism, and generated
quantities unchanged). Fit at q=0.05, 800/400 iterations.

| | v4.9 standard (with serology) | v4.9 serology-ablation |
|---|---|---|
| HMC gate | PASS (div 0, max Rhat 1.0008, min ESS 570) | PASS (div 0, max Rhat 1.0079, min ESS 460) |
| 2017 total (median) | 88,710 (83% of 107,421) | 89,400 (83% of 107,421) |
| 2017 peak (median) | 5,336 (52% of 10,301) | 5,432 (53% of 10,301) |
| 2022 total (median, post-seed) | 50,628 (104%) | 50,692 (104%) |
| 2022 peak (median) | 3,364 (96%) | 3,332 (95%) |
| 2018 immune fraction | 0.2628 | 0.2638 |
| phi_obs (median) | 7.614 | 7.669 |
| A_year[2017] (median [95%]) | 1.509 [1.127, 2.022] | 1.489 [1.098, 1.964] |

Removing the serology likelihood entirely changes essentially nothing:
2017 total and peak ratios are identical to two significant figures, and
the model-implied 2018 immune fraction (0.264) stays close to the
observed Juazeiro value (0.255) *even without the serology term forcing
it there* — the case data and transmission dynamics alone already imply
that level of cumulative infection.

**This is CASE B: the serology anchor is NOT responsible for the 2017
underprediction.** The remaining hypothesis — that the smooth
two-harmonic R0(t) structure (even with year-specific amplitude A_year)
cannot represent the unusually sharp 2017 transmission episode — is the
one carried forward into Phases 2a/2b below. (The arithmetic check that
prompted this ablation — that ~365,000 additional infections at q=0.05
would be needed to close the 2017 gap, ~4 percentage points of the Ceará
population — is consistent with the ablation finding no serology
pressure: closing that gap would require the transmission dynamics
themselves to generate more infections, not a change in how the single
serology data point is weighted.)

## 4–5. v4.9 → v5.0 code changes and dynamic-R priors

### v5.0a (weekly AR(1) — ARCHIVED, computationally infeasible)

`renewal_ceara_v5_0a_dynamic_R_direct_serology.stan`: v4.9's full harmonic
+ A_year backbone renamed `mu_R[t]` (unchanged formula); `log_R0[t] =
mu_R[t] + delta_R[t]` where `delta_R` is a non-centred stationary AR(1)
process over all N=573 weeks (`rho_R ~ beta(8,2)` [mean 0.8], `sigma_R_dynamic
~ normal(0,0.15)` half-normal, both pre-registered), re-centred to zero
mean within each calendar year. Direct Juazeiro serology anchor retained
unchanged (justified by the Phase 1 CASE B result — no transport-offset
variant was needed).

**Result:** after ~20 minutes only ~50/1200 planned iterations had
completed; 83–86% of transitions hit max treedepth (16,383 leapfrog
steps/iteration) and 4–7 real divergences appeared per chain even during
early warmup. Projected full-run time: many hours, with no assurance of a
valid posterior at the end. **Archived as a computationally infeasible
parameterisation** (`outputs/FAILED_v5_0a_weekly_AR1_computationally_infeasible/FAILURE_SUMMARY.md`)
— explicitly NOT treated as evidence against the underlying scientific
hypothesis (limited short-timescale R0 deviation). Likely cause: the
prescribed within-year mean-centring removes `delta_R`'s level but not
its *shape*, and a highly persistent AR(1) (`rho_R` prior mean 0.8) can
still trace slow within-year curves overlapping with the harmonic+A_year
term's own shape, producing a funnel-like joint posterior across 573
sequentially correlated latent states.

### v5.0b (low-rank K=8 cyclic basis — PILOT ONLY, STOPPED)

`renewal_ceara_v5_0b_lowrank_dynamic_R.stan`: same `mu_R[t]` backbone;
`delta_R[t] = sum_k B_dynamic[t,k] * theta[year_id[t],k]`, `theta[y,k] =
tau_dynamic * z_dynamic[y,k]` (non-centred), with only Y×K = 11×8 = 88
latent deviation parameters plus one global `tau_dynamic` — no weekly
latent state at all. `B_dynamic` is a periodic piecewise-linear ("hat")
basis over K=8 equally spaced knots (≈6.5-week spacing) across the annual
cycle, built in R
(`04_build_lowrank_basis_and_prior_predictive.R`), each row summing to 1,
and **pre-centred within each calendar year in R** (each column's
within-year mean subtracted) so that *any* choice of `theta[y,]` yields a
within-year-zero-mean deviation by construction — no Stan-side centring
needed, and no possibility of duplicating `year_effect`'s role.

**Prior, pre-registered via prior-predictive simulation before any
fitting** (`tau_dynamic ~ normal(0, 0.15)` half-normal): at a knot, the
implied multiplicative R0 deviation has median 1.00, 95% interval
[0.75, 1.34]; P(|deviation| > 50%) = 1.8%; P(deviation doubles or halves
R0) = 0.19%. This favours small deviations and makes large ones rare, as
required.

**Computational pilot** (150 warmup / 100 sampling, 4 chains, q=0.05,
same HMC config as v4.9):

| | Warmup (150 draws/chain) | Sampling (100 draws/chain) |
|---|---|---|
| Chain 1 | 93% at max treedepth, 8 divergent | **100%** at max treedepth, 0 divergent |
| Chain 2 | 82% at max treedepth, 15 divergent | **100%** at max treedepth, 0 divergent |
| Chain 3 | 89% at max treedepth, 11 divergent | **100%** at max treedepth, 0 divergent |
| Chain 4 | 91% at max treedepth, 8 divergent | 98% at max treedepth, 2 divergent |

BFMI: 0.88–1.29 (healthy). Max Rhat: **NA** ("chains have not mixed").
Bulk/tail ESS: reported unreliable. Wall time: 2621 s (43.7 min) for 250
total iterations → 2.62 s/iteration; extrapolated full run (1200
iterations): **≈3.5 hours**.

98–100% of post-warmup transitions hitting the maximum treedepth is
**pervasive max-treedepth behaviour**, exactly the condition the
pre-registered stop rule names explicitly. Divergences were suppressed
during warmup only by adapting to an extremely small step size — itself
the signature of a poorly conditioned, stiff posterior, not a transient
warmup artefact that more iterations would resolve. **Per the
pre-registered rule, no full production run was launched.** This pilot
(and its diagnostics) is preserved at
`outputs/v5_0b_q0.05_pilot/` (`PILOT_RESULT_STOP_DECISION.md`).

## 6–9. 2017 / 2022 / serology / susceptibility comparison

No valid v5.0 posterior exists (v5.0a computationally infeasible, v5.0b
pilot geometrically pathological), so the only reliable comparison is
v4.9 standard vs. v4.9 serology-ablation, already reported in full in
Section 2–3 above. Both give: 2017 total ≈83% of observed, 2017 peak
≈52–53% of observed (peak week within 1 week of observed: predicted
2017-04-23 vs. observed 2017-04-30); 2022 total ≈104% of observed, 2022
peak ≈95–96% of observed (predicted peak week 2022-05-08 vs. observed
2022-05-15); 2018 immune fraction ≈26.3% vs. observed Juazeiro 25.5%.
Susceptible proportion trajectories (v4.9, q=0.05) drop from 100% to
≈74% by end of 2017 and to ≈63–64% after the 2022 recurrence — ample
remaining susceptibility at every stage, consistent with earlier findings
that susceptible exhaustion is not the limiting mechanism.

## 10. HMC diagnostics summary

| model | stage | divergences | treedepth hits | max Rhat | min ESS | verdict |
|---|---|---|---|---|---|---|
| v4.9 (q=0.05) | canary (400/200) | 0 | 0 | 1.019 | 291 | marginal, resolved by more iterations |
| v4.9 (q=0.05) | longer (800/400) | 0 | 0 | 1.001 | 570 | **PASS, clean** |
| v4.9 serology-ablation (q=0.05) | 800/400 | 0 | 0 | 1.008 | 460 | **PASS, clean** |
| v5.0a weekly AR(1) (q=0.05) | partial (~50/1200) | 4–7/chain (warmup) | 83–86% | n/a (incomplete) | n/a | **computationally infeasible** |
| v5.0b low-rank K=8 (q=0.05) | pilot (150/100) | 0–2/chain (sampling) | 98–100% (sampling) | NA | too low | **pathological geometry, STOP** |

## 11. Decision: FREEZE v4.9

Per the pre-registered final stop rule ("If v5.0b: does not materially
improve 2017, OR produces pathological HMC, OR requires effectively
unconstrained dynamic deviations, then STOP long-horizon model
development... freeze the stable v4.9 fixed-q model"): **v5.0b produced
pathological HMC in its computational pilot (criterion 2 of the stop
rule), independent of whether it would have improved the 2017 fit. No
v5.1 will be created.**

**v4.9 (`renewal_ceara_v4_9_hierarchical_seasonality.stan`) is frozen as
the final parsimonious Ceará long-horizon renewal reconstruction.**

Documented limitation, to be carried forward explicitly: *the
parsimonious model reproduces cumulative epidemic burden, susceptibility,
and the 2022 recurrence reasonably well, but smooths the exceptionally
sharp 2017 epidemic peak (posterior median weekly peak ≈52–53% of
observed, while the annual total is calibrated at ≈83% of observed).*
This was tested and is not attributable to the Juazeiro serology anchor
(Phase 1); it was not resolved by either a weekly (v5.0a) or a low-rank
(v5.0b) short-timescale R0 deviation, both of which failed on
computational/geometric grounds before a scientific verdict on 2017 could
even be reached.

## 12. Should the full fixed-q sensitivity sweep now be rerun?

**Yes.** The model structure is now frozen (v4.9), which is the
pre-registered condition for running the full sweep ("Only after the
final model structure is frozen should the full fixed-q sensitivity
sweep be rerun"). Note that q=0.05, 0.10, 0.15, and 0.20 already have
valid, HMC-passing v4.9 fits from earlier in this development sequence
(`28_v4_9_hierarchical_seasonality/outputs/`); only q=0.25 and q=0.30
remain to complete the pre-registered 0.05–0.30 grid under the frozen
v4.9 structure. This sweep should NOT assign posterior probability to q
(per standing instruction) — it exists to characterise how the same
qualitative 2017/2022 pattern (2022 well reproduced, 2017 peak
systematically ~50–55% of observed) behaves across the full plausible q
range, as was done for v4.7/v4.8.

## Files produced

- `02_Script/stan/renewal_ceara_v4_9_serology_ablation.stan` (diagnostic only)
- `02_Script/stan/renewal_ceara_v5_0a_dynamic_R_direct_serology.stan` (archived, infeasible)
- `02_Script/stan/renewal_ceara_v5_0b_lowrank_dynamic_R.stan` (pilot only, stopped)
- `28_v4_9_hierarchical_seasonality/scripts/04_fit_v4_9_serology_ablation.R`
- `29_v5_0_dynamic_R/scripts/01_fit_v5_0a.R`
- `29_v5_0_dynamic_R/scripts/02_plot_v5_0a_six_panel.R`, `03_plot_v4_9_vs_v5_0a_diagnostic.R` (unused — no valid v5.0a fit)
- `29_v5_0_dynamic_R/scripts/04_build_lowrank_basis_and_prior_predictive.R`
- `29_v5_0_dynamic_R/scripts/05_fit_v5_0b.R`
- `29_v5_0_dynamic_R/outputs/FAILED_v5_0a_weekly_AR1_computationally_infeasible/FAILURE_SUMMARY.md`
- `29_v5_0_dynamic_R/outputs/v5_0b_q0.05_pilot/PILOT_RESULT_STOP_DECISION.md`
- This document.

All v1–v4.9 model versions, scripts, and outputs are untouched. No v5.1
will be created without review.
