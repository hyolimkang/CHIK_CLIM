# v5.0b computational pilot result — STOP decision

Pilot settings: 150 warmup / 100 sampling, 4 chains, q=0.05, otherwise
identical HMC configuration to v4.9 (dense_e, adapt_delta=0.95,
max_treedepth=14). Purpose: diagnose sampling geometry before committing
to a full production run, per explicit pre-registration.

## Result

| | Warmup (150 draws) | Sampling (100 draws) |
|---|---|---|
| Chain 1 | 93% at max treedepth, 8 divergent | **100%** at max treedepth, 0 divergent |
| Chain 2 | 82% at max treedepth, 15 divergent | **100%** at max treedepth, 0 divergent |
| Chain 3 | 89% at max treedepth, 11 divergent | **100%** at max treedepth, 0 divergent |
| Chain 4 | 91% at max treedepth, 8 divergent | 98% at max treedepth, 2 divergent |

- Wall time: 2621 s (43.7 min) for 250 total iterations -> **2.62
  s/iteration** (per chain, averaged across the 4 parallel chains).
- Extrapolated full run (800 warmup + 400 sampling = 1200 iterations):
  **~12,580 s ≈ 3.5 hours**.
- BFMI: 0.88–1.29 (healthy) across all 4 chains.
- Bulk/tail ESS: reported "too low to be reliable" by rstan.
- Max Rhat: **NA** ("largest R-hat is NA, indicating chains have not
  mixed").

## Decision

**STOP, per the pre-registered rule**: "If the pilot again shows
pervasive max-treedepth behaviour or clearly pathological geometry: STOP.
Do not create another parameterisation. Freeze v4.9 as the final Ceara
model."

98–100% of post-warmup sampling transitions hitting the maximum treedepth
is pervasive, not marginal (contrast with v4.9's occasional
borderline-Rhat cases, which resolved cleanly with more iterations and
never showed treedepth saturation). Divergences were largely eliminated
during warmup only by the step-size adapting to an extremely small value
(consistent with 98–100% treedepth saturation), which is itself the
signature of a very poorly conditioned / stiff posterior geometry, not a
transient warmup artefact that more iterations would resolve. Extrapolating
to a ~3.5-hour full run would very likely still yield an unreliable
posterior (NA/undefined Rhat, degenerate ESS) at the end, matching the
explicit instruction not to chase this further.

**No full production run of v5.0b was launched.** This pilot result is
preserved as-is (chain CSVs under `chains/`) as evidence supporting the
decision to freeze v4.9. See
`29_v5_0_dynamic_R/V5_0_DYNAMIC_R_ASSESSMENT.md` for the full write-up.
