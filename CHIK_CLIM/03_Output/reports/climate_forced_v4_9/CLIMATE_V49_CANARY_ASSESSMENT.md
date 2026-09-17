# Climate-Forced v4.9 Canary Assessment — CE @ q=0.05

## Model

`renewal_ceara_v4_9_climate_forced_canary.stan` = frozen
`renewal_ceara_v4_9_hierarchical_seasonality.stan` + exactly one change:
`log_R0[t]` gains `+ beta_T*z_T_anom[t] + beta_P*z_P_anom[t]`, with
`beta_T, beta_P ~ Normal(0, 0.15)`. Everything else (q=0.05 fixed,
observation model, GI, demographic accounting, seed mechanism, Juazeiro
serology likelihood, all other priors) is byte-identical to the frozen
file. Fit: 4 chains, 800 warmup + 400 sampling, adapt_delta=0.95,
max_treedepth=14 — same settings that passed the baseline.

## A. HMC

| | baseline v4.9 | climate-forced canary |
|---|---|---|
| divergences | 0 | 0 |
| max treedepth hits | 0 | 0 |
| max Rhat (classic) | 1.001 | 1.010 |
| max Rhat (rank-normalised) | 1.003 | 1.010 |
| min ESS bulk (rank-norm.) | 571 | 470 |
| min ESS tail (rank-norm.) | 581 | 378 |

Both pass the repo's strict gate (Rhat ≤ 1.01) under **both** Rhat
conventions. The climate model sits right at the 1.01 boundary rather
than comfortably under it (unlike the baseline) — a real but small
cost in mixing quality from adding 2 parameters, not a pathology
(0 divergences, 0 treedepth hits, healthy ESS under either convention).

## B. Reconstruction stability

**Adding climate does NOT materially change the long-term
reconstruction.** At every one of 5 checkpoints (end-2016/18/20/22/25),
S_prop differs by at most **0.27%** (relative) between baseline and
climate-forced fits — visually indistinguishable
(`CE_canary_S_prop_comparison.png`). Key parameters move by ≤2.2%:
`alpha_R` −0.45%, `beta_sin1` +2.20%, `beta_cos1` +0.47%,
`sigma_season_year` −1.07%, `phi_obs` +0.41%, cumulative infection 2025
−0.43%. This is the central Phase-3 question ("does climate radically
change inferred susceptibility?") and the answer is clearly **no**.

## C. Climate identifiability / confounding

`cor(beta_T, beta_P) = 0.377` (moderate — expected, given the underlying
climate anomalies themselves are correlated at −0.389). Correlations of
beta_T/beta_P with the pre-existing structural parameters are all small:
`alpha_R` (0.01, −0.09), `beta_sin1` (0.03, 0.09), `beta_cos1` (0.05,
0.06), `sigma_season_year` (−0.14, −0.21), `phi_obs` (0.03, 0.04),
min-S/final-immune draws (≤0.08 in magnitude). **No strong confounding
with the existing harmonic/year structure** — consistent with the
Phase-1 orthogonality check (anomalies vs. harmonics, max \|cor\|=0.019).
The one real cross-term is beta_T–beta_P themselves, not
climate-vs-existing-structure.

## D. PPC

Essentially unchanged, marginally better: weekly 95% coverage 96.3% ->
96.9%, weekly RMSE 458.3 -> 455.5. Annual totals track within the same
ratios as baseline for every year 2015-2025 (e.g. 2017: both 0.826;
2022: 0.038 vs 0.030 above observed). PPC did not need to improve and
did not deteriorate.

## E. Climate effect

| parameter | median | 95% CrI | P(>0) |
|---|---|---|---|
| beta_T | −0.016 | [−0.071, 0.032] | 26% |
| beta_P | 0.024 | [−0.004, 0.052] | 96% |

Temperature: **no detectable effect** (CI straddles 0 both sides,
consistent with the first-epidemic climate GAM's finding of a
non-significant temperature term). Precipitation: **weak, suggestively
positive** — 96% of posterior mass is positive and the CI's lower bound
sits just barely below 0 — directionally consistent with the
first-epidemic GAM's significant positive precipitation association,
but not decisive on its own at this single-state, single-lag-spec scale.
`climate_multiplier(t)` stays in a narrow, plausible range
[0.930, 1.099] over the whole 2015-2025 series — a modest forcing, not
an extreme one.

## Recommendation: **PASS TO MULTI-STATE**

Rationale: the canary is computationally clean (0 divergences under a
standard, previously-validated iteration budget), does not destabilise
the existing long-term susceptibility reconstruction (<0.3% change at
every checkpoint), shows no material confounding with the pre-existing
seasonal/year structure, and does not degrade PPC. The climate signal
itself is modest and only partially decisive (precipitation suggestive,
temperature null) — this is an acceptable, expected outcome for a
single-state canary with a deliberately conservative prior; it is not a
disqualifying finding, since Phase 2/3's purpose was to establish a
*stable, interpretable* forcing mechanism, not to prove a strong climate
effect. The one thing to carry forward explicitly into Phase 4: the
climate model's Rhat sits at the strict gate's boundary (1.010, both
conventions) rather than comfortably under it — pooling across states in
Phase 4 should watch this metric per state and not assume it will stay
this close to the boundary as N or state count grows.

## Outputs

- `02_Script/stan/renewal_ceara_v4_9_climate_forced_canary.stan`
- `outputs/ce_canary/ce_climate_forced_canary_q0.05.rds`
- Tables: `CE_canary_HMC_comparison.csv`,
  `CE_canary_checkpoint_comparison.csv`,
  `CE_canary_parameter_comparison.csv`,
  `CE_canary_confounding_correlations.csv`,
  `CE_canary_annual_PPC_comparison.csv`, `CE_canary_PPC_summary.csv`,
  `CE_canary_climate_effect_summary.csv`,
  `CE_canary_climate_multiplier_posterior.csv`
- Figures: `CE_canary_beta_posteriors.png`,
  `CE_canary_climate_multiplier_posterior.png`,
  `CE_canary_S_prop_comparison.png`

## STOP RULE

Per instruction, this assessment stops here for review before Phase 4
(shared climate coefficients pooled across BA/RJ/MT/CE). No multi-state
joint model has been built yet.
