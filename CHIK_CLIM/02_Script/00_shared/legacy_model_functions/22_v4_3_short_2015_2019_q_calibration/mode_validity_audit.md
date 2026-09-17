# Mode Validity Audit — v4.3 short-period (2015–2019) q calibration

Scope: audits whether the low-q (~0.05) and high-q (~0.83) modes found in the
v4.3 weekly Fit A/Fit B posteriors are both genuine features of the
data-generating process, or artifacts of a misspecified weekly observation
likelihood. Nothing in the transmission model, q prior, kappa_sero, or
sampler was modified during this audit.

## 1. Bridge-sampling methodology check

Two bridge-sampling computations were run:

1. **Informal ("trapped chains")** — treated the 3 chains that stayed in the
   low-q basin (`modeB_lowclean`) and the 4 chains that stayed in the
   high-q basin (`modeB_high`) as if each were a sample from a properly
   truncated posterior, without declaring truncation in Stan.
2. **Rigorous (explicit partition)** — two dedicated Stan models,
   `renewal_ceara_v4_3_weak_q_truncated.stan` (`eta_q` bounded
   `<upper=logit(0.25)>`) and `..._truncated_high.stan`
   (`<lower=logit(0.25)>`), with `q_cut = 0.25` chosen from the empirically
   empty valley in the pooled Fit B canary draws (q ∈ [0.09, 0.63] had ~0
   posterior mass) **before** any mass calculation. Because a Stan
   distribution statement's `target +=` contribution is unaffected by a
   parameter's declared bounds (only the unconstraining Jacobian differs),
   `bridge_sampler()` on each truncated fit yields a correctly normalized
   estimate of ∫ over that region of the **original, untruncated** posterior.

Result: `mode_audit_partition_bridge_sampling_CORRECTED.csv`

| region | log marginal likelihood | relative posterior mass |
|---|---|---|
| q < 0.25 (low) | −1475.29 | 0.65% |
| q ≥ 0.25 (high) | −1470.26 | 99.35% |

The rigorous partition estimate (99.35% high-q) agrees closely with the
informal trapped-chains estimate (~99%), so the mismatch between these two
methods is **not** the source of the high-q mode's dominance — both
calculations agree that, taken at face value under the weekly likelihood,
the model's own posterior overwhelmingly favours the high-q region. This
number should therefore **not** be read as "high-q is correct because
bridge sampling favoured it": the remainder of this audit shows why the
model producing that posterior is itself suspect.

## 2. Epidemiological consistency of the two weekly modes

From `mode_audit_epidemiological_consistency.csv` and
`mode_audit_serology_contrast.csv`:

| | low-q mode | high-q mode |
|---|---|---|
| q (median) | 0.050 | 0.83 |
| exp(alpha_R) | 1.97 | 1.06 |
| R0/Reff growth median | Reff=1.01 | Reff=1.05 |
| frac. growth-weeks Reff>1 | 52.1% | 37.5% |
| immune fraction, 2018 (model) | 28.5% [21.3, 34.7] | 1.52% [1.2, 2.1] |
| Juazeiro observed (103/404) | 25.5% — **inside** interval | 25.5% — **far outside** interval |
| sero_pred (posterior pred., /404) | 113 [58, 179] | 4 [0, 29] |

Per the explicit instruction not to reject the high-q mode merely because
exp(alpha_R) < 1: the high-q mode is **not** rejected on that basis — it
does still produce Reff > 1 during 37.5% of growth weeks, i.e. real (if
weaker) epidemic growth is mechanically possible. The high-q mode is
instead flagged because its posterior-predictive serology distribution
(median 4/404, 97.5th percentile 29/404) makes the observed 103/404 result
**essentially impossible** under that mode, even though the Juazeiro
likelihood term is a genuine part of Fit B's posterior (i.e. Fit B's own
generative model, in the high-q region, cannot explain its own serology
data point).

## 3. Weekly observation-likelihood residual autocorrelation (Section 3)

Randomized quantile (PIT) residuals were computed for the weekly NB case
likelihood in both modes (`weekly_residual_acf.csv`,
`weekly_residual_acf.png`):

| lag | low-q ACF | high-q ACF |
|---|---|---|
| 1 | 0.858 | 0.806 |
| 2 | 0.810 | 0.759 |
| 4 | 0.696 | 0.612 |
| 8 | 0.491 | 0.464 |

Lag-1 autocorrelation of ~0.81–0.86 in **both** modes means the 313
nominally independent weekly NB likelihood terms are, in information terms,
equivalent to a far smaller number of truly independent observations
(rough AR(1)-implied effective N on the order of 25–30). This is a
pseudo-replication problem in the observation model, not a feature of a
particular mode.

## 4. Spatial plausibility of the high-q mode (independent of serology)

Using Ceará confirmed-case incidence (same `CLASSI_FIN==13` case
definition as the model's `C[t]`, entirely independent of the serosurvey
likelihood), cumulative 2015-01-04–2018-09-15 (`juazeiro_spatial_plausibility.csv`):

- Juazeiro do Norte incidence per capita: 7.29 × 10⁻⁴
- Ceará state average incidence per capita: 1.554 × 10⁻²
- Juazeiro / Ceará ratio: **0.047×** (Juazeiro reports *far fewer* cases
  per capita than the state average)
- Juazeiro's percentile rank among 184 Ceará municipalities: 23rd
  percentile (i.e. most municipalities report *more* per-capita cases than
  Juazeiro)
- 90th/95th/max percentile ratio-to-state among all municipalities: 1.34× /
  1.72× / 3.67× (max)

The high-q mode requires Juazeiro's true infection prevalence to be **~17×**
the state average (to reconcile a ~1.5% model-implied state immune fraction
with the observed 25.5% Juazeiro serology). But Juazeiro's own reported
case incidence is 0.047× the state average — the **opposite direction**,
and the required 17× ratio exceeds even the most extreme observed
municipality ratio (3.67×) by more than a factor of 4.

Quantified reporting-contrast requirement: implied relative ascertainment
of Juazeiro vs. the state average would need to be ≈ 0.05 / 17 = **0.00276**,
i.e. Juazeiro's case reporting would need to be **~362.5-fold** less
complete than the state average to reconcile a true 17× prevalence excess
with a 0.047× observed case-incidence ratio. This is treated as a
plausibility diagnostic, not exact proof (municipal reporting completeness
genuinely does vary spatially) — but a >300-fold differential ascertainment
gap, in the direction required, has no support in the case data.

## 5. 4-week block-aggregated observation likelihood (v4.4 diagnostic)

A diagnostic model (`renewal_ceara_v4_3_4week_diagnostic.stan`,
`23_v4_4_four_week_observation/`) held the entire latent weekly process
(renewal recursion, S/U bookkeeping, transmission, GI, seasonality, annual
effects, imports, q definition, Juazeiro likelihood, q prior) identical to
Fit A/Fit B, and changed **only** the case observation likelihood: 313
weekly NB terms → 78 non-overlapping 4-week block NB terms with a freshly
estimated `phi_obs_block` (not copied from any weekly `phi_obs`). The
trailing 313th week was excluded from the likelihood (documented rule,
`include_in_block=0`) while the latent recursion still ran through it.
Both fits used multimodal initialization: chains 1–2 in q≈0.03–0.08,
chains 3–4 in q≈0.70–0.90, specifically to test — not select — whether
both regions would still lead to persistent modes.

Result (`four_week_sensitivity_comparison.csv`):

| | Fit A4 (no serology) | Fit B4 (+ Juazeiro) |
|---|---|---|
| q, all 4 chains | **merged**: 0.31–0.37 (chain means) | **merged**: 0.10–0.13 (chain means) |
| q median [95% CrI] | 0.329 [0.156, 0.675] | 0.113 [0.043, 0.226] |
| exp(alpha_R) | 0.928 [0.816, 1.035] | 1.138 [0.983, 1.281] |
| immune fraction, 2018 | 5.1% [1.9, 16.8] (below observed 25.5%) | 21.6% [10.2, 32.5] (**covers** observed 25.5%) |
| sero_pred median [95% CrI] | 20 [1, 76] | 84 [26, 158] (**covers** observed 103) |
| block residual ACF, lag 1 | 0.790 | 0.817 |
| HMC gate | **FAIL** (19 divergences, max Rhat 1.057, min ESS 53.1) | FAIL, marginal (0 divergences, 0 TD hits, max Rhat 1.034, min ESS 164.7 — passes ESS, fails only the Rhat threshold) |

**The chains initialized in the high-q region (0.70–0.90) did not stay
there in either fit** — they moved down and merged with the low-q-initialized
chains into a single, much narrower region in both A4 and B4. This is
strong, direct evidence that the extreme separation between the weekly
low-q and high-q modes was substantially an artifact of the weekly
likelihood's pseudo-replication (Section 3): with only 78 nominal
observation units instead of 313, the case likelihood's ability to pull
the posterior into an extreme, serology-defying region is greatly reduced,
and with serology included (B4) the resulting single region is
well-calibrated against both the Juazeiro serosurvey and its own
posterior-predictive case totals (96.9% of blocks covered at 95%).

However, two qualifications are essential:

- **The HMC gate still fails in both A4 and B4.** A4 has genuine
  divergences (19), indicating real residual sampling difficulty, not
  merely a slow-mixing but otherwise healthy geometry. B4 is close to
  passing (0 divergences/treedepth hits) but still exceeds the Rhat
  threshold — chains 1–4 have not fully equilibrated even though they no
  longer occupy visibly separate basins.
- **Block-level residual autocorrelation is essentially unchanged**
  (lag-1 ACF 0.79–0.82) from the weekly analysis (0.81–0.86). Aggregating
  into 4-week blocks reduced the *number* of pseudo-independent
  observation units (313→78) and thereby diluted the case likelihood's
  leverage over q, but it did **not** fix the underlying temporal
  dependence in the observation error. The observation model remains
  materially misspecified with respect to conditional independence.

## Summary of evidence

- Rigorous bridge sampling (methodologically confirmed correct) says the
  weekly model's raw posterior favours high-q at ~99.35%.
- But: weekly residual ACF ~0.81–0.86 at lag 1 in both modes shows the
  weekly likelihood substantially overstates its own information content.
- But: independent spatial case-incidence data (not the serology
  likelihood) shows Juazeiro reports *fewer*, not far more, cases per
  capita than the Ceará average — directly contradicting what the high-q
  mode requires, by roughly two orders of magnitude beyond the most
  extreme municipality observed.
- Under 4-week block aggregation (diagnostic, latent process otherwise
  identical), the extreme high-q/low-q separation **disappears** — both
  multimodal-initialized fits converge to a single, narrower, intermediate
  region, and with serology included that region reproduces the observed
  Juazeiro result and case totals reasonably well.
- Neither the weekly nor the 4-week model achieves a clean HMC gate pass,
  and block-level residual autocorrelation remains as severe as the
  weekly analysis — so this 4-week result should be read as confirmation
  of the pseudo-replication mechanism, not as a validated production
  posterior for q.

**Correct framing (per explicit instruction):** Do NOT say "low-q is
correct because serology agrees" or "high-q is correct because bridge
sampling favored it." The evidence instead shows: the weekly likelihood is
seriously misspecified with respect to temporal residual dependence,
independent spatial epidemiology makes the high-q implied susceptibility
history difficult to defend, and the extreme bimodality itself substantially
dissolves once that pseudo-replication is reduced — but the resulting
4-week model is not yet a clean, well-specified replacement, since its own
residual autocorrelation remains just as severe and its HMC gate does not
pass. Mode support must be reassessed under a properly-specified
observation model, not adjudicated by comparing raw posterior mass between
the two current mis-specified fits.
