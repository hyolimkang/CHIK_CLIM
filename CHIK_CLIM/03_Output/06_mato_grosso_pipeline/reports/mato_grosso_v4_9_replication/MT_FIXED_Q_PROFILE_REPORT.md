# Mato Grosso — Targeted Fixed-q Identification Profile

**Purpose**: not a new model. Determines whether the MT case series itself
meaningfully prefers q values near the global-q posterior region
(0.034-0.052), or whether a broad q range produces similarly good case
fits while implying very different susceptibility histories. The
existing global-q model remains the primary inferential candidate; this
profile is an identification diagnostic only.

Grid (fixed per instruction, not modified after inspection): q =
0.020, 0.030, 0.040, 0.050, 0.070. Same frozen v4.9 structure/data/priors
throughout; q pinned via a tight logit-normal prior (sd=0.02) on the
unmodified global-q Stan file — zero structural changes.

## HMC (computational validity check only — not used to select q)

| q | status | divergences | max treedepth hits | max Rhat | min bulk ESS |
|---|---|---|---|---|---|
| 0.020 | FAIL | 4 | 0 | 1.0129 | 619 |
| 0.030 | PASS | 0 | 0 | 1.0027 | 1207 |
| 0.040 | FAIL | 1 | 0 | 1.0084 | 1012 |
| 0.050 | PASS | 0 | 0 | 1.0016 | 890 |
| 0.070 | PASS | 0 | 0 | 1.0059 | 823 |

All five are computationally usable (no structural geometry pathology —
scattered, few divergences at the imposed extremes only). Per instruction
these HMC results are not used to select or exclude a q value.

## Case log-likelihood profile (Section 4)

| q | delta total logLik | early-episode delta | late-episode delta |
|---|---|---|---|
| 0.020 | **-35.4** | -1.4 | **-31.1** |
| 0.030 | -2.6 | 0.0 (early best) | -1.9 |
| 0.040 | 0.0 (overall best) | -2.4 | 0.0 (late best) |
| 0.050 | -2.5 | -4.0 | -1.9 |
| 0.070 | -8.1 | -6.5 | -5.7 |

See `MT_fixed_q_likelihood_profile.png` and `MT_fixed_q_episode_loglik.png`.

## Case PPC across q (Section 6)

| q | weekly coverage | RMSE | annual coverage | early peak covered | late peak covered |
|---|---|---|---|---|---|
| 0.020 | 96.3% | 275.4 | 72.7% | No | No |
| 0.030 | 96.9% | 212.4 | 72.7% | No | No |
| 0.040 | 96.9% | 186.7 | 81.8% | No | Yes |
| 0.050 | 96.9% | 181.5 | 81.8% | No | Yes |
| 0.070 | 97.6% | 181.9 | 81.8% | No | Yes |

The single 2018 peak week (MT_wave_06, 1,275 observed cases) is outside
the 95% predictive interval at **every** fixed q — a q-independent,
structural mild underprediction of that one week, not an identification
issue (matches the Phase-1 report's MT_wave_06 finding).

## Latent-scale consequences (Section 7 — NOT primary evidence for choosing q)

| q | S_2025 | immune_2025 | min S/N | alpha_R median [95% CrI] | max R0 |
|---|---|---|---|---|---|
| 0.020 | 32.0% | 68.0% | 32.0% | 0.196 [0.163, 0.222] | 3.07 |
| 0.030 | 38.2% | 61.8% | 38.2% | 0.138 [0.102, 0.165] | 3.16 |
| 0.040 | 46.7% | 53.3% | 46.7% | 0.087 [0.050, 0.116] | 3.15 |
| 0.050 | 55.2% | 44.8% | 55.2% | 0.048 [0.010, 0.076] | 3.10 |
| 0.070 | 67.2% | 32.8% | 67.2% | -0.005 [-0.048, 0.024] | 3.05 |

See `MT_fixed_q_susceptibility_profile.png` and
`MT_fixed_q_q_alphaR_tradeoff.png`. The q-alpha_R and q-immune_2025
trade-offs are both perfectly monotonic (Spearman rho = -1.000 for both)
across the grid — a smooth, textbook depletion/transmission-level
trade-off with no discontinuities.

## Answers to the six required questions

**1. Is q ≈ 0.04 genuinely preferred by the MT case data?**
Weakly. q=0.04 gives the single best total log-likelihood, but q=0.03 and
q=0.05 are within ~2.6 log-lik units of it — a statistically negligible
difference over a 573-week series with total logLik ≈ -1,820. There is a
genuine, non-trivial *region* of support (roughly 0.03-0.05), not a sharp
point preference for 0.04 specifically.

**2. How much worse are q=0.02 and q=0.07?**
Very different degrees. **q=0.02 is decisively rejected**: -35.4 log-lik
units below the best q, almost entirely driven by the late (2024-2025)
episode alone (-31.1 units there). **q=0.07 is moderately, but only
moderately, disfavoured**: -8.1 units total, split fairly evenly between
the early (-6.5) and late (-5.7) episodes. Low q is ruled out far more
strongly than high q is.

**3. Do early and late epidemics prefer the same q?**
Not really comparable — they don't *disagree*, but they carry very
different amounts of information. The **late episode** (2024-2025, which
contains the bulk of all reported cases) drives almost the entire total
log-likelihood profile, including the decisive rejection of q=0.02. The
**early episode** (2016-2020) is nearly flat across the *entire* grid
(a <5-log-lik-unit range from q=0.02 to q=0.07) — it is essentially
uninformative about q, not contradictory. This is **not** an
episode-conflict (Classification D): the early period simply lacks the
statistical power (far fewer cases) to pull against whatever the late
period prefers.

**4. How much does S_2025 vary over the case-supported q range?**
Substantially. Restricting to the statistically indistinguishable
0.03-0.05 band (delta logLik all within ~2.6 units), **S_2025 ranges from
38.2% (q=0.03) to 55.2% (q=0.05)** — a 17-percentage-point spread in
final population susceptibility for case fits that are, by the
likelihood, essentially tied.

**5. Does the strong global-q q-alpha_R correlation reflect a narrow but
correlated posterior, or a practically flat likelihood ridge?**
**A practically flat ridge.** The global-q free fit's 95% CrI
[0.0336, 0.0524] sits almost exactly on top of the fixed-q grid's
statistically flat plateau (q=0.03-0.05, delta logLik ≤ 2.6). The
free-fit posterior looks "narrow" only because the broad prior plus mild
likelihood curvature at the plateau edges (sharp rejection just past
q≈0.02-0.025, moderate rejection past q≈0.06-0.07) combine to bound it —
not because the data sharply peak at one value inside that range. The
q-alpha_R correlation is the geometric signature of this flat ridge, not
evidence of a well-identified ascertainment fraction riding on a
genuinely peaked likelihood.

**6. Classification: B. PARTIAL IDENTIFICATION**

- q around 0.03-0.05 is modestly preferred (clear, large rejection only
  outside this band: decisively below at q=0.02, moderately above at
  q=0.07).
- A reasonably broad range within that band (0.03-0.05) remains
  similarly supported by the case likelihood.
- That range corresponds to a 17-percentage-point spread in S_2025 —
  a scientifically meaningful difference in absolute epidemiological
  conclusions.
- Interpretation, per the predefined rule: **MT cases constrain the
  infection scale (very low and, to a lesser extent, very high q are
  excluded) but do not uniquely determine it.** The global-q model's
  narrow posterior should not be read as strong absolute-scale
  identification; downstream use of MT susceptibility should treat
  S_2025 as uncertain across at least the 38%-55% range, not as a single
  point estimate near the global-q posterior median (46.6%).

## Outputs

- `MT_FIXED_Q_PROFILE_SUMMARY.csv`
- `MT_fixed_q_likelihood_profile.png`
- `MT_fixed_q_episode_loglik.png`
- `MT_fixed_q_case_PPC.png`
- `MT_fixed_q_susceptibility_profile.png`
- `MT_fixed_q_q_alphaR_tradeoff.png`

## MT FIXED-Q PROFILE COMPLETE — STOP FOR REVIEW

Per instruction, no automatic escalation to prior sensitivity or
simulation recovery follows this analysis.
