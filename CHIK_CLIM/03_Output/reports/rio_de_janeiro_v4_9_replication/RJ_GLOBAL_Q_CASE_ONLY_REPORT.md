# Rio de Janeiro Phase 1 -- Global-q Case-Only Report

## Section 10: HMC gate

| Chain | Divergent (post-warmup) | Max-treedepth hits | Mean stepsize | Mean leapfrog | Mean accept | BFMI |
|---|---|---|---|---|---|---|
| 1 | 1 | 0 | 0.0168 | 253 | 0.974 | 0.839 |
| 2 | 1 | 0 | 0.0235 | 164 | 0.942 | 0.870 |
| 3 | 1 | 0 | 0.0164 | 252 | 0.972 | 0.916 |
| 4 | 1 | 0 | 0.0169 | 253 | 0.973 | 0.898 |

Overall: **4 post-warmup divergences** (0.13% of 3000 total post-warmup
draws), **0** treedepth-14 hits, **max Rhat = 1.003**, **min bulk ESS =
997**, **min tail ESS = 1050**, all BFMI healthy (0.84-0.92).

**Classification (initial run, adapt_delta=0.95): HMC NEAR-PASS.**
(`hmc_pass=FALSE` only because of the strict `divergences==0` requirement;
every other criterion is excellent.) All 4 chains completed normally with
no chain-specific collapse and no cross-chain branch separation (no
Pernambuco-like pathology).

**Rescue attempt (adapt_delta 0.95 -> 0.98, SAME 1200/750 iterations,
SAME model/priors/inits):** justified because the 4 divergences were
sparse (0.13% of post-warmup draws), with zero treedepth-14 hits and
otherwise excellent Rhat/ESS/BFMI -- the established "benign noise,
adaptation-resolvable" signature in this project, not a structural
geometry failure. Result: **0 divergences, 0 treedepth-14 hits, max
Rhat=1.006, min bulk ESS=752, BFMI 0.82-1.00 -- `hmc_pass=TRUE`, full HMC
PASS.** q posterior essentially unchanged (median 0.0150, 95% CrI
[0.0130, 0.0179] vs 0.0149 [0.0130, 0.0181] at adapt_delta=0.95) --
confirms the divergences were sampler noise, not a sign of a different
underlying posterior. Saved:
`outputs/caseonly_ad098/rj_global_q_case_only_ad098.rds`.

## Section 11: q / S identification diagnostics

| Quantity | Value |
|---|---|
| q prior median [95% CrI] | 0.105 [0.017, 0.458] |
| **q posterior median [95% CrI]** | **0.0149 [0.0130, 0.0181]** |
| Posterior/prior 95% width ratio | 0.012 (posterior 98.8% narrower) |
| cor(logit_q, alpha_R) | **-0.855** |
| cor(q, immune_2025) | **-0.864** |
| cor(q, min S/N) | **+0.863** |
| cor(q, max R0) | **-0.711** |

**The posterior is dramatically narrower than the prior, but ALL FOUR
correlation diagnostics exceed |0.7|** -- per instruction, this is NOT
interpreted as clean identification. This is the same q-alpha_R-depletion
ridge signature seen throughout this project (Ceara, Bahia, Pernambuco):
the case data pin down a strong RELATIONSHIP between q and transmission
scale/depletion, not q independently. Unlike Pernambuco, however, the
ridge here is NOT causing HMC pathology (near-clean sampling) and the
resulting q range is narrow in absolute terms -- but narrowness under a
strong ridge reflects how tightly the depletion dynamics of one large,
well-resolved epidemic (2018-2020) constrain the (q, alpha_R) combination
GIVEN the model structure, not necessarily genuine external identification
of q's absolute value.

## Checkpoints (S/N, cumulative immune fraction, R0)

| Checkpoint | S/N (median [95% CrI]) | Immune fraction (median [95% CrI]) | R0 (median) |
|---|---|---|---|
| End-2016 | 0.994 [0.993, 0.995] | 0.006 [0.005, 0.007] | 1.38 |
| End-2018 | 0.856 [0.830, 0.883] | 0.144 [0.117, 0.170] | 1.83 |
| End-2020 | 0.559 [0.510, 0.624] | 0.441 [0.376, 0.490] | 2.31 |
| End-2022 | 0.565 [0.516, 0.629] | 0.435 [0.371, 0.484] | 2.23 |
| End-2025 | 0.553 [0.502, 0.619] | 0.447 [0.381, 0.498] | 2.00 |

Nearly all cumulative depletion happens during the single massive
2018-2020 epidemic (S/N falls from 0.99 to 0.56), consistent with RJ's
case history (this one episode = 64% of the decade's total cases).
S/N and immune fraction are essentially flat 2020-2025, consistent with
the small size of the 2023-24/2024-25 recurrences relative to the now
partly-depleted pool.

## Section 12-13: case PPC and latent trajectories

Weekly and annual PPC (`RJ_global_q_case_PPC.png`, `RJ_global_q_annual_PPC.png`)
show an EXCELLENT fit -- every major wave (2016, 2017-18, 2018-19-20,
2023-24, 2024-25) is well captured in both magnitude and timing, with no
epidemic systematically sacrificed to fit another:

| Wave | Observed total | Predicted total (median) | Observed peak week | Predicted peak week (median) |
|---|---|---|---|---|
| RJ_wave_04 (2016) | 1,314 | 1,301 | 2016-04-10 | 2016-04-24 |
| RJ_wave_05 (2017) | 3,725 | 4,014 | 2017-01-22 | 2017-04-23 |
| RJ_wave_06 (2018) | 27,994 | 26,287 | 2018-04-29 | 2018-05-06 |
| RJ_wave_07 (2019-20) | 81,290 | 78,082 | 2019-05-05 | 2019-04-21 |
| RJ_wave_10 (2023-24) | 4,476 | 4,546 | 2024-04-07 | 2024-04-14 |
| RJ_wave_11 (2024-25) | 1,833 | 1,756 | 2025-02-16 | 2025-02-23 |

All predicted totals within ~8% of observed; all predicted peak weeks
within 1-13 weeks of observed (the wave_05 gap, ~13 weeks, is the largest,
still small relative to that episode's ~10-month duration). This is
markedly cleaner than any homogeneous Pernambuco case-only fit produced in
this project.

Latent trajectory figure: `RJ_global_q_trajectories.png` (labelled
explicitly as case-only model-implied, NOT validated absolute
susceptibility, per instruction).

## Section 14: identification escalation rule -- TRIGGERED

Per the pre-registered escalation criteria, a targeted fixed-q profile is
warranted because: **"q posterior is strongly correlated with alpha_R or
S"** is clearly met (|cor| 0.71-0.86 across all four checks), even though
the posterior itself is narrow and HMC is near-clean. This does NOT
warrant a large default fixed-q sweep (0.02/0.05/0.10/0.20/0.30) -- those
values are far outside the region the case-only posterior already
supports (q 95% CrI: 0.013-0.018) and would mostly explore territory the
data have already strongly disfavoured. A SMALL targeted grid centred on
the revealed region is the appropriate next step (proceeding per explicit
authorisation).

## Section 15: Phase-1 classification

**B. COMPUTATIONALLY CLEAN BUT q/S IDENTIFICATION UNCERTAIN.**

- HMC is now fully clean after the adapt_delta=0.98 rescue (0 divergences,
  0 treedepth hits, Rhat<=1.006, BFMI healthy; `cor(logit_q, alpha_R) =
  -0.839`, consistent with the original run's -0.855) -- not
  classification C (computational failure).
- Case fit is excellent across every wave -- not classification D
  (scientific case-fit failure).
- But the q posterior's narrowness coexists with strong (|cor|>0.7)
  correlation to alpha_R, immune_2025, min S/N, and max R0 -- this rules
  out classification A (clean identifiable candidate) per the explicit
  instruction not to equate a narrow posterior with identification.

**No final RJ susceptibility claim is made at this stage.**

## Outputs

- `RJ_MODEL_AUDIT.md`, `RJ_v49_input_audit.md`, `RJ_SEROLOGY_AUDIT.md`,
  `RJ_GLOBAL_Q_CASE_ONLY_REPORT.md` (this file)
- `RJ_weekly_cases_clean.csv`, `RJ_case_audit.csv`, `RJ_episode_audit.csv`,
  `rio_de_janeiro_v49_input_audit.csv`, `RJ_serology_inventory.csv`,
  `RJ_global_q_HMC.csv`, `RJ_global_q_posterior_summary.csv`,
  `RJ_global_q_checkpoints.csv`, `RJ_global_q_annual_PPC.csv`,
  `RJ_global_q_wave_peak_check.csv`
- `RJ_weekly_cases_2015_2025.png`, `RJ_episode_plot.png`,
  `RJ_global_q_case_PPC.png`, `RJ_global_q_annual_PPC.png`,
  `RJ_global_q_trajectories.png`

**"RJ GLOBAL-q CASE-ONLY PHASE COMPLETE -- STOP FOR REVIEW"**

Per explicit authorisation already given, proceeding next to a SMALL
targeted fixed-q profile (centred on q~0.010-0.020) as the escalation
rule (Section 14) permits -- not a large default sweep, and not serology
yet.
